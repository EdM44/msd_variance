#####MOVE WARMING LEVELS TO INPUT_DATA DIRECTORY, SET THAT AS A VARIABLE, USE FILE.PATH
#######################################################################################
#' ---
#' title: "Code to perform uncertainty partitioning"
#' author: "Ed Maurer"
#' date: "December 20, 2024"
#' ---
# 
# Adapted from Michalek, A. (2023). Disentangling the sources of uncertainties 
# in the projection of floods risks in Iowa, HydroShare
# https://doi.org/10.4211/hs.62102d9b9bc64b5e8efd3cde8192cd18

################################################################################
# Initial Section
#       --- setup common variables, directories 
#
# Loop through each cell, doing the following tasks:
# Task 1. Apply a smoothing window to the data. (Moving Average)
#       --- 11 year window.
#
# Task 2. Set up data for use with the rest of the code. Make it in anomaly form.
#       --- Reference period set with hist_styr nyrs in initial section
#
# Task 3. Fit Loess Curves to  anomalies.
#       --- Fit per model per downscaling method per experiment
#
# Task 4. Calculate Internal Variability, V(x)
#       --- For each warming level,x, calculate variability
#       --- Requires treating each model separately, since different years are used
#
# Task 5. Calculate Model Uncertainty M(x).
#
# Task 6. Calculate Downscaling Uncertainty D(x).
#
# Task 7. Calculate Experiemnt (MSD definition) Uncertainty S(x).
#
# Task 8. Calculate total variance T(x) = V(x) + M(t) + D(x) + S(x)
#
# Task 9. Create plots of the percentages of uncertainty, if specified
#
# Task 10. Add uncertainties for cell to final list
# 
# Write output to a csv file for later analysis and plotting
#
###########Initial setup and common variables##################################
rm(list=ls(all=TRUE))
library(tidyverse)

hist_styr <- 1970
nyrs <- 30           # number of years to define climatology

# Define impact variable of interest and input files
variable <- "duration"     # intensity or duration

save_plots <- TRUE
plot_dir <- file.path(paste0("Z:/zraid_data2/nicaragua/variance_partition/plots_",variable))
#plot_dir <- file.path(paste0("./plots_",variable))
if (!dir.exists(plot_dir)) dir.create(plot_dir)

save_dir <- file.path("./data_out")
if (!dir.exists(save_dir)) dir.create(save_dir)
out_fname <- paste0(variable,"_uncertainty_all.csv")

# level of warming (Celsius above preindustrial) at which to evaluate GCM projections
# 0 indicates historic avg, not preindustrial.
warming_levels <- c("0", "1.5", "2.0", "3.0") 

scen <- 'ssp585'  #GCM scenario used to calculate warming levels
# File includes all 10 GCMs shared in both NASA-NEX and CIL ensembles
fnwarming <- paste0("CIL_gcm_yearsofwarming_",scen,".csv")
warming_level_years <- read_csv(fnwarming, col_names = TRUE, col_types = "cddd")
warming_yrs <- pivot_longer(warming_level_years, !Model, names_to = "warming_level", values_to = "midyear")

cmip6_future <- read_csv(paste0("./input_data/",variable,"_ssp585_all.csv"))
cmip6_hist <- read_csv(paste0("./input_data/",variable,"_historical_all.csv"))
cells <- unique(cmip6_hist$cell)

### Function to accumulate values - returns NA if all are NA, otherwise ignores NA
sumna <- function(x) {
  if(all(is.na(x))) NA else sum(x, na.rm = TRUE)
}

# BEGIN LOOP FOR ALL CELLS #######################################################

all_cells_list <- list()

### Run for all grid cells
for(ii in 1:length(cells)){
  ### Set up Model data for a cell
  hist_sub <- subset(cmip6_hist,cmip6_hist$cell==cells[ii])
  fut_sub <- subset(cmip6_future,cmip6_future$cell==cells[ii])
  model_data <- bind_rows(hist_sub,fut_sub) |> select(-c(cell)) |> rename(Year = time)
  
  # type <- "Flood Peak" XXX only used for plotting
  
  ###Clean variables
  rm(hist_sub,fut_sub)
  
  ### Get these to use in create of dataframes and for loops later
  models <- unique(model_data$model)
  expers <- unique(model_data$exper) 
  ds_methods <- unique(model_data$ds_method) 
  years <- unique(model_data$Year) 
  
  # Task 1 SMOOTH TIME SERIES #####################################################
  
  window <- 11 #years 

  dflist <- list()
  i <- 0
  for(d in ds_methods){
    for(s in expers){
      for(m in models){
        i <- i+1
        dflist[[i]] <- filter(model_data, ds_method == d, exper == s, model == m) |>
          mutate(smooth = ifelse(is.na(values), NA, 
                                 zoo::rollapply(values, width = window, FUN=mean, 
                                         na.rm=TRUE, fill="extend"))) |>
          select(-values) |>
          rename(values = smooth)
      }
    }
  }
  average_df <- bind_rows(dflist)

  ###Clean up variables
  rm(dflist,m,s,d)
  
  # Task 2 CALCULATE HISTORIC PERIOD ANOMALIES #####################################
  
  anom_df <- average_df
  for(d in ds_methods){
    for(s in expers){
      for(m in models){
        hist_mean <- filter(average_df, ds_method == d, exper == s, model == m) |>
          filter(Year >= hist_styr & Year <= (hist_styr + nyrs - 1)) |>
          summarise(mean(values, na.rm = TRUE)) |> as.numeric()

        idx <- which(average_df$model==m & average_df$exper==s & average_df$ds_method==d)
        anom_df$values[idx] <- average_df$values[idx] - hist_mean
        #Keep variables clean
        rm(hist_mean, idx)
      }
    }
  }

  ###Clean up variables
  rm(average_df,m,s)
  
  # Task 3 FIT LOESS CURVE TO ANOMALIES #################################################
  
  i <- 0
  for(d in ds_methods){
    for(s in expers){
      for(m in models){
        i <- i+1
       
        ###Temporary dataframe for loess fitting
        temp_df <- filter(anom_df, ds_method == d, exper == s, model == m)
        
        # Check for grid cells outside the domain or with no data (or more than 90% NA)
        invalid_loess <- FALSE
        if(all(is.na(temp_df$values)) | sum(is.na(temp_df$values)) > length(temp_df$values)*0.9){
          dflist[[i]] <- data.frame(Modeled = NA, Predicted = NA, Residuals = NA, Year = years,
                                    ds_method = d, model = m, exper = s)
          invalid_loess <- TRUE
        } else {
          ###Fit Loess Curve for each scenario and model. Default is getOption("na.action"): "na.omit"
          loess.model <- loess(values ~ Year, data = temp_df, span=0.75,
                               control = loess.control(surface = "direct")) #Use default setting
          ###Add to dataframe
          dflist[[i]] <- data.frame(Modeled = loess.model$y,
                                    Predicted = loess.model$fitted,
                                    Residuals = loess.model$residuals,
                                    Year = loess.model$x,
                                    ds_method = d,
                                    model = m,
                                    exper = s)
        }

        if(save_plots & !(invalid_loess)) {
          ##Save loess curve fit to value anomalies for diagnostics
          pdir <- file.path(plot_dir,"loess_plots")
          if (!dir.exists(pdir)) {dir.create(pdir)}
          fname <- paste0("loess_cell_", ii,"_",m,"_",s,"_",d, "_",variable,".png")
          png(filename = file.path(pdir,fname),width=3*469,height = 1.5*335)
          par(cex=1.3)
          plot(loess.model$x,loess.model$fitted,col='black',type="l",lty=1,lwd=3,ylab = variable, xlab='Year',
               ylim = c(min(temp_df$values, na.rm=TRUE),max(temp_df$values, na.rm=TRUE)))
          lines(loess.model$x,loess.model$fitted,col='black',lty=1,lwd=3)
          points(loess.model$x,loess.model$y,pch=21,col="black",bg="red")
          dev.off()
        }
        
        ###Clean up variables
        suppressWarnings(rm(temp_df,loess.model))
        
      }
    }
  }
  loess_df <- bind_rows(dflist)
  ####Clean up variables
  rm(dflist,d,s,m)
  
  # Task 4 INTERNAL VARIABILITY ####################################################
  
  # Use variance for entire nyrs period at each warming level
  # Note 0 warming is used for the historic period as a placeholder

  dflist <- list()
  i <- 0

  # Get variance of residuals for each model across all downscaling methods and experiments
  for(warm_lvl in warming_levels){
    i <- i + 1
    
    Internal.variability <- NA
    validmodels <- 0
    for(m in models){
      #get year for warming level for GCM
      if (warm_lvl == "0") {
        ywarm1 <- hist_styr
        ywarm2 <- ywarm1 + nyrs - 1
      } else {
        ymid <- filter(warming_yrs, Model == m & warming_level == warm_lvl)$midyear
        ywarm1 <- ymid - (floor(nyrs/2))
        ywarm2 <- ywarm1 + nyrs - 1
      }
      temp_var <-  var(filter(loess_df, model == m & Year >= ywarm1 & Year <= ywarm2)$Residuals, 
                       na.rm = TRUE)

      if (!is.na(temp_var)) validmodels <- validmodels + 1
      Internal.variability <- sumna(c(Internal.variability, temp_var))
      rm(temp_var)
    }

    # divide by number of models with valid (non-NA) data to apply equal weighting
    dflist[[i]] <- data.frame(Level=warm_lvl, Internal = Internal.variability/validmodels)
  }
  df.internal <- bind_rows(dflist)

  # Clean up variables
  rm(m, warm_lvl, ywarm1, ywarm2, dflist)
  
  # Task 5 MODEL VARIABILITY ##########################################################
  
  # Use specific time windows for warming levels for each model
  dflist <- list()
  loess_df_warm_lvl_list <- list()
  
  i <- 0
  for(warm_lvl in warming_levels) {
    i <- i + 1
    
    dflist_loess_subset = list()
    
    j <- 0
    for(m in models){
      j <- j + 1
      # Create a subset with all models with appropriate warming years
      if (warm_lvl == "0") {
        ywarm1 <- hist_styr
        ywarm2 <- ywarm1 + nyrs - 1
      } else {
        ymid <- filter(warming_yrs, Model == m & warming_level == warm_lvl)$midyear
        ywarm1 <- ymid - (floor(nyrs/2))
        ywarm2 <- ywarm1 + nyrs - 1
      }
      dflist_loess_subset[[j]] <- filter(loess_df, model == m & Year >= ywarm1 & Year <= ywarm2)
    }
    # Create new data frame for specified warming level i
    loess_df_warm_lvl <- bind_rows(dflist_loess_subset)
    # save the subsets by warming level years for next steps
    loess_df_warm_lvl_list[[i]] <- loess_df_warm_lvl
    
    ###Clean Variables
    rm(dflist_loess_subset)
    nvalid <- 0
    Model.variability <- NA
    
    for(d in ds_methods) {
      for(s in expers) {
        temp_var <- var(filter(loess_df_warm_lvl, ds_method == d & exper==s)$Predicted, 
                        na.rm = TRUE)
        if (!is.na(temp_var)) nvalid <- nvalid + 1
        Model.variability <- sumna(c(Model.variability, temp_var))
      }
    }
    dflist[[i]] <- data.frame(Level=warm_lvl, Model.Uncertainty = Model.variability/nvalid)
    
    ###Clean Variables
    rm(Model.variability, temp_var)
  }
  df.model.uncertainty <- bind_rows(dflist)
  
  ###Clean Variables
  rm(warm_lvl, s, d, dflist)
  
  # Task 6 DOWNSCALING VARIABILITY ######################################################

  dflist <- list()
  
  i <- 0
  for(warm_lvl in warming_levels) {
    i <- i + 1

    loess_df_warm_lvl <- loess_df_warm_lvl_list[[i]]
    
    nvalid <- 0
    Downscaling.variability <- NA
    
    for(m in models) {
      for(s in expers) {
        temp_var <- var(filter(loess_df_warm_lvl, model == m & exper==s)$Predicted, 
                        na.rm = TRUE)
        if (!is.na(temp_var)) nvalid <- nvalid + 1
        Downscaling.variability <- sumna(c(Downscaling.variability, temp_var))
      }
    }
    dflist[[i]] <- data.frame(Level=warm_lvl, Downscaling.Uncertainty = Downscaling.variability/nvalid)
    
    ###Clean Variables
    rm(Downscaling.variability, temp_var)
  }
  df.downscaling.uncertainty <- bind_rows(dflist)
  
  ###Clean Variables
  rm(warm_lvl, s, m, dflist)

  # Task 7 EXPERIMENT (IMPACT DEFINITION) UNCERTAINTY #######################################
  dflist <- list()
  
  i <- 0
  for(warm_lvl in warming_levels) {
    i <- i + 1
    
    loess_df_warm_lvl <- loess_df_warm_lvl_list[[i]]
    
    nvalid <- 0
    Experiment.variability <- NA
    
    for(m in models) {
      for(d in ds_methods) {
        temp_var <- var(filter(loess_df_warm_lvl, model == m & ds_method == d)$Predicted, 
                        na.rm = TRUE)
        if (!is.na(temp_var)) nvalid <- nvalid + 1
        Experiment.variability <- sumna(c(Experiment.variability, temp_var))
      }
    }
    dflist[[i]] <- data.frame(Level=warm_lvl, Experiment.Uncertainty = Experiment.variability/nvalid)
    
    ###Clean Variables
    rm(Experiment.variability, temp_var)
  }
  df.experiment.uncertainty <- bind_rows(dflist)
  
  ###Clean Variables
  rm(warm_lvl, d, m, dflist)
  
  # Task 8 ASSEMBLE TOTAL VARIANCE DATA FRAME ############################################
  
  uncertainty.list <- list(df.internal, df.model.uncertainty, df.downscaling.uncertainty, df.experiment.uncertainty)
  Total.Variance <- uncertainty.list |> reduce(full_join, by='Level')
  Total.Variance$Level <- as.numeric(as.character(Total.Variance$Level))

  colnames(Total.Variance) <- c("Warming_Level","Internal", "Model", "Downscaling", "Experiment")
  Total.Variance$Total <- Total.Variance$Internal + Total.Variance$Model + 
    Total.Variance$Downscaling + Total.Variance$Experiment
  Total.Variance$ID <- cells[ii]
  
  ###Clean up Variables
  rm(df.internal, df.model.uncertainty, df.downscaling.uncertainty, df.experiment.uncertainty)
  
  # Task 9 CREATE PLOTS FOR CURRENT CELL ##################################################
  
  if(save_plots) {
    
    ### Plot of Variances
    png(filename=file.path(plot_dir,paste0("var_lines_cell_",cells[ii],".png")),res=600,units = "in",width = 5,height = 4)
    p <- pivot_longer(Total.Variance, cols = 2:6, names_to = "Uncertainty", values_to = "Variance") |>
      ggplot(aes(x=Warming_Level, y=Variance, color = Uncertainty)) +
      geom_line() +
      ylab("Variance") +
      xlab("Warming Level") +
      annotate("text",  x=-Inf, y = Inf, label = paste0("Cell ",ii), vjust=1, hjust=0) +
      theme_bw()
    plot(p)
    dev.off()
    
    ###Plot a break down of total Variance 
    plot_df <- data.frame(Warming_Level=Total.Variance$Warming_Level,
                        Internal = (Total.Variance$Total-Total.Variance$Internal)/Total.Variance$Total,
                        Downscaling = (Total.Variance$Total-Total.Variance$Internal-Total.Variance$Downscaling)/Total.Variance$Total,
                        Experiment = (Total.Variance$Total-Total.Variance$Internal-Total.Variance$Downscaling-Total.Variance$Experiment)/Total.Variance$Total)

    png(filename=file.path(plot_dir,paste0("fraction_cell_",cells[ii],".png")),res=600,units = "in",width = 5,height = 4)
    p <- ggplot(plot_df,aes(x=Warming_Level)) +
      geom_ribbon(aes(ymin = Internal, ymax = 1,fill = "Internal"),color="black") +
      geom_ribbon(aes(ymin = Downscaling, ymax = Internal,fill = "Downscaling"),color="black",linewidth=1.05) +
      geom_ribbon(aes(ymin = Experiment, ymax = Downscaling,fill = "Experiment"),color="black",linewidth=1.05) +
      geom_ribbon(aes(ymin = 0, ymax = Experiment,fill = "Model"),color="black",linewidth=1.05) +
      geom_hline(yintercept = seq(0,1,0.1),color="white",alpha=0.5) +
      geom_vline(xintercept = c(0, 1.5, 2.0, 3.0),color="white",alpha=0.5) +
      scale_fill_manual(values = c("Internal"="orange3","Downscaling"="palegreen3","Model"="blue3", "Experiment"="darkred")) +
      ylab("Fraction of Total Variance") +
      xlab("Warming Level") +
      scale_x_continuous(expand=c(0,0),limits = c(0,3),breaks=c(0, 1.5, 2.0, 3.0))+
      scale_y_continuous(expand=c(0,0),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))+
      theme_bw()+
      #annotate("text",  x=-Inf, y = Inf, label = paste0("Cell ",ii), vjust=1, hjust=0) +
      theme(
        #panel.grid.minor = element_blank(),
        #legend.position="none",
        legend.position.inside=c(.87,.12),
        legend.text=element_text(size=8),
        #legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.15, "in"),
        plot.title = element_text(hjust = 0.5,size = 10),
        strip.text = element_text(size=12))
      #facet_wrap(vars(Type))
    plot(p)
    dev.off()
  }
  # Task 10 ADD CURRENT CELL TO OUTPUT LIST #############################################

  all_cells_list[[ii]] <- Total.Variance
  
  rm(Total.Variance)
  
  print(paste0("Completed cell: ",ii," of ",length(cells)))
}

df.total.uncertainty <- bind_rows(all_cells_list)
write_csv(df.total.uncertainty, file.path(save_dir,out_fname))
