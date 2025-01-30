# Assemble data frames to be compatible with approach of Michalek et al., GRL, 2023
# Output moved to 'input_data' for next step

library(terra)
library(tidyverse)

statistic <- "intensity"  #intensity or duration

ds_methods = c("cil", "nasa")

scens <- c("historical", "ssp585")

#ind <- paste0("C:/Users/EdMaurer/backups/projects/nicaragua/cmip6/data/out")
#ind <- paste0("Z:/zraid_data2/nicaragua/cmip6/out")

maskfile <- "../nasa_mask.nc"

models <- c('BCC-CSM2-MR', 'CMCC-ESM2', 'EC-Earth3-Veg-LR', 'GFDL-ESM4',
            'INM-CM5-0',  'MIROC6', 'MPI-ESM1-2-HR', 'NorESM2-MM')

# First get cell numbers that are on land using the NASA mask
rmask <- rast(file.path(maskfile))
df.mask <- as.data.frame(rmask, cells=TRUE, xy=TRUE, wide=FALSE, na.rm=FALSE) |>
  filter(values == 1) |>
  select(-c(layer, values))

for(scen in scens) {
  i <- 0
  alldatalist <- list()
  df.all <- NULL
  for (dateshift in c("default", "early", "late")) {
    for (ds_method in ds_methods) {
    
    if(dateshift == "default" | dateshift == "") {
      ind <- paste0("Z:/zraid_data2/nicaragua/cmip6/out")
    } else if (dateshift == "early") {
      ind <- paste0("Z:/zraid_data2/nicaragua/cmip6/out_",dateshift)
    } else if (dateshift == "late") {
      ind <- paste0("Z:/zraid_data2/nicaragua/cmip6/out_",dateshift)
    } else {
      stop("dateshift must be default, early, or late")
    }
    
      for(mod in models) {
        df.mod <- NULL
        i <- i+1
        
        r <- rast(file.path(ind,paste0(ds_method,"_",mod,"_",scen,"_",statistic,".nc")))
        
        df.mod <- as.data.frame(r, cells=TRUE, xy=TRUE, time=TRUE, wide=FALSE, na.rm=FALSE) |>
          filter(cell %in% df.mask$cell)
        df.mod$model <- mod
        df.mod$ds_method <- ds_method
        df.mod$exper <- dateshift
        
        #first time through save cell number and coordinate info
        if(i == 1) {
          # just use year 1
          yr1 <- df.mod$time[1]
          df.mod1 <- df.mod |> filter(time == yr1)
          df.cells <- data.frame(cell = df.mod1$cell, lon = df.mod1$x, lat = df.mod1$y)
          write_csv(df.cells, "cell_lon_lat.csv")
          rm(df.mod1)
        }

        #  change occasional NaN values to NA
        df.mod <- df.mod |> mutate_all(~ifelse(is.nan(.), NA, .))
        
        # Do one rudimentary check
        pctnull <- sum(is.na(df.mod$values))/length(df.mod$values)
        cat(paste(ds_method, scen, dateshift, mod, pctnull, "\n"))

        # Save to new data frame
        df.mod <- select(df.mod, -c(x, y, layer))
        alldatalist[[i]] <- df.mod
      }
    }
  }
  df.all <- bind_rows(alldatalist)
  write_csv(df.all, paste0(statistic,"_",scen,"_all.csv"))
}

