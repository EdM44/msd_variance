#calculate the years to analyze for each GCM for a given warming level

library(tidyverse)

#cmip and experiment
cmip <- 'cmip6'
#scen <- 'ssp245' 
scen <- 'ssp585' 
#warming level in C relative to pre-industrial: 1850-1900
levels <- c("1.0", "1.5", "2.0", "3.0")

#GCMs for each ensemble
if(cmip == 'cmip6') {
  GCMs <- c('BCC-CSM2-MR', 'CMCC-ESM2', 'EC-Earth3-Veg-LR', 'FGOALS-g3', 'GFDL-ESM4',
            'INM-CM5-0', 'MIROC6', 'MPI-ESM1-2-HR', 'NorESM2-MM', 'UKESM1-0-LL')
} else {
  GCMs <- c('ACCESS-1.0', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CMCC-CMS',
            'CNRM-CM5', 'GFDL-CM3', 'HadGEM2-CC', 'HadGEM2-ES', 'MIROC5')
}

#path to NASA-NEX list of downscaled runs
fpath <- paste0("../")
fn <- file.path(fpath,"CIL_cmip6_runs.csv")
#fn <- file.path(fpath_nasa,"NASA-NEX_cmip6_runs.csv")
dfn <- readr::read_delim(fn, comment = '#', trim_ws = FALSE)

#path to warming input files
fpath <- paste0("../warming_levels/",cmip,"_all_ens/csv")

#name of input file
if(cmip == 'cmip6') {
  fname <- paste0(cmip,"_warming_levels_all_ens_1850_1900_grid.csv")
} else {
  fname <- paste0(cmip,"_warming_levels_all_ens_1850_1900.csv")
}

f <- file.path(fpath,fname)

df <- readr::read_delim(f, comment = '#', trim_ws = FALSE)
#remove leading spaces in column names
names(df) <- gsub(" ", "", names(df))

df <- rename(df,c('exper'='exp'))
#trim white space throughout data frame
df <- data.frame(lapply(df, trimws), stringsAsFactors = FALSE)

df.out <- data.frame()
for (GCM in GCMs) {
  #find variant and grid from nasa list
  model_row <- subset(dfn, dfn$Model==GCM)
  variant <- model_row$Variant
  #grid <- model_row$Grid
  #calculations of start and end used: 
  #start_year = int(central_year - 20 / 2)
  #end_year = int(central_year + (20 / 2 - 1))
  d <- subset(df, model == GCM & exper == scen & ensemble == variant)
  years <- d[d$warming_level %in% c(levels), "start_year"]
  if (length(years) < length(levels)) years <- c(years,rep(NA,length(levels)-length(years)))
  central_years <- as.numeric(years) + 10
  df.out <- rbind(df.out,c(GCM,central_years))
}
colnames(df.out) <- c("Model",levels)
write_csv(df.out,paste0("CIL_gcm_yearsofwarming_",scen,".csv"))
