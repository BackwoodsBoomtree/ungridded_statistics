library(terra)
library(ncdf4)

## OCO Flags ##
# cloud_flag_abp:  Clouds: 0 - Classified clear, 1 - Classified cloudy, 2 - Not classified, all other values undefined; not used in SIF processing
# qc: Quality Flag: 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated

input_dirs  <- list.dirs("G:/OCO3/B10/extracted/africa", full.names = TRUE, recursive = FALSE)
out_tag    <- "_2019-2021_sifd"
out_dir     <- "G:/Africa/csv/ecoregions/mask_Dans/OCO3/"
years       <- c(2019:2022)
time        <- "month"
variable    <- "Daily_SIF_740nm"
filters     <- c("qc", "qc", "cloud_flag_abp", "Mode")
threshs     <- c(-1, 2, 2, 0)
direct      <- c("gt", "lt", "lt", "eq")
annual      <- FALSE # Compresses output to singular annual values; ie monthly means over many years

file_df <- function(input_dir, year, time) {
  file_list <- list.files(input_dir, pattern = "*.nc", full.names = TRUE, recursive = TRUE)
  
  if (time == "8-day") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0((year),"-12-31")), by="days")
    
    # Create data frame with column for each 8-day file list
    for (i in 1:46) {
      
      sub_dates <- dates[(i * 8 - 7):(i * 8)]
      
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  if (time == "16-day") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0((year + 1),"-12-31")), by="days")
    
    # Create data frame with column for each 16-day file list
    for (i in 1:23) {
      
      sub_dates <- dates[(i * 16 - 15):(i * 16)]
      
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  if (time == "month") {
    dates <- seq(as.Date(paste0(year,"-01-01")), as.Date(paste0(year,"-12-31")), by="days")
    
    df <- data.frame(matrix(ncol = 12, nrow = 31))
    # Create data frame with column for each month
    for (i in 1:12){
      if (i < 10) {
        m <- paste0("0", i)
      } else {
        m <- as.character(i)
      }
      
      sub_dates <- subset(dates, format.Date(dates, "%m") == m)
      sub_files <- c()
      
      for (j in 1:length(sub_dates)) {
        
        check_file <- file_list[grepl(sub_dates[j], file_list)]
        
        if (length(check_file) != 0) {
          sub_files <- c(sub_files, check_file)
        } else {
          sub_files <- c(sub_files, NA)
        }
      }
      
      # Force length to 31
      if (length(sub_files) < 31) {
        sub_files <- sub_files[1:31]
      }
      
      if (i == 1) {
        df <- cbind(sub_files)
      } else {
        df <- cbind(df, sub_files)
      }
    }
  }
  
  return(df)
}
get_ts  <- function(df_f, variable, time, filters, threshs, direct) {
  
  annual_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(annual_df) <- c("Mean", "SD", "SEM", "n", "filename")
  
  if (time == "8-day") {
    t <- 46
  } else if (time == "16-day") {
    t <- 23
  } else if (time == "month") {
    t <- 12
  }
  
  for (i in 1:t) {
    
    df_t <- df_f[, i]
    df_t <- df_t[!is.na(df_t)]
    
    if (length(df_t) != 0) {
      for (j in 1:length(df_t)) {
        nc <- nc_open(df_t[j])
        
        # Get data for this time step
        data <- data.frame(var = ncvar_get(nc, variable))
        colnames(data)[1] <- variable
        
        # Get filters for this time step
        if (!is.null(filters)) {
          for (f in 1:length(unique(filters))) {
            data <- cbind(data, f = ncvar_get(nc, unique(filters)[f]))
            colnames(data)[(f + 1)] <- unique(filters)[f]
          }
        }
        
        nc_close(nc)
        
        if (j == 1){
          ts_data <- data
        } else {
          ts_data <- rbind(ts_data, data)
        }
      }
      
      # filter the data
      if (!is.null(filters)) {
        for (f in 1:length(filters)){
          
          loc <- match(filters[f], names(ts_data)) # locates column of filter
          
          if (direct[f] == "lt"){
            ts_data <- ts_data[ts_data[, (loc)] < threshs[f],]
            message(paste0("Keeping ", filters[f], " values < ", threshs[f]))
          } else if (direct[f] == "gt"){
            ts_data <- ts_data[ts_data[, (loc)] > threshs[f],]
            message(paste0("Keeping ", filters[f], " values > ", threshs[f]))
          } else if (direct[f] == "eq"){
            ts_data <- ts_data[ts_data[, (loc)] == threshs[f],]
            message(paste0("Keeping ", filters[f], " values == ", threshs[f]))
          } else if (direct[f] == "neq"){
            ts_data <- ts_data[ts_data[, (loc)] != threshs[f],]
            message(paste0("Keeping ", filters[f], " values != ", threshs[f]))
          }
        }
      }
      
      annual_df[nrow(annual_df) + 1,] <- c(mean(ts_data[, 1], na.rm = TRUE),
                                           sd(ts_data[, 1], na.rm = TRUE), 
                                           sd(ts_data[, 1], na.rm = TRUE) / (sqrt(length(ts_data[, 1]))),
                                           length(ts_data[, 1]),
                                           basename(df_t[j]))
      
    } else {
      annual_df[nrow(annual_df) + 1,] <- c(NA, NA, NA, NA, NA)
    }
  }
  
  return(annual_df)
}


# if (annual == FALSE) {
#   # Append values from each year to each other and output
#   for (j in 1:length(years)) {
#     
#     files   <- file_df(input_files, years[j], time)
#     ts_data <- get_ts(files, variable, time, filters, threshs, direct)
#     
#     if (j == 1) {
#       ts_out <- ts_data
#     } else {
#       ts_out <- rbind(ts_out, ts_data)
#     }
#   }
#   write.csv(ts_out, paste0(out_dir, out_name, ".csv"), row.names = FALSE)
#   
# } else if (annual == TRUE) {
#   
#   # Combine annual dfs so we can calculate across all
#   for (j in 1:length(years)) {
#     year_files   <- file_df(input_files, years[j], time)
#     if (j == 1) {
#       all_files <- year_files
#     } else {
#       all_files <- rbind(all_files, year_files)
#     }
#   }
#   
#   ts_data <- get_ts(all_files, variable, time, filters, threshs, direct)
#   write.csv(ts_data, paste0(out_dir, out_name, ".csv"), row.names = FALSE)
#   
# }




# Write time series csv for each folder
for (i in 1:length(input_dirs)){
  for (j in 1:length(years)) {
    
    files   <- file_df(input_dirs[i], years[j], time)
    
    if (all(is.na(files)) == FALSE) {
      ts_data <- get_ts(files, variable, time, filters, threshs, direct)
      
      if (j == 1) {
        ts_out <- ts_data
      } else {
        ts_out <- rbind(ts_out, ts_data)
      }
    }
  }
  
  if (all(is.na(files)) == FALSE) {
    out_name   <- basename(files[!is.na(files)][1])
    out_name   <- paste0(substr(out_name, 1, nchar(out_name)-14), out_tag, ".csv")
    write.csv(ts_out, paste0(out_dir, out_name), row.names = FALSE)
    message(paste0("Saved ", out_name))
  } else {
    message(paste0("Skipped due to no files: ", input_dirs[i]))
  }
}