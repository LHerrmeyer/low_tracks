library(pacman)
pacman::p_load(tidyverse, raster, ncdf4, docstring,
               fs, maps, USAboundaries, sf, rnaturalearth)

fmt_lat_lon <- function(lats, lons){
  latlon_info <- sprintf("%d:%d:%.2f %d:%d:%.2f",
                         lons[1],
                         as.integer((lons[2]-lons[1])/lons[3]),
                         lats[3],
                         lats[1],
                         as.integer((lats[2]-lats[1])/lats[3]),
                         lats[3]
  )
  return(latlon_info)
}

download_and_clean <- function(url, lats, lons, fname,
                               loader, var_regex){
  cache_dir <- file.path(path_home_r(),".track_data")
  dir.create(cache_dir, showWarnings=FALSE)
  nc_file_path <- file.path(cache_dir, fname)
  latlon_info <- fmt_lat_lon(lats, lons)
  
  if(file.exists(nc_file_path)){
    return(loader(nc_file_path))
  }
  
  # Download model data
  grib_temp <- tempfile()
  reproj_grib_temp <- tempfile()
  download.file(url, grib_temp, mode="wb")
  
  # Filter to wanted variables and reproject to lat-lon
  txt <- system(sprintf("wgrib2 %s",grib_temp),intern=TRUE)
  var_str <- txt[lapply(txt, function(x)grepl(
    var_regex,x)) == TRUE]
  reproj_cmd <- sprintf("wgrib2 -i %s -new_grid latlon %s %s",
                        grib_temp, latlon_info, reproj_grib_temp)
  system(reproj_cmd,input=var_str,intern=TRUE)
  
  # Convert from GRIB to NetCDF, clean up temp files, and return
  system(sprintf("wgrib2 %s -netcdf %s",
                 reproj_grib_temp, nc_file_path),intern=TRUE)
  unlink(c(grib_temp, reproj_grib_temp))
  return(loader(nc_file_path))
}

get_raster.gfs <- function(my_date, my_hour){
  # https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-004-0.5-degree/analysis/201008/20100809/gfs_4_20100809_0000_000.grb2
  base_url <- "https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-004-0.5-degree/analysis/"
  fmt_str <- paste0(base_url,
    "%s/%s/",
    "gfs_4_%s_%02d00_000.grb2"
  )
  date_str <- format(my_date, "%Y%m%d")
  url <- sprintf(fmt_str,
                 format(my_date, "%Y%m"), 
                        date_str, date_str, my_hour)
  
  file_fmt_str <- "gfs_%s_%d.nc"
  fname <- sprintf(file_fmt_str, date_str, my_hour)
  
  var_regex = "GRD:500|GRD:10 m above|PRMSL"
  lons <- c(-140, -60, 0.5)
  lats <- c(20, 70, 0.5)
  
  load_raster <- function(fpath){
    #PRMSL_meansealevel, UGRD_500mb, VGRD_500mb, UGRD_10maboveground, VGRD_10maboveground 
    r <- brick(
      raster(fpath, varname="PRMSL_meansealevel"),
      raster(fpath, varname="UGRD_10maboveground"),
      raster(fpath, varname="VGRD_10maboveground"),
      raster(fpath, varname="UGRD_500mb"),
      raster(fpath, varname="VGRD_500mb")
    )
    raster_df <- as.data.frame(rasterToPoints(r))
    colnames(raster_df) <- 
      c("x","y","MSLP","U10","V10","U500","V500")
    #print(raster_df %>% head())
    raster_df <- raster_df %>%
      mutate(Wind = sqrt(U10^2 + V10^2)) %>%
      dplyr::select(-U10,-V10)
    return(raster_df)
  }
  
  # Load the data and return it
  return(download_and_clean(url, lats=lats, lons=lons,
                            fname=fname, loader=load_raster,
                            var_regex=var_regex))
}

get_raster.hrrr <- function(my_date, my_hour){
  #' Get raster from HRRR
  #' 
  #' @description Get raster data in a data.frame from
  #' the HRRR model at a specificed date an hour
  #' from the NOAA archive on AWS.
  #'
  #' @param my_date Date to retrieve model from.
  #' @param my_hour Hour to retrieve model frm.
  #'
  #' @return A data.frame of the raster grabbed from the model
  #' with MSLP, 500mb wind U and V directions, and 10m wind speed.
  
  # Create initial variables
  fmt_str <- "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.%s/conus/hrrr.t%02dz.wrfsfcf00.grib2"
  file_fmt_str <- "hrrr_%s_%d.nc"
  var_regex = "GRD:500|WIND:10|MSLMA"
  date_str <- format(my_date, "%Y%m%d")
  lons <- c(-130, -60, 0.25)
  lats <- c(20, 60, 0.25)
  url <- sprintf(fmt_str, date_str, my_hour)
  fname <- sprintf(file_fmt_str, date_str, my_hour)
  #"-130:130:0.5 25:50:0.5"

  # Create loader callback function
  load_raster <- function(fpath){
    # vars UGRD_500mb, VGRD_500mb, MSLMA_meansealevel, WIND_10maboveground
    r <- brick(
      raster(fpath, varname="MSLMA_meansealevel"),
      raster(fpath, varname="WIND_10maboveground"),
      raster(fpath, varname="UGRD_500mb"),
      raster(fpath, varname="VGRD_500mb")
    )
    raster_df <- as.data.frame(rasterToPoints(r))
    colnames(raster_df) <- c("x","y","MSLP","Wind","U500","V500")
    return(raster_df)
  }
  # Load the data and return it
  return(download_and_clean(url, lats=lats, lons=lons,
                            fname=fname, loader=load_raster,
                            var_regex=var_regex))
}

get_raster <- function(my_date, my_hour, model="hrrr"){
  if(model == "gfs"){
    return(get_raster.gfs(my_date, my_hour))
  }
  else if(model == "hrrr"){
    return(get_raster.hrrr(my_date, my_hour))
  }
}

get_track <- function(min_date, max_date, min_time = 0,
                      model = "hrrr",
                      center = NULL,
                      center_radius = 2, wind_radius = 6,
                      min_dot = -1, max_dist = 10,
                      filter_callback = NULL,
                      timestep = 6){
  #' Get track of low-pressure area
  #'
  #' @description This function gets 
  #' the track of a low-pressure system from HRRR
  #' at a specified date range.
  #'
  #' @param min_date Date. First date to get tracking data from.
  #' @param max_date Date. Last date to get tracking data from.
  #' @param min_time Numeric. Time (integer hours) to start tracking
  #' at on min_date
  #' @param model String. Which model ("hrrr" or "gfs") to use.
  #' Defaults to 'hrrr'.
  #' @param center c(lat,lon). Location that low pressure center starts at.
  #' Defaults to NULL, meaning to just find the point of lowest pressure.
  #' @param center_radius Numeric. How many degrees away from the given center
  #' to look for a low pressure minimum. Does not matter if center=NULL.
  #' @param wind_radius Numeric. How many degrees away from the cyclone
  #' to look for maximum wind speeds (Radius of maximum winds).
  #' Defaults to 6 (414 mi at the equator), which is normally sufficient
  #' for extratropical cyclones.
  #' @param min_dot Numeric (-1 to 1). The minimum dot product between the 
  #' normalized vector of the new low center location minus the old,
  #' and the normalize vector of the 500mb wind at the old location.
  #' These dot products correspond to the the cosine of the angle between
  #' these vectors. 0 means a 90 degree angle, -1 a 180 degree angle,
  #' and 1 meaning a 0 degree angle (same direction). This parameter is
  #' used to make sure a low track doesn't move the opposite
  #' direction of the 500mb wind. This isn't really very useful,
  #' just leave it at the defaults (-1), as the cyclone in early
  #' stages can move a bit in the opposite direction.
  #' @param max_dist Numeric. How many degs can a low be from last low
  #' position. Defaults to 10.
  #' @param filter_callback Function(df, lat, lon). A function to take
  #' to filter data based on raster data frame, last position lat, and last
  #' position lon. Return TRUE to keep a point, FALSE to delete it
  #' @param timestep Numeric. How many hours apart to get each model run.
  #' Defaults to 6.
  #'
  #' @return A data.frame of the locations of the low at
  #' timestep intervals, with central pressure, central
  #' pressure on a pressure scale (the original Saffir-Simpson
  #' scale), the date, time, the maximum wind speed, and
  #' that wind speed on the Saffir-Simpson scale.  
  
  # Create initial variables
  date_dist <- as.numeric(max_date-min_date)
  track_df <- data.frame(Date=as.Date(character()),
                         hour=as.numeric(),
                         lat=as.numeric(),
                         lon=as.numeric(),
                         pres=as.numeric(),
                         wind_max=as.numeric())
  last_pos <- NULL
  last_wind_vec <- NULL
  # Grab and process HRRR files for each day
  for(i in 0:date_dist){
    cur_date <- min_date + i
    hours <- seq(from=0,to=18,by=timestep)
    # If we are on the first day, use only the hours after min_time
    if(i == 0){
      hours = hours[hours >= min_time]
    }
    # Add raster data to data frame for each hour in the day
    for(hour in hours){
      raster_df <- get_raster(cur_date, hour, model)
      # Check if we are using a center, then
      # find low within acceptable range of center
      if(!is.null(center)){
        raster_df <- raster_df %>% filter(
          abs(y-center[1]) < center_radius &
          abs(x-center[2]) < center_radius
        )
        center <- NULL
      }
      # If we have a last position, use sanity checks
      # to make sure our next low position isn't too far away
      # (we don't want a global minimum on the other side
      # of the map).
      if(!is.null(last_pos)){
        normalize <- function(x){x/sqrt(sum(x^2))}
        lat_last <- as.numeric(last_pos[2])
        lon_last <- as.numeric(last_pos[1])
        # Filter low position using callback
        if(!is.null(filter_callback)){
          raster_df <- raster_df %>%
            filter(filter_callback(., lat_last, lon_last))
          #print(raster_df %>% head())
        }
        # Check if low position is within
        # an acceptable range according to degree
        # distance.
        raster_df <- raster_df %>% filter(
          abs(x - lon_last) < max_dist &
          abs(y - lat_last) < max_dist
        ) %>%
          # Check if low is within acceptable range
          # according to dot product.
          mutate(diff_x = x - lon_last, diff_y = y - lat_last) %>%
          rowwise() %>%
          mutate(dot = normalize(c(diff_x,diff_y)) %*%
                                   normalize(last_wind_vec)) %>%
          ungroup() %>%
          filter(dot > min_dot)
      }
      # Find low pressure center
      min_entry <- raster_df %>%
        filter(rank(MSLP, ties.method="first") == 1)
      xy <- as.numeric(min_entry[c(1,2)])
      pres <- min_entry[3]
      # Find max wind within radius of maximum winds
      wind_raster <- raster_df %>%
        filter((abs(x - xy[1]) < wind_radius) &
               (abs(y - xy[2]) < wind_radius))
      max_wind <- wind_raster %>% 
        mutate(Wind2 = -1 * Wind) %>%
        filter(rank(Wind2, ties.method="first") == 1) %>% 
        pull(Wind)
      max_wind <- as.numeric(max_wind)
      # Add tracking information to data frame
      track_df[nrow(track_df) + 1, ] = c(as.character(cur_date),
                                         hour,
                                         as.numeric(xy[2]),
                                         as.numeric(xy[1]),
                                         pres,
                                         max_wind)
      # Save last position and last wind vector
      # for sanity checks of next low position.
      last_pos <- xy
      last_wind_vec <- c(min_entry$U500, min_entry$V500)
    }
  }
  # We now have all the tracking data.
  # Make sure all tracking information is numeric
  track_df$hour <- as.numeric(track_df$hour)
  track_df$lat <- as.numeric(track_df$lat)
  track_df$lon <- as.numeric(track_df$lon)
  track_df$pres <- as.numeric(track_df$pres)/100
  track_df$wind_max <- as.numeric(track_df$wind_max)
  # Find low intensity using Saffir-Simpson pressure scale.
  track_df <- track_df %>% mutate(p_scale = case_when(
    pres >= 1005 ~ "TD",
    pres < 1005 & pres >= 995 ~ "TS",
    pres < 995 & pres >= 980 ~ "TC1",
    pres < 980 & pres >= 965 ~ "TC2",
    pres < 965 & pres >= 945 ~ "TC3",
    pres < 945 & pres >= 920 ~ "TC4",
    pres < 920 ~ "TC5"
  ))
  # Find intensity using Saffir-Simpson wind scale.
  track_df <- track_df %>% mutate(saffir_scale = case_when(
    wind_max >=70 ~ "TC5",
    wind_max < 32 & wind_max >= 18 ~ "TS",
    wind_max < 42 & wind_max >= 32 ~ "TC1",
    wind_max < 49 & wind_max >= 42 ~ "TC2",
    wind_max < 58 & wind_max >= 49 ~ "TC3",
    wind_max < 70 & wind_max >= 58 ~ "TC4",
    wind_max > 0 ~ "TD"
  ),
  wind_max_mph = wind_max * 2.237)
  # Add a day column
  track_df$day = as.numeric(format(track_df$Date, format="%d"))
  
  # Return the tracking data
  return(track_df)
}

# Plot the track from data frame
plot_track <- function(track_df, scale='p'){
  states <- us_states() %>%
    filter(!(state_abbr %in% c("AK","HI")))
  world <- ne_countries(returnclass='sf')
  my_plot <- ggplot() + 
    geom_sf(data=world) +
    geom_sf(data=states) +
    geom_path(data=track_df, aes(x=lon,y=lat),group=1) +
    coord_sf(xlim=c(-140,-60),ylim=c(25,70))
  if(scale == 'p'){
    my_plot <- my_plot + geom_point(data=track_df, aes(x=lon,y=lat,color=as.factor(p_scale)),size=2) 
  }
  else if(scale == 'ss'){
    my_plot <- my_plot + geom_point(data=track_df, aes(x=lon,y=lat,color=as.factor(saffir_scale)),size=2)
  }
  else{
    my_plot <- my_plot + geom_point(data=track_df, aes(x=lon,y=lat,color=as.factor(day)),size=2)
  }
  return(my_plot)
}

# tdf <- get_track(as.Date("2021-12-13"),as.Date("2021-12-18"),min_time=12, center=c(43.97,-120.59),max_dist=8,filter_callback=function(df,y,x)(df$x>(x-0.3)),center_radius=2,timestep=6)
# plot_track(tdf,scale='day')

if(FALSE){
tdf <- get_track(as.Date("2021-12-13"),as.Date("2021-12-18"),
model='gfs',min_time = 12, center=c(44, -120),
center_radius=3, max_dist = 8,
filter_callback=function(df,y,x){
  (df$x-x)>-1 &
    !(df$y<47 & between(df$x,-118,-114))
})
plot_track(tdf, scale='day')
}