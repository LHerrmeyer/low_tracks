library(pacman)
pacman::p_load(tidyverse, raster, ncdf4, docstring,
               fs, maps, USAboundaries, sf, rnaturalearth,
               ggrepel)

fmt_lat_lon <- function(lats, lons){
  #' Format latitude and longitude data.
  #'
  #' @description Format latitude and longitude data to use
  #' in `wgrib2` reprojection
  #' 
  #' @param lats Vector c(min, max, res). A vector of minimum
  #' and maximum latitudes, along with resolution. 
  #' @param lons Vector c(min, max, res). A vector of minimum
  #' and maximum longitudes, along with resolution.
  #' 
  #' @return A string of the lat-lon projection to use in
  #' `wgrib2`
  
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
  #' Download and format data from a specific url.
  #' 
  #' @description Download and format GRIB data from a 
  #' specific url into a data.frame containing raster data.
  #' 
  #' @param url String. A url to download GRIB data from.
  #' @param lats Vector c(min, max, res). A vector of minimum
  #' and maximum latitudes, along with resolution. 
  #' @param lons Vector c(min, max, res). A vector of minimum
  #' and maximum longitudes, along with resolution.
  #' @param fname String. The filename to save the NetCDF file as.
  #' @param loader Function(path_to_file). A function that 
  #' takes a path to a NetCDF files and returns a data.frame of
  #' the the raster data in that file
  #' @param var_regex String. A regex to filter variables
  #' in the grib file, by filtering the output of `wgrib2`
  #' 
  #' @return A data.frame of raster data downloaded
  #' and converted from the URL.
  
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
  #' Get raster from GFS
  #' 
  #' @description Get raster data in a data.frame from
  #' the GFS model at a specificed date an hour
  #' from the NOAA archive.
  #'
  #' @param my_date Date to retrieve model from.
  #' @param my_hour Hour to retrieve model frm.
  #'
  #' @return A data.frame of the raster grabbed from the model
  #' with MSLP, 500mb wind U and V directions, 10m wind speed,
  #' and geopotential heights for the 900mb, 600mb, and 300mb levels.
  
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

  # deltaZ = max Z in 500 km radius - min Z in 500 km radius
  # d(deltaZ)/d(ln p) 600-900 hPa
  # d(deltaZ)/d(ln p) 300-600 hPa
  # Variable Z represents geopotential height
  #"HGT:900 mb|HGT:600 mb|HGT:300mb"
  var_regex = "GRD:500|GRD:10 m above|PRMSL|HGT:900 mb|HGT:600 mb|HGT:300 mb"
  lons <- c(-180, 180, 0.5)
  lats <- c(-90, 90, 0.5)
  
  load_raster <- function(fpath){
    #PRMSL_meansealevel, UGRD_500mb, VGRD_500mb, UGRD_10maboveground, VGRD_10maboveground 
    r <- brick(
      raster(fpath, varname="PRMSL_meansealevel"),
      raster(fpath, varname="UGRD_10maboveground"),
      raster(fpath, varname="VGRD_10maboveground"),
      raster(fpath, varname="UGRD_500mb"),
      raster(fpath, varname="VGRD_500mb"),
      raster(fpath, varname="HGT_900mb"),
      raster(fpath, varname="HGT_600mb"),
      raster(fpath, varname="HGT_300mb")
    )
    raster_df <- as.data.frame(rasterToPoints(r))
    colnames(raster_df) <- 
      c("x","y","MSLP","U10m","V10m","U500","V500","Z900","Z600","Z300")
    raster_df <- raster_df %>%
      mutate(Wind = sqrt(U10m^2 + V10m^2)) %>%
      dplyr::select(-U10m,-V10m)
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
  #' Get raster data from weather models.
  #' 
  #' @description Get raster for a specific date and time from
  #' a specficic weather model.
  #' 
  #' @param my_date Date. Date to get raster data from.
  #' @param my_hour Integer. Hour to get raster data from.
  #' @param model String. Model to get raster data from.
  #' Options are "hrrr" and "gfs". Defaults to "hrrr".
  #' 
  #' @returns A data.frame of raster data from the model at
  #' the specified date and time.

  if(model == "gfs"){
    return(get_raster.gfs(my_date, my_hour))
  }
  else if(model == "hrrr"){
    return(get_raster.hrrr(my_date, my_hour))
  }
}
normalize <- function(x){x/sqrt(sum(x^2))}

get_track <- function(min_date, max_date, min_time = 0,
                      model = "hrrr",
                      center = NULL,
                      center_radius = 2, wind_radius = 6,
                      min_dot = -1, max_dist = 10,
                      filter_callback = NULL,
                      timestep = 6){
  #' Get track of low-pressure area.
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
  #' 
  #' @examples
  #' # Get plot of Dec 2021 North American winter storm (Winter Storm Bankston),
  #' # starting in Washington, from dates Dec 13 to Dec 18,
  #' # with the low track moving a maximum distance of 8 degrees,
  #' # using the GFS model, starting at 12Z,
  #' # and making sure that the track does not move westward
  #' # and making sure that the position is not further north than
  #' # 47 N when the latitude is between 118 W and 114 W.
  #' bankston.gfs <- get_track(as.Date("2021-12-13"),as.Date("2021-12-18"),
  #'        model='gfs',min_time = 12, center=c(44, -120),
  #'        center_radius=3, max_dist = 8,
  #'        filter_callback=function(df,y,x){
  #'          (df$x-x)>-1 &
  #'          !(df$y<47 & between(df$x,-118,-114))
  #'        })
  #'        
  
  # Create initial variables
  date_dist <- as.numeric(max_date-min_date)
  track_df <- data.frame(Date=as.Date(character()),
                         hour=as.numeric(),
                         lat=as.numeric(),
                         lon=as.numeric(),
                         pres=as.numeric(),
                         wind_max=as.numeric(),
                         VTL=as.numeric(),
                         VTU=as.numeric(),
                         B=as.numeric())
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
        lat_last <- as.numeric(last_pos[2])
        lon_last <- as.numeric(last_pos[1])
        # Filter low position using callback
        if(!is.null(filter_callback)){
          raster_df <- raster_df %>%
            filter(filter_callback(., lat_last, lon_last))
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
      # Find neg_VTL, neg_VTU, and B if applicable
      VTL <- NA
      VTU <- NA
      B <- NA
      # Find movement of storm vector to calculate B.
      mov_vec <- NULL
      if(!is.null(last_pos)){
        lat_last <- as.numeric(last_pos[2])
        lon_last <- as.numeric(last_pos[1])
        # Remember, x is lon, y is lat.
        mov_vec <- normalize(c(xy[1] - lon_last, xy[2] - lat_last))
      }
      if("Z600" %in% colnames(min_entry)){
        # Find nearby points within 500 km (~4 degs)
        nearby <- raster_df %>%
          filter(x - xy[1] < 10, y - xy[2] < 10) %>%
          filter(sqrt((x - xy[1])^2 + (y - xy[2])^2) <= 4.5) %>%
          mutate(delta_x = x - xy[1], delta_y = y - xy[2])
        # Find delta Z (zmax-zmin within the radius)
        deltaZ_300 <- max(nearby$Z300) - min(nearby$Z300)
        deltaZ_600 <- max(nearby$Z600) - min(nearby$Z600)
        deltaZ_900 <- max(nearby$Z900) - min(nearby$Z900)
        
        # Calculate B (asymmetry)
        # The right side of an ET cyclone is the warm sector
        #tdf <- get_track(as.Date("2021-10-24"),as.Date("2021-11-02"),
        #model='gfs',min_time = 0, center=c(35.02, 13.39),
        #center_radius=10, max_dist = 8, wind_radius=3)
        # Right means right of current storm motion
        # U wind is X-axis, Positive U wind is from the west
        
        # If we do not know the storm direction,
        # assume the storm direction is the same as the 500mb wind
        if(is.null(mov_vec)){
          mov_vec <- normalize(c(-min_entry$U500,-min_entry$V500))
        }

        # Get right-hand (clockwise) perpendicular vector: (x,y) -> (y, -x)
        # Eg (2, 1) -> (1, -2)
        # I probable have something backwards somewhere
        # https://stackoverflow.com/questions/4780119/2d-euclidean-vector-rotations
        right_vector <- c(mov_vec[2], -mov_vec[1])
        # Where the dot product with the right hand vector is positive, the
        # vector must also be on the right side
        nearby_directions <- nearby %>% 
          rowwise() %>%
          mutate(direction_dot = normalize(c(delta_x,delta_y)) %*%
                   right_vector) %>%
          ungroup()
        left <- nearby_directions %>% filter(direction_dot < 0)
        right <- nearby_directions %>% filter(direction_dot > 0)
        
        # Calculate the mean thicknesses on each side. A larger thickness
        # means a larger temperature.
        right_600_900_thik_mean <- right %>%
          mutate(thik = Z600 - Z900) %>% pull(thik) %>% mean()
        left_600_900_thik_mean <- left %>%
          mutate(thik = Z600 - Z900) %>% pull(thik) %>% mean()
        h <- ifelse(xy[2] >= 0, 1, -1)
        B <- h * (right_600_900_thik_mean - left_600_900_thik_mean)
        
        # Calculate Thermal Wind
        # Thermal Wind is derivative of depth intensity wrt pressure
        # Bigger deltaZ meanas bigger differnce between highest Z and lowest,
        # meaning more intense low. This intensity is also known as
        # height perturbation.
        # If intensity increases wrt pressure (meaning higher pressure = more intensity),
        # this is warm core. Lower pressures are at higher altitudes (1000mb is near surface, 300mb is near stratosphere, etc)
        # Using a central finite difference
        # Find vertical wind in lower levels
        VTL <- (deltaZ_600 - deltaZ_900) / (log(600) - log(900))
        # Find vertical wind in upper levels
        VTU <- (deltaZ_300 - deltaZ_600) / (log(300) - log(600))
      }
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
                                         max_wind,
                                         VTL,
                                         VTU,
                                         B)
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
  
  # Determine if storm is extratropical, subtropical, or tropical
  num <- 4
  track_df <-
    track_df %>%
    mutate(B = as.numeric(B),
           VTL = as.numeric(VTL),
           VTU = as.numeric(VTU)) %>%
    mutate(
      Bs = zoo::rollmean(B, k=num, fill=NA),
      VTLs = zoo::rollmean(VTL, k=num, fill=NA),
      VTUs = zoo::rollmean(VTU, k=num, fill=NA)
    ) %>%
    mutate(cyclone_type = case_when(
      VTLs < 0 & VTUs < 0 ~ "ET",
      VTLs > 0 & VTUs < 0 ~ "SU",
      VTLs > 0 & VTUs > 0 ~ "TC",
      VTLs < 0 & VTUs > 0 ~ "ET",
      TRUE ~ "ET"
    )) %>%
    mutate(cyclone_type = factor(cyclone_type, levels = c("ET","SU","TC")))
  # Make sure we don't have any false positive subtropical cyclones
  # (transitioning from tropical)
  track_df <- track_df %>% mutate(cyclone_type = as.character(cyclone_type))
  for(i in 2:nrow(track_df)){
    row <- track_df[i, ]
    last_row <- track_df[i-1, ]
    row$cyclone_type = ifelse(row$cyclone_type == "SU" & last_row$cyclone_type == "TC",
                              "TC",row$cyclone_type)
    track_df[i, ] = row
  }
  track_df <- track_df %>%
    mutate(cyclone_type = factor(cyclone_type, levels = c("ET","SU","TC")))
  # Set a column for HURDAT format of cyclone type/intensity
  # https://www.aoml.noaa.gov/hrd/data_sub/newHURDAT.html
  track_df <- track_df %>%
    mutate(hurdat_scale = case_when(
      saffir_scale %in% c("TC1","TC2","TC3","TC4","TC5") ~ "HU",
      cyclone_type == "ET" ~ "EX",
      cyclone_type == "SU" & saffir_scale == "TD" ~ "SD",
      cyclone_type == "SU" & saffir_scale == "TS" ~ "SS",
      cyclone_type == "TC" & saffir_scale == "TD" ~ "TD",
      cyclone_type == "TC" & saffir_scale == "TS" ~ "TS",
      TRUE ~ as.character(NA)
    ))
  # Return the tracking data
  return(track_df)
}

# Plot the track from data frame
plot_track <- function(track_df, scale='p', region="us"){
  #' Plot the low pressure track
  #'
  #' @description Get a ggplot object of the low track on a map
  #' of North America.
  #'
  #' @param track_df Data.frame. A data frame with low positions,
  #' generated by get_track()
  #' @param scale String. A string with which scale to use on the plot.
  #' Options are 'p' (pressure scale, see Arthur 2013), 
  #' 'ss' (Saffir-Simpson wind scale),
  #' and 'day' to show the day of the track.
  #' @param region String. Region of the world to plot in.
  #' Options are "world", "eu", and "us" (default).
  #'
  #' @return A ggplot object of the graph.
  #' 
  #' @references
  #' Arthur, W. C. & Woolf, H. M. 2013. 
  #' Assessment of tropical cyclone risk in the Pacific region: 
  #' analysis of changes in key tropical cyclone parameters. 
  #' Record 2013/23. Geoscience Australia: Canberra
  #' 
  #' @examples 
  #' # Plot path of Cyclone Apollo (2021)
  #' apollo <- get_track(as.Date("2021-10-24"),as.Date("2021-11-02"),
  #'                  model='gfs',min_time = 0, center=c(35.02, 13.39),
  #'                  center_radius=10, max_dist = 8, wind_radius=3)
  #' plot_track(apollo,region="eu",scale="ss")
  #'
  #' # Plot path of Hurricane Ida (2021)
  #' ida <- get_track(as.Date("2021-08-26"),as.Date("2021-09-05"),
  #'   model='gfs',min_time = 0, center=c(15.8, -74.8),
  #'   center_radius=3, max_dist = 6)
  #' plot_track(ida,region="us",scale="p")

  # Create the initial plot
  states <- us_states() %>%
    filter(!(state_abbr %in% c("AK","HI")))
  world <- ne_countries(returnclass='sf')
  # We set inherit.aes to false or else the geom_sf will try to find
  # the cyclone_type column where it does not exist.
  my_plot <- ggplot(track_df, aes(shape=cyclone_type)) + 
    geom_sf(data=world, inherit.aes = FALSE) +
    geom_sf(data=states, inherit.aes = FALSE) +
    geom_path(aes(x=lon,y=lat),group=1) +
    scale_shape_manual(name="", values = c(17,15,16), drop=FALSE)
  
  # Set plot limits based on region
  if(region == "us"){
    my_plot <- my_plot + coord_sf(xlim=c(-140,-60),ylim=c(25,70))
  }
  else if(region == "eu"){
    my_plot <- my_plot + coord_sf(xlim=c(-15,55),ylim=c(30,65))
  }
  else if(region == "world"){
    my_plot <- my_plot + coord_sf(xlim=c(-180,180),ylim=c(-90,90))
  }
  
  # Set points based on scale
  if(scale == 'p'){
    my_plot <- my_plot + geom_point(aes(x=lon,y=lat,color=as.factor(p_scale)),size=2) 
  }
  else if(scale == 'ss'){
    my_plot <- my_plot + geom_point(aes(x=lon,y=lat,color=as.factor(saffir_scale)),size=2)
  }
  else if(scale == 'day'){
    my_plot <- my_plot + geom_point(aes(x=lon,y=lat,color=as.factor(day)),size=2)
  }
  return(my_plot)
}


plot_phase <- function(track_df, B=FALSE, smooth=0){
  #' A function to plot a phase diagram of the cyclone.
  #' 
  #' @description Plots a phase diagram of a cyclone. The possible phase diagrams
  #' plot either a) Asymmetry (B) vs. lower level (900 hPa - 600 hPa) thermal wind 
  #' (low level warm vs cold core) (-VTL) or b) Upper level (600 hPa - 300 hPa) 
  #' thermal wind (upper level warm vs cold core) (-VTU) vs lower level
  #' (900 hPa - 600 hPa) thermal wind (lower level warm vs cold core) (-VTL).
  #' It is important to know that currently only the GFS model supports
  #' 
  #' @param track_df Data.frame. Data.frame of cyclone data, downloaded from
  #' get_track. Currently, only the GFS model can be used for phase information.
  #' @param B Logical. Whether or not to plot a B vs -VTL plot (TRUE) or a -VTU vs 
  #' -VTL plot (FALSE). Defaults to FALSE.
  #' @param smooth Numeric. How many observations (in 6 hour periods for GFS and
  #' HRRR models) to smooth the phase diagram with a moving average. Defaults to
  #' 0 (meaning to not smooth).
  #' 
  #' @return A ggplot object of the phase graph of the cyclone.
  #' 
  #' @references Hart, R. E. (2003).
  #' A Cyclone Phase Space Derived from Thermal Wind and Thermal Asymmetry, 
  #' Monthly Weather Review, 131(4), 585-616. 
  #' Retrieved Jan 23, 2022, from 
  #' https://journals.ametsoc.org/view/journals/mwre/131/4/1520-0493_2003_131_0585_acpsdf_2.0.co_2.xml
  #' 
  #' @examples
  #' # Get phase diagrams for Hurricane Ida, and smooth it with
  #' # a 24-hr moving average (or 4 observations every 6 hours)
  #' # Get the track data for Hurricane Ida
  #' ida <- get_track(as.Date("2021-08-26"),as.Date("2021-09-05"),
  #'   model='gfs',min_time = 0, center=c(15.8, -74.8),
  #'   center_radius=3, max_dist = 6)
  #' # Plot a VTU vs VTL chart.
  #' plot_phase(ida, smooth = 4)
  #' # Plot a B vs VTL chart.
  #' plot_phase(ida, B = TRUE, smooth = 4)

  if(smooth > 0){
    track_df <- track_df %>% 
      mutate(B = zoo::rollmean(B, k=smooth, fill=NA),
             VTL = zoo::rollmean(VTL, k=smooth, fill=NA),
             VTU = zoo::rollmean(VTU, k=smooth, fill=NA))
  }
  my_plot <- NULL
  if(!B){
    my_plot <- ggplot(data=track_df, aes(x=VTL,y=VTU)) +
      geom_hline(yintercept=0) +
      ylim(c(-600, 300)) +
      ylab("-VTU (600 hPa-300 hPa Thermal Wind)")
  }
  else{
    my_plot <- ggplot(data=track_df, aes(x=VTL, y=B)) +
      geom_hline(yintercept=10) +
      ylim(c(-25, 125)) +
      ylab("B (900 hPa-600 hPa Storm Relative Thickness Symmetry)")
  }
  my_plot <- my_plot +
    theme_bw() +
    geom_text_repel(data=track_df %>% filter(hour == 12),
              aes(label=day),size=4, color="black") +
    xlim(c(-600,300)) +
    geom_vline(xintercept=0) +
    geom_point(data=track_df, aes(color=as.factor(day)),size=3) +
    geom_path(data=track_df, color="red") +
    xlab("-VTL (900 hPa-600 hPa Thermal Wind)")
  return(my_plot)
}

# $ track --format hurdat2 --input /mnt/c/Users/Logan/Desktop/stat/low_tracks/apollo.txt --output ./apollo.png --extra 1 --res 2700
# to_hurdat(apollo2, "apollo.txt")
to_hurdat <- function(track_df, filename){
  # Function to format coordinates
  # Example: 45.5N or 83.2W
  fmt_crd <- function(crd, type){
    hemisphere <- "N"
    if(type == "lat"){
      hemisphere <- ifelse(crd >= 0, "N", "S")
    }
    else if(type == "lon"){
      hemisphere <- ifelse(crd >= 0, "E", "W")
    }
    return(str_c(
      format(abs(round(crd,1)), nsmall=1),
      hemisphere
    ))
  }
  len <- str_pad(nrow(track_df), 2, "left", 0)
  # Create the header
  rows <- c(sprintf("AL012021,%s,     %02d,",
                  str_pad("UNNAMED", 19, "left", " "),
                  nrow(track_df)))
  # Create the rows
  for(i in 1:nrow(track_df)){
    row <- track_df[i,]
    # Create the initial row text with the data
    row_txt <- str_c(
      format(row$Date,"%Y%m%d"), # Date
      sprintf(" %02d00",row$hour), # Time (e.g. 1200)
      "  ", # Code for special entries (landfall, intensity peak)
      str_pad(row$hurdat_scale,3,"left"), # Status on HURDAT2 Scale,
      str_pad(fmt_crd(row$lat, "lat"),6,"left"), # Latitude
      str_pad(fmt_crd(row$lon, "lon"),7,"left"), # Longitude
      str_pad(  # Wind speed in knots (converted from m/s),
        round(row$wind_max * 1.94)
        ,4,"left"),
      str_pad(round(row$pres),5,"left"), # Minimum Sea Level Pressure (hPa)
    sep=",")
    # Make sure to add the missing fields
    # These fields are mainly wind speeds in each sector
    row_txt <- str_c(
      row_txt,
      ",",
      strrep(" -999,",12)
    )
    rows <- c(rows, row_txt)
  }
  # This extra line is needed for wptc-track to recognize it
  rows <- c(rows, "AL022021,            UNNAMED,      1,")
  file_conn <- file(filename, "wb") # Open in bin mode so we can use LF not CRLF
  writeLines(rows, file_conn, sep="\n")
  close(file_conn)
}