pacman::p_load(tidyverse, raster, ncdf4, rNOMADS)
min_date <- as.Date("2021-12-13")
max_date <- as.Date("2021-12-16")
date_dist <- as.numeric(max_date-min_date)
track_df <- data.frame(Date=as.Date(character()),
                       time=as.numeric(),
                       lat=as.numeric(),
                       lon=as.numeric())
get_raster <- function(my_date){
  
}
for(i in 0:date_dist){
  cur_date <- min_date + i
  baseurl <- "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl?dir=%2Fnam."
  searchurl <- paste0(baseurl,format(cur_date,"%Y%m%d"))
  model.parameters <- ParseModelPage(searchurl)
  levels <- c("mean_sea_level")
  vars <- c("PRMSL")
  mypredlist <- model.parameters$pred[grep("f00", model.parameters$pred)]
  hours <- c(0:3)*6
  print("Here!")
  for(hour in hours){
    print(sprintf("%s:%s",cur_date,hour))
    #print(mypredlist)
    grep_str <- sprintf("t%02dz",hour)
    my.pred <- mypredlist[grep(grep_str, mypredlist)]
    model.info <- GribGrab(searchurl, my.pred, levels, vars)
    file_path <- model.info[[1]]$file.name
    latlon_info <- "-130:130:0.5 25:50:0.5"
    print("Reprojecting...")
    reproj_file <- tempfile()
    reproj_cmd <- sprintf("wgrib2 %s -new_grid latlon %s %s",
                        file_path, latlon_info, reproj_file)
    system(reproj_cmd, intern=TRUE)
    print("Converting to NetCDF...")
    netcdf_file <- tempfile()
    conv_cmd <- sprintf("wgrib2 %s -netcdf %s",
                        reproj_file, netcdf_file)
    system(conv_cmd, intern=TRUE)
    print("Loading raster...")
    myraster <- raster(netcdf_file)
    print("Filling NA values...")
    my_mean <- 101500
    myraster <- reclassify(myraster, cbind(NA,my_mean))
    xy <- xyFromCell(myraster, which.min(myraster))
    print(xy)
    print(cur_date)
    print(as.character(cur_date))
    track_df[nrow(track_df) + 1, ] = c(as.character(cur_date),
                                       hour,
                                       as.numeric(xy[2]),
                                       as.numeric(xy[1]))
    print("Removing old files...")
    unlink(reproj_file)
    unlink(netcdf_file)
    unlink(file_path)
  }
}
track_df$hour <- as.numeric(track_df$hour)
track_df$lat <- as.numeric(track_df$lat)
track_df$lon <- as.numeric(track_df$lon)
my_plot <- ggplot() + 
  geom_sf(data=states) + 
  geom_point(data=track_df, aes(x=lon,y=lat),col="red",size=2) + 
  geom_path(data=track_df, aes(x=lon,y=lat),group=1) +
  coord_sf(xlim=c(-130,-65),ylim=c(25,50))