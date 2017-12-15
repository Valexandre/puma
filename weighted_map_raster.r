
getMatrixElementsByLinearIndex <- function(mat, row_idx, col_idx){
  nr = nrow(mat)
  return(mat[row_idx + nr*(col_idx-1)])
}
createWeightedGrid <- function(data, variable=NULL, kernel_lat_sd=1/69, tile_latitude_size=0.005,
                               max_influence_sd = 2){
  data = data[(!is.na(data$longitude)) & !is.na(data$latitude),]
  if (!is.null(variable)){
    data = data[!is.na(data[,variable]),]
  }
  longitude = data$longitude
  latitude = data$latitude
  range_lon = range(longitude, na.rm=TRUE)
  range_lat = range(latitude, na.rm=TRUE)
  med_lon = mean(range_lon)
  med_lat = mean(range_lat)
  lon_scale_factor = cos(med_lon * pi/180)
  tile_longitude_size = tile_latitude_size * lon_scale_factor
  
  lat_tile_values = seq(range_lat[1], range_lat[2], tile_latitude_size)
  lon_tile_values = seq(range_lon[1], range_lon[2], tile_longitude_size)
  
  n_lat_tiles = length(lat_tile_values)
  n_lon_tiles = length(lon_tile_values)
  
  weight_matrix = matrix(0, ncol=n_lon_tiles, nrow=n_lat_tiles)
  
  lat_to_index <- function(latitude){
    round(n_lat_tiles * (latitude - range_lat[1])/diff(range_lat)) + 1
  }
  
  lon_to_index <- function(longitude){
    round(n_lon_tiles * (longitude - range_lon[1])/diff(range_lon))
  }
  
  lat_index_radius = round(n_lat_tiles * max_influence_sd * tile_latitude_size/diff(range_lat))
  lon_index_radius = round(n_lon_tiles * max_influence_sd * tile_longitude_size/diff(range_lon))
  
  calculate_latitude_bounds <- function(latitude){
    lat_index = lat_to_index(latitude)
    bounds = c(pmax(lat_index - lat_index_radius, 1), 
               pmin(lat_index+lat_index_radius, n_lat_tiles))
    return(bounds[1]:bounds[2])
  }
  
  calculate_longitude_bounds <- function(longitude){
    lon_index = lon_to_index(longitude)
    bounds = c(pmax(lon_index - lon_index_radius, 1), 
               pmin(lon_index+lon_index_radius, n_lon_tiles))
    return(bounds[1]:bounds[2])
  }
  
  if (!is.null(variable)){
    x = data[,variable]
  }
  else{
    x = rep(1, nrow(data))
  }
  
  calculate_distance <- function(latitude, longitude, grid){
    sqrt(((grid$longitude-longitude)*lon_scale_factor)^2 + 
           (grid$latitude-latitude)^2)
  }
  for (i in 1:nrow(data)){
    x_latitude = latitude[i]
    x_longitude = longitude[i]
    valid_lat_indexes = calculate_latitude_bounds(x_latitude)
    valid_lon_indexes = calculate_longitude_bounds(x_longitude)
    latlon_grid = expand.grid(lat_index = valid_lat_indexes, 
                              lon_index = valid_lon_indexes)
    latlon_grid$latitude = lat_tile_values[latlon_grid$lat_index]
    latlon_grid$longitude = lon_tile_values[latlon_grid$lon_index]
    
    new_weights = 1/(pi*kernel_lat_sd) * exp(-((calculate_distance(x_latitude, 
                                                     x_longitude,
                                                     latlon_grid)))^2/
                                               (2*kernel_lat_sd))
    weight_matrix[latlon_grid$lat_index + 
                    n_lat_tiles * (latlon_grid$lon_index-1)] = 
      weight_matrix[latlon_grid$lat_index + 
                      n_lat_tiles * (latlon_grid$lon_index-1)] + 
      new_weights * x[i]
  }
  final_grid = expand.grid(lat_index = 1:n_lat_tiles, lon_index=1:n_lon_tiles)
  final_grid$longitude = lon_tile_values[final_grid$lon_index]
  final_grid$latitude = lat_tile_values[final_grid$lat_index]
  final_grid$weight = getMatrixElementsByLinearIndex(weight_matrix,
                                                     final_grid$lat_index,
                                                     final_grid$lon_index)
  return(final_grid)
}
