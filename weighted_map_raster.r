#there is a bug in the ggtern code, so I will fix it
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  #CHANGE 1: make n longer if n is length 1
  if (length(n)==1)
    n=rep(n, 2)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n[1]) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n[2]) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  #CHANGE 2: select the appropriate element of `n` for each matrix`
  z <- (matrix(rep(w,n[1]), nrow=n[1], ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n[1], nx)) %*% t(matrix(dnorm(ay), n[2], nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

#still need to debug kde2d.weighted()
createWeightedGrid.kde <- function(data, variable=NULL, kernel_lat_sd=1./69, tile_latitude_size=0.005, padding_lat=0.005){
  data = data[(!is.na(data$longitude)) & !is.na(data$latitude),]
  if (!is.null(variable)){
    data = data[!is.na(data[,variable]),]
  }
  longitude = data$longitude
  latitude = data$latitude
  med_lon = mean(range(longitude))
  med_lat = mean(range(latitude))
  lon_scale_factor = cos(med_lon * pi/180)
  padding_lon = lon_scale_factor * padding_lat
  range_lon = range(longitude, na.rm=TRUE) + c(-1, 1) * padding_lon
  range_lat = range(latitude, na.rm=TRUE) + c(-1, 1) * padding_lat
  
  
  tile_longitude_size = tile_latitude_size * lon_scale_factor
  padding_lon = lon_scale_factor * padding_lat
  
  lat_tile_values = seq(range_lat[1], range_lat[2], tile_latitude_size)
  lon_tile_values = seq(range_lon[1], range_lon[2], tile_longitude_size)
  
  n_lat_tiles = length(lat_tile_values)
  n_lon_tiles = length(lon_tile_values)
  if (!is.null(variable)){
    w = data[,variable]
  }
  else{
    w = rep(1, nrow(data))
  }
  
  kde_list = kde2d.weighted(longitude, latitude, 
                        h=kernel_lat_sd * c(1/lon_scale_factor, 1),
                        n=100,#,n=as.integer(c(n_lon_tiles, n_lat_tiles)), 
                        lims = c(range_lon, range_lat), 
                        w=w )
  egrid = expand.grid(longitude=kde_list$x, latitude=kde_list$y)
  egrid$weight = as.numeric(kde_list$z)
  return(egrid)
}


getMatrixElementsByLinearIndex <- function(mat, row_idx, col_idx){
  nr = nrow(mat)
  return(mat[row_idx + nr*(col_idx-1)])
}


#about 2.5 times slower than kde2d from MASS, but allows usage of a variable to produce a weighted density

## Function creates a kernel density grid, either based on density of points (default) 
##          or a weighted density
# variable - character of a variable that is used to weight the density
# kernel_lat_sd - the standard deviation/bandwidth of the kernel in latitude units; default 
#                 is approximately 1 mile
# tile_latitude_size - the length of each tile's latitude dimension; longitude will scale 
#                      with this to be roughly square
# max_influence_sd - the maximum number of standard deviations in either axis that a point
#                    can create influence
# padding_lat - the amount of latitude added in each direction where the density will still
#               be calculated; padding for longitude is calculated from this
# normalize_weights - if TRUE, then the kernel density will be normalized over each affected 
#                     region for each data point
createWeightedGrid <- function(data, variable=NULL, kernel_lat_sd=1/69, tile_latitude_size=0.005,
                               max_influence_sd = 2, padding_lat=0.005, normalize_weights=TRUE){
  data = data[(!is.na(data$longitude)) & !is.na(data$latitude),]
  if (!is.null(variable)){
    data = data[!is.na(data[,variable]),]
  }
  longitude = data$longitude
  latitude = data$latitude
  med_lon = mean(range(longitude))
  med_lat = mean(range(latitude))
  lon_scale_factor = cos(med_lon * pi/180)
  padding_lon = lon_scale_factor * padding_lat
  range_lon = range(longitude, na.rm=TRUE) + c(-1, 1) * padding_lon
  range_lat = range(latitude, na.rm=TRUE) + c(-1, 1) * padding_lat
 
  
  tile_longitude_size = tile_latitude_size * lon_scale_factor
  padding_lon = lon_scale_factor * padding_lat
  
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
    if (normalize_weights){
      new_weights = new_weights/(sum(new_weights) + 1e-8)
    }
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




#similar to above, but with different method of calculation
#slower, needs correction
createWeightedGrid2 <- function(data, variable=NULL, kernel_lat_sd=1/69, tile_latitude_size=0.005,
                               max_influence_sd = 2, padding_lat=0.005){
  data = data[(!is.na(data$longitude)) & !is.na(data$latitude),]
  if (!is.null(variable)){
    data = data[!is.na(data[,variable]),]
  }
  longitude = data$longitude
  latitude = data$latitude
  simple_data = data.frame(latitude=data$latitude, longitude=data$longitude)
  med_lon = mean(range(longitude))
  med_lat = mean(range(latitude))
  lon_scale_factor = cos(med_lon * pi/180)
  padding_lon = lon_scale_factor * padding_lat
  range_lon = range(longitude, na.rm=TRUE) + c(-1, 1) * padding_lon
  range_lat = range(latitude, na.rm=TRUE) + c(-1, 1) * padding_lat
  
  
  tile_longitude_size = tile_latitude_size * lon_scale_factor
  padding_lon = lon_scale_factor * padding_lat
  
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
  for (i in 1:n_lat_tiles){
    for (j in 1:n_lon_tiles){
      weight_matrix[i, j] = sum(x * 1/(pi*kernel_lat_sd) * exp(-((calculate_distance(
        lat_tile_values[i],
        lon_tile_values[j],
        simple_data
        
      )))^2/(2*kernel_lat_sd)))
    }
  }
  final_grid = expand.grid(lat_index = 1:n_lat_tiles, lon_index=1:n_lon_tiles)
  final_grid$longitude = lon_tile_values[final_grid$lon_index]
  final_grid$latitude = lat_tile_values[final_grid$lat_index]
  final_grid$weight = getMatrixElementsByLinearIndex(weight_matrix,
                                                     final_grid$lat_index,
                                                     final_grid$lon_index)
  return(final_grid)
}

