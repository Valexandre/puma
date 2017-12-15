library(ggplot2)
library(ggmap)
library(cetcolor)

source('weighted_map_raster.r')

party_ready = read.csv('data/party_ready.csv')
party_sideways = read.csv('data/party_sideways.csv')
party_ready$duration = as.numeric(as.POSIXct(party_ready$closed_date) - 
                                    as.POSIXct(party_ready$created_date))

regular_weights = createWeightedGrid(party_ready)
duration_weights = createWeightedGrid(party_ready, variable='duration')
duration_weights$normalized = duration_weights$weight/(1e-10 + regular_weights$weight)

nyc <- get_map('new york city, new york', maptype='toner', source='stamen')

ggmap(nyc) + geom_tile(data=duration_weights, aes(x=longitude, y=latitude, 
                                               fill=normalized,
                                               alpha=normalized)) + 
  scale_alpha_continuous() + 
  scale_color_gradientn(colors=cet_pal(5))

antilog_formatter = function(x){
  round(exp(x)-0.1)
}

ggmap(nyc) + coord_cartesian() + 
  geom_raster(data=duration_weights[duration_weights$normalized > exp(0.101),], aes(x=longitude, y=latitude, 
                                                  fill=log(0.1 + normalized),
                                                  alpha=log(0.1 + normalized))) + 
  scale_alpha_continuous() + 
  scale_fill_gradientn(colors=cet_pal(5), label=antilog_formatter) 
