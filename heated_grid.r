library(ggplot2)
library(ggmap)
library(cetcolor)
library(microbenchmark)
library(MASS)

source('weighted_map_raster.r')

party_ready = read.csv('data/party_ready.csv')
party_sideways = read.csv('data/party_sideways.csv')
party_ready$duration = as.numeric(as.POSIXct(party_ready$closed_date) - 
                                    as.POSIXct(party_ready$created_date))

party_ready_realistic = with(party_ready, party_ready[duration > 0 & duration < 3000 & !is.na(longitude),])

#let's compare createWeightedGrid performance to the one based on kde2d.weighted
microbenchmark(createWeightedGrid(party_ready), times=1)
microbenchmark(createWeightedGrid.kde(party_ready), times=1)


default_density = createWeightedGrid(party_ready_realistic)
weighted_density = createWeightedGrid(party_ready_realistic, variable='duration')
weighted_density$normalized = weighted_density$weight/(1e-10 + default_density$weight)
weighted_density$default_density = default_density$weight

nyc <- get_map('new york city, new york', maptype='toner', source='stamen', zoom=11)


antilog_formatter = function(x){
  round(exp(x))
}

r5 = sqrt(5)

#units of density is 
ggmap(nyc) + coord_cartesian() + 
  geom_raster(data=weighted_density[weighted_density$normalized > exp(0.101),], 
              aes(x=longitude, y=latitude, 
                  fill=log( default_density),
                  alpha=log(default_density))) + 
  scale_alpha_continuous() + 
  scale_fill_gradientn('Density of Complaints', 
                       colors=cet_pal(5), label=antilog_formatter, breaks = log(c(1/r5, 1,r5, 4, 4*r5, 20, 20*r5, 100, 100*r5 ))) +
  ggtitle('Density of Noise Complaints in NYC (per year)') +
  guides(alpha='none')

ggmap(nyc) + coord_cartesian() + 
  geom_raster(data=weighted_density[weighted_density$normalized > exp(0.101),], 
              aes(x=longitude, y=latitude, 
                  fill=log( normalized),
                  alpha=log(normalized))) + 
  scale_alpha_continuous() + 
  scale_fill_gradientn('Average Response Time (minutes)', 
                       colors=cet_pal(5), label=antilog_formatter) +
  ggtitle('Average Response Time to Noise Complaints in NYC') +
  guides(alpha='none')
