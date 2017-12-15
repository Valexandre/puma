# puma
visualization based off R tutorial from Chicago RLadies meetup

See [http://maxcandocia.com/article/2017/Dec/15/overlaying-heatmaps-on-geo/](http://maxcandocia.com/article/2017/Dec/15/overlaying-heatmaps-on-geo/) for visualizations and explanations.

The most interesting thing to someone just looking for copy/paste code is probably the `createWeightedGrid()` function in "weighted_map_raster.r". It computs weighted densities across a data frame with longitude/latitude coordinates. Note that it works relatively well for smaller regions, such as cities, or possibly states, but does not generalize as well due to the non-cartesian nature of those coordinates. 
