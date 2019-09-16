## Macro-design app -- generates Halton points over a defined study area, with
## stratification

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(mapview)
library(leaflet)
library(BalancedSampling)
library(readr)

library(mgcv)

source("../02_macro-design/generate_Halton_pts.R")

##### User-specified parameters

nSim <- 1 # number of simulated animal populations to use
use_halton <- TRUE # use Halton cameras or cameras read in from file (non-random)
npts <- 12 # number of Halton points to use
expD <- 0.5 # expected density per 100km2
study_area_cutoff <- 0 # estimation area includes all sites with Pr(occ) > this cutoff (0 to use all cameras)
pred_area_cutoff <- 0 # prediction area includes all sites with Pr(occ) > this cutoff (0 to extrapolate to entire area)

##### end parameters

## get occupancy results
occupancies <- st_read("data/Mongolia_SL_Occ.shp")
occupancies %>% st_crs()
occupancies <- occupancies %>% st_set_crs("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")

# merge study area covariates into occupancy df
mong_sf <- occupancies

## going to simulate animal counts in each cell from a "true" density surface, which is calculated from
## a smoothed version of the estimated occupancy surface provided in the 'occupancies' df.

# smooth the occupancies with a gam
mong_pts <- st_centroid(mong_sf) 
mong_pts$x <- st_coordinates(mong_pts)[,1]
mong_pts$y <- st_coordinates(mong_pts)[,2]

k <- 5 # amount of smoothing
gam0 <- gam(Occ_1200 ~ te(x, y, k = k, bs=c("cs","cs")), data = mong_pts, family = "binomial")
mong_pts$smoothOcc_1200 <- predict(gam0, type = "response")

# merge the predictions back into the original dataset (NB! ST_JOIN DOESN'T WORK HERE NOT SURE WHY)
mong_sf$smoothOcc_1200 <- mong_pts$smoothOcc_1200

###### generate camera trap locations from a lat-long file 

## get camera locations
cams <- read_csv("data/Mongolia PAWS Cams2019.csv") %>% filter(!is.na(Longitude)) %>% filter(!is.na(Latitude))
cams <- st_as_sf(cams, coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")

# add buffer around each camera location (how far it sees)
buff_cams <- st_buffer(cams, dist = 12000)

###### generate camera trap locations from halton points

# now same with halton points

draw <- 10000
# set.seed(initial_seed)
my_seeds <- round(c(runif(1, 1, 1e6), runif(1, 1, 1e6)))
pts <- generate_Halton_pts(n = draw, seeds = my_seeds, bases = c(2,3))

# scale and shift Halton points to fit the bounding box of the study region
bb <- st_bbox(mong_sf)
scale <- c(bb$xmax, bb$ymax) - c(bb$xmin, bb$ymin)
shift <- c(bb$xmin, bb$ymin)

pts <- t(t(pts) * scale + shift)
pts <- data.frame(pts)
colnames(pts) <- c("X", "Y")

# make the data frame an sf object
pts_sf <- st_as_sf(pts, coords = c("X", "Y"), crs = "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")

# remove any points that fall outside the study area
pts_sf <- pts_sf %>% filter(apply(st_intersects(pts_sf, mong_sf, sparse = FALSE), 1, any)) 

# add variable indicating whether camera is within the estimation area
pts_sf <- pts_sf %>% mutate(in_est_area = apply(st_intersects(pts_sf, mong_sf, sparse = FALSE), 1, any)) 

# add buffer around each camera location (how far it sees)
buff_pts <- st_buffer(pts_sf, dist = 40000)

# use the occupancy probabilities to generate a number of animals per cell

# need total number of animals on the mask
expD_per_cell <- expD * (15 * 15) / 100
N <- expD_per_cell * nrow(mong_sf)

# cells on density surface have intensity lambda = -ln(1 - smoothOcc_1200) [from P(X=0) for Poisson]
# scaled so that sum of lambdas over all cells = N
mong_sf <- mong_sf %>% mutate(lambda = -log(1 - smoothOcc_1200)) %>%
  mutate(lambda = lambda / sum(lambda) * N)

bias_table <- data.frame(`True N` = as.integer(),
                         `Estimated N` = as.integer(),
                         `% bias` = as.numeric())
survey_table <- data.frame(n = as.integer(),
                             n_seen = as.integer(),
                             n_per_cell = as.numeric())
for(i in 1:nSim){
  
  # now generate animals on the mask
  mong_sf <- mong_sf %>% mutate(n_true = rpois(nrow(mong_sf), lambda))
  sum(mong_sf$n_true)
  
  # select whether to use cams input from file (buff_cams) or halton cams (buff_pts)
  # if halton cams indicate how many as well e.g. buff_pts[1:npts, ]
  seed <- round(runif(1, min = 1, max = nrow(buff_pts) - npts))
  halton_cams <- buff_pts[seed:(seed + npts - 1), ]
  
  if(use_halton){ cams <- halton_cams } else { cams <- buff_cams }
  
  ## going to estimate potential biases by:
  ## 1) defining an estimation area using some cut-off on occupancy probability
  ## 2) "estimate" D by counting up how many number of animals seen / number of cells in estimation area
  ## 3) define an area to extrapolate the results to, again using a cut-off on occupancy probability
  ## 4) work out predicted Nhat and, by comparing to the true N, the % bias
  
  ## 1 & 3) defining study and extrapolation areas using some cut-offs on occupancy probability
  estimation_area <- mong_sf %>% filter(smoothOcc_1200 > study_area_cutoff) 
  prediction_area <- mong_sf %>% filter(smoothOcc_1200 > pred_area_cutoff) 
  
  # add variable indicating whether camera is within the estimation area
  cams <- cams %>% mutate(in_est_area = apply(st_intersects(cams, estimation_area, sparse = FALSE), 1, any)) 
  
  # add variable indicating whether each cell intersections with an cameras field of view 
  # ... that means animals in that cell are detected
  estimation_area <- estimation_area %>% 
    mutate(detected = apply(st_intersects(estimation_area, cams, sparse = FALSE), 1, any))
  
  ## 2) "estimate" D by counting up how many number of animals seen / number of cells in study area
  
  survey_results <- estimation_area %>% data.frame() %>%
    filter(detected == TRUE) %>% summarize(n = sum(n_true),
                                           n_seen_cells = n(),
                                           n_per_cell = n / n_seen_cells)
  
  ## 4) work out predicted Nhat and, by comparing to the true N, the % bias
  prediction_area <- prediction_area %>% mutate(pred_n = survey_results$n_per_cell,
                                                pred_error = pred_n - lambda) 
  
  bias_table <- rbind(bias_table, 
                      prediction_area %>% data.frame() %>% 
                        summarize(`True N` = sum(n_true),
                                  `Estimated N` = sum(pred_n),
                                  `% bias` = 100 * (sum(pred_n) - sum(n_true)) / sum(n_true)))
  
  survey_table <- rbind(survey_table,
                        estimation_area %>% data.frame() %>%
                          filter(detected == TRUE) %>% summarize(n = sum(n_true),
                                                                 n_seen_cells = n(),
                                                                 n_per_cell = n / n_seen_cells))
  
}

bias_table_means <- bias_table %>% summarize_all(mean)
survey_table_means <- survey_table %>% summarize_all(mean)
bias_table_means
survey_table_means

# plot expected densities over the estimation area
pal <- colorNumeric("YlOrRd", domain = mong_sf$lambda)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(estimation_area, crs = 4326), fillOpacity = 0.7, weight = 1, color = ~pal(lambda)) %>%
  addPolygons(data = st_transform(cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
  addLegend(pal = pal, values = mong_sf$lambda, opacity = 0.7, position = "topright", title = "Expected N") 

# plot true counts over the estimation area
pal <- colorFactor("YlOrRd", domain = mong_sf$n_true)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(estimation_area, crs = 4326), fillOpacity = 0.7, weight = 1, color = ~pal(n_true)) %>%
  addPolygons(data = st_transform(cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
  addLegend(pal = pal, values = mong_sf$n_true, opacity = 0.7, position = "topright", title = "Realised N") 

# plot errors in prediction area
pal <- colorNumeric("viridis", domain = c(-1.05,1.5))
#pal <- colorNumeric("viridis", domain = prediction_area$pred_error)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(prediction_area, crs = 4326), fillOpacity = 0.4, weight = 1, opacity = 0.4, color = ~pal(pred_error)) %>%
  addPolygons(data = st_transform(cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
  addLegend(pal = pal, values = c(-1.05,1.5), opacity = 0.7, position = "topright", title = "Error") 
 
# library(webshot)
# mapshot(m, file = "biased_error.png")
