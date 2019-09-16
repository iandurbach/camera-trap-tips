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


# ## get study area shape file
# study_area <- st_read("data/Camera Trap Activities/CameraTrapAllOrganization_2.shp")
# study_area %>% st_crs()
# study_area <- study_area %>% st_set_crs("+proj=utm +zone=46 +datum=WGS84 +units=m +no_defs")
# study_area <- study_area %>% st_transform(crs = "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")
# study_area_planned <- study_area %>% filter(Progress == "Planned") %>% 
#   filter(CameraTrap == "Munkhtogtokh")
# saveRDS(study_area_planned, file = "mesh-size/data/area_munkhtogtokh.Rds")

# load the micro design Rds 
cams <- readRDS("mesh-size/data/traps_pbar_munkhtogtokh.Rds")
cams <- st_as_sf(cams, coords = c("x", "y"), crs = "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs")

# add buffer around each camera location (how far it sees)
buff_cams <- st_buffer(cams, dist = 1000)
# add variable indicating whether camera is within the estimation area
buff_cams <- buff_cams %>% mutate(in_est_area = apply(st_intersects(cams, estimation_area, sparse = FALSE), 1, any))

# load the mesh Rds 
mesh <- readRDS("mesh-size/data/mesh_munkhtogtokh.Rds")

# load the study area
mong_sf <- readRDS("mesh-size/data/mong_sf.Rds")

# load the study area
sap_1 <- readRDS("mesh-size/data/area_munkhtogtokh.Rds")

# add buffer for an area bigger than the biggest possible mesh you might want
sap_1_mega_buffered <- st_buffer(sap_1, dist = 50000)

# place a grid over the area, with cells 3km x 3km (adjust as needed)
my_mega_grid <- st_make_grid(sap_1_mega_buffered, cellsize = c(3000, 3000))
# add smoothed occupancies in each grid cell
my_mega_grid <- st_join(st_sf(my_mega_grid), mong_sf %>% select(smoothOcc_1200), largest = TRUE)
# replace any NA occupancies with zeros
my_mega_grid <- my_mega_grid %>% mutate(smoothOcc_1200 = replace_na(smoothOcc_1200, 0))

# # remove any cells outside survey area
# my_grid <- st_intersection(my_grid, sap_1) 
# # add smoothed occupancies in each grid cell
# my_grid <- st_join(st_sf(my_grid), mong_sf %>% select(smoothOcc_1200), largest = TRUE)

# plot the cells and smoothed occupancy probabilities
pal <- colorNumeric("YlOrRd", domain = my_mega_grid$smoothOcc_1200)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = st_transform(sap_1, crs = 4326), fill = FALSE, weight = 3) %>%
  addPolygons(data = st_transform(sap_1_mega_buffered, crs = 4326), color = "blue", fill = FALSE, weight = 3) %>%
  addPolygons(data = st_transform(my_mega_grid, crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal(smoothOcc_1200)) %>%
  addPolygons(data = st_transform(buff_cams, crs = 4326), color = "black", weight = 3) 

# use the occupancy probabilities to generate a number of animals per cell

# need total number of animals on the mask
expD <- 0.5 # expected density per 100km2
expD_per_cell <- expD * (3 * 3) / 100
N <- expD_per_cell * nrow(my_mega_grid)

# cells on density surface have intensity lambda = -ln(1 - smoothOcc_1200) [from P(X=0) for Poisson]
# scaled so that sum of lambdas over all cells = N
my_mega_grid <- my_mega_grid %>% mutate(lambda = -log(1 - smoothOcc_1200)) %>%
  mutate(lambda = lambda / sum(lambda) * N)

# now generate animals on the mask
my_mega_grid <- my_mega_grid %>% mutate(n_true = rpois(nrow(my_mega_grid), lambda))
sum(my_mega_grid$n_true)

# plot the cells and actual animal locs 
pal <- colorFactor("YlOrRd", domain = my_mega_grid$n_true)
leaflet() %>%
  addTiles() %>%
  addPolygons(data = st_transform(my_mega_grid, crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal(n_true)) %>%
  addPolygons(data = st_transform(buff_cams, crs = 4326), color = "black", weight = 3) 

my_mega_grid_centroids <- st_centroid(my_mega_grid)

g0 <- 1
sigma <- 8000
dists <- st_distance(my_mega_grid_centroids, cams) 
dists <- g0 * exp(-dists^2 / (2 * sigma^2))
my_mega_grid <- my_mega_grid %>% 
  mutate(det_prob = 1 - apply(1-dists, 1, prod),
         detected = runif(nrow(my_mega_grid)) < det_prob)

my_mega_grid_centroids <- st_centroid(my_mega_grid)
all_animals <- my_mega_grid_centroids %>% filter(n_true > 0) %>% st_buffer(dist = 1000)

# plot the cells and prob of detection
pal <- colorNumeric("YlOrRd", domain = my_mega_grid$det_prob)
pal2 <- colorFactor("Set1", domain = all_animals$detected)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(my_mega_grid, crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal(det_prob)) %>%
  addPolygons(data = st_transform(buff_cams, crs = 4326),  fillOpacity = 1, color = "black", weight = 1) %>%
  addPolygons(data = st_transform(all_animals, crs = 4326), fillOpacity = 1, color = ~pal2(detected), weight = 1) 

# add a mesh with some user-specified buffer 
sap_1_mesh <- st_buffer(sap_1, dist = 5000)

# plot the cells and prob of detection
pal <- colorNumeric("YlOrRd", domain = my_mega_grid$det_prob)
pal2 <- colorFactor("Set1", domain = all_animals$detected)
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(sap_1_mesh, crs = 4326), fill = FALSE, weight = 3) %>%
  addPolygons(data = st_transform(my_mega_grid, crs = 4326), fillOpacity = 0.4, weight = 1, color = ~pal(det_prob)) %>%
  addPolygons(data = st_transform(buff_cams, crs = 4326),  fillOpacity = 1, color = "black", weight = 1) %>%
  addPolygons(data = st_transform(all_animals, crs = 4326), fillOpacity = 1, color = ~pal2(detected), weight = 1) 

# remove any cells outside mesh
my_mega_grid_on_mesh <- st_intersection(my_mega_grid, sap_1_mesh) 
# estimate effective survey area is sum of detection probabilities over all cells on the mesh
esa <- sum(my_mega_grid_on_mesh$det_prob) 
# true effective survey area is sum of detection probabilities over all cells everywhere
true_esa <- sum(my_mega_grid$det_prob) 
# bias?
(esa - true_esa) / true_esa

# number of animals seen
Nhat <- all_animals %>% filter(detected) %>% nrow()
# estimate of density per cell
Nhat / esa
# estimate of density per 100km2
Dhat <- (Nhat / esa) * 100 / (3 * 3)
# % error
(Dhat - expD) / expD 

# true number of animals on mesh
N_on_mesh <- st_intersection(all_animals, sap_1_mesh)  %>% filter(detected) %>% nrow()
# estimate of density per 100km2
Dhat2 <- (N_on_mesh / esa) * 100 / (3 * 3)
# % error
(Dhat2 - expD) / expD 











## get occupancy results
occupancies <- st_read("mesh-size/data/Mongolia_SL_Occ.shp")
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

# use the occupancy probabilities to generate a number of animals per cell

# need total number of animals on the mask
expD <- 0.5 # expected density per 100km2
expD_per_cell <- expD / (15 * 15) * 100
N <- expD_per_cell * nrow(mong_sf)

# cells on density surface have intensity lambda = -ln(1 - smoothOcc_1200) [from P(X=0) for Poisson]
# scaled so that sum of lambdas over all cells = N
mong_sf <- mong_sf %>% mutate(lambda = -log(1 - smoothOcc_1200)) %>%
  mutate(lambda = lambda / sum(lambda) * N)

# now generate animals on the mask
mong_sf <- mong_sf %>% mutate(n_true = rpois(nrow(mong_sf), lambda))
sum(mong_sf$n_true)

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

###### bias demo starts here

# select whether to use cams input from file (buff_cams) or halton cams (buff_pts)
# if halton cams indicate how many as well e.g. buff_pts[1:10, ]
npts <- 10
seed <- round(runif(1, min = 1, max = nrow(buff_pts) - npts))
cams <- buff_pts[seed:(seed + npts - 1), ]

cams <- buff_cams

## going to estimate potential biases by:
## 1) defining an estimation area using some cut-off on occupancy probability
## 2) "estimate" D by counting up how many number of animals seen / number of cells in estimation area
## 3) define an area to extrapolate the results to, again using a cut-off on occupancy probability
## 4) work out predicted Nhat and, by comparing to the true N, the % bias

## 1 & 3) defining study and extrapolation areas using some cut-offs on occupancy probability
study_area_cutoff <- 0 # 0 to include all cameras
pred_area_cutoff <- 0 # 0 to extrapolate to entire area

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
                                              pred_error = lambda - pred_n) 

prediction_area %>% data.frame() %>% summarize(N = sum(n_true),
                                               Nhat = sum(pred_n),
                                               perc_bias = (Nhat - N) / N)

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
pal <- colorNumeric("viridis", domain = c(-0.35,0.15))
leaflet() %>%
  addTiles() %>%
  addProviderTiles("Esri.WorldStreetMap") %>%
  addPolygons(data = st_transform(prediction_area, crs = 4326), fillOpacity = 0.4, weight = 1, opacity = 0.4, color = ~pal(pred_error)) %>%
  addPolygons(data = st_transform(cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
  addLegend(pal = pal, values = prediction_area$pred_error, opacity = 0.7, position = "topright", title = "Error") 
