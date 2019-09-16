library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(mapview)
library(leaflet)
library(BalancedSampling)
library(readr)
library(gridExtra)
library(mgcv)

source("generate_Halton_pts.R")

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

ui <- fluidPage(
  
  titlePanel("Sampling favourable habitats biases abundance estimates"),
  # titlePanel("No extrapolation from biased sampling"), 
  
  column(8,
         column(12, leafletOutput("map1"), 
                absolutePanel(bottom = 20, left = 20, 
                              fluidRow(
                                column(6, actionButton(inputId = "calc", label = "Calculate bias"))),
                              column(6, checkboxInput(inputId = "halton", label = "Random sampling")))),
         column(12, leafletOutput("map2"),
                absolutePanel(top = 10, right = 20, 
                              fluidRow(
                                column(6, tableOutput("bias")))))),
  
  column(4,
         p("This app shows that abundance estimates are biased if sampling is done in mostly 
           favourable habitats"),
         p("Click the `Calculate bias` button to run a survey with cameras mostly in high density areas. Cameras
         are indicated in black. Animal locations are randomly generated each run (set the plot type to 
         `Realised N` to see where they are)."),
         p("To see what would happen with a spatially balanced sample, select the `Random sampling` box and 
         `Calculate bias` again. The bias in the balanced sample will usually be 
           substantially lower! Increase the number of simulations to get an idea of longer run patterns -- the 
           numbers in the table are means."),
         p("Use the sliders to limit the extent of the survey area (for example, sampling very favourable areas)
           and/or only predict abundance in a subregion of the study area."), 
         p(),
         sliderInput("cutoff_exparea",
                     "Survey area cut-off",
                     min = 0,
                     max = floor(10 * max(mong_pts$smoothOcc_1200)) / 10,
                     step = 0.1,
                     value = 0),
         sliderInput(inputId = "cutoff_predarea",
                     label = "Prediction area cut-off",
                     min = 0,
                     max = floor(10 * max(mong_pts$smoothOcc_1200)) / 10,
                     step = 0.1,
                     value = 0),
         numericInput(inputId = "D",
                      label = "Expected D/100km2",
                      min = 0.01,
                      step = 0.05,
                      value = c(0.7)),
         selectInput(inputId = "plot_type",
                     label = "Data to plot",
                     choices = c("Expected N" = "exp_n",
                                 "Realised N" = "real_n")),
         numericInput(inputId = "Nsim",
                      label = "Number of simulations",
                      min = 1,
                      step = 1,
                      value = c(1)))
  
)

server <- function(input, output) {
  
  bias_table <- data.frame(`True N` = as.integer(),
                           `Estimated N` = as.integer(),
                           `% bias` = as.numeric())
  
  bias_calcs <- eventReactive(input$calc, {
    
    ## going to estimate potential biases by:
    ## 1) defining an estimation area using some cut-off on occupancy probability
    ## 2) "estimate" D by counting up how many number of animals seen / number of cells in estimation area
    ## 3) define an area to extrapolate the results to, again using a cut-off on occupancy probability
    ## 4) work out predicted Nhat and, by comparing to the true N, the % bias
    
    # use the occupancy probabilities to generate a number of animals per cell
    
    # need total number of animals on the mask
    expD <- input$D # expected density per 100km2
    expD_per_cell <- expD / (15 * 15) * 100
    N <- expD_per_cell * nrow(mong_sf)
    
    # cells on density surface have intensity lambda = -ln(1 - smoothOcc_1200) [from P(X=0) for Poisson]
    # scaled so that sum of lambdas over all cells = N
    mong_sf <- mong_sf %>% mutate(lambda = -log(1 - smoothOcc_1200)) %>%
      mutate(lambda = lambda / sum(lambda) * N)
    
    for(i in 1:input$Nsim){
      
      # now generate animals on the mask
      mong_sf <- mong_sf %>% mutate(n_true = rpois(nrow(mong_sf), lambda))
      sum(mong_sf$n_true)
      
      # select whether to use cams input from file (buff_cams) or halton cams (buff_pts)
      # if halton cams indicate how many as well e.g. buff_pts[1:10, ]
      npts <- 10
      seed <- round(runif(1, min = 1, max = nrow(buff_pts) - npts))
      
      if(input$halton){
        cams <- buff_pts[seed:(seed + npts - 1), ]
      } else { cams <- buff_cams }
      
      ## 1 & 3) defining study and extrapolation areas using some cut-offs on occupancy probability
      study_area_cutoff <- input$cutoff_exparea # 0 to include all cameras
      pred_area_cutoff <- input$cutoff_predarea # 0 to extrapolate to entire area
      
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
      
      bias_table <- rbind(bias_table, prediction_area %>% data.frame() %>% summarize(`True N` = sum(n_true),
                                                                                     `Estimated N` = sum(pred_n),
                                                                                     `% bias` = 100 * (sum(pred_n) - sum(n_true)) / sum(n_true)))
      if(i == 1){
        prediction_area$cellID <- 1:nrow(prediction_area)
        prediction_area_all <- prediction_area
      } else {
        prediction_area$cellID <- 1:nrow(prediction_area)
        prediction_area_all <- rbind(prediction_area_all, prediction_area)
      }
      
    }
    
    bias_table_means <- bias_table %>% summarize_all(mean)
    prediction_area <- prediction_area_all %>% group_by(cellID) %>% summarize(pred_n = mean(pred_n),
                                                                              pred_error = mean(pred_error))
    
    return(list(estimation_area = estimation_area, prediction_area = prediction_area, 
                cams = cams, bias_table = bias_table_means, all_area = mong_sf))
    
  })
  
  # only include aspects of the map that won't need to change dynamically
  output$map1 <- renderLeaflet({
    p1 <- leaflet() %>% addTiles() %>%
      addProviderTiles("Esri.WorldStreetMap") %>%
      flyToBounds(95, 42, 101, 52)
    
  })
  
  output$map2 <- renderLeaflet({
    leaflet() %>% addTiles() %>%
      addProviderTiles("Esri.WorldStreetMap") %>%
      flyToBounds(95, 42, 101, 52)
  })
  
  # plot expected densities over the estimation area
  observe({
    if(input$plot_type == "exp_n"){
      pal <- colorNumeric("YlOrRd", domain = bias_calcs()$all_area$lambda)
      leafletProxy("map1") %>% clearControls() %>% 
        clearGroup("lambdas") %>% clearGroup("cams") %>% clearGroup("ns") %>% clearGroup("haltons") %>%
        addPolygons(group = "lambdas", data = st_transform(bias_calcs()$estimation_area, crs = 4326), fillOpacity = 0.7, weight = 1, color = ~pal(lambda)) %>%
        addPolygons(group = "cams", data = st_transform(bias_calcs()$cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
        addLegend(pal = pal, values = bias_calcs()$all_area$lambda, opacity = 0.7, position = "bottomright", title = "Expected N")
    } else {
      pal <- colorFactor("YlOrRd", domain = bias_calcs()$all_area$n_true)
      leafletProxy("map1") %>% clearControls() %>% 
        clearGroup("lambdas") %>% clearGroup("cams") %>% clearGroup("ns") %>% clearGroup("haltons") %>%
        addPolygons(group = "ns", data = st_transform(bias_calcs()$estimation_area, crs = 4326), fillOpacity = 0.7, weight = 1, color = ~pal(n_true)) %>%
        addPolygons(group = "haltons", data = st_transform(bias_calcs()$cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
        addLegend(pal = pal, values = bias_calcs()$all_area$n_true, opacity = 0.7, position = "bottomright", title = "Realised N") 
    }
  })
  
  
  # plot errors in prediction area
  pal <- colorNumeric("viridis", domain = c(0.5,-0.3))
  observe({
    leafletProxy("map2") %>% clearControls() %>% 
      clearGroup("errors") %>% clearGroup("haltons") %>%
      addPolygons(group = "errors", data = st_transform(bias_calcs()$prediction_area, crs = 4326), fillOpacity = 0.4, weight = 1, opacity = 0.4, color = ~pal(pred_error)) %>%
      addPolygons(group = "haltons", data = st_transform(bias_calcs()$cams %>% filter(in_est_area), crs = 4326), color = "black", opacity = 0.5, weight = 1) %>%
      addLegend(pal = pal, values = bias_calcs()$prediction_area$pred_error, opacity = 0.7, position = "bottomright", title = "Error") 
  })
  
  output$bias <- renderTable({
    bias_calcs()$bias_table
  }, digits = 0, width = "300px", align = "c")
  
  
}

shinyApp(ui = ui, server = server)
