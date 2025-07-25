
# IUCN Biodiversity analyst practical excercise
# 2025 July 25
# Laura Berman
# To be considered, completed exercises must be returned via email to nicholas.macfarlane@iucn.org by 11:59pm Washington DC time (EDT) on Thursday July 31st, 2025.

# Clean code is above.
# Data exploration is at the end of the file.

library(readxl)
library(writexl)
library(sf)
library(mapview)
library(raster)
library(terra)
library(ggplot2)

# setwd ------------------------------------------------------------------------
setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Job Application Portfolio/IUCN biodiversity analyst/IUCN Practical")


################################################################################
# Part 1: STAR Analysis
################################################################################

# Species Threat Abatement and Restoration (STAR) metric
# estimate the potential reduction in species extinction risk
# AOH = area of habitat. % of total habitat held within a cell


# load data --------------------------------------------------------------------
RL_cat <- read_excel("Part 1 - STAR/species_RLcategory.xlsx") 
BraulioCarrillo <- read_sf("Part 1 - STAR/Braulio Carillo shp/braulio_carrillo.shp")

areas_path <- "Part 1 - STAR/Species Area of Habitat"
tif_files <- list.files(areas_path, pattern = "\\.tif$", full.names = TRUE)
sp_rasters <- lapply(tif_files, rast)


# match projections ------------------------------------------------------------
BraulioCarrillo_proj <- st_transform(BraulioCarrillo, crs(sp_rasters[[1]]))


# 1. Calculate species habitat proportion within the park ----------------------

# Convert sf to terra vect
park_vect <- vect(BraulioCarrillo_proj)

# Initialize output list
results <- lapply(seq_along(sp_rasters), function(i) {
  r <- sp_rasters[[i]]
  name <- names(r)
  
  # find the total global habitat
  global_sum <- global(r, "sum", na.rm = TRUE)[1, 1]
  
  # find total habitat within the park (the sum of intersecting pixels, weighted by how much of the pixel intersects with the park)
  park_sum <- extract(r, park_vect, weights = TRUE, fun = sum, na.rm = TRUE)[1, 2]
  
  data.frame(
    Raster = name,
    GlobalHabitat = global_sum,
    HabitatInPark = park_sum
  )
})

# Combine results into a single data frame
habitat_summary <- do.call(rbind, results)

# get proportion in park
habitat_summary$PropInPark <- habitat_summary$HabitatInPark / habitat_summary$GlobalHabitat


# 2. Calculate STAR score within the park --------------------------------------

# add weighted risk column
RL_cat$weightedRisk <- c(300, 300, 400, 200)

# merge habitat_summary and risk level
names(habitat_summary)[names(habitat_summary) == "Raster"] <- "Species"
sp_STAR_score <- merge(habitat_summary, RL_cat)

# star score of a species at a site is the % of its population at the site * extinction risk
sp_STAR_score$STARscore <- sp_STAR_score$PropInPark * sp_STAR_score$weightedRisk

# save it
write_xlsx(sp_STAR_score, "speciesSTARscore_BraulioCarrillo.xlsx")




################################################################################
# Part 2: remote sensing analysis
################################################################################

# load data --------------------------------------------------------------------
AGBD <- rast("Part 2 - Remote Sensing/above_ground_biomass_density.tif")
treeCover <- rast("Part 2 - Remote Sensing/tree_cover.tif")
field <- st_read("Part 2 - Remote Sensing/field.geojson")


# 1. calculate the average tree cover of the field -----------------------------
extract(treeCover, field, weights = TRUE, fun = mean, na.rm = TRUE) # 31.606 % tree cover

# 2. calculate the total biomass of the field ----------------------------------

# how big is one pixel?
AGBD_proj <- project(AGBD, "ESRI:54009")
pixel_area_ha <- prod(res(AGBD_proj)) / 10000 # 0.015 hectares

# sum up the pixels
sum_of_pixels <- extract(AGBD, field, weights = TRUE, fun = sum, na.rm = TRUE)[2] # 17001.06 

# divide sum by pixel size to get tons/hectare
sum_of_pixels * pixel_area_ha # 266.8 tons of biomass




################################################################################
# Data Exploration
################################################################################

# This is code I used to explore the dataset, aid in finding the best clean code, and double check the accuracy of results.

# PART 1 -----------------------------------------------------------------------

A.annae_area <- rast("Part 1 - STAR/Species Area of Habitat/Agalychnis annae.tif")
E.fimb_area <- rast("Part 1 - STAR/Species Area of Habitat/Ecnomiohyla fimbrimembra.tif")
I.hold_area <- rast("Part 1 - STAR/Species Area of Habitat/Incilius holdridgei.tif")
L.vib_area <- rast("Part 1 - STAR/Species Area of Habitat/Lithobates vibicarius.tif")

mapview(BraulioCarrillo) # hello Costa Rica
BraulioCarrillo
mapview(A.annae_area) # a local species
A.annae_area # WGS84
summary(A.annae_area)
mapview(E.fimb_area) # also local
E.fimb_area
mapview(I.hold_area)
I.hold_area
summary(I.hold_area) # only NAs. Is this species extinct?
# Incilius holdridgei is critically endangered, known from one location. https://www.iucnredlist.org/species/54664/54358615
mapview(L.vib_area) # another local species

mapview(BraulioCarrillo) + mapview(A.annae_area)

BraulioCarrillo_proj <- st_transform(BraulioCarrillo, crs(A.annae_area))
BraulioCarrillo_sp <- as(BraulioCarrillo_proj, "Spatial")
r_crop <- terra::crop(A.annae_area, BraulioCarrillo_sp)
mapview(BraulioCarrillo) + mapview(r_crop) # eyeballing it, ~11 "present" pixels fit inside BraulioCarrillo
x <- as.data.frame(r_crop) # pixel values range from 0-1. 
x <- global(A.annae_area, "sum", na.rm = TRUE)

# The sum of all A.annae_area pixels is 151, not 1
# individual pixel values range from 0-1
# so I am assuming that pixel values are the proportion of that pixel which is suitable habitat
# rather than the percent of the total global habitat held within the pixel (which is what I thought the instructions said)

x <- raster::extract(A.annae_area, BraulioCarrillo_proj, weights=TRUE, na.rm=TRUE, fun=sum) # 11.08. That is close to my eyeball estimate. 
x$`Agalychnis annae`


# PART 2 -----------------------------------------------------------------------

mapview(AGBD) # enormous
mapview(field) # in cambodia
mapview(treeCover)
mapview(field) + mapview(treeCover) # only need a fraction of the raster
AGBD
treeCover # same CRS
field

# 1. calculate the average tree cover of the field
extract(treeCover, field, weights = TRUE, fun = mean, na.rm = TRUE)

# 2. calculate the total biomass of the field
extract(AGBD, field, weights = TRUE, fun = sum, na.rm = TRUE)

# is one pixel one hectare?
AGBD_proj <- project(AGBD, "ESRI:54009")
pixel_area_ha <- prod(res(AGBD_proj)) / 10000
