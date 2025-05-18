# Libraries needed for spatial distribution analysis
library(ineq)    # For Gini index calculation
library(ggplot2) # For plotting
library(dplyr)   # For data manipulation
library(tidyr)   # For reshaping data
library(sf)
library(viridis)
library(maps)
library(rnaturalearth)
library(patchwork)

################################################################################
# Load world map data
world_map=ne_countries(scale = "medium", returnclass = "sf")

# Calculate Gini coefficient and D95 (simplified versions)
calcGini=function(x) {
  x=sort(x)
  n=length(x)
  return(sum((2*seq(1:n)-n-1)*x)/(n*sum(x)))}

calcD95=function(biomass) {
  sorted_biomass=sort(biomass,decreasing=TRUE)
  cumsum_biomass=cumsum(sorted_biomass)
  total_biomass=sum(biomass)
  threshold_index=min(which(cumsum_biomass>=0.95*total_biomass))
  return(threshold_index/length(biomass))}

# Set random seed for reproducibility
set.seed(123)

minCell=25

# Example 1: Even, Wide Distribution
grid=expand.grid(lon=seq(-80,0,by=1), lat=seq(-20,40,by=1))
# Increase mean and decrease SD for more even pattern
grid$biomass=abs(rnorm(nrow(grid), mean=15, sd=1.5))  

gini1=round(calcGini(grid$biomass), 2)
d95_1=round(calcD95( grid$biomass), 2)

evenPlot=ggplot() +
  # Use tiles instead of points for continuous appearance
  geom_tile(data=grid, aes(x=lon, y=lat, fill=biomass)) +
  geom_polygon(data=map_data("world"), aes(x=long, y=lat, group=group),
               fill="darkgreen", color="gray70", size=0.2) +
  # Adjust contrast with direction parameter and limit values
  scale_fill_viridis(option="inferno", 
                     limits=c(min(grid$biomass), max(grid$biomass))) +
  coord_quickmap(xlim=c(-80, 0), ylim=c(-20, 40)) +
  labs(title="Even, Wide Distribution",
       subtitle=paste("Gini =", gini1, " D95 =", d95_1),
       x="Longitude", y="Latitude", fill="Biomass") +
  theme_minimal()


# Create example distribution data for tuna (more clustered) and swordfish (more diffuse)
set.seed(123)
grid_points=expand.grid(
  lon = seq(-80, 0, by = 1),
  lat = seq(-20, 40, by = 1)
)

# Create distribution patterns
# Tuna: concentrated in hotspots (higher Gini)
tuna_centers=data.frame(
  lon = c(-70, -45, -30),
  lat = c(15, 5, 25)
)

tuna_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(tuna_centers)) {
    dist=sqrt((grid_points$lon[i] - tuna_centers$lon[j])^2 + 
                   (grid_points$lat[i] - tuna_centers$lat[j])^2)
    tuna_biomass[i]=tuna_biomass[i] + 100 * exp(-0.1 * dist)
  }
}
tuna_biomass=tuna_biomass + rnorm(length(tuna_biomass), 0, 5)
tuna_biomass[tuna_biomass <minCell]=0

# Swordfish: more widely distributed (lower Gini)
swordfish_centers=data.frame(
  lon = c(-65, -30, -50, -40),
  lat = c(10, 30, -10, 15)
)

swordfish_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(swordfish_centers)) {
    dist=sqrt((grid_points$lon[i] - swordfish_centers$lon[j])^2 + 
                   (grid_points$lat[i] - swordfish_centers$lat[j])^2)
    swordfish_biomass[i]=swordfish_biomass[i] + 80 * exp(-0.05 * dist)}}

swordfish_biomass=swordfish_biomass + rnorm(length(swordfish_biomass), 0, 5)
swordfish_biomass[swordfish_biomass <minCell]=0

tuna_gini=round(calcGini(tuna_biomass), 2)
tuna_d95=round(calcD95(tuna_biomass), 2)
swordfish_gini=round(calcGini(swordfish_biomass), 2)
swordfish_d95=round(calcD95(swordfish_biomass), 2)

# Combine data for plotting
grid_points$tuna=tuna_biomass
grid_points$swordfish=swordfish_biomass

# Plot tuna distribution
tunaPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = tuna)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = paste0("Schooling"),
    subtitle = paste0("Gini = ", tuna_gini, ", D95 = ", tuna_d95),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

# Plot swordfish distribution
swordfishPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = swordfish)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = paste0("Solitary"),
    subtitle = paste0("Gini = ", swordfish_gini, ", D95 = ", swordfish_d95),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

# Simulate range contraction for bigeye tuna
# Create historical (wide) distribution
historical_centers=data.frame(
  lon = c(-70, -55, -40, -30, -65, -45),
  lat = c(0, 10, 20, 5, -10, -5)
)

historical_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(historical_centers)) {
    dist=sqrt((grid_points$lon[i] - historical_centers$lon[j])^2 + 
                   (grid_points$lat[i] - historical_centers$lat[j])^2)
    historical_biomass[i]=historical_biomass[i] + 90 * exp(-0.08 * dist)
  }
}
historical_biomass=historical_biomass + rnorm(length(historical_biomass), 0, 5)
historical_biomass[historical_biomass <minCell]=0

# Create recent (contracted) distribution
recent_centers=data.frame(
  lon = c(-70, -55, -45),
  lat = c(0, 10, -5)
)

recent_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(recent_centers)) {
    dist=sqrt((grid_points$lon[i] - recent_centers$lon[j])^2 + 
                   (grid_points$lat[i] - recent_centers$lat[j])^2)
    recent_biomass[i]=recent_biomass[i] + 100 * exp(-0.07 * dist)
  }
}
recent_biomass=recent_biomass + rnorm(length(recent_biomass), 0, 5)
recent_biomass[recent_biomass <minCell]=0

# Calculate spatial indicators
hist_gini=round(calcGini(historical_biomass), 2)
hist_d95=round(calcD95(historical_biomass), 2)
recent_gini=round(calcGini(recent_biomass), 2)
recent_d95=round(calcD95(recent_biomass), 2)

# Combine for plotting
grid_points$historical=historical_biomass
grid_points$recent=recent_biomass

# Plot range contraction
historicalPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = historical)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = "Historical Distribution (1970s)",
    subtitle = paste0("Gini = ", hist_gini, ", D95 = ", hist_d95),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

recentPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = recent)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = "Recent Distribution (2020s)",
    subtitle = paste0("Gini = ", recent_gini, ", D95 = ", recent_d95, 
                      " - Range Contraction Example"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

# Simulate reduction in number of schools
# Many schools initially
many_schools=data.frame(
  lon = c(-75, -65, -60, -45, -35, -70, -50, -40),
  lat = c(15, 5, 25, 10, 30, -5, -15, 20)
)

many_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(many_schools)) {
    dist=sqrt((grid_points$lon[i] - many_schools$lon[j])^2 + 
                (grid_points$lat[i] - many_schools$lat[j])^2)
    many_biomass[i]=many_biomass[i] + 90 * exp(-0.15 * dist)
  }
}
many_biomass=many_biomass + rnorm(length(many_biomass), 0, 3)
many_biomass[many_biomass <minCell]=0

# Fewer schools (reduced number but similar size)
few_schools=data.frame(
  lon = c(-75, -60, -45),
  lat = c(15, 25, 10)
)

few_biomass=numeric(nrow(grid_points))
for(i in 1:nrow(grid_points)) {
  for(j in 1:nrow(few_schools)) {
    dist=sqrt((grid_points$lon[i] - few_schools$lon[j])^2 + 
                (grid_points$lat[i] - few_schools$lat[j])^2)
    few_biomass[i]=few_biomass[i] + 90 * exp(-0.15 * dist)
  }
}
few_biomass=few_biomass + rnorm(length(few_biomass), 0, 3)
few_biomass[few_biomass <minCell]=0

# Calculate indicators
many_gini=round(calcGini(many_biomass), 2)
many_d95=round(calcD95(many_biomass), 2)
few_gini=round(calcGini(few_biomass), 2)
few_d95=round(calcD95(few_biomass), 2)

# Add to data
grid_points$many_schools=many_biomass
grid_points$few_schools=few_biomass

# Plot school reduction
manyPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = many_schools)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = "Many Schools",
    subtitle = paste0("Gini = ", many_gini, ", D95 = ", many_d95),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

fewPlot=ggplot() +
  geom_tile(data = grid_points, aes(x = lon, y = lat, fill = few_schools)) +
  geom_sf(data = world_map, fill = "darkgreen", color = "gray70", size = 0.2) +
  scale_fill_viridis(option = "plasma", name = "Biomass") +
  coord_sf(xlim = c(-80, 0), ylim = c(-20, 40)) +
  labs(
    title = "Fewer Schools",
    subtitle = paste0("Gini = ", few_gini, ", D95 = ", few_d95,
                     ""),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )


evenPlot



# Distribution patterns
tuna_swordfish=swordfishPlot+tunaPlot

# Range contraction
range_contraction=historicalPlot+recentPlot

# School reduction
school_reduction=manyPlot+fewPlot

# Display plots separately
tuna_swordfish
range_contraction
school_reduction


