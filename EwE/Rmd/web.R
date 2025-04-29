library(igraph)
library(NetIndices)

# Read the CSV (skip empty trailing columns)
df=read.csv("../data/Table 1.csv",check.names=FALSE,
                                                   stringsAsFactors=FALSE)

df[is.na(df)]=0

dmns=list(predator=df[1:42,2],prey=dimnames(df)[[2]][-(1:2)])
df  =df[1:42,-c(1:2)]
mat =as.matrix(df,nrow=42,ncol=42,dimnames=dmns)
rownames(mat)=dmns[[1]]
names(dimnames(mat))=c("Predator","Prey")

rownames(mat)=dmns[[1]]
colnames(mat)=dmns[[1]]

g = graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)

# Assign node names (igraph uses the row names by default, but you can set explicitly)
V(g)$name = dmns[[1]]

V(g)$name

plot(
  g,
  vertex.label = V(g)$name,  # Use names as labels
  vertex.size = 10,
  edge.arrow.size = 0.3,
  main = "Trophic Web")


yellowfin = "Yellowfin"
connected_nodes = neighbors(g, v = yellowfin, mode = "all")$name  # Prey and predators
all_highlighted = c(yellowfin, connected_nodes)

# 3. Find all edges connected to Yellowfin
connected_edges = E(g)[.inc(yellowfin)]  # Updated .inc() syntax


# 4. Set colors
V(g)$color = "gray90"  # Default color
V(g)$color[V(g)$name == yellowfin] = "gold"  # Yellowfin itself
V(g)$color[V(g)$name %in% connected_nodes] = "lightblue"  # Connected nodes

E(g)$color = "gray80"  # Default edge color
E(g)$color[connected_edges] = "firebrick"  # Highlighted edges
E(g)$width = ifelse(E(g) %in% connected_edges, 2, 0.8)  # Thicker highlighted edges

# 5. Create trophic level layout
trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)

# 6. Plot
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Yellowfin Network Connections")


################################################################################
albacore = "Albacore"
connected_nodes = neighbors(g, v = albacore, mode = "all")$name  # Prey and predators
all_highlighted = c(albacore, connected_nodes)

# 3. Find all edges connected to Albacore
connected_edges = E(g)[.inc(albacore)]  # Updated .inc() syntax


# 4. Set colors
V(g)$color = "gray90"  # Default color
V(g)$color[V(g)$name == albacore] = "gold"  # Albacore itself
V(g)$color[V(g)$name %in% connected_nodes] = "lightblue"  # Connected nodes

E(g)$color = "gray80"  # Default edge color
E(g)$color[connected_edges] = "firebrick"  # Highlighted edges
E(g)$width = ifelse(E(g) %in% connected_edges, 2, 0.8)  # Thicker highlighted edges

# 5. Create trophic level layout
trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)

# 6. Plot
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Albacore Network Connections")


################################################################################
bluefin = "Bluefin"
connected_nodes = neighbors(g, v = bluefin, mode = "all")$name  # Prey and predators
all_highlighted = c(bluefin, connected_nodes)

# 3. Find all edges connected to Bluefin
connected_edges = E(g)[.inc(bluefin)]  # Updated .inc() syntax


# 4. Set colors
V(g)$color = "gray90"  # Default color
V(g)$color[V(g)$name == bluefin] = "gold"  # Bluefin itself
V(g)$color[V(g)$name %in% connected_nodes] = "lightblue"  # Connected nodes

E(g)$color = "gray80"  # Default edge color
E(g)$color[connected_edges] = "firebrick"  # Highlighted edges
E(g)$width = ifelse(E(g) %in% connected_edges, 2, 0.8)  # Thicker highlighted edges

# 5. Create trophic level layout
trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)

# 6. Plot
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Bluefin Network Connections")


################################################################################
bigeye = "Bigeye"
connected_nodes = neighbors(g, v = bigeye, mode = "all")$name  # Prey and predators
all_highlighted = c(bigeye, connected_nodes)

# 3. Find all edges connected to Bigeye
connected_edges = E(g)[.inc(bigeye)]  # Updated .inc() syntax


# 4. Set colors
V(g)$color = "gray90"  # Default color
V(g)$color[V(g)$name == bigeye] = "gold"  # Bigeye itself
V(g)$color[V(g)$name %in% connected_nodes] = "lightblue"  # Connected nodes

E(g)$color = "gray80"  # Default edge color
E(g)$color[connected_edges] = "firebrick"  # Highlighted edges
E(g)$width = ifelse(E(g) %in% connected_edges, 2, 0.8)  # Thicker highlighted edges

# 5. Create trophic level layout
trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)

# 6. Plot
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Bigeye Network Connections")



################################################################################
skipjack = "Skipjack"
connected_nodes = neighbors(g, v = skipjack, mode = "all")$name  # Prey and predators
all_highlighted = c(skipjack, connected_nodes)

# 3. Find all edges connected to Skipjack
connected_edges = E(g)[.inc(skipjack)]  # Updated .inc() syntax


# 4. Set colors
V(g)$color = "gray90"  # Default color
V(g)$color[V(g)$name == skipjack] = "gold"  # Skipjack itself
V(g)$color[V(g)$name %in% connected_nodes] = "lightblue"  # Connected nodes

E(g)$color = "gray80"  # Default edge color
E(g)$color[connected_edges] = "firebrick"  # Highlighted edges
E(g)$width = ifelse(E(g) %in% connected_edges, 2, 0.8)  # Thicker highlighted edges

# 5. Create trophic level layout
trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)

# 6. Plot
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Skipjack Network Connections")



################################################################################


# List your focal species
focal_species = c("Yellowfin", "Skipjack", "Albacore", "Bigeye")


# Find all nodes connected to any focal species (prey and predators)
connected_nodes = unique(unlist(
  lapply(focal_species, function(sp) neighbors(g, sp, mode = "all")$name)
))
all_highlighted_nodes = unique(c(focal_species, connected_nodes))

# Find all edges connected to any focal species
connected_edges = unlist(
  lapply(focal_species, function(sp) E(g)[.inc(sp)]), use.names = FALSE
)

# Default colors
V(g)$color = "gray90"
E(g)$color = "gray80"
E(g)$width = 0.8

# Highlight focal species
V(g)$color[V(g)$name %in% focal_species] = "gold"

# Highlight directly connected nodes
V(g)$color[V(g)$name %in% setdiff(connected_nodes, focal_species)] = "lightblue"

# Highlight connected edges
E(g)$color[connected_edges] = "firebrick"
E(g)$width[connected_edges] = 2


trophic_levels = TrophInd(mat)$TL
layout = cbind(runif(nrow(mat)), trophic_levels)


plot(
  g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.size = 12,
  vertex.label.cex = 0.7,
  edge.arrow.size = 0.4,
  main = "Highlighted Tuna and Connections"
)

################################################################################

library(igraph)
library(NetIndices)

# 2. Create weighted graph
g = graph_from_adjacency_matrix(
  mat,
  mode = "directed",
  weighted = TRUE
)

# 3. Set visual properties
# Base styles
V(g)$size = 8
V(g)$label.cex = 0.7
V(g)$frame.color = NA
E(g)$arrow.size = 0.3

# Highlight tuna species and connections
focal_species = c("Yellowfin", "Skipjack", "Albacore", "Bigeye")
connected_edges = unlist(lapply(focal_species, function(sp) E(g)[.inc(sp)]))
connected_nodes = unique(unlist(neighbors(g, focal_species, mode = "all")))

# Color coding
V(g)$color = "gray90"
V(g)$color[V(g)$name %in% focal_species] = "#FFD700"  # Gold
V(g)$color[V(g)$name %in% connected_nodes] = "#87CEEB"  # Sky blue

E(g)$color = "gray80"
E(g)$color[connected_edges] = "firebrick"

# Edge widths scaled to weights
E(g)$width = sqrt(E(g)$weight) * 3  # Adjust scaling factor as needed

# 4. Trophic level layout
trophic_levels = TrophInd(mat)$TL
set.seed(123)  # For reproducible layout
layout = cbind(
  x = runif(vcount(g), min = -1, max = 1),
  y = trophic_levels
)

# 5. Plot
par(mar = c(0,0,1,0))
plot(g, 
     layout = layout,
     main = "North Atlantic Trophic Web with Tuna Connections",
     vertex.label.color = "black",
     edge.curved = 0.1)

################################################################################
# 1. Define your species groups
tunas = c("Yellowfin", "Skipjack", "Albacore", "Bigeye", "Bluefin")
sharks = c("Pelagic sharks")
billfish = c("Billfishes")


# 2. Highlight settings for each group
group_colors = list(
  tunas = list(node = "#FFD700", edge = "#FF4500"),  # Gold nodes, orange edges
            sharks = list(node = "#2E8B57", edge = "#228B22"),  # SeaGreen nodes, ForestGreen edges
            billfish = list(node = "brown", edge = "brown")  # SeaGreen nodes, ForestGreen edges
)

# 3. Find connections for both groups
highlight_data = lapply(list(tunas, sharks,billfish), function(group) {
  edges = unlist(lapply(group, function(sp) E(g)[.inc(sp)]))
  nodes = unique(c(group, unlist(neighbors(g, group, mode = "all"))))
  list(edges = edges, nodes = nodes)
})

# 4. Assign colors and styles
V(g)$color = "gray90"
E(g)$color = "gray80"
E(g)$width = sqrt(E(g)$weight) * 3  # Keep weighted edges

# Tunas group
V(g)$color[V(g)$name %in% highlight_data[[1]]$nodes] = group_colors$tunas$node
E(g)$color[highlight_data[[1]]$edges] = group_colors$tunas$edge

# Sharks group
V(g)$color[V(g)$name %in% highlight_data[[2]]$nodes] = group_colors$sharks$node
E(g)$color[highlight_data[[2]]$edges] = group_colors$sharks$edge

# Billfish group
V(g)$color[V(g)$name %in% highlight_data[[2]]$nodes] = group_colors$billfish$node
E(g)$color[highlight_data[[2]]$edges] = group_colors$billfish$edge


# 5. Enhanced layout with trophic levels
trophic_levels = TrophInd(mat)$TL
set.seed(123)
layout = cbind(runif(vcount(g)), trophic_levels)

# 6. Plot with legend
plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Trophic Web: Tunas vs Sharks/Whales")

# Add legend
legend("topleft",
       legend = c("Tunas", "Sharks", "Billfishes", "Connected", "Other"),
       pch = 21,
       pt.bg = c(group_colors$tunas$node, group_colors$sharks$node, group_colors$billfish$node, "lightblue", "gray90","pink"),
       col = "black",
       pt.cex = 2)

################################################################################

# Define your groups
tunas <- c("Yellowfin", "Skipjack", "Albacore", "Bigeye", "Bluefin")
sharks <- c("Pelagic sharks")
billfish <- c("Billfishes")

# Okabe-Ito colorblind-safe palette
group_colors <- list(
  tunas = list(node = "#0072B2", edge = "#005683"),      # Blue
#  sharks = list(node = "#E69F00", edge = "#B06F00"),     # Orange
  billfish = list(node = "#D55E00", edge = "#A13B00")    # Vermillion
)

# Find connections for all groups
highlight_data <- lapply(list(tunas, billfish, sharks)[1:2], function(group) {
  edges <- unlist(lapply(group, function(sp) E(g)[.inc(sp)]))
  nodes <- unique(c(group, unlist(neighbors(g, group, mode = "all"))))
  list(edges = edges, nodes = nodes)
})

# Assign colors and styles
V(g)$color <- "#DDDDDD"   # Light grey for others
E(g)$color <- "#BBBBBB"
E(g)$width <- sqrt(E(g)$weight) * 3

# Tunas 
V(g)$color[V(g)$name %in% highlight_data[[1]]$nodes] <- group_colors$tunas$node
E(g)$color[highlight_data[[1]]$edges] <- group_colors$tunas$edge

# Billfishes 
V(g)$color[V(g)$name %in% highlight_data[[2]]$nodes] <- group_colors$billfish$node
E(g)$color[highlight_data[[2]]$edges] <- group_colors$billfish$edge

# Sharks
#V(g)$color[V(g)$name %in% highlight_data[[2]]$nodes] <- group_colors$sharks$node
#E(g)$color[highlight_data[[3]]$edges] <- group_colors$sharks$edge

# Layout and plot
trophic_levels <- TrophInd(mat)$TL
set.seed(123)
layout <- cbind(runif(vcount(g)), trophic_levels)


# 1. Define your mesopelagic taxa (edit as needed)
mesopelagic <- c("Lg. Meso fish", "Sm. Meso fish")

# 2. Set all edge colors to default (e.g., gray)
E(g)$color <- "gray80"

# 3. Find edges where either end is a mesopelagic node
edge_list <- ends(g, E(g), names = TRUE)
meso_edges <- which(edge_list[,1] %in% mesopelagic | edge_list[,2] %in% mesopelagic)

# 4. Color those edges red
E(g)$color[meso_edges] <- "red"

plot(g, 
     layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 12,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.4,
     main = "Trophic Web: Tunas, Billfishes")

legend("topleft",
       legend = c("Tunas", "Sharks", "Billfishes", "Other"),
       pch = 21,
       pt.bg = c(group_colors$tunas$node, group_colors$sharks$node, group_colors$billfish$node, "#DDDDDD"),
       col = "black",
       pt.cex = 2)


