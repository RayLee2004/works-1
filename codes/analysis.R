# Editor: Li Ruihong (Department of Environmental Science, Chongqing University)
# this file includes: (1-2) richness, diversity, synchrony and stability of taxonomic groups 
#                     (3) meta-web metrics (4) pSEM (5) model compensating effects (6) LLM

# Load packages
packages <- c(
  'openxlsx', 'tidyr', 'vegan', 'dplyr', 'codyn', 'tibble', 
  'igraph', 'bipartite', 'purrr', 'Hmisc', 'piecewiseSEM',
  'tidyverse','rstatix','car','lme4','lmerTest','multilevelTools',
  'extraoperators','JWileymisc','effectsize','influence.ME','GGally',
  'MuMIn','sjstats'
)
lapply(packages, library, character.only = TRUE)

### Part.1 richness, diversity, synchrony and stability of each taxonomic group

# Read datasets
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(2:8)]
all_groups <- lapply(split(data, data$group), function(df) df[, -1])
rm(list = setdiff(ls(), 'all_groups'))

# Convert absolute abundance to relative abundance
standardize <- function(x) {
  species <- as.data.frame(x$species)
  colnames(species)[1] <- 'species'
  x <- x[, -1]
  col_sums <- colSums(x)
  x <- sweep(x, 2, col_sums, FUN = '/')
  x <- cbind(species, x)
  return(x)
}
all_groups_standardized <- lapply(all_groups, standardize)

# Transpose data tables
transpose <- function(x) {
  current_df <- as.data.frame(x)
  col_names  <- current_df[, 1]
  current_df <- as.data.frame(t(current_df[, -1]))
  colnames(current_df) <- col_names
  return(current_df)
}
all_groups_transposed <- lapply(all_groups, transpose)

# Calculate species richness
richness <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    richness_values <- vegan::specnumber(x[[y]])
    temp_df <- data.frame(richness = as.numeric(richness_values))
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}
richness_single_groups <- richness(all_groups_transposed)

# Calculate Shannon diversity
diversity <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    diversity_values <- vegan::diversity(x[[y]])
    temp_df <- data.frame(diversity = as.numeric(diversity_values))
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}
diversity_single_groups <- diversity(all_groups_transposed)

# Convert data table to long format
convert_to_longdata <- function(df) {
  df %>%
    pivot_longer(
      cols      = -1,
      names_to  = c('sampling_site', 'sampling_order'),
      names_sep = '_',
      values_to = 'abundance'
    ) %>%
    rename(species = 1) %>%
    mutate(sampling_order = as.integer(sampling_order))
}
all_groups_longdata <- lapply(all_groups, convert_to_longdata)

# Calculate community's synchrony
synchrony <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    synchrony_values <- codyn::synchrony(
      x[[y]],
      species.var   = 'species',
      time.var      = 'sampling_order',
      abundance.var = 'abundance',
      replicate.var = 'sampling_site'
    )
    colnames(synchrony_values)[2] <- substitute(y)
    result_df <- cbind(result_df, synchrony_values[2])
  }
  return(result_df)
}
synchrony_single_groups <- synchrony(all_groups_longdata)

# Calculate community's stability
stability <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    stability_values <- codyn::community_stability(
      x[[y]],
      time.var      = 'sampling_order',
      abundance.var = 'abundance',
      replicate.var = 'sampling_site'
    )
    colnames(stability_values)[2] <- substitute(y)
    result_df <- cbind(result_df, stability_values[2])
  }
  return(result_df)
}
stability_single_groups <- stability(all_groups_longdata)

### Part.2 richness, diversity, synchrony and stability of multiple groups

multi_groups <- convert_to_longdata(transpose(bind_rows(all_groups)))
richness_multi_groups  <- data.frame(multi_groups = vegan::specnumber(multi_groups_transposed))
diversity_multi_groups <- data.frame(multi_groups = vegan::diversity(multi_groups_transposed))
synchrony_multi_groups <- synchrony(list(multi_groups = multi_groups_longdata))
stability_multi_groups <- stability(list(multi_groups = multi_groups_longdata))

# Combine all results and export
community_metrics_collection <- list(
  richness  = cbind(richness_single_groups, richness_multi_groups),
  diversity = cbind(diversity_single_groups, diversity_multi_groups),
  synchrony = cbind(synchrony_single_groups, multi_groups = synchrony_multi_groups[, 2]),
  stability = cbind(stability_single_groups, multi_groups = stability_multi_groups[, 2])
)
write.xlsx(
  community_metrics_collection, 
  'C:/Users/23926/Desktop/works/#1 datasets and codes/community_metrics_collection.xlsx'
)

### Part.3 meta-web metrics

# Read datasets
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(1:8)]
species_col <- 'species'
for (colname in names(data)[-1]) {
  tmp <- data[data[[colname]] != 0, species_col, drop = FALSE]
  assign(colname, tmp, envir = .GlobalEnv)
}
rm(data, tmp, colname, species_col)

# Export occurrence list from loaded objects
all_objects <- ls()
occurrences <- list()
for (x in all_objects) {
  if (is.data.frame(get(x))) {
    occurrences[[x]] <- get(x)
  }
}
occurrences <- lapply(occurrences, as.data.frame)
rm(list = setdiff(ls(), 'occurrences'))

# Read original network
original_network <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/trophic interaction 0-1 adjacent matrix.xlsx'
)[, -(1:2)]
diag(original_network) <- 0
original_network <- as.data.frame(original_network)
rownames(original_network) <- colnames(original_network)

# Initialize lists to store graphs and matrices of webs
dfs      <- list()
matrices <- list()
graphs   <- list()

# Extract subnetwork data tables
get_df <- function(x) {
  occurrence_list <- x[, ncol(x)]
  current_df <- original_network[occurrence_list, occurrence_list]
  rownames(current_df) <- colnames(current_df)
  dfs <<- c(dfs, list(current_df))
}
lapply(occurrences, get_df)

# Extract matrices and graphs
get_matrix_graph <- function(x) {
  matrix <- as.matrix(x)
  graph  <- igraph::graph_from_adjacency_matrix(
    matrix, 
    mode = 'directed', 
    weighted = TRUE, 
    diag = FALSE
  )
  # Retain only connected nodes
  none_isolated_nodes <- igraph::V(graph)[igraph::degree(graph) > 0]
  graph2 <- igraph::induced_subgraph(graph, vids = none_isolated_nodes)
  matrices <<- c(matrices, list(matrix))
  graphs   <<- c(graphs, list(graph2))
}
lapply(dfs, get_matrix_graph)

# Name the network objects
names(dfs)      <- names(occurrences)
names(matrices) <- names(occurrences)
names(graphs)   <- names(occurrences)

# Calculate path length, no.nodes, connectance, modularity, nestedness, robustness and vulnerability of each network
cal_chain_metrics <- function(x) {
  result_df <- data.frame(
    sampling_site_order = character(0),
    connectance         = numeric(0),
    max_path_length     = numeric(0),
    mean_path_length    = numeric(0),
    no_nodes            = numeric(0)
  )
  for (y in names(x)) {
    current_graph <- x[[y]]
    metrics <- data.frame(
      sampling_site_order = y,
      connectance         = edge_density(current_graph),
      max_path_length     = diameter(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      mean_path_length    = mean_distance(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      no_nodes            = sum(igraph::degree(current_graph) > 0)
    )
    result_df <- rbind(result_df, metrics)
  }
  return(result_df)
}

cal_structure_metrics <- function(x, y) {
  result_df1 <- data.frame(sampling_site_order = character(0), modularity = numeric(0))
  for (z in names(y)) {
    current_community <- cluster_walktrap(y[[z]])
    mod_df <- data.frame(
      sampling_site_order = z,
      modularity          = modularity(current_community, weights = E(current_community)$weight)
    )
    result_df1 <- rbind(result_df1, mod_df)
  }
  
  result_df2 <- data.frame(nestedness = numeric(0), vulnerability = numeric(0), robustness = numeric(0))
  
  vulnerability <- function(n) {
    gen_vul <- networklevel(n, index = 'vulnerability', weighted = TRUE)
    gen_vul[!grepl('generality', names(gen_vul))]
  }
  
  for (w in names(x)) {
    current_matrix <- x[[w]]
    net_metrics <- data.frame(
      nestedness    = nested(current_matrix, 'wine'),
      vulnerability = vulnerability(current_matrix),
      robustness    = robustness(
        second.extinct(current_matrix, participant = 'higher', method = 'random', nrep = 10, details = FALSE)
      )
    )
    result_df2 <- rbind(result_df2, net_metrics)
  }
  
  return(cbind(result_df1, result_df2))
}

# Combine all results and export
chain_metrics <- cal_chain_metrics(graphs)
network_indices <- cal_structure_metrics(matrices, graphs)
network_metrics_collection <- cbind(chain_metrics, network_indices[, -1])
write.xlsx(
  network_metrics_collection, 
  'C:/Users/23926/Desktop/works/#1 datasets and codes/network_metrics_collection.xlsx'
)

### Part.4 piecewise structural equation models

# Read datasets
rm(list = ls())
data_average <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_average.xlsx')
data_original <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_original.xlsx')

# Model 1: relation of richness, diversity, synchrony and stability
model1 <- psem(
  lm(invertebrate_richness ~ fish_richness + fish_diversity, data_average),
  lm(insect_richness ~ invertebrate_richness + fish_diversity, data_average),
  lm(algae_richness ~ zooplankton_richness + insect_richness, data_average),
  lm(zooplankton_richness ~ fish_richness + invertebrate_richness + insect_richness, data_average),
  lm(fungi_richness ~ zooplankton_richness, data_average),
  lm(bacteria_richness ~ zooplankton_richness + insect_richness, data_average),
  lm(fish_diversity ~ fish_richness, data_average),
  lm(invertebrate_diversity ~ invertebrate_richness + fish_richness, data_average),
  lm(insect_diversity ~ insect_richness + invertebrate_richness + fish_richness, data_average),
  lm(algae_diversity ~ algae_richness, data_average),
  lm(zooplankton_diversity ~ zooplankton_richness + invertebrate_richness + insect_richness, data_average),
  lm(fungi_diversity ~ fungi_richness + zooplankton_richness, data_average),
  lm(bacteria_diversity ~ bacteria_richness + zooplankton_richness, data_average),
  lm(fish_synchrony ~ fish_diversity + algae_synchrony, data_average),
  lm(invertebrate_synchrony ~ invertebrate_diversity + fish_richness, data_average),
  lm(insect_synchrony ~ insect_diversity + invertebrate_synchrony, data_average),
  lm(algae_synchrony ~ algae_diversity + zooplankton_synchrony + insect_synchrony, data_average),
  lm(zooplankton_synchrony ~ zooplankton_diversity + invertebrate_synchrony + insect_synchrony, data_average),
  lm(fungi_synchrony ~ fungi_diversity + zooplankton_synchrony, data_average),
  lm(bacteria_synchrony ~ bacteria_diversity + zooplankton_synchrony, data_average),
  lm(fish_stability ~ fish_synchrony + zooplankton_synchrony + insect_synchrony, data_average),
  lm(invertebrate_stability ~ invertebrate_synchrony + zooplankton_synchrony, data_average),
  lm(insect_stability ~ insect_synchrony + algae_synchrony, data_average),
  lm(zooplankton_stability ~ zooplankton_synchrony, data_average),
  lm(algae_stability ~ algae_synchrony, data_average),
  lm(fungi_stability ~ fungi_synchrony, data_average),
  lm(bacteria_stability ~ bacteria_synchrony, data_average),
  lm(multi_richness ~ fish_richness + zooplankton_richness + algae_richness + bacteria_richness + invertebrate_richness, data_average),
  lm(multi_diversity ~ multi_richness, data_average),
  lm(multi_synchrony ~ multi_diversity + fish_synchrony, data_average),
  lm(multi_stability ~ multi_synchrony + algae_stability, data_average)
)
# Model1 results
summary(model1)

# Model 2: richness impact on network
model2 <- psem(
  lm(mean_path_length ~ fish_richness + zooplankton_richness + invertebrate_richness, data_original),
  lm(connectance ~ mean_path_length + fish_richness + zooplankton_richness, data_original),
  lm(nestedness ~ mean_path_length + connectance + fish_richness + zooplankton_richness + invertebrate_richness, data_original),
  lm(modularity ~ mean_path_length + connectance + fish_richness + zooplankton_richness + nestedness, data_original),
  lm(robustness ~ connectance + mean_path_length + nestedness + modularity + fish_richness + zooplankton_richness + invertebrate_richness, data_original),
  lm(vulnerability ~ modularity + robustness + connectance + nestedness + modularity + mean_path_length + fish_richness + invertebrate_richness + zooplankton_richness, data_original)
)
# Model2 results
summary(model2)

# Model 3: diversity impact on network
model3 <- psem(
  lm(mean_path_length ~ fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
  lm(connectance ~ mean_path_length + fish_diversity + zooplankton_diversity, data_original),
  lm(nestedness ~ mean_path_length + connectance + fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
  lm(modularity ~ mean_path_length + connectance + fish_diversity + zooplankton_diversity + invertebrate_diversity + nestedness, data_original),
  lm(robustness ~ connectance + mean_path_length + nestedness + modularity + fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
  lm(vulnerability ~ modularity + robustness + connectance + nestedness + modularity + mean_path_length + fish_diversity + invertebrate_diversity + zooplankton_diversity, data_original)
)
# Model3 results
summary(model3)

### Part.5 model compensating effects (model 'yellow' as an example)

# Read datasets (replace module name here.)
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(1:8)]

# Delete species belonging to module nodes from dataset
deleted_species <- read.table(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/wgcna results/wgcna_network/yellow.network.nodes.txt'
)[, 1]
data <- data[!data[, 1] %in% deleted_species, ]

# Repeat meta-web metrics calculation (as in Part.3)
species_col <- 'species'
for (colname in names(data)[-1]) {
  tmp <- data[data[[colname]] != 0, species_col, drop = FALSE]
  assign(colname, tmp, envir = .GlobalEnv)
}
rm(data, tmp, colname, species_col)

all_objects <- ls()
occurrences <- list()
for (x in all_objects) {
  if (is.data.frame(get(x))) {
    occurrences[[x]] <- get(x)
  }
}
occurrences <- lapply(occurrences, as.data.frame)
rm(list = setdiff(ls(), 'occurrences'))

original_network <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/trophic interaction 0-1 adjacent matrix.xlsx'
)[, -(1:2)]
diag(original_network) <- 0
original_network <- as.data.frame(original_network)
rownames(original_network) <- colnames(original_network)

dfs      <- list()
matrices <- list()
graphs   <- list()

get_df <- function(x) {
  occurrence_list <- x[, ncol(x)]
  current_df <- original_network[occurrence_list, occurrence_list]
  rownames(current_df) <- colnames(current_df)
  dfs <<- c(dfs, list(current_df))
}
lapply(occurrences, get_df)

get_matrix_graph <- function(x) {
  matrix <- as.matrix(x)
  graph  <- igraph::graph_from_adjacency_matrix(
    matrix, 
    mode = 'directed', 
    weighted = TRUE, 
    diag = FALSE
  )
  none_isolated_nodes <- igraph::V(graph)[igraph::degree(graph) > 0]
  graph2 <- igraph::induced_subgraph(graph, vids = none_isolated_nodes)
  matrices <<- c(matrices, list(matrix))
  graphs   <<- c(graphs, list(graph2))
}
lapply(dfs, get_matrix_graph)

names(dfs)      <- names(occurrences)
names(matrices) <- names(occurrences)
names(graphs)   <- names(occurrences)

cal_chain_metrics <- function(x) {
  result_df <- data.frame(
    sampling_site_order = character(0),
    connectance         = numeric(0),
    max_path_length     = numeric(0),
    mean_path_length    = numeric(0),
    no_nodes            = numeric(0)
  )
  for (y in names(x)) {
    current_graph <- x[[y]]
    metrics <- data.frame(
      sampling_site_order = y,
      connectance         = edge_density(current_graph),
      max_path_length     = diameter(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      mean_path_length    = mean_distance(current_graph, directed = TRUE, weights = E(current_graph)$weight),
      no_nodes            = sum(igraph::degree(current_graph) > 0)
    )
    result_df <- rbind(result_df, metrics)
  }
  return(result_df)
}

cal_structure_metrics <- function(x, y) {
  result_df1 <- data.frame(sampling_site_order = character(0), modularity = numeric(0))
  for (z in names(y)) {
    current_community <- cluster_walktrap(y[[z]])
    mod_df <- data.frame(
      sampling_site_order = z,
      modularity          = modularity(current_community, weights = E(current_community)$weight)
    )
    result_df1 <- rbind(result_df1, mod_df)
  }
  
  result_df2 <- data.frame(nestedness = numeric(0), vulnerability = numeric(0), robustness = numeric(0))
  
  vulnerability <- function(n) {
    gen_vul <- networklevel(n, index = 'vulnerability', weighted = TRUE)
    gen_vul[!grepl('generality', names(gen_vul))]
  }
  
  for (w in names(x)) {
    current_matrix <- x[[w]]
    net_metrics <- data.frame(
      nestedness    = nested(current_matrix, 'wine'),
      vulnerability = vulnerability(current_matrix),
      robustness    = robustness(
        second.extinct(current_matrix, participant = 'higher', method = 'random', nrep = 10, details = FALSE)
      )
    )
    result_df2 <- rbind(result_df2, net_metrics)
  }
  
  return(cbind(result_df1, result_df2))
}

chain_metrics <- cal_chain_metrics(graphs)
network_indices <- cal_structure_metrics(matrices, graphs)
network_metrics_collection <- cbind(chain_metrics, network_indices[, -1])

write.xlsx(network_metrics_collection, 'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/module:yellow_network_metrics.xlsx')

### Part.6 fit linear mixed-effects models (LLM) of metrics in study and environmental factors

# Read datasets 
rm(list = ls())
env_data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/environmental quality.xlsx')
metrics_data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_original.xlsx')
data <- cbind(env_data,metrics_data[,-1])
# Standardize numeric vectors
data[, 3:ncol(data)] <- scale(data[, 3:ncol(data)])

## Environment - Algae_richness
model <- lmer(algae_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Algae_diversity
model <- lmer(algae_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Bacteria_richness
model <- lmer(bacteria_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Bacteria_diversity
model <- lmer(bacteria_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Fungi_richness
model <- lmer(fungi_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Fungi_diversity
model <- lmer(fungi_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Zooplankton_richness
model <- lmer(zooplankton_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Zooplankton_diversity
model <- lmer(zooplankton_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Insect_richness
model <- lmer(insect_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Insect_diversity
model <- lmer(insect_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Benthic invertebrate_richness
model <- lmer(invertebrate_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Benthic invertebrate_diversity
model <- lmer(invertebrate_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Fish_richness
model <- lmer(fish_richness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Fish_diversity
model <- lmer(fish_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Multi_richness
model <- lmer(multi_richness~Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Multi_diversity
model <- lmer(multi_diversity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Mean Path
model <- lmer(mean_path_length ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Connectance
model <- lmer(connectance ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Modularity
model <- lmer(modularity ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Nestedness
model <- lmer(nestedness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Robustness
model <- lmer(robustness ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

## Environment - Vulnerability
model <- lmer(vulnerability ~ Temperature + pH + DO + COD + NH3N + TP + TN + EC+ TUB + (1|year),data)
summary(model)
car::vif(model) # Collinearity of the independent variables
performance::r2(model)# R2 statistic value

