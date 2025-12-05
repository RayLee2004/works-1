#     CONTENT OF THIS R CODE SCRIPTS
#  1. Richness, Shannon diversity, Synchrony and Stability Calculating 
#  2. Meta-Network Metrics Calculating
#  3. Piecewise Structural Equation Models (pSEM) Fitting
#  4. Beta Diversity Decomposition of Taxonomic Groups
#  5. Model Compensating Effects Analysis
#  6. Model Compensating Effects Plots
#  7. Linear Mixed Effects Models (LLM) Fitting of Env/Eco Factors
#  8. RDA Analysis and Plot of Taxonomic Groups - Environment
#  9. PCoA & Ecotraj Analysis and Plot of Taxonomic Groups
# 10. PCoA & Path Length of Multiple Ecological Indices 
# 11. Procrustes Analysis for Pairings of Taxonomic Groups


# Editor: Li Ruihong (Department of Environmental Science, Chongqing University)

###### ================================================================= ######
######         Load the packages required for this script to run         ######
###### ================================================================= ######

pre.packages <- c(
  'openxlsx', 'tidyr', 'vegan', 'dplyr', 'codyn', 'tibble', 'igraph', 
  'bipartite', 'purrr', 'Hmisc', 'piecewiseSEM','tidyverse','rstatix',
  'car','lme4','lmerTest','multilevelTools','extraoperators','JWileymisc',
  'effectsize','influence.ME','GGally','MuMIn','sjstats','labdsv','betapart',
  'reshape2','ggtern','ggplot2', 'pairwiseAdonis', 'ecotraj','stringr', 
  'smacof','ggrepel') # package list

# check installed packages
installed <- installed.packages()[, "Package"] 
# download uninstalled packages
install.packages(pre.packages[!pre.packages %in% installed], dependencies = TRUE) 
# load all packages
lapply(pre.packages, library, character.only = TRUE) 

###### ================================================================= ######
######   1. Richness, Shannon diversity, Synchrony and Stability Calculating 
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

######   1-1 calculation for single taxonomic group
{
  # data pre-processing
  data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(2:8)] # read dataset
  all_groups <- lapply(split(data, data$group), function(df) df[, -1]) # separate different groups
  rm(list = setdiff(ls(), 'all_groups'))
  
  # convert absolute abundance to relative abundance
  standardize <- function(x) {
    species <- as.data.frame(x$species)
    colnames(species)[1] <- 'species'
    x <- x[, -1]
    x <- cbind(species, sweep(x, 2, colSums(x), FUN = '/'))
    return(x)
  }
  all_groups_standardized <- lapply(all_groups, standardize) # list of standardized data sets
  
  # transpose data tables
  transpose <- function(x) {
    current_df <- as.data.frame(x)
    col_names  <- current_df[, 1]
    current_df <- as.data.frame(t(current_df[, -1]))
    colnames(current_df) <- col_names
    return(current_df)
  }
  all_groups_transposed <- lapply(all_groups, transpose) # list of transposed data sets
  
  # calculate species richness
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
  richness_single_groups <- richness(all_groups_transposed) # list of richness results
  
  # calculate Shannon diversity
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
  diversity_single_groups <- diversity(all_groups_transposed) # list of diversity results
  
  # convert data table to long format table
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
  all_groups_longdata <- lapply(all_groups, convert_to_longdata) # list of data sets in long data format
  
  # calculate community's synchrony
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
  synchrony_single_groups <- synchrony(all_groups_longdata) # list of synchrony results
  
  # calculate community's stability
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
  stability_single_groups <- stability(all_groups_longdata) # list of synchrony results
}

######   1-2 calculation for multiple taxonomic group
{
  # Combine all data sets and conduct all analyses
  multi_groups <- convert_to_longdata(transpose(bind_rows(all_groups)))
  richness_multi_groups  <- data.frame(multi_groups = vegan::specnumber(multi_groups_transposed))
  diversity_multi_groups <- data.frame(multi_groups = vegan::diversity(multi_groups_transposed))
  synchrony_multi_groups <- synchrony(list(multi_groups = multi_groups_longdata))
  stability_multi_groups <- stability(list(multi_groups = multi_groups_longdata))
}

######   1-3 combine all results and export
{
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
}

###### ================================================================= ######
######   2. Meta-Network Metrics Calculating
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

######   2-1 data pre-processing: get list of species existence
{ 
  # separate different groups
  data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(1:8)]
  species_col <- 'species'
  for (colname in names(data)[-1]) {
    tmp <- data[data[[colname]] != 0, species_col, drop = FALSE]
    assign(colname, tmp, envir = .GlobalEnv)
  }
  rm(data, tmp, colname, species_col) # clear redundant objects
  
  # get existence list from current objects in the R-environment
  all_objects <- ls()
  occurrences <- list()
  for (x in all_objects) {
    if (is.data.frame(get(x))) {
      occurrences[[x]] <- get(x)
    }
  }
  occurrences <- lapply(occurrences, as.data.frame) # save the list
  rm(list = setdiff(ls(), 'occurrences')) # clear redundant objects
}

######   2-2 data pre-processing: get list of data.frame, matrix and graph of network
{
  # read original network
  original_network <- read.xlsx(
    'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/trophic interaction 0-1 adjacent matrix.xlsx'
  )[, -(1:2)]
  diag(original_network) <- 0 # ignore the cannibalism
  rownames(original_network) <- colnames(original_network)
  
  # initialize the lists
  dfs      <- list()
  matrices <- list()
  graphs   <- list()
  
  # extract dfs
  get_df <- function(x) {
    occurrence_list <- x[, ncol(x)]
    current_df <- original_network[occurrence_list, occurrence_list]
    rownames(current_df) <- colnames(current_df)
    dfs <<- c(dfs, list(current_df)) # save dfs
  }
  lapply(occurrences, get_df)
  
  # extract matrices and graphs
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
    matrices <<- c(matrices, list(matrix)) # save matrices
    graphs   <<- c(graphs, list(graph2)) # save graphs
  }
  lapply(dfs, get_matrix_graph)
  
  # rename the network objects
  names(dfs)      <- names(occurrences)
  names(matrices) <- names(occurrences)
  names(graphs)   <- names(occurrences)
}

######   2-3 calculation: path length, no.nodes, connectance,
######       modularity, nestedness, robustness and vulnerability
{
  # function to calculate connectance / path_length / no_nodes
  cal_network_1 <- function(x) {
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
  network_indices1 <- cal_network_1(graphs) # execute the analysis
  
  # function to calculate modularity / nestedness / vulnerability / robustness
  cal_network_2 <- function(x, y) {
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
  network_indices2 <- cal_network_2(matrices, graphs) # execute the analysis
}

######   2-4 combine all results and export
{
  network_metrics_collection <- cbind(network_indices1, network_indices2[, -1])
  write.xlsx(
    network_metrics_collection, 
    'C:/Users/23926/Desktop/works/#1 datasets and codes/network_metrics_collection.xlsx'
  )
}

###### ================================================================= ######
######   3. Piecewise Structural Equation Models (pSEM) Fitting
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# read data sets
data_average <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_average.xlsx') 
data_original <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_original.xlsx') 

# Model 1: relation of richness, diversity, synchrony and stability
{
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
}
summary(model1) # Model 1: summary

# Model 2: richness impact on network
{
  model2 <- psem(
    lm(mean_path_length ~ fish_richness + zooplankton_richness + invertebrate_richness, data_original),
    lm(connectance ~ mean_path_length + fish_richness + zooplankton_richness, data_original),
    lm(nestedness ~ mean_path_length + connectance + fish_richness + zooplankton_richness + invertebrate_richness, data_original),
    lm(modularity ~ mean_path_length + connectance + fish_richness + zooplankton_richness + nestedness, data_original),
    lm(robustness ~ connectance + mean_path_length + nestedness + modularity + fish_richness + zooplankton_richness + invertebrate_richness, data_original),
    lm(vulnerability ~ modularity + robustness + connectance + nestedness + modularity + mean_path_length + fish_richness + invertebrate_richness + zooplankton_richness, data_original)
  )
}
summary(model2) # Model 2: summary

# Model 3: diversity impact on network
{
  model3 <- psem(
    lm(mean_path_length ~ fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
    lm(connectance ~ mean_path_length + fish_diversity + zooplankton_diversity, data_original),
    lm(nestedness ~ mean_path_length + connectance + fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
    lm(modularity ~ mean_path_length + connectance + fish_diversity + zooplankton_diversity + invertebrate_diversity + nestedness, data_original),
    lm(robustness ~ connectance + mean_path_length + nestedness + modularity + fish_diversity + zooplankton_diversity + invertebrate_diversity, data_original),
    lm(vulnerability ~ modularity + robustness + connectance + nestedness + modularity + mean_path_length + fish_diversity + invertebrate_diversity + zooplankton_diversity, data_original)
  )
}
summary(model3) # Model 3: summary

###### ================================================================= ######
######    4. Beta Diversity Decomposition of Taxonomic Groups
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# define the list and order of the data sets
datasets <- c('fish', 'invertebrate', 'insect', 'zooplankton', 
              'fungi', 'bacteria', 'algae', 'multigroups')

# define the functions for batch analysis and plotting
analyze_beta_diversity <- function(dataset_name) {
  # read data set
  file_path <- paste0('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/', 
                      dataset_name, '.xlsx')
  data <- read.xlsx(file_path)
  rownames(data) <- data[, 1]
  data <- as.data.frame(t(data[, -1]))
  data[data > 0] <- 1
  
  # calculate repl, richdiff and total diversity
  beta_com <- beta.div.comp(data, coef = "J", quant = FALSE, save.abc = FALSE)
  # repl
  repl <- as.matrix(beta_com$repl)
  repl_df <- melt(repl)
  repl_df <- repl_df[repl_df$Var1 != repl_df$Var2, ]
  repl_df <- repl_df[as.numeric(repl_df$Var1) < as.numeric(repl_df$Var2), ]
  rownames(repl_df) <- NULL
  colnames(repl_df)[3] <- "Repl"
  # richdiff
  rich <- as.matrix(beta_com$rich)
  rich_df <- melt(rich)
  rich_df <- rich_df[rich_df$Var1 != rich_df$Var2, ]
  rich_df <- rich_df[as.numeric(rich_df$Var1) < as.numeric(rich_df$Var2), ]
  rownames(rich_df) <- NULL
  colnames(rich_df)[3] <- "RichDiff"
  # total diversity
  data_df <- merge(repl_df, rich_df)
  data_df$BDtotal <- data_df$Repl + data_df$RichDiff
  
  # remove duplicate pairings
  data_df <- data_df[
    substr(data_df$Var1, 1, regexpr("_", data_df$Var1) - 1) == 
      substr(data_df$Var2, 1, regexpr("_", data_df$Var2) - 1), 
  ]
  rownames(data_df) <- NULL
  
  # output beta-decomposition summary
  summary_stats <- data.frame(
    Dataset = dataset_name,
    `BDtotal (total diversity)` = mean(data_df$BDtotal, na.rm = TRUE),  
    `Repl (turnover)` = mean(data_df$Repl, na.rm = TRUE),              
    `RichDiff (nestedness)` = mean(data_df$RichDiff, na.rm = TRUE),
    `Repl/BDtotal` = mean(data_df$Repl / data_df$BDtotal, na.rm = TRUE),
    `RichDiff/BDtotal` = mean(data_df$RichDiff / data_df$BDtotal, na.rm = TRUE),
    `Number of sample pairs` = nrow(data_df),
    check.names = FALSE
  )
  
  # plot
  p <- ggtern(
    data = data_df,  
    aes(x = 1 - BDtotal,  
        y = Repl,         
        z = RichDiff)    
  ) +
    stat_density_tern(
      geom = "polygon",   
      aes(fill = ..level.., alpha = ..level..),
      base = "ilr",       
      n = 100,            
      bins = 20,           
      alpha = 0.225
    ) +
    scale_fill_viridis_c(
      guide = guide_colorbar(
        barwidth = 1,   
        barheight = 5, 
        label.size = 0.8,
        title.position = "left", 
        title.hjust = 0.5
      )
    ) +
    geom_point(
      size = 3.75, 
      color = "grey", 
      fill = "grey",   
      alpha = 0.8     
    ) +
    scale_alpha(range = c(0.2, 0.6), guide = "none") +
    theme_bw() +
    theme_showarrows() +
    labs(
      title = paste0("Beta Diversity Decomposition - ", 
                     tools::toTitleCase(dataset_name), " (Jaccard)"),
      x = "Similarity",
      y = "Repl",
      z = "RichDiff"
    ) +
    theme(
      tern.axis.arrow = element_line(size = 0.7, color = "black", alpha = 0.5),
      tern.axis.arrow.sep = 0.09,
      tern.axis.title.T = element_blank(),
      tern.axis.title.L = element_blank(),
      tern.axis.title.R = element_blank(),
      tern.axis.arrow.text.T = element_text(size = 26.5, vjust = 0.3),
      tern.axis.arrow.text.L = element_text(size = 26.5, vjust = 0.3),
      tern.axis.arrow.text.R = element_text(size = 26.5, vjust = 0.6),
      tern.axis.text = element_text(
        size = 19,          
        color = "black",    
        family = "Arial"    
      ),
      plot.title = element_blank(),
      legend.position = c(0.81, 0.81), 
      legend.box.just = "right",       
      legend.background = element_blank(), 
      legend.text = element_text(size = 16),  
      legend.title = element_text(size = 18)
    )
  
  # save plot
  ggsave(
    paste0("C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/jac_", 
           dataset_name, ".png"),
    p,
    width = 10,
    height = 10,
    dpi = 100
  )
  
  # return summary
  return(summary_stats)
}

# conduct batch analysis
all_results <- list()
for (dataset in datasets) {
  tryCatch({
    result <- analyze_beta_diversity(dataset)
    all_results[[length(all_results) + 1]] <- result
  })
}

# batch output summary
print(do.call(rbind, all_results))

###### ================================================================= ######
######   5. Model Compensating Effects Analysis (model 'yellow' as an example)
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# read data sets
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_dataset.xlsx')[, -(1:8)]

# delete species belonging to the specific module
deleted_species <- read.table(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/wgcna results/wgcna_network/yellow.network.nodes.txt'
)[, 1] # replace module's name (yellow) here
data <- data[!data[, 1] %in% deleted_species, ]

# repeat meta-web metrics calculation (as in Part 2-3)
{
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
  
  cal_network_1 <- function(x) {
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
  network_indices1 <- cal_network_1(graphs)
  
  cal_network_2 <- function(x, y) {
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
  network_indices2 <- cal_network_2(matrices, graphs)
}

# combine all results and export
yellow.module.comp <- cbind(network_indices1, network_indices2[, -1])
write.xlsx(yellow.module.comp, 'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/module:yellow.comp.effect.xlsx')

###### ================================================================= ######
######   6. Model Compensating Effects Plots
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# save horizontal line and CV values of each index
{
  metrics <- list(
    list(
      name = "Mean Path Length",
      baseline_values = c(1.162488958, 1.233940422, 1.468541456, 1.294841298),
      legend_labels = c(
        "2022" = "2022  CV=0.005", 
        "2023" = "2023  CV=0.007",
        "2024" = "2024  CV=0.019",  
        "2025" = "2025  CV=0.007"
      ),
      y_limits = c(1.14, 1.5)
    ),
    list(
      name = "Connectance",
      baseline_values = c(0.1373143333, 0.1222612264, 0.06741378756, 0.09476242467),
      legend_labels = c(
        "2022" = "2022  CV=0.019", 
        "2023" = "2023  CV=0.021",
        "2024" = "2024  CV=0.170",  
        "2025" = "2025  CV=0.028"
      ),
      y_limits = c(0.05, 0.16)
    ),
    list(
      name = "Modularity",
      baseline_values = c(0.05526768333, 0.08118164028, 0.124806275, 0.05914249222),
      legend_labels = c(
        "2022" = "2022  CV=0.038", 
        "2023" = "2023  CV=0.048",
        "2024" = "2024  CV=0.071",  
        "2025" = "2025  CV=0.057"
      ),
      y_limits = c(0.04, 0.145)
    ),
    list(
      name = "Nestedness",
      baseline_values = c(84.87469389, 85.62822939, 82.36244733),
      legend_labels = c(
        "2022" = "2022  CV=0.009", 
        "2023" = "2023  CV=0.008",
        "2024" = "2024  CV=0.007",  
        "2025" = "2025  CV=0.007"
      ),
      y_limits = c(80, 88)
    ),
    list(
      name = "Vulnerability",
      baseline_values = c(214.3640539, 119.9470125, 279.8278948),
      legend_labels = c(
        "2022" = "2022  CV=0.043", 
        "2023" = "2023  CV=0.054",
        "2024" = "2024  CV=0.043",  
        "2025" = "2025  CV=0.048"
      ),
      y_limits = c(105, 300)
    ),
    list(
      name = "Robustness",
      baseline_values = c(0.4424948867, 0.4356060338, 0.2837738379, 0.3604388324),
      legend_labels = c(
        "2022" = "2022  CV=0.009", 
        "2023" = "2023  CV=0.015",
        "2024" = "2024  CV=0.059",  
        "2025" = "2025  CV=0.016"
      ),
      y_limits = c(0.2, 0.455)
    )
  )
}

colors <- c("2022" = "#647ADD",
            "2023" = "#C0C0C0",
            "2024" = "#FFDF69",
            "2025" = "#F5664D")

# define the functions for batch analysis and plotting
plot_metric <- function(metric_info) {
  # clear redundant objects
  rm(list = setdiff(ls(), c("metric_info", "colors")))
  
  # read data sets
  file_path <- paste0(
    'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_',
    tolower(gsub(" ", "_", metric_info$name)), 
    '.xlsx'
  )
  data_ave <- read.xlsx(file_path, 1)
  data_range <- read.xlsx(file_path, 2)
  module_levels <- unique(data_ave$module)
  data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
  data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))
  
  # plot
  p <- ggplot() +
    geom_ribbon(
      data = data_range,
      aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
      alpha = 0.2 
    ) +
    geom_line(
      data = data_ave,
      aes(x = module_id, y = value, color = factor(year)),
      size = 1.25
    ) +
    geom_point(
      data = data_ave,
      aes(x = module_id, y = value, color = factor(year)),
      size = 3
    ) +
    scale_x_continuous(
      breaks = 1:length(module_levels), 
      labels = module_levels
    ) +
    scale_fill_manual(
      values = colors,
      labels = metric_info$legend_labels
    ) +
    scale_color_manual(
      values = colors,
      labels = metric_info$legend_labels
    ) +
    scale_y_continuous(
      breaks = metric_info$baseline_values, 
      labels = function(x) {  
        sprintf("%.3f", x)
      }
    ) +
    coord_cartesian(ylim = metric_info$y_limits) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
      axis.text.y = element_text(size = 15),                    
      axis.title = element_blank(),                                 
      plot.title = element_text(size = 34, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.grid = element_blank(),                       
      strip.text = element_text(size = 12),         
      legend.position = c(0.08, 0.05), 
      legend.justification = c(0.35, 0.15),  
      legend.background = element_rect(fill = NA, color = NA),  
      legend.key = element_blank(),     
      legend.title = element_text(size = 13), 
      legend.text = element_text(size = 13)   
    )
  
  # add horizontal line
  for (value in metric_info$baseline_values) {
    p <- p + geom_hline(
      yintercept = value, 
      color = "black", 
      linetype = "dashed", 
      size = 0.7, 
      alpha = 0.35
    )
  }
  
  # save plot
  output_path <- paste0(
    'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/', 
    metric_info$name, 
    '.png'
  )
  ggsave(output_path, plot = p, width = 9, height = 6, dpi = 300)
}

# conduct batch analysis
results <- list()
for (i in seq_along(metrics)) {
  metric_info <- metrics[[i]]
  tryCatch({
    result <- plot_metric(metric_info)
    results[[length(results) + 1]] <- result
  })
}

###### ================================================================= ######
######   7. Linear Mixed Effects Models (LLM) Fitting of Env/Eco Factors
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# read data sets 
env_data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/environmental quality.xlsx')
metrics_data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/indices_original.xlsx')
data <- cbind(env_data,metrics_data[,-1])

# standardize vectors
data[, 3:ncol(data)] <- scale(data[, 3:ncol(data)])

# define the function: batch execute the lmer model and output results
run_lmer_models <- function(response_vars, data) {
  model_results <- list()
  for (resp in response_vars) {
    cat("\n=====================================\n")
    cat(paste("## Environment -", resp), "\n")
    cat("=====================================\n")
    formula_str <- paste(resp, "~ Temperature + pH + DO + COD + NH3N + TP + TN + EC + TUB + (1|year)")
    model <- lmer(as.formula(formula_str), data = data)
    cat("\n--- Model Summary ---\n")
    print(summary(model)) 
    cat("\n--- VIF (Collinearity) ---\n")
    print(car::vif(model))
    cat("\n--- R² Value ---\n")
    print(performance::r2(model))
    model_results[[resp]] <- model
    cat("\n-------------------------------------\n\n")
  }
  return(model_results)
}

# define the response variables
{ 
  response_variables <- c(
    "algae_richness", "algae_diversity","bacteria_richness", 
    "bacteria_diversity","fungi_richness", "fungi_diversity",
    "zooplankton_richness", "zooplankton_diversity","insect_richness", 
    "insect_diversity","invertebrate_richness", "invertebrate_diversity",
    "fish_richness", "fish_diversity","multi_richness", 
    "multi_diversity","mean_path_length", "connectance",
    "modularity", "nestedness","robustness", "vulnerability"
  ) 
}

# conduct all analyses
all_models <- run_lmer_models(response_vars = response_variables, data = data)
# optional: call a certain model alone:
# summary(all_models[["fish_richness"]])

###### ================================================================= ######
######   8. RDA Analysis and Plot of Taxonomic Groups - Environment
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# define the list and order of the data sets
datasets <- c('fish', 'invertebrate', 'insect', 'zooplankton', 
              'fungi', 'bacteria', 'algae', 'multigroups')

mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')

# define the functions for batch analysis and plotting
analyze_rda <- function(dataset_name) {
  
  # read data sets
  species_file <- paste0('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/',
                         dataset_name, '.xlsx')
  data <- read.xlsx(species_file)
  rownames(data) <- data[, 1]
  data <- hellinger(as.data.frame(t(data[, -1])))
  groups_file <- 'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx'
  groups <- read.xlsx(groups_file)
  env_file <- 'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/environmental quality.xlsx'
  env <- read.xlsx(env_file)
  rownames(env) <- env[, 1]
  env <- as.data.frame(scale(as.data.frame(env[, -(1:2)])))
  
  # DCA analysis
  dca_result <- decorana(data)
  print(dca_result)
  
  # RDA analysis
  rda_result <- rda(data ~ ., data = env)  
  print(anova(rda_result))
  
  # print RDA results
  set.seed(123)
  print(anova(rda_result, by = "term"))
  print(anova(rda_result, by = "axis"))
  print(RsquareAdj(rda_result))
  
  # extract RDA coords
  rda_sites <- as.data.frame(rda_result$CCA$u)
  rda_sites$sample <- rownames(rda_sites)
  colnames(rda_sites)[1:2] <- c("RDA1", "RDA2")
  coords <- merge(rda_sites, groups[, c("sample", "group")], by = "sample")
  coords$group <- as.factor(coords$group)
  
  # extract env.factors' arrows
  env_arrows <- as.data.frame(rda_result$CCA$biplot[, 1:2])
  env_arrows$var <- rownames(env_arrows)
  colnames(env_arrows) <- c("RDA1", "RDA2", "var")
  
  # extract the explanation proportions of RDA1 and RDA2 
  rda1 <- round(rda_result$CCA$eig[1] / sum(rda_result$CCA$eig) * 100, 2)
  rda2 <- round(rda_result$CCA$eig[2] / sum(rda_result$CCA$eig) * 100, 2)
  
  # extract centroids
  centroids <- coords %>% 
    group_by(group) %>% 
    summarise(
      RDA1 = mean(RDA1),
      RDA2 = mean(RDA2)
    ) %>%
    mutate(label = 1:n())
  
  # plot
  p <- ggplot(coords, aes(x = RDA1, y = RDA2, color = group)) +
    stat_ellipse(aes(fill = group),
                 geom = "polygon",
                 alpha = 0.22,
                 linewidth = 0,
                 color = NA) +
    geom_point(shape = 16, fill = "white", size = 3.8, stroke = 1.1) +
    geom_segment(data = env_arrows,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 linewidth = 0.7,
                 color = "black") +
    geom_text_repel(data = env_arrows,
                    aes(x = RDA1, y = RDA2, label = var),
                    color = "red",
                    size = 6.5,
                    min.segment.length = 0,
                    box.padding = 0.6,
                    nudge_y = 0.04, 
                    segment.color = NA) +
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    xlab(paste0("RDA1 (", rda1, "%)")) +
    ylab(paste0("RDA2 (", rda2, "%)")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.1, color = "black"),
      axis.title = element_text(size = 30, color = "black"),
      axis.text  = element_text(size = 27.5, color = "black"),
      plot.title = element_text(size = 27.5, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA), 
      strip.text = element_text(size = 12),         
      scale_x_discrete(expand = expansion(mult = c(0, 0))),
      legend.position = 'none',
      plot.margin = unit(c(0.55, 0.55, 0.55, 0.55), 'cm')
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.6, linewidth = 0.9)
  
  # save plot
  output_file <- paste0('C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/rda_',
                        dataset_name, '.png')
  ggsave(output_file, plot = p, width = 8, height = 8, dpi = 600)
  
  # return rda summary
  summary_stats <- list(
    dataset = dataset_name,
    n_samples = nrow(data),
    n_vars = ncol(data),
    rda1_explained = rda1,
    rda2_explained = rda2,
    total_explained = round(sum(rda_result$CCA$eig) / sum(rda_result$CCA$eig) * 100, 2)
  )
  
  return(summary_stats)
}

# conduct batch analysis
all_results <- list()
for (dataset in datasets) {
  tryCatch({
    result <- analyze_rda(dataset)
    all_results[[length(all_results) + 1]] <- result
  })
}

# batch output summary
print(do.call(rbind, lapply(all_results, as.data.frame)))

###### ================================================================= ######
######   9. PCoA & Ecotraj Analysis and Plot of Taxonomic Groups
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

######   9-1 PCoA & Ecotraj Analysis and Plot:  algae
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/algae.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_algae.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-2 PCoA & Ecotraj Analysis and Plot:  bacteria
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/bacteria.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_bacteria.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-3 PCoA & Ecotraj Analysis and Plot:  fungi
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/fungi.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_fungi.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-4 PCoA & Ecotraj Analysis and Plot:  zooplankton
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/zooplankton.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_zooplankton.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-5 PCoA & Ecotraj Analysis and Plot:  insect
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/insect.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_insect.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-6 PCoA & Ecotraj Analysis and Plot:  benthic invertebrate
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/invertebrate.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_benthic invertebrate.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-7 PCoA & Ecotraj Analysis and Plot:  fish
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/fish.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_fish.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   9-8 PCoA & Ecotraj Analysis and Plot:  multigroups
{
  ######   analysis of PCoA & Ecotraj 
  {
    # data pre-processing
    data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/multigroups.xlsx')
    groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
    rownames(data) <- data[, 1]
    data <- hellinger(as.data.frame(t(data[, -1]))) # standardize method: Hellinger 
    d <- as.matrix(vegdist(data, method = 'bray')) # distance method: Bray-Curtis 
    
    # extract 'entities'
    entities <- sub('_.*', '', groups$sample)
    
    # extract 'surveys'
    surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
    
    # define trajectory and execute pcoa analysis
    trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
    
    # extract pcoa coordinates
    trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
    coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
    coords$Entity <- entities
    coords$Survey <- surveys
    coords$Entity <- as.factor(coords$Entity)
    coords$Survey <- as.factor(coords$Survey)
    colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
    
    # extract the explanation proportions of PC1 and PC2 
    pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
    pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  }
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_multigroups.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

###### ================================================================= ######
######   10. PCoA & Path Length of Multiple Ecological Indices 
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# read data sets 
{
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/multi_indices.xlsx',2)
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]
data <- scale(t(data[, -1])) # standardize method: Z-Score (to ensure the comparative rationality) 
d <- as.matrix(vegdist(data, method = 'euclidean',scale.=FALSE)) # distance method: Euclidean
}

######   10-1 PCoA & Ecotraj analysing and plotting
# (repeat analysing and plotting process (as in Part 6))
{
  # extract 'entities'
  entities <- sub('_.*', '', groups$sample)
  
  # extract 'surveys'
  surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18)) 
  
  # define trajectory and execute pcoa analysis
  trajectory_pcoa <- trajectoryPCoA(defineTrajectories(d, entities, surveys))
  
  # extract pcoa coordinates
  trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
  coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
  coords$Entity <- entities
  coords$Survey <- surveys
  coords$Entity <- as.factor(coords$Entity)
  coords$Survey <- as.factor(coords$Survey)
  colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
  
  # extract the explanation proportions of PC1 and PC2 
  pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
  pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]
  
  ######   plot pre-processing 
  {
    mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
    group_labels <- c('2022', '2023', '2024', '2025')
    
    # extract plotting data of hulls
    hulls <- coords %>%
      group_by(group) %>%
      slice(chull(dim1, dim2))
    
    # extract plotting data of centroids
    centroids <- coords %>%
      group_by(group) %>%
      summarise(
        dim1 = mean(dim1),
        dim2 = mean(dim2)
      ) %>%
      arrange(factor(group, levels = group_labels)) %>%
      mutate(label = 1:n())
    
    # extract plotting data of traj's segments
    traj_segments <- data.frame(
      x = centroids$dim1[-nrow(centroids)],  # start x
      y = centroids$dim2[-nrow(centroids)],  # start y
      xend = centroids$dim1[-1],             # end x
      yend = centroids$dim2[-1],             # end y
      group = group_labels[-length(group_labels)]
    )
  }
  
  ######   plot 
  {
    p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
      geom_polygon(
        data = hulls,
        aes(fill = as.factor(group), color = as.factor(group)),
        alpha = 0.4,
        linewidth = 0.4
      ) +
      geom_point(size = 0, shape = 21, stroke = 2.5, alpha = 0.5) +
      geom_segment(
        data = traj_segments,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = 'black',
        linewidth = 1,
        alpha = 0.8,
        arrow = arrow(length = unit(0.5, 'cm'), type = 'closed')
      ) +
      geom_point(
        data = centroids,
        aes(x = dim1, y = dim2),
        color = 'black',
        fill = mycol,
        shape = 21,
        size = 5.5,
        stroke = 1.15
      ) +
      geom_text(
        data = centroids,
        aes(x = dim1, y = dim2, label = label),
        color = 'black',
        size = 13,
        hjust = -0.2,
        vjust = -0.2
      ) +
      xlab(paste0('PCoA 1 (', pc1, '%)')) +
      ylab(paste0('PCoA 2 (', pc2, '%)')) +
      scale_color_manual(values = mycol, labels = group_labels, name = 'Year') +
      scale_fill_manual(values = mycol, labels = group_labels, name = 'Year') +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 33.5, color = 'black'),
        axis.text = element_text(size = 27.5, color = 'black'),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
      geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9)
  }
  
  ######   save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_multigroups.png', 
      plot = p, width = 8, height = 8, dpi = 600
    )
  }
}

######   10-2 loadings plot of PCoA Axis
{
  # plot pre-processing   
  coords_raw <- coords[,(1:2)]
  loadings <- data.frame(
    dim1 = apply(as.data.frame(data), 2, function(z) cor(z, coords_raw$dim1)),
    dim2 = apply(as.data.frame(data), 2, function(z) cor(z, coords_raw$dim2))
  )
  loadings$var <- rownames(loadings)
  loadings$label_clean <- gsub("_PC1$|_diversity$", "", loadings$var)
  loadings$category <- ifelse(grepl("_PC1$", loadings$var), "PCoA (PC1)",
                              ifelse(grepl("_diversity$", loadings$var), "Diversity",
                                     "Network Index"))
  loadings <- loadings %>%
    dplyr::mutate(
      hjust_label = case_when(
        category == "PCoA (PC1)" ~ 0.5,
        category == "Diversity" ~ -0.05,
        category == "Network Index" ~ 0.5
      ),
      vjust_label = case_when(
        category == "PCoA (PC1)" ~ 2.15,
        category == "Diversity" ~ 2.15,
        category == "Network Index" ~ -1.65
      )
    )
  fill_colors <- c("PCoA (PC1)" = "#F5664D", "Diversity" = "#FFDF67", "Network Index" = "#647ADD")
  
  p_load <- ggplot(loadings, aes(x = dim1, y = dim2)) +
    geom_segment(
      aes(x = 0, y = 0, xend = dim1, yend = dim2),
      linewidth = 0.5,
      color = "black",
      alpha = 0.5
    ) +
    geom_point(
      aes(fill = category),
      shape = 22,
      size = 8.5,
      color = "black",
      stroke = 1.025
    ) +
    geom_text_repel(
      aes(label = label_clean),
      size = 5.15,
      max.overlaps = 10,
      segment.size = 0.6,
      segment.color = "black",
      min.segment.length = 1,
      segment.alpha = 1,
      box.padding = unit(0.8, "lines"),
      force = 1.5,
      show.legend = FALSE
    ) +
    annotate("point",
             x = 0, y = 0,
             size = 3,
             color = "black",
             shape = 16,
             alpha = 0.8
    ) +
    xlab(paste0("PCoA1 (", pc1, "%)")) +
    ylab(paste0("PCoA2 (", pc2, "%)")) +
    scale_fill_manual(values = fill_colors) +
    theme_bw() +
    theme(
      panel.border = element_rect(linewidth = 1.2, color = "black"),
      panel.grid = element_blank(),
      axis.title = element_text(size = 33.5, color = "black"),
      axis.text = element_text(size = 27.5, color = "black"),
      legend.position = c(0.025, 0.975),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12.5),
      plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
    ) +
    geom_vline(xintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
    geom_hline(yintercept = 0, color = 'black', alpha = 0.6, linetype = 'dashed', linewidth = 0.9) +
    xlim(min(loadings$dim1) - 0.1, max(loadings$dim1) + 0.1) +
    ylim(min(loadings$dim2) - 0.1, max(loadings$dim2) + 0.1)
  
  # save plot
  ggsave(
    'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Multiindices_loadings.png',
    plot = p_load, width = 8, height = 8, dpi = 600
  )
}

######   10-3 scatter plotting of path length / distance
{
  # plot pre-processing
  {
    traj_length <- as.data.frame(trajectoryLengths(defineTrajectories(d, entities, surveys))[,(1:3)])
    colnames(traj_length) <- c("1-2", "2-3", "3-4")
    entity_list <- unique(entities)
  }
  
  # calculate distance between point 1/2 and point 4
  {
    dist_1_4 <- numeric(length(entity_list))
    dist_2_4 <- numeric(length(entity_list))
    for (i in 1:length(entity_list)) {
      entity <- entity_list[i]
      entity_indices <- which(entities == entity)
      entity_indices <- entity_indices[order(surveys[entity_indices])]
      dist_1_4[i] <- d[entity_indices[1], entity_indices[4]]
      dist_2_4[i] <- d[entity_indices[2], entity_indices[4]]
    }
    traj_length$"2-4" <- dist_2_4
    traj_length$"1-4" <- dist_1_4
  }
  
  # convert data table to long format table
  {
    data_long <- pivot_longer(
      traj_length,
      cols = c("1-2", "2-3", "3-4", "2-4", "1-4"),
      names_to = "Segment",
      values_to = "Length"
    )
  }
  
  # calculate mean values
  {
    mean_values <- data_long %>%
      group_by(Segment) %>%
      dplyr::summarize(Mean = mean(Length))
    data_long$Segment <- factor(data_long$Segment, levels = c("1-2", "2-3", "3-4", "2-4", "1-4"))
    mean_values$Segment <- factor(mean_values$Segment, levels = c("1-2", "2-3", "3-4", "2-4", "1-4"))
  }
  
  # anova test of the path lengths
  {
    anova_result <- aov(Length ~ Segment, data = data_long)
    summary(anova_result)
    TukeyHSD(anova_result)
  }
  
  # plot
  {
    segment_colors <- c(
      "1-2" = "#F5664D",
      "2-3" = "#FFDF67",
      "3-4" = "#647ADD",
      "2-4" = "#9C27B0",
      "1-4" = "#4CAF50"
    )
    p <- ggplot(data_long, aes(x = Segment, y = Length, fill = Segment)) +
      geom_jitter(
        position = position_jitter(width = 0.105, height = 0),
        color = "black",
        stroke = 0,
        size = 9,
        alpha = 0.43,
        shape = 21
      ) +
      geom_point(
        data = mean_values,
        aes(x = Segment, y = Mean, fill = Segment),
        color = "black",
        shape = 21,
        size = 16,
        stroke = 1.2
      ) +
      geom_text(
        data = mean_values,
        aes(x = Segment, y = Mean + 0.8, label = round(Mean, 2)),
        size = 9,
        fontface = 'bold',
        color = "black"
      ) +
      annotate(
        "text",
        label = "F=33.29  P=0.001",
        size = 9,
        color = "black",
        x = -Inf,
        y = -Inf,
        hjust = -0.04,
        vjust = -0.8
      ) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.45, size = 0.4, color = "black", show.legend = FALSE) +
      scale_color_manual(values = segment_colors) +
      scale_fill_manual(values = segment_colors) +
      labs(x = "Segment", y = "Path Length / Distance") +
      theme_bw() +
      theme(
        panel.border = element_rect(linewidth = 1.2, color = 'black'),
        panel.grid = element_blank(),
        axis.title = element_text(size = 35, color = 'black'),
        axis.text = element_text(size = 29.5, color = 'black'),
        axis.ticks.length = unit(0.4, "cm"),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = 'none',
        plot.margin = unit(c(0.15, 0.55, 0.15, 0.55), 'cm')
      ) +
      scale_y_continuous(limits = c(2.5, 12.5), breaks = seq(2.5, 12.5, by = 2.5))
  }
  
  # save plot
  {
    ggsave(
      'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/scatter_Multiindices.png',
      plot = p,
      width = 8,
      height = 8,
      dpi = 600
    )
  }
}

###### ================================================================= ######
######   11. Procrustes Analysis for Pairings of Taxonomic Groups
###### ================================================================= ######
rm(list = ls()) # clear the R-environment

# read data sets
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_datasets.Rdata')

# define the list and order of the data sets
{
  datasets_list <- list(
    fish = fish,
    invertebrate = invertebrate,
    insect = insect,
    zooplankton = zooplankton,
    fungi = fungi,
    bacteria = bacteria,
    algae = algae
  )
}

# define the functions for batch analysis and plotting
perform_procrustes <- function(name1, name2, data1, data2) {
  # standardize the data and extract the PCoA coordinates
  transpose <- function(y) {
    current_df <- as.data.frame(y)
    species    <- current_df[, 1]
    mat        <- t(current_df[, -1])
    colnames(mat) <- species
    return(as.data.frame(mat))
  }
  
  pcoa1 <- cmdscale(vegdist(decostand(transpose(data1), method = 'hellinger'), method = 'bray'))
  pcoa2 <- cmdscale(vegdist(decostand(transpose(data2), method = 'hellinger'), method = 'bray'))
  
  # extract Procrustes coords
  set.seed(123)
  pro_test <- protest(pcoa1, pcoa2, permutations = 999)
  Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
  Pro_rotation <- data.frame(pro_test$rotation)
  colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
  
  # plot
  plotdf_group <- data.frame(group = c(name1, name2))
  legend_data <- data.frame(
    group = c(name1, name2),
    color = c('#E41A1C', '#377EB8')
  )
  
  p <- ggplot(data = Pro_coords) +
    geom_point(aes(X1, X2, color = name1), size = 7, shape = 16) +
    geom_point(aes(Dim1, Dim2, color = name2), size = 7, shape = 16) +
    geom_segment(
      aes(x = X1, y = X2, xend = (X1 + Dim1)/2, yend = (X2 + Dim2)/2),
      arrow = arrow(length = unit(0, 'cm')),
      color = '#E41A1C',
      size = 0.9
    ) +
    geom_segment(
      aes(x = (X1 + Dim1)/2, y = (X2 + Dim2)/2, xend = Dim1, yend = Dim2),
      arrow = arrow(length = unit(0.2, 'cm')),
      color = '#377EB8',
      size = 0.9
    ) +
    scale_color_manual(
      values = setNames(c('#E41A1C', '#377EB8'), c(name1, name2)),
      name = 'Taxonomic Group',
      limits = c(name1, name2)
    ) +
    labs(x = 'Dimension 1', y = 'Dimension 2') +
    scale_x_continuous(labels = function(x) sprintf('%.2f', x)) +
    scale_y_continuous(labels = function(x) sprintf('%.2f', x)) +
    geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.95) +
    geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.95) +
    geom_abline(intercept = 0, slope = Pro_rotation[1,2]/Pro_rotation[1,1], size = 0.85) +
    geom_abline(intercept = 0, slope = Pro_rotation[2,2]/Pro_rotation[2,1], size = 0.85) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      axis.ticks.length = unit(0.4, 'lines'),
      axis.ticks = element_line(color = 'black'),
      axis.line = element_line(colour = 'black'),
      axis.title.x = element_text(colour = 'black', size = 37),
      axis.title.y = element_text(colour = 'black', size = 37),
      axis.text = element_text(colour = 'black', size = 26),
      legend.position = c(0.16, 0.9),
      legend.key.size = unit(1.5, 'cm'),
      legend.title = element_text(size = 23.5, face = 'bold'),
      legend.text = element_text(size = 23.5),
      legend.background = element_rect(color = 'transparent', fill = alpha('white', 0.1)),
      legend.key = element_rect(fill = 'transparent', color = 'transparent')
    ) +
    annotate(
      'text',
      label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
      x = Inf,
      y = Inf,
      hjust = 1.2,
      vjust = 1.35,
      size = 11.75
    )
  
  # save plot
  ggsave(
    sprintf('C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/%s-%s.png', 
            name1, name2),
    plot = p,
    width = 11,
    height = 11,
    dpi = 200
  )
  
  # return the summary of the analysis results
  return(list(
    pair = paste(name1, "-", name2),
    M2 = pro_test$ss,
    p_value = pro_test$signif,
    residuals = residuals(procrustes(pcoa1, pcoa2, symmetric = TRUE))
  ))
}

# conduct batch pairwise analysis
results <- list()
dataset_names <- names(datasets_list)

# pair in the required order
for (i in 1:(length(dataset_names) - 1)) {
  name1 <- dataset_names[i]
  for (j in (i + 1):length(dataset_names)) {
    name2 <- dataset_names[j]
    result <- perform_procrustes(name1, name2, datasets_list[[name1]], datasets_list[[name2]])
    results[[length(results) + 1]] <- result
  }
}

# output results
for (i in 1:length(results)) {
  result <- results[[i]]
  cat(sprintf("%s: M2 = %.3f, p = %.3f\n", 
              result$pair, result$M2, result$p_value))
}
