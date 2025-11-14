# Editor: Li Ruihong (Department of Environmental Science, Chonqing University)
# this file includes: (1) pcoa and ecotrajectories (2) procrustes analysis
#                     (3) model compensating effects 

# Load packages
packages <- c(
  'vegan', 'ggplot2', 'openxlsx', 'pairwiseAdonis', 'ecotraj', 
  'dplyr', 'stringr', 'smacof', 'labdsv', 'tidyr'
)
lapply(packages, library, character.only = TRUE)

### Part.1 pcoa and ecotrajectories

## Algae
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/algae.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Algae.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

## Bacteria
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/bacteria.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Bacteria.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

## Fungi
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/fungi.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Fungi.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

## Insect
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/insect.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Insect.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

## Zooplankton
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/zooplankton.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Zooplankton.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

## Fish
# Read data
rm(list = ls())
data <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/fish.xlsx')
groups <- read.xlsx('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/datasets_by taxonomic groups/# pcoa_groups.xlsx')
rownames(data) <- data[, 1]

# Standardize datasets
data <- hellinger(as.data.frame(t(data[, -1])))

# Calculate bray distance
d <- as.matrix(vegdist(data, method = 'bray'))

# Extract entities
entities <- sub('_.*', '', groups$sample)
surveys <- year_vec <- as.numeric(rep(c('1', '2', '3', '4'), times = 18))  # Extract surveys

# Analyze trajectory and extract pcoa coordinates
trajectory <- defineTrajectories(d, entities, surveys)
trajectory_pcoa <- trajectoryPCoA(trajectory)
trajectory_pcoa_point <- as.data.frame(trajectory_pcoa$points)
coords <- as.data.frame(trajectory_pcoa_point)[, 1:2]
coords$Entity <- entities
coords$Survey <- surveys
coords$Entity <- as.factor(coords$Entity)
coords$Survey <- as.factor(coords$Survey)
colnames(coords)[1] <- 'dim1'
colnames(coords)[2] <- 'dim2'
pc1 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[1]
pc2 <- round((trajectory_pcoa$eig / sum(trajectory_pcoa$eig)) * 100, 2)[2]

# Plot
colnames(coords) <- c('dim1', 'dim2', 'site', 'group')
mycol <- c('#647ADD', '#C0C0C0', '#FFDF67', '#F5664D')
group_labels <- c('2022', '2023', '2024', '2025')

hulls <- coords %>%
  group_by(group) %>%
  slice(chull(dim1, dim2))

centroids <- coords %>%
  group_by(group) %>%
  summarise(
    dim1 = mean(dim1), 
    dim2 = mean(dim2)
  ) %>%
  arrange(factor(group, levels = group_labels)) %>%
  mutate(label = 1:n())

traj_segments <- data.frame(
  x = centroids$dim1[-nrow(centroids)],  # Start x
  y = centroids$dim2[-nrow(centroids)],  # Start y
  xend = centroids$dim1[-1],             # End x
  yend = centroids$dim2[-1],             # End y
  group = group_labels[-length(group_labels)]
)

p <- ggplot(coords, aes(x = dim1, y = dim2, color = as.factor(group))) +
  geom_polygon(
    data = hulls, 
    aes(fill = as.factor(group), color = as.factor(group)), 
    alpha = 0.58, linewidth = 0.4
  ) +
  geom_point(size = 0, shape = 21, stroke = 2.5, fill = 'white', alpha = 0.75) +
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
    color = 'black', fill = mycol, shape = 21, size = 7.5
  ) +
  geom_text(
    data = centroids,
    aes(x = dim1, y = dim2, label = label),
    color = 'black',
    size = 13,
    hjust = -0.2, vjust = -0.2
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

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/pcoa_Fish.png', 
  plot = p, width = 8, height = 8, dpi = 600
)

### Part.2 procrustes analysis 

## Fish - Algae
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged_datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
algae_hel <- decostand(transpose(algae), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_algae <- cmdscale(vegdist(algae_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_algae, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_algae, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Algae')) 
legend_data <- data.frame(
  group = c('Fish', 'Algae'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Algae'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Algae' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Algae')
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
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Algae.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

## Fish - Bacteria
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
bacteria_hel <- decostand(transpose(bacteria), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_bacteria <- cmdscale(vegdist(bacteria_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_bacteria, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_bacteria, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Bacteria')) 
legend_data <- data.frame(
  group = c('Fish', 'Bacteria'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Bacteria'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Bacteria' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Bacteria')
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
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Bacteria.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

## Fish - Fungi
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
fungi_hel <- decostand(transpose(fungi), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_fungi <- cmdscale(vegdist(fungi_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_fungi, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_fungi, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Fungi')) 
legend_data <- data.frame(
  group = c('Fish', 'Fungi'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Fungi'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Fungi' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Fungi')
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
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Fungi.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

## Fish - Benthic invertebrate
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
invertebrate_hel <- decostand(transpose(invertebrate), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_invertebrate <- cmdscale(vegdist(invertebrate_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_invertebrate, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_invertebrate, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Benthic invertebrate')) 
legend_data <- data.frame(
  group = c('Fish', 'Benthic invertebrate'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Benthic invertebrate'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Benthic invertebrate' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Benthic invertebrate')
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
    legend.position = c(0.21, 0.9),  
    legend.key.size = unit(1.5, 'cm'),
    legend.title = element_text(size = 23.5, face = 'bold'),
    legend.text = element_text(size = 23.5),
    legend.background = element_rect(color = 'transparent', fill = alpha('white', 0.1)),  
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Benthic invertebrate.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

## Fish - Insect
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
insect_hel <- decostand(transpose(insect), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_insect <- cmdscale(vegdist(insect_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_insect, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_insect, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Insect')) 
legend_data <- data.frame(
  group = c('Fish', 'Insect'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Insect'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Insect' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Insect')
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
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Insect.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

## Fish - Zooplankton
# Read data
rm(list = ls())
load('C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/merged datasets.Rdata')

# Standardize datasets and extract pcoa coordinates
transpose <- function(y) {
  current_df <- as.data.frame(y)
  species    <- current_df[, 1]
  mat        <- t(current_df[, -1])
  colnames(mat) <- species
  return(as.data.frame(mat))
}

fish_hel <- decostand(transpose(fish), method = 'hellinger')
zooplankton_hel <- decostand(transpose(zooplankton), method = 'hellinger')
pcoa_fish <- cmdscale(vegdist(fish_hel, method = 'bray'))
pcoa_insect <- cmdscale(vegdist(zooplankton_hel, method = 'bray'))

# Conduct procrustes analysis
pro <- procrustes(pcoa_fish, pcoa_insect, symmetric = TRUE)
plot(pro, kind = 2)
residuals(pro) 

set.seed(123)
pro_test <- protest(pcoa_fish, pcoa_insect, permutations = 999) 
pro_test$ss # M2 statistic value
pro_test$signif # M2 significance

# Extract protes coordinates
Pro_coords <- cbind(data.frame(pro_test$Yrot), data.frame(pro_test$X))
Pro_rotation <- data.frame(pro_test$rotation)
colnames(Pro_coords) <- c('X1', 'X2', 'Dim1', 'Dim2')
rm(list = setdiff(ls(), c('Pro_coords', 'Pro_rotation', 'pro_test')))

# Plot
plotdf_group <- data.frame(group = c('Fish', 'Zooplankton')) 
legend_data <- data.frame(
  group = c('Fish', 'Zooplankton'),
  color = c('#E41A1C', '#377EB8') 
)

p <- ggplot(data = Pro_coords) +
  geom_point(aes(X1, X2, color = 'Fish'), size = 7, shape = 16) + 
  geom_point(aes(Dim1, Dim2, color = 'Zooplankton'), size = 7, shape = 16) + 
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
    values = c('Fish' = '#E41A1C', 'Zooplankton' = '#377EB8'),
    name = 'Taxonomic Group',  
    limits = c('Fish', 'Zooplankton')
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
    legend.key = element_rect(fill = 'transparent', color = 'transparent'),  
  ) +
  annotate(
    'text', 
    label = sprintf('M2 = %.3f\np = %.3f', pro_test$ss, pro_test$signif),
    x = Inf, y = Inf, hjust = 1.2, vjust = 1.35,
    size = 11.75
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Fish-Zooplankton.png', 
  plot = p, width = 11, height = 11, dpi = 200
)

### Part.3 model compensating effects

## Mean Path Length
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_mean_path_length.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_mean_path_length.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(1.162488958, 1.233940422, 1.468541456, 1.294841298)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.005", 
  "2023" = "2023  CV=0.007",
  "2024" = "2024  CV=0.019",  
  "2025" = "2025  CV=0.007"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(1.14, 1.5)) +
  geom_hline(
    yintercept = 1.162488958, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 1.233940422, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 1.468541456, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 1.294841298, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Mean Path Length.png', 
  plot = p, width = 9, height = 6, dpi = 300
)

## Connectance
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_connectance.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_connectance.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(0.1373143333, 0.1222612264, 0.06741378756, 0.09476242467)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.019", 
  "2023" = "2023  CV=0.021",
  "2024" = "2024  CV=0.170",  
  "2025" = "2025  CV=0.028"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(0.05, 0.16)) +
  geom_hline(
    yintercept = 0.1373143333, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.1222612264, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.06741378756, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.09476242467, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Connectance.png', 
  plot = p, width = 9, height = 6, dpi = 300
)

## Modularity
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_modularity.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_modularity.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(0.05526768333, 0.08118164028, 0.124806275, 0.05914249222)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.038", 
  "2023" = "2023  CV=0.048",
  "2024" = "2024  CV=0.071",  
  "2025" = "2025  CV=0.057"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(0.04, 0.145)) +
  geom_hline(
    yintercept = 0.05526768333, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.08118164028, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.124806275, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.05914249222, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Modularity.png', 
  plot = p, width = 9, height = 6, dpi = 300
)

## Nestedness
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_nestedness.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_nestedness.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(85.61900944, 84.87469389,82.36244733)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.009", 
  "2023" = "2023  CV=0.008",
  "2024" = "2024  CV=0.007",  
  "2025" = "2025  CV=0.007"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(80, 88)) +
  geom_hline(
    yintercept = 85.62, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 84.87469389, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 85.62822939, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 82.36244733, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Nestedness.png', 
  plot = p, width = 9, height = 6, dpi = 300
)

## Vulnerability
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_vulnerability.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_vulnerability.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(213.3328018, 119.9470125, 279.8278948)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.043", 
  "2023" = "2023  CV=0.054",
  "2024" = "2024  CV=0.043",  
  "2025" = "2025  CV=0.048"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(105, 300)) +
  geom_hline(
    yintercept = 213.3328018, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 214.3640539, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 119.9470125, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 279.8278948, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Vulnerability.png', 
  plot = p, width = 9, height = 6, dpi = 300
)

## Robustness
# Read data
rm(list = ls())
data_ave <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_robustness.xlsx', 
  1
)
data_range <- read.xlsx(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/datasets/module compensating effects/plotting files/plot_robustness.xlsx', 
  2
)
module_levels <- unique(data_ave$module)
data_ave$module_id <- as.numeric(factor(data_ave$module, levels = module_levels))
data_range$module_id <- as.numeric(factor(data_range$module, levels = module_levels))

# Define baseline values 
baseline_values <- c(0.4424948867, 0.4356060338,0.2837738379, 0.3604388324)

# Define legend texts 
legend_labels <- c(
  "2022" = "2022  CV=0.009", 
  "2023" = "2023  CV=0.015",
  "2024" = "2024  CV=0.059",  
  "2025" = "2025  CV=0.016"
)

# Plot
p <- ggplot() +
  geom_ribbon(
    data = data_range,
    aes(x = module_id, ymin = MIN, ymax = MAX, fill = factor(year)),
    alpha = 0.12 
  ) +
  geom_line(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 1
  ) +
  geom_point(
    data = data_ave,
    aes(x = module_id, y = value, color = factor(year)),
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:length(module_levels), 
    labels = module_levels
  ) +
  scale_fill_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_color_manual(
    values = c(
      "2022" = "#647ADD", 
      "2023" = "#C0C0C0", 
      "2024" = "#FFDF69", 
      "2025" = "#F5664D"
    ),
    labels = legend_labels
  ) +
  scale_y_continuous(
    breaks = baseline_values, 
    labels = function(x) {  
      sprintf("%.3f", x)
    }
  ) +
  coord_cartesian(ylim = c(0.2, 0.455)) +
  geom_hline(
    yintercept = 0.4424948867, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.4356060338, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.2837738379, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  geom_hline(
    yintercept = 0.3604388324, color = "black", linetype = "dashed", 
    size = 0.7, alpha = 0.35
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 15),                    
    axis.title = element_blank(),                                 
    plot.title = element_text(size = 34, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid = element_blank(),                       
    strip.text = element_text(size = 12),         
    scale_x_discrete(expand = expansion(mult = c(0, 0))),
    legend.position = c(0.08, 0.05), 
    legend.justification = c(0.35, 0.15),  
    legend.background = element_rect(fill = NA, color = NA),  
    legend.key = element_blank(),     
    legend.title = element_text(size = 13), 
    legend.text = element_text(size = 13)   
  )

# Save plots
ggsave(
  'C:/Users/23926/Desktop/works/#1 datasets and codes/codes/figures/Robustness.png', 
  plot = p, width = 9, height = 6, dpi = 300

)
