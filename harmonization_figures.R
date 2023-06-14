
library(feather)
library(tidyverse)
library(Rtsne)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(rstatix)
library(dplyr)
library(irlba)
library(sjPlot)
library(mmrm)
library(margins)
library(emmeans)
library(ggeffects)
library(ggpmisc)
library(mosaic)
library(ggpmisc)
library(GGally)

#### FIGURE HARMONIZATION ####

base <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/neurocombat_wide.csv', sep = ",", header = TRUE)
no_combat <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/unharm_wide.csv', sep = ",", header = TRUE)
neuroharmonize <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/neuroharm_wide.csv', sep = ",", header = TRUE)
opcombat <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/opcombat_wide.csv', sep = ",", header = TRUE)


neuroharmonize <- neuroharmonize %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- neuroharmonize %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(neuroharmonize, by="ID")


neuroharmonized = ggscatter(fit_spectral_df, x = "tSNE1",
                       y = "tSNE2",
                       color = "center", ellipse = T, 
                       mean.point = T, 
                       star.plot = F, alpha=0.4, ggtheme = theme_pubr(border = T, base_family = "mono"), title="neuroHarmonize", 
                       subtitle = "",   ellipse.alpha = 0.04, size =1.5, shape=20, mean.point.size=4)
neuroharmonized = ggpar(p = neuroharmonized, xlim=c(-30,30),ylim=c(-30,30))  + geom_hline(yintercept=0, linetype="dashed", color = "gray")  +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  scale_fill_discrete(name="",labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  scale_color_discrete(name="",
                       labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))

neuroharmonized


# CALCULATE CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)


# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_neuroharmonize <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_neuroharmonize


combat <- base
combat <- combat %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 30
set.seed(42)
fit_spectral <- combat %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(combat, by="ID")


normalized = ggscatter(fit_spectral_df, x = "tSNE1",
                       y = "tSNE2",
                       color = "center", ellipse = T, 
                       mean.point = T, 
                       star.plot = F, alpha=0.4, ggtheme = theme_pubr(border = T, base_family = "mono"), title="neuroCombat", 
                       subtitle = "",   ellipse.alpha = 0.04, size =1.5, shape=20, mean.point.size=4 )
normalized = ggpar(p = normalized, xlim=c(-30,30),ylim=c(-30,30))  + geom_hline(yintercept=0, linetype="dashed", color = "gray")  +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  scale_fill_discrete(name="",labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  scale_color_discrete(name="",
                       labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))

normalized


# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))
# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)


# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_combat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_combat











base_wide <- no_combat
base_wide <- base_wide %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- base_wide %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(base_wide, by="ID")


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(base_wide, by="ID")


no_normalized = ggscatter(fit_spectral_df, x = "tSNE1",
                          y = "tSNE2",
                          color = "center", ellipse = T, 
                          mean.point = T,
                          star.plot = F, alpha = 0.4, ggtheme = theme_pubr(border = T, base_family = "mono"), title = "Not Harmonized", 
                          subtitle = "",   ellipse.alpha = 0.04, size =1.5, shape=20, mean.point.size=4)
no_normalized = ggpar(p = no_normalized, xlim=c(-30,30),ylim=c(-30,30)) + geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray")+
  scale_fill_discrete(name="",labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  scale_color_discrete(name="",
                       labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))

no_normalized




# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]



# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)



# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_no_combat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
theme(
  plot.title = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  axis.text = element_text(face = "bold", size = 12))

heatmap_no_combat







opcombat <- opcombat %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- opcombat %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(opcombat, by="ID")


opcombated = ggscatter(fit_spectral_df, x = "tSNE1",
                            y = "tSNE2",
                            color = "center", ellipse = T, 
                            mean.point = T, 
                            star.plot = F, alpha=0.4, ggtheme = theme_pubr(border = T, base_family = "mono"), title="OPNCombat-GMM", 
                            subtitle = "",   ellipse.alpha = 0.04, size =1.5, shape=20, mean.point.size=4)
opcombated = ggpar(p = opcombated, xlim=c(-30,30),ylim=c(-30,30))  + geom_hline(yintercept=0, linetype="dashed", color = "gray")  +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  scale_fill_discrete(name="",labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  scale_color_discrete(name="",
                       labels=c("California","Finland", "Iowa","LEMON","Oslo")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))

opcombated


# GET CENTROIDS

centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)


# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_opcombat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_opcombat


figure_normalization   <- ggarrange(no_normalized, normalized,neuroharmonized,opcombated,
                                    ncol = 4, nrow = 1, common.legend = T,  font.label = list(face = "bold", family = "mono"))

figure_normalization

figure_centroids   <- ggarrange(heatmap_no_combat, heatmap_combat, heatmap_neuroharmonize,heatmap_opcombat,
                                ncol = 4, nrow = 1, font.label = list(face = "bold", family = "mono"))

figure_centroids   <- annotate_figure(figure_centroids, top = text_grob("", 
                                      color = "black", family = 'mono',face = "bold", size = 16))




final_fig <- ggarrange(figure_normalization, figure_centroids, labels = c("A", "B"),
                       ncol = 1, nrow = 2, font.label = list(face = "bold", family = "mono"))


tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/tsne_dist.tiff", units="in", width=16, height=8, res=500)
final_fig
dev.off()












#### FIGURE NO ELLIPSE ####
base <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/neurocombat_wide.csv', sep = ",", header = TRUE)
no_combat <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/unharm_wide.csv', sep = ",", header = TRUE)
neuroharmonize <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/neuroharm_wide.csv', sep = ",", header = TRUE)
opcombat <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/opcombat_wide.csv', sep = ",", header = TRUE)
harmonizer_nan <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/tsne/harmonizer_wide.csv', sep = ",", header = TRUE)
harmonizer_nan[is.na(harmonizer_nan)] <- 0

neuroharmonize <- neuroharmonize %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- neuroharmonize %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(neuroharmonize, by="ID")



# Define the colors from the Set1 seaborn palette
set1_palette <- c("#4C72B0", "#8172B3","#55A868", "#DD8452",  "#C44E52")


# Scatter plot
neuroharmonized <- ggscatter(fit_spectral_df, x = "tSNE1", y = "tSNE2", 
                             color = "center", ellipse = FALSE, 
                             mean.point = TRUE, star.plot = FALSE, alpha = 0.4,
                             ggtheme = theme_pubr(border = TRUE, base_family = "mono"), 
                             title = "neuroHarmonize", subtitle = "", ellipse.alpha = 0.04,
                             size = 1.5, shape = 20, mean.point.size = 4)

# Set the color and fill scales to the Set1 seaborn palette
neuroharmonized <- neuroharmonized +
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

# Set plot limits and add dashed lines
neuroharmonized <- neuroharmonized +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

# Customize the plot titles and axis labels
neuroharmonized <- neuroharmonized +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

neuroharmonized









# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)



# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_neuroharmonize <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_neuroharmonize


combat <- base
combat <- combat %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- combat %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(combat, by="ID")



# Scatter plot
normalized <- ggscatter(fit_spectral_df, x = "tSNE1", y = "tSNE2", 
                             color = "center", ellipse = FALSE, 
                             mean.point = TRUE, star.plot = FALSE, alpha = 0.4,
                             ggtheme = theme_pubr(border = TRUE, base_family = "mono"), 
                             title = "neuroCombat", subtitle = "", ellipse.alpha = 0.04,
                             size = 1.5, shape = 20, mean.point.size = 4)

# Set the color and fill scales to the Set1 seaborn palette
normalized <- normalized +
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

# Set plot limits and add dashed lines
normalized <- normalized +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

# Customize the plot titles and axis labels
normalized <- normalized +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

normalized






# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)

# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_combat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_combat











base_wide <- no_combat
base_wide <- base_wide %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- base_wide %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(base_wide, by="ID")


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(base_wide, by="ID")



# Scatter plot
no_normalized <- ggscatter(fit_spectral_df, x = "tSNE1", y = "tSNE2", 
                        color = "center", ellipse = FALSE, 
                        mean.point = TRUE, star.plot = FALSE, alpha = 0.4,
                        ggtheme = theme_pubr(border = TRUE, base_family = "mono"), 
                        title = "Unharm. Zero", subtitle = "", ellipse.alpha = 0.04,
                        size = 1.5, shape = 20, mean.point.size = 4)

# Set the color and fill scales to the Set1 seaborn palette
no_normalized <- no_normalized +
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

# Set plot limits and add dashed lines
no_normalized <- no_normalized +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

# Customize the plot titles and axis labels
no_normalized <- no_normalized +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

no_normalized








# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]



# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)



# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_no_combat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_no_combat








opcombat <- opcombat %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- opcombat %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(opcombat, by="ID")





# Scatter plot
opcombated <- ggscatter(fit_spectral_df, x = "tSNE1", y = "tSNE2", 
                           color = "center", ellipse = FALSE, 
                           mean.point = TRUE, star.plot = FALSE, alpha = 0.4,
                           ggtheme = theme_pubr(border = TRUE, base_family = "mono"), 
                           title = "OPNComBat-GMM", subtitle = "", ellipse.alpha = 0.04,
                           size = 1.5, shape = 20, mean.point.size = 4)

# Set the color and fill scales to the Set1 seaborn palette
opcombated <- opcombated +
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

# Set plot limits and add dashed lines
opcombated <- opcombated +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

# Customize the plot titles and axis labels
opcombated <- opcombated +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

opcombated







# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)



# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)


heatmap_opcombat <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
    labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_opcombat



harmonizer_nan <- harmonizer_nan %>%
  mutate(ID=row_number())

#### FIT TSNE MODEL FOR SPECTRAL METRICS perplexity 20
set.seed(42)
fit_spectral <- harmonizer_nan %>%
  dplyr::select(!c(participant_id, age, gender, center)) %>%
  column_to_rownames(var = "ID") %>%
  scale() %>%
  Rtsne(perplexity = 30, check_duplicates = TRUE)


fit_spectral_df <- fit_spectral$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(ID=row_number())

fit_spectral_df <- fit_spectral_df %>%
  inner_join(harmonizer_nan, by="ID")





# Scatter plot
harmonizered <- ggscatter(fit_spectral_df, x = "tSNE1", y = "tSNE2", 
                        color = "center", ellipse = FALSE, 
                        mean.point = TRUE, star.plot = FALSE, alpha = 0.4,
                        ggtheme = theme_pubr(border = TRUE, base_family = "mono"), 
                        title = "HarmonizR", subtitle = "", ellipse.alpha = 0.04,
                        size = 1.5, shape = 20, mean.point.size = 4)

# Set the color and fill scales to the Set1 seaborn palette
harmonizered <- harmonizered +
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

# Set plot limits and add dashed lines
harmonizered <- harmonizered +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-30, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

# Customize the plot titles and axis labels
harmonizered <- harmonizered +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

harmonizered







# GET CENTROIDS
centroids <- fit_spectral_df %>%
  group_by(center) %>%
  summarise(tSNE1_centroid = mean(tSNE1), tSNE2_centroid = mean(tSNE2))

# CALCULATE THE EUCLIDEAN DISTANCE BETWEEN CENTROIDS
distances <- as.matrix(dist(centroids[,2:3]))

# GET DISTANCE DIFF BETWEEN CENTERS FROM THE DIAGONAL OF THE MATRIX
distances_between_centroids <- distances[upper.tri(distances)]


# ADD LABELS OF EACH RESEARCH CENTER TO THE MATRIX
rownames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")
colnames(distances) <- c("California","Finland", "Iowa","LEMON","Oslo")

# RESHAPE MATRIX INTO A DATAFRAME
df <- reshape2::melt(distances)

# locations of lower triangle
lt <- lower.tri(distances, diag = FALSE)
ltdf <- data.frame(row = rownames(distances)[row(distances)[lt]],
                   col = colnames(distances)[col(distances)[lt]],
                   values = distances[lt],
                   stringsAsFactors = FALSE)



heatmap_harmonizer <- ggplot(ltdf, aes(x = row, y = col, fill = values)) + geom_raster() +
  geom_text(aes(label = round(values,2), size = 2), family = 'mono') +
  scale_fill_gradient(limits = c(0, 16), n.breaks = 30, low = "white", high = "#980043")+
  theme_minimal(base_family = "mono") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = "", y = "", title = paste0( "Av. Dist. ", round(mean(ltdf$values), 2)))+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", size = 12))

heatmap_harmonizer





figure_normalization   <- ggarrange(no_normalized, normalized,neuroharmonized,opcombated,harmonizered,
                                    ncol = 5, nrow = 1, common.legend = T,  font.label = list(face = "bold", family = "mono"))

figure_normalization

figure_centroids   <- ggarrange(heatmap_no_combat, heatmap_combat, heatmap_neuroharmonize,heatmap_opcombat,heatmap_harmonizer,
                                ncol = 5, nrow = 1, font.label = list(face = "bold", family = "mono"))

figure_centroids   <- annotate_figure(figure_centroids, top = text_grob("", 
                                                                        color = "black", family = 'mono',face = "bold", size = 16))




final_fig <- ggarrange(figure_normalization, figure_centroids, labels = c("A", "B"),
                       ncol = 1, nrow = 2, font.label = list(face = "bold", family = "mono"))


tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/tsne_noelipse_dist.tiff", units="in", width=19, height=8, res=500)
final_fig
dev.off()






final_fig




#### CORRELATION PLOTS AGE ####


delta_lm = lm(value ~ age + gender + center, data = delta)
plot_model(delta_lm, type = "pred", terms = c("age", "center"),
           show.data = T, show.values = T, show.p = T) +
           labs(title = paste("Adj R2 = ",signif(summary(delta_lm)$adj.r.squared, 2),
                     " p-value =",signif(summary(delta_lm)$coef[2,4], 1)))


ggplotRegression <- function(fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 =", signif(summary(fit)$adj.r.squared, 2),
                       "p-value =", signif(summary(fit)$coef[2, 4], 1)))
}

ggplotRegression(lm(value ~ age + gender + center, data = delta))

fit$model = lm(value ~ age + gender + center, data = delta)



#### BANDPOWERS ####

unh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/unharmonized_hc.csv', sep = ",", header = TRUE)
combat_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neurocombat.csv', sep = ",", header = TRUE)
neuroh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neuroharmonize.csv', sep = ",", header = TRUE)
opcomb_roi <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/opcombat.csv', sep = ",", header = TRUE)
harmonizer_roi <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/harmonizer_long.csv', sep = ",", header = TRUE)


covars <- unh_roi %>%
  dplyr::select(participant_id, age, gender, center)%>%
  distinct(participant_id, .keep_all = TRUE)

# Crear vector con los nombres de las ROIs
ROI_names <- c("anterior", "central", "posterior")

# Crear una lista con los canales de cada ROI
ROI_channels <- list(anterior = c("Fp1", "Fp2", "AF3", "AF4", "F3", "F4", "F7", "F8", "Fz"),
                     central  = c("FC1", "FC2", "FC5", "FC6", "T7", "T8", "C3", "C4", "Cz", "CP5", "CP6", "CP1", "CP2"),
                     posterior = c("P3", "P4", "P7", "P8", "Pz", "PO3", "PO4", "Oz", "O1", "O2"))

# Crear una función que asigne cada canal a su respectiva ROI
get_ROI <- function(channel) {
  for (ROI in ROI_names) {
    if (channel %in% ROI_channels[[ROI]]) {
      return(ROI)
    }
  }
}

# Crear variable ROI en el data frame original
unh_roi$ROI <- sapply(unh_roi$channel, get_ROI)
combat_roi$ROI <- sapply(combat_roi$channel, get_ROI)
neuroh_roi$ROI <- sapply(neuroh_roi$channel, get_ROI)
opcomb_roi$ROI <- sapply(opcomb_roi$channel, get_ROI)
harmonizer_roi$ROI <- sapply(harmonizer_roi$channel, get_ROI)


# Load r_squared dataset and filter data with r_squared < 0.8
rsq_data <- read.table('D:/data_analysis/papers_alberto/combat_eeg/unharmonized_error.csv', sep = ",", header = TRUE)

# Seleccionar las variables deseadas del nuevo dataset
rsq_data <- rsq_data[, c("participant_id", "channel", "r_squared")]

# Unir los datasets basado en las variables "participant_id" y "channel"
unh_roi <- merge(unh_roi, rsq_data, by = c("participant_id", "channel"))
combat_roi <- merge(combat_roi, rsq_data, by = c("participant_id", "channel"))
neuroh_roi <- merge(neuroh_roi, rsq_data, by = c("participant_id", "channel"))
opcomb_roi <- merge(opcomb_roi, rsq_data, by = c("participant_id", "channel"))
harmonizer_roi <- merge(harmonizer_roi, rsq_data, by = c("participant_id", "channel"))


##### ####






###### NAN ####
unh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/unh_nan_long.csv', sep = ",", header = TRUE)
harmonizer_roi <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/harmonizer_long.csv', sep = ",", header = TRUE)


covars <- unh_roi %>%
  dplyr::select(participant_id, age, gender, center)%>%
  distinct(participant_id, .keep_all = TRUE)

# Crear vector con los nombres de las ROIs
ROI_names <- c("anterior", "central", "posterior")

# Crear una lista con los canales de cada ROI
ROI_channels <- list(anterior = c("Fp1", "Fp2", "AF3", "AF4", "F3", "F4", "F7", "F8", "Fz"),
                     central  = c("FC1", "FC2", "FC5", "FC6", "T7", "T8", "C3", "C4", "Cz", "CP5", "CP6", "CP1", "CP2"),
                     posterior = c("P3", "P4", "P7", "P8", "Pz", "PO3", "PO4", "Oz", "O1", "O2"))

# Crear una función que asigne cada canal a su respectiva ROI
get_ROI <- function(channel) {
  for (ROI in ROI_names) {
    if (channel %in% ROI_channels[[ROI]]) {
      return(ROI)
    }
  }
}

# Crear variable ROI en el data frame original
unh_roi$ROI <- sapply(unh_roi$channel, get_ROI)
harmonizer_roi$ROI <- sapply(harmonizer_roi$channel, get_ROI)


# Load r_squared dataset and filter data with r_squared < 0.8
rsq_data <- read.table('D:/data_analysis/papers_alberto/combat_eeg/unharmonized_error.csv', sep = ",", header = TRUE)

# Seleccionar las variables deseadas del nuevo dataset
rsq_data <- rsq_data[, c("participant_id", "channel", "r_squared")]

# Unir los datasets basado en las variables "participant_id" y "channel"
unh_roi <- merge(unh_roi, rsq_data, by = c("participant_id", "channel"))
harmonizer_roi <- merge(harmonizer_roi, rsq_data, by = c("participant_id", "channel"))

# #### median by roi ####
#


unh_roi <-  unh_roi  %>%  filter(r_squared >=0.8)
harmonizer_roi <-  harmonizer_roi  %>%  filter(r_squared >=0.8)



unh_roi <- unh_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))
harmonizer_roi <- harmonizer_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))


unh_roi <- unh_roi %>%
  inner_join(covars, by="participant_id")

harmonizer_roi <- harmonizer_roi %>%
  inner_join(covars, by="participant_id")



##### ####



delta <-  harmonizer_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")

theta <-  harmonizer_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")

alpha <-  harmonizer_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")

beta <-  harmonizer_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")

harmonizer_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 4,
                          common.legend = T, legend = "none")




#### Correlations bands #####
# Define the colors from the Set1 seaborn palette
set1_palette <- c("#4C72B0", "#8172B3","#55A868", "#DD8452",  "#C44E52")


delta <-  unh_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")


delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))


theta <-  unh_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  unh_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  unh_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

unh_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 1,
                   common.legend = T, legend = "none")


delta <-  harmonizer_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  harmonizer_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  harmonizer_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  harmonizer_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

harmonizer_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 1,
                     common.legend = T, legend = "none")



figure_power = ggarrange(unh_pw, harmonizer_pw, ncol = 1,
                         common.legend = T, labels = c('A', 'B'),
                         legend = "none")


#### #####


tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/medianpost_powers_age_nan.tiff", 
    units="in", width=15, height=10, res=500)
figure_power
dev.off()

figure_power

#### Correlations #####
interc <-  unh_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  unh_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  unh_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  unh_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")

dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

unh_p = ggarrange(dominant_p, interc_p, slope_p, nrow = 1,
                  common.legend = T, legend = "none")

interc <-  harmonizer_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "offset (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  harmonizer_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  harmonizer_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  harmonizer_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")
dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

harmonizer_p = ggarrange(dominant_p, interc_p, slope_p, ncol = 3,
                    common.legend = T, legend = "none")


figure_aperiodic = ggarrange(unh_p, harmonizer_p, ncol = 1,
                             common.legend = T, labels = c('A', 'B'),
                             legend = "none")

#### #####

tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/medianpost_fooof_age_nan.tiff", 
    units="in", width=15, height=10, res=500)
figure_aperiodic
dev.off()

figure_aperiodic









#### NOT RUN ####
ggplotRegression <- function(fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(x= 'Age (years)', y = 'exponent (a.u)',
         title = paste("Adj R2 =", signif(summary(fit)$adj.r.squared, 2),
                       "p-value =", signif(summary(fit)$coef[2, 4], 1)))
}
slope <-  unh_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")
a = ggplotRegression(lm(value ~ age + gender, data = slope))

slope <-  combat_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")
b = ggplotRegression(lm(value ~ age + gender, data = slope))

slope <-  neuroh_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")
c = ggplotRegression(lm(value ~ age + gender, data = slope))

slope <-  opcomb_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")
d = ggplotRegression(lm(value ~ age + gender, data = slope))

ggarrange(a,b,c,d, nrow=2, ncol=2)







ggplotRegression <- function(fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(x= 'Age (years)', y = 'offset (a.u)',
         title = paste("Adj R2 =", signif(summary(fit)$adj.r.squared, 2),
                       "p-value =", signif(summary(fit)$coef[2, 4], 1)))
}
interc <-  unh_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")
a = ggplotRegression(lm(value ~ age + gender, data = interc))

interc <-  combat_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")
b = ggplotRegression(lm(value ~ age + gender, data = interc))

interc <-  neuroh_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")
c = ggplotRegression(lm(value ~ age + gender, data = interc))

interc <-  opcomb_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")
d = ggplotRegression(lm(value ~ age + gender, data = interc))

ggarrange(a,b,c,d, nrow=2, ncol=2)


































ggscatter(delta, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, delta = 0.3, color = 'center',ylab = "Not-harmonized delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(delta, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, delta = 0.3,ylab = "Not-harmonized delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))

theta <-  raw  %>%  filter(feature == "theta_pw")

ggscatter(theta, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, theta = 0.3, color = 'center',ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(theta, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, theta = 0.3,ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))



alpha <-  raw  %>%  filter(feature == "alpha_pw")

ggscatter(alpha, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, alpha = 0.3, color = 'center',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(alpha, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, alpha = 0.3,ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))


beta <-  raw  %>%  filter(feature == "beta_pw")

ggscatter(beta, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, beta = 0.3, color = 'center',ylab = "Not-harmonized beta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(beta, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, beta = 0.3,ylab = "Not-harmonized beta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))





exponent <-  raw  %>%  filter(feature == "exponent")

ggscatter(exponent, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, exponent = 0.3, color = 'center',ylab = "Not-harmonized exponent Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(exponent, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, exponent = 0.3,ylab = "Not-harmonized Aperiodic exponent (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))




offset <-  raw  %>%  filter(feature == "offset")

ggscatter(offset, x = "age", y = "value",
          add = "reg.line",                         # Add regression line
          conf.int = T,                          # Add confidence interval
          size=0.5, offset = 0.3, color = 'center',ylab = "Not-harmonized offset Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+          # Color by groups "cyl")+
  stat_cor(method = "spearman", cor.coef.name = "rho", label.y = c(-0.6),  label.x = 0, size = 5) + # Add correlation coefficient 
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"))



ggscatter(offset, x = "age", y = "value",
          add = "loess",                         # Add regression line
          conf.int = T,# Add confidence interval
          size=0.5, offset = 0.3,ylab = "Not-harmonized Aperiodic offset (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))

















#### BANDPOWERS ####

unh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/unharmonized_hc.csv', sep = ",", header = TRUE)
combat_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neurocombat.csv', sep = ",", header = TRUE)
neuroh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neuroharmonize.csv', sep = ",", header = TRUE)
opcomb_roi <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/opcombat.csv', sep = ",", header = TRUE)


covars <- unh_roi %>%
  dplyr::select(participant_id, age, gender, center)%>%
  distinct(participant_id, .keep_all = TRUE)

# Crear vector con los nombres de las ROIs
ROI_names <- c("anterior", "central", "posterior")

# Crear una lista con los canales de cada ROI
ROI_channels <- list(anterior = c("Fp1", "Fp2", "AF3", "AF4", "F3", "F4", "F7", "F8", "Fz"),
                     central  = c("FC1", "FC2", "FC5", "FC6", "T7", "T8", "C3", "C4", "Cz", "CP5", "CP6", "CP1", "CP2"),
                     posterior = c("P3", "P4", "P7", "P8", "Pz", "PO3", "PO4", "Oz", "O1", "O2"))

# Crear una función que asigne cada canal a su respectiva ROI
get_ROI <- function(channel) {
  for (ROI in ROI_names) {
    if (channel %in% ROI_channels[[ROI]]) {
      return(ROI)
    }
  }
}

# Crear variable ROI en el data frame original
unh_roi$ROI <- sapply(unh_roi$channel, get_ROI)
combat_roi$ROI <- sapply(combat_roi$channel, get_ROI)
neuroh_roi$ROI <- sapply(neuroh_roi$channel, get_ROI)
opcomb_roi$ROI <- sapply(opcomb_roi$channel, get_ROI)


# Load r_squared dataset and filter data with r_squared < 0.8
rsq_data <- read.table('D:/data_analysis/papers_alberto/combat_eeg/unharmonized_error.csv', sep = ",", header = TRUE)

# Seleccionar las variables deseadas del nuevo dataset
rsq_data <- rsq_data[, c("participant_id", "channel", "r_squared")]

# Unir los datasets basado en las variables "participant_id" y "channel"
unh_roi <- merge(unh_roi, rsq_data, by = c("participant_id", "channel"))
combat_roi <- merge(combat_roi, rsq_data, by = c("participant_id", "channel"))
neuroh_roi <- merge(neuroh_roi, rsq_data, by = c("participant_id", "channel"))
opcomb_roi <- merge(opcomb_roi, rsq_data, by = c("participant_id", "channel"))


# #### median by roi ####
#


unh_roi <-  unh_roi  %>%  filter(r_squared >=0.8)
combat_roi <-  combat_roi  %>%  filter(r_squared >=0.8)
neuroh_roi <-  neuroh_roi  %>%  filter(r_squared >=0.8)
opcomb_roi <-  opcomb_roi  %>%  filter(r_squared >=0.8)



unh_roi <- unh_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))
combat_roi <- combat_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))
neuroh_roi <- neuroh_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))
opcomb_roi <- opcomb_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value = mean(value))


unh_roi <- unh_roi %>%
  inner_join(covars, by="participant_id")

combat_roi <- combat_roi %>%
  inner_join(covars, by="participant_id")

neuroh_roi <- neuroh_roi %>%
  inner_join(covars, by="participant_id")

opcomb_roi <- opcomb_roi %>%
  inner_join(covars, by="participant_id")

##### ####
#### Correlations bands #####
delta <-  unh_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")


delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  unh_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  unh_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  unh_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

unh_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 4,
                   common.legend = T, legend = "none")


delta <-  combat_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  combat_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  combat_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  combat_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

ncomb_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 4,
                     common.legend = T, legend = "none")





delta <-  neuroh_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  neuroh_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  neuroh_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  neuroh_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

nharm_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 4,
                     common.legend = T, legend = "none")



delta <-  opcomb_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  opcomb_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  opcomb_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  opcomb_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

opcombat_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, nrow = 4,
                        common.legend = T, legend = "none")




delta <-  harmonizer_roi  %>%  filter(feature == "delta_pw")
delta <-  delta  %>%  filter(ROI == "posterior")

delta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Delta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
delta_p = ggscatter(delta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 1.5),
                    delta = 0.3,size = 1, alpha = 0.25, 
                    ylab = "Delta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Delta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

theta <-  harmonizer_roi  %>%  filter(feature == "theta_pw")
theta <-  theta  %>%  filter(ROI == "posterior")

theta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Theta (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

theta_p = ggscatter(theta, y="value", x="age",
                    color="center",ylim=c(0,2),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 1.5),
                    ylab = "Theta (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Theta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

alpha <-  harmonizer_roi  %>%  filter(feature == "alpha_pw")
alpha <-  alpha  %>%  filter(ROI == "posterior")

alpha %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,3),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Alpha (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
alpha_p =   ggscatter(alpha, y="value", x="age", ylim=c(0,3),
                      color="center",
                      cor.coef = TRUE, cor.method = "pearson",
                      cor.coef.coord = c(20, 2.5),
                      size = 1, delta = 0.3, alpha = 0.25, 
                      ylab = "Alpha (a.u)",
                      xlab = "Age (years)", 
                      ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Alpha (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

beta <-  harmonizer_roi  %>%  filter(feature == "beta_pw")
beta <-  beta  %>%  filter(ROI == "posterior")

beta %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,2),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "Beta (a.u)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
beta_p =   ggscatter(beta, y="value", x="age", ylim=c(0,2),
                     color="center",
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, 1.5),
                     size = 1, delta = 0.3, alpha = 0.25, 
                     ylab = "Beta (a.u)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("Beta (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

harmonizer_pw = ggarrange(delta_p, theta_p, alpha_p, beta_p, ncol = 1,
                          common.legend = T, legend = "none")


figure_power = ggarrange(unh_pw, ncomb_pw, nharm_pw, opcombat_pw,harmonizer_pw,
                         ncol = 5, common.legend = T, 
                         labels = c('A', 'B','C', 'D', 'E'),
                         legend = "none")


#### #####


tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/medianpost_powers_age_all.tiff", 
    units="in", width=15, height=10, res=500)
figure_power
dev.off()

figure_power
#### Correlations #####
interc <-  unh_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  unh_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  unh_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  unh_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")

dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

unh_p = ggarrange(dominant_p, interc_p, slope_p, nrow = 3,
                  common.legend = T, legend = "none")

interc <-  combat_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "offset (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  combat_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  combat_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  combat_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")
dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

ncomb_p = ggarrange(dominant_p, interc_p, slope_p, nrow = 3,
                    common.legend = T, legend = "none")


interc <-  neuroh_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "offset (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  neuroh_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  neuroh_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")


dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  neuroh_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")
dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

nharm_p = ggarrange(dominant_p, interc_p, slope_p, nrow = 3,
                    common.legend = T,  legend = "none")


interc <-  opcomb_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "offset (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  opcomb_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))



dominant <-  opcomb_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dfv <-  opcomb_roi  %>%  filter(feature == "extalphas_bw")
dfv <-  dfv  %>%  filter(ROI == "posterior")
dfv %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(0,8),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "DF bw (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dfv_p =   ggscatter(dfv, y="value", x="age", ylim=c(0,8),
                    color="center",
                    cor.coef = TRUE, cor.method = "pearson",
                    cor.coef.coord = c(20, 8),
                    size = 1, delta = 0.3, alpha = 0.25, 
                    ylab = "DF bw (Hz)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

opcombat_p = ggarrange(dominant_p, interc_p, slope_p, nrow = 3,
                       common.legend = T, legend = "none")




interc <-  harmonizer_roi  %>%  filter(feature == "offset")
interc <-  interc  %>%  filter(ROI == "posterior")

interc %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center", 
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "offset (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
interc_p = ggscatter(interc, y="value", x="age",
                     color="center",ylim=c(-12,-6),
                     cor.coef = TRUE, cor.method = "pearson",
                     cor.coef.coord = c(20, -6.5),
                     delta = 0.3,size = 1, alpha = 0.25, 
                     ylab = "offset (a.u.)",
                     xlab = "Age (years)", 
                     ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("offset (", mu, "V"^2,")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

slope <-  harmonizer_roi  %>%  filter(feature == "exponent")
slope <-  slope  %>%  filter(ROI == "posterior")

slope %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age",
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "exponent (a.u.)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")

slope_p = ggscatter(slope, y="value", x="age",
                    color="center",ylim=c(-1.5,2.5),
                    cor.coef = TRUE, cor.method = "pearson",
                    delta = 0.3, size = 1, alpha = 0.25, 
                    cor.coef.coord = c(20, 2),
                    ylab = "exponent (a.u.)",
                    xlab = "Age (years)", 
                    ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  labs(y = expression(paste("exponent (", mu, "V"^2, " Hz"^-1, ")")),
       x = "Age (years)")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

dominant <-  harmonizer_roi  %>%  filter(feature == "extalphas_cf")
dominant <-  dominant  %>%  filter(ROI == "posterior")

dominant %>%
  mutate(center = factor(center)) %>%
  ggscatter(., y="value", x="age", ylim=c(4,15),
            color="center",
            cor.coef = TRUE, cor.method = "pearson",
            delta = 0.3, alpha = 0.25, 
            ylab = "iAPF (Hz)",
            xlab = "Age (years)", 
            ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")
dominant_p =   ggscatter(dominant, y="value", x="age", ylim=c(4,15),
                         color="center",
                         cor.coef = TRUE, cor.method = "pearson",
                         cor.coef.coord = c(20, 14),
                         size = 1, delta = 0.3, alpha = 0.25, 
                         ylab = "iAPF (Hz)",
                         xlab = "Age (years)", 
                         ggtheme = theme_pubr(border = T, base_family = "mono")) +
  geom_smooth(method = "lm", se = T, color = "black")+
  scale_color_manual(values = set1_palette, name = "", 
                     labels = c("California", "Finland", "Iowa", "LEMON", "Oslo")) +
  scale_fill_manual(values = set1_palette, name = "", 
                    labels = c("California", "Finland", "Iowa", "LEMON", "Oslo"))

harmonizer_p = ggarrange(dominant_p, interc_p, slope_p, ncol = 1,
                         common.legend = T, legend = "none")


figure_aperiodic = ggarrange(unh_p, ncomb_p, nharm_p, opcombat_p, harmonizer_p, 
                             ncol = 5,
                             common.legend = T, labels = c('A', 'B','C', 'D', 'E'),
                             legend = "none")

#### #####
figure_aperiodic

tiff("D:/data_analysis/papers_alberto/combat_eeg/figures/medianpost_fooof_age.tiff", 
    units="in", width=15, height=10, res=500)
figure_aperiodic
dev.off()


































#### LMM ####
unh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/unharmonized_hc.csv', sep = ",", header = TRUE)
combat_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neurocombat.csv', sep = ",", header = TRUE)
neuroh_roi <- read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/neuroharmonize.csv', sep = ",", header = TRUE)
opcomb_roi <-read.table(file = 'D:/data_analysis/papers_alberto/combat_eeg/opcombat.csv', sep = ",", header = TRUE)


covars <- unh_roi %>%
  dplyr::select(participant_id, age, gender, center)%>%
  distinct(participant_id, .keep_all = TRUE)

# Crear vector con los nombres de las ROIs
ROI_names <- c("anterior", "central", "posterior")

# Crear una lista con los canales de cada ROI
ROI_channels <- list(anterior = c("Fp1", "Fp2", "AF3", "AF4", "F3", "F4", "F7", "F8", "Fz"),
                     central  = c("FC1", "FC2", "FC5", "FC6", "T7", "T8", "C3", "C4", "Cz", "CP5", "CP6", "CP1", "CP2"),
                     posterior = c("P3", "P4", "P7", "P8", "Pz", "PO3", "PO4", "Oz", "O1", "O2"))

# Crear una función que asigne cada canal a su respectiva ROI
get_ROI <- function(channel) {
  for (ROI in ROI_names) {
    if (channel %in% ROI_channels[[ROI]]) {
      return(ROI)
    }
  }
}

# Crear variable ROI en el data frame original
unh_roi$ROI <- sapply(unh_roi$channel, get_ROI)
combat_roi$ROI <- sapply(combat_roi$channel, get_ROI)
neuroh_roi$ROI <- sapply(neuroh_roi$channel, get_ROI)
opcomb_roi$ROI <- sapply(opcomb_roi$channel, get_ROI)


unh_roi <- unh_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value_roi = mean(value))
combat_roi <- combat_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value_roi = mean(value))
neuroh_roi <- neuroh_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value_roi = mean(value))
opcomb_roi <- opcomb_roi %>% group_by(participant_id, feature, ROI) %>% summarise (value_roi = mean(value))


unh_roi <- unh_roi %>%
  inner_join(covars, by="participant_id")

combat_roi <- combat_roi %>%
  inner_join(covars, by="participant_id")

neuroh_roi <- neuroh_roi %>%
  inner_join(covars, by="participant_id")

opcomb_roi <- opcomb_roi %>%
  inner_join(covars, by="participant_id")

ggscatter(neuroh_roi, x = "age", y = "value",
          color = 'center',size=0.25, alpha = 1, facet.by='feature',
          ylab = "Not-harmonized DF Bandwidth (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))

unh_roi$ROI <- as.factor(unh_roi$ROI)
unh_roi$feature <- as.factor(unh_roi$feature)
unh_roi$gender <- factor(unh_roi$gender)
unh_roi$participant_id <- as.factor(unh_roi$participant_id)
unh_roi$feat_ROI <- as.factor(paste0(unh_roi$feature,'_',unh_roi$ROI))

combat_roi$ROI <- as.factor(combat_roi$ROI)
combat_roi$feature <- as.factor(combat_roi$feature)
combat_roi$gender <- factor(combat_roi$gender)
combat_roi$participant_id <- as.factor(combat_roi$participant_id)
combat_roi$feat_ROI <- as.factor(paste0(combat_roi$feature,'_',combat_roi$ROI))

neuroh_roi$ROI <- as.factor(neuroh_roi$ROI)
neuroh_roi$feature <- as.factor(neuroh_roi$feature)
neuroh_roi$gender <- factor(neuroh_roi$gender)
neuroh_roi$participant_id <- as.factor(neuroh_roi$participant_id)
neuroh_roi$feat_ROI <- as.factor(paste0(neuroh_roi$feature,'_',neuroh_roi$ROI))

opcomb_roi$ROI <- as.factor(opcomb_roi$ROI)
opcomb_roi$feature <- as.factor(opcomb_roi$feature)
opcomb_roi$gender <- factor(opcomb_roi$gender)
opcomb_roi$participant_id <- as.factor(opcomb_roi$participant_id)
opcomb_roi$feat_ROI <- as.factor(paste0(opcomb_roi$feature,'_',opcomb_roi$ROI))



##### NOT HARMONIZED & LMM ####
unh_roi <- data.frame(unh_roi)
unh_roi <- unh_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mmrm2 <- mmrm(value_roi_z ~ age_z + gender + 
              feature*ROI + us(feat_ROI|participant_id), 
              data=unh_roi, reml=T)
tab_model(mmrm2, digits = 4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)

# Delta
delta_roi <-  unh_roi  %>%
  filter(feature == "delta_pw")
delta_roi <- data.frame(delta_roi)
delta_roi <- delta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_delta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=delta_roi, reml=T)
tab_model(mm_delta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_delta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Estimated delta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# theta
theta_roi <-  unh_roi  %>%  filter(feature == "theta_pw")
theta_roi <- data.frame(theta_roi)
theta_roi <- theta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_theta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=theta_roi, reml=T)
tab_model(mm_theta, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_theta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Est. theta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# alpha
alpha_roi <-  unh_roi  %>%  filter(feature == "alpha_pw")
alpha_roi <- data.frame(alpha_roi)
alpha_roi <- alpha_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_alpha <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=alpha_roi, reml=T)
tab_model(mm_alpha, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_alpha, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated alpha Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')



# beta

beta_roi <-  unh_roi  %>%  filter(feature == "beta_pw")
beta_roi <- data.frame(beta_roi)
beta_roi <- beta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_beta <- mmrm(value_roi_z ~ age_z +gender + 
                  ROI + us(ROI|participant_id), 
                data=beta_roi, reml=T)
tab_model(mm_beta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_beta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated beta Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# extalpha_cf
extalpha_cf_roi <-  unh_roi  %>%  filter(feature == "extalphas_cf")
extalpha_cf_roi <- data.frame(extalpha_cf_roi)
extalpha_cf_roi <- extalpha_cf_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalpha_cf <- mmrm(value_roi_z ~ age_z + gender + 
                         ROI + us(ROI|participant_id), 
                       data=extalpha_cf_roi, reml=T, na.action(na.omit))
tab_model(mm_extalpha_cf, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalpha_cf, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Dominant Frequency (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# extalphas_bw
extalphas_bw_roi <-  unh_roi  %>%  filter(feature == "extalphas_bw")
extalphas_bw_roi <- data.frame(extalphas_bw_roi)
extalphas_bw_roi <- extalphas_bw_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalphas_bw <- mmrm(value_roi_z ~ age_z + gender + 
                          ROI + us(ROI|participant_id), 
                        data=extalphas_bw_roi, reml=T)
tab_model(mm_extalphas_bw, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalphas_bw, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Peak Freq. Bandwidth (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# exponent
exponent_roi <-  unh_roi  %>%  filter(feature == "exponent")
exponent_roi <- data.frame(exponent_roi)
exponent_roi <- exponent_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_exponent <- mmrm(value_roi_z ~ age_z + gender + 
                      ROI + us(ROI|participant_id), 
                    data=exponent_roi, reml=T)
tab_model(mm_exponent,  digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_exponent, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Exponent (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# offset
offset_roi <-  unh_roi  %>%  filter(feature == "offset")
offset_roi <- data.frame(offset_roi)
offset_roi <- offset_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_offset <- mmrm(value_roi_z ~ age_z +  gender + 
                    ROI + us(ROI|participant_id), 
                  data=offset_roi, reml=T)
tab_model(mm_offset, digits=4,   string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_offset, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Offset (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

tab_model(mm_delta, mm_theta, mm_alpha, mm_beta, 
          mm_extalpha_cf, mm_extalphas_bw, mm_exponent, mm_offset,
          df.method = "satterthwaite",
          dv.labels = c("Delta Pw", "Theta Pw", "Alpha Pw", "Beta Pw",
                        "Dom. Freq", "DFV", "Exponent", "Offset"),
          string.p = "FDR",  digits=3,  show.df=T,  p.adjust='fdr', collapse.ci =T,
          show.aic = T)




# OBSERVED DATA PLOTS

ggscatter(delta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(theta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(alpha_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(beta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalphas_bw_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized DF Bandwidth (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")



ggscatter(exponent_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Exponent (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")


ggscatter(offset_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Offset (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")



#### NEUROCOMBAT - LMM ####

# WORKING PLOT & LMM
combat_roi <- data.frame(combat_roi)
combat_roi <- combat_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mmrm2 <- mmrm(value_roi_z ~ age_z + gender + 
                feature*ROI + us(feat_ROI|participant_id), 
              data=combat_roi, reml=T)
tab_model(mmrm2, digits = 4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)

# Delta
delta_roi <-  combat_roi  %>%
  filter(feature == "delta_pw")
delta_roi <- data.frame(delta_roi)
delta_roi <- delta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_delta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=delta_roi, reml=T)
tab_model(mm_delta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_delta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Estimated delta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# theta
theta_roi <-  combat_roi  %>%  filter(feature == "theta_pw")
theta_roi <- data.frame(theta_roi)
theta_roi <- theta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_theta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=theta_roi, reml=T)
tab_model(mm_theta, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_theta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Est. theta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# alpha
alpha_roi <-  combat_roi  %>%  filter(feature == "alpha_pw")
alpha_roi <- data.frame(alpha_roi)
alpha_roi <- alpha_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_alpha <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=alpha_roi, reml=T)
tab_model(mm_alpha, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_alpha, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated alpha Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')



# beta

beta_roi <-  combat_roi  %>%  filter(feature == "beta_pw")
beta_roi <- data.frame(beta_roi)
beta_roi <- beta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_beta <- mmrm(value_roi_z ~ age_z +gender + 
                  ROI + us(ROI|participant_id), 
                data=beta_roi, reml=T)
tab_model(mm_beta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_beta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated beta Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# extalpha_cf
extalpha_cf_roi <-  combat_roi  %>%  filter(feature == "extalphas_cf")
extalpha_cf_roi <- data.frame(extalpha_cf_roi)
extalpha_cf_roi <- extalpha_cf_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalpha_cf <- mmrm(value_roi_z ~ age_z + gender + 
                         ROI + us(ROI|participant_id), 
                       data=extalpha_cf_roi, reml=T, na.action(na.omit))
tab_model(mm_extalpha_cf, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalpha_cf, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Dominant Frequency (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# extalphas_bw
extalphas_bw_roi <-  combat_roi  %>%  filter(feature == "extalphas_bw")
extalphas_bw_roi <- data.frame(extalphas_bw_roi)
extalphas_bw_roi <- extalphas_bw_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalphas_bw <- mmrm(value_roi_z ~ age_z + gender + 
                          ROI + us(ROI|participant_id), 
                        data=extalphas_bw_roi, reml=T)
tab_model(mm_extalphas_bw, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalphas_bw, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Peak Freq. Bandwidth (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# exponent
exponent_roi <-  combat_roi  %>%  filter(feature == "exponent")
exponent_roi <- data.frame(exponent_roi)
exponent_roi <- exponent_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_exponent <- mmrm(value_roi_z ~ age_z + gender + 
                      ROI + us(ROI|participant_id), 
                    data=exponent_roi, reml=T)
tab_model(mm_exponent,  digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_exponent, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Exponent (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# offset
offset_roi <-  combat_roi  %>%  filter(feature == "offset")
offset_roi <- data.frame(offset_roi)
offset_roi <- offset_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_offset <- mmrm(value_roi_z ~ age_z +  gender + 
                    ROI + us(ROI|participant_id), 
                  data=offset_roi, reml=T)
tab_model(mm_offset, digits=4,   string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_offset, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Offset (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

tab_model(mm_delta, mm_theta, mm_alpha, mm_beta, 
          mm_extalpha_cf, mm_extalphas_bw, mm_exponent, mm_offset,
          df.method = "satterthwaite",
          dv.labels = c("Delta Pw", "Theta Pw", "Alpha Pw", "Beta Pw",
                        "Dom. Freq", "DFV", "Exponent", "Offset"),
          string.p = "FDR", title = "neuroCombat",  digits=3,  show.df=T,  p.adjust='fdr', collapse.ci =T,
          show.aic = T)





# HARMONIZED DATA PLOTS

ggscatter(delta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(theta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(alpha_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(beta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalphas_bw_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized DF Bandwidth (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")



ggscatter(exponent_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Exponent (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")


ggscatter(offset_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Offset (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")

#### NEUROHARMONIZE- LMM ####

# WORKING PLOT & LMM
neuroh_roi<- data.frame(combat_roi)
neuroh_roi<- neuroh_roi%>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mmrm2 <- mmrm(value_roi_z ~ age_z + gender + 
                feature*ROI + us(feat_ROI|participant_id), 
              data=combat_roi, reml=T)
tab_model(mmrm2, digits = 4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)

# Delta
delta_roi <-  neuroh_roi %>%
  filter(feature == "delta_pw")
delta_roi <- data.frame(delta_roi)
delta_roi <- delta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_delta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=delta_roi, reml=T)
tab_model(mm_delta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_delta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Estimated delta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# theta
theta_roi <-  neuroh_roi %>%  filter(feature == "theta_pw")
theta_roi <- data.frame(theta_roi)
theta_roi <- theta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_theta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=theta_roi, reml=T)
tab_model(mm_theta, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_theta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Est. theta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# alpha
alpha_roi <-  neuroh_roi %>%  filter(feature == "alpha_pw")
alpha_roi <- data.frame(alpha_roi)
alpha_roi <- alpha_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_alpha <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=alpha_roi, reml=T)
tab_model(mm_alpha, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_alpha, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated alpha Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')



# beta

beta_roi <-  neuroh_roi %>%  filter(feature == "beta_pw")
beta_roi <- data.frame(beta_roi)
beta_roi <- beta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_beta <- mmrm(value_roi_z ~ age_z +gender + 
                  ROI + us(ROI|participant_id), 
                data=beta_roi, reml=T)
tab_model(mm_beta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_beta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated beta Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# extalpha_cf
extalpha_cf_roi <-  neuroh_roi %>%  filter(feature == "extalphas_cf")
extalpha_cf_roi <- data.frame(extalpha_cf_roi)
extalpha_cf_roi <- extalpha_cf_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalpha_cf <- mmrm(value_roi_z ~ age_z + gender + 
                         ROI + us(ROI|participant_id), 
                       data=extalpha_cf_roi, reml=T, na.action(na.omit))
tab_model(mm_extalpha_cf, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalpha_cf, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Dominant Frequency (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# extalphas_bw
extalphas_bw_roi <-  neuroh_roi %>%  filter(feature == "extalphas_bw")
extalphas_bw_roi <- data.frame(extalphas_bw_roi)
extalphas_bw_roi <- extalphas_bw_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalphas_bw <- mmrm(value_roi_z ~ age_z + gender + 
                          ROI + us(ROI|participant_id), 
                        data=extalphas_bw_roi, reml=T)
tab_model(mm_extalphas_bw, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalphas_bw, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Peak Freq. Bandwidth (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# exponent
exponent_roi <-  neuroh_roi %>%  filter(feature == "exponent")
exponent_roi <- data.frame(exponent_roi)
exponent_roi <- exponent_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_exponent <- mmrm(value_roi_z ~ age_z + gender + 
                      ROI + us(ROI|participant_id), 
                    data=exponent_roi, reml=T)
tab_model(mm_exponent,  digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_exponent, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Exponent (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# offset
offset_roi <-  neuroh_roi %>%  filter(feature == "offset")
offset_roi <- data.frame(offset_roi)
offset_roi <- offset_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_offset <- mmrm(value_roi_z ~ age_z +  gender + 
                    ROI + us(ROI|participant_id), 
                  data=offset_roi, reml=T)
tab_model(mm_offset, digits=4,   string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_offset, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Offset (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

tab_model(mm_delta, mm_theta, mm_alpha, mm_beta, 
          mm_extalpha_cf, mm_extalphas_bw, mm_exponent, mm_offset,
          df.method = "satterthwaite",
          dv.labels = c("Delta Pw", "Theta Pw", "Alpha Pw", "Beta Pw",
                        "Dom. Freq", "DFV", "Exponent", "Offset"),
          string.p = "FDR", title = "neuroHarmonize",  digits=3,  show.df=T,  p.adjust='fdr', collapse.ci =T,
          show.aic = T)
# HARMONIZED DATA PLOTS

ggscatter(delta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(theta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(alpha_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(beta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Beta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1,
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")

ggscatter(extalphas_bw_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized DF Bandwidth (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")



ggscatter(exponent_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Exponent (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")


ggscatter(offset_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Offset (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")


#### OPN COMBAT_GMM - LMM ####

# WORKING PLOT & LMM
opcomb_roi<- data.frame(combat_roi)
opcomb_roi<- opcomb_roi%>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mmrm2 <- mmrm(value_roi_z ~ age_z + gender + 
                feature*ROI + us(feat_ROI|participant_id), 
              data=combat_roi, reml=T)
tab_model(mmrm2, digits = 4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)

# Delta
delta_roi <-  opcomb_roi %>%
  filter(feature == "delta_pw")
delta_roi <- data.frame(delta_roi)
delta_roi <- delta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_delta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=delta_roi, reml=T)
tab_model(mm_delta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_delta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Estimated delta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# theta
theta_roi <-  opcomb_roi %>%  filter(feature == "theta_pw")
theta_roi <- data.frame(theta_roi)
theta_roi <- theta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_theta <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=theta_roi, reml=T)
tab_model(mm_theta, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_theta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Z-Est. theta Power (a.u.)', xlab = 'Z-Age (years)',
      tlab = 'Gender')

# alpha
alpha_roi <-  opcomb_roi %>%  filter(feature == "alpha_pw")
alpha_roi <- data.frame(alpha_roi)
alpha_roi <- alpha_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_alpha <- mmrm(value_roi_z ~ age_z + gender + 
                   ROI + us(ROI|participant_id), 
                 data=alpha_roi, reml=T)
tab_model(mm_alpha, digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_alpha, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated alpha Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')



# beta

beta_roi <-  opcomb_roi %>%  filter(feature == "beta_pw")
beta_roi <- data.frame(beta_roi)
beta_roi <- beta_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_beta <- mmrm(value_roi_z ~ age_z +gender + 
                  ROI + us(ROI|participant_id), 
                data=beta_roi, reml=T)
tab_model(mm_beta, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_beta, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Estimated beta Power (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# extalpha_cf
extalpha_cf_roi <-  opcomb_roi %>%  filter(feature == "extalphas_cf")
extalpha_cf_roi <- data.frame(extalpha_cf_roi)
extalpha_cf_roi <- extalpha_cf_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalpha_cf <- mmrm(value_roi_z ~ age_z + gender + 
                         ROI + us(ROI|participant_id), 
                       data=extalpha_cf_roi, reml=T, na.action(na.omit))
tab_model(mm_extalpha_cf, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalpha_cf, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Dominant Frequency (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# extalphas_bw
extalphas_bw_roi <-  opcomb_roi %>%  filter(feature == "extalphas_bw")
extalphas_bw_roi <- data.frame(extalphas_bw_roi)
extalphas_bw_roi <- extalphas_bw_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))


mm_extalphas_bw <- mmrm(value_roi_z ~ age_z + gender + 
                          ROI + us(ROI|participant_id), 
                        data=extalphas_bw_roi, reml=T)
tab_model(mm_extalphas_bw, digits=4, string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_extalphas_bw, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Peak Freq. Bandwidth (Hz)', xlab = 'Age (years)',
      tlab = 'Gender')


# exponent
exponent_roi <-  opcomb_roi %>%  filter(feature == "exponent")
exponent_roi <- data.frame(exponent_roi)
exponent_roi <- exponent_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))
mm_exponent <- mmrm(value_roi_z ~ age_z + gender + 
                      ROI + us(ROI|participant_id), 
                    data=exponent_roi, reml=T)
tab_model(mm_exponent,  digits=4,  string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_exponent, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Exponent (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

# offset
offset_roi <-  opcomb_roi %>%  filter(feature == "offset")
offset_roi <- data.frame(offset_roi)
offset_roi <- offset_roi %>% 
  mutate(value_roi_z = scale(value_roi),
         age_z= scale(age))

mm_offset <- mmrm(value_roi_z ~ age_z +  gender + 
                    ROI + us(ROI|participant_id), 
                  data=offset_roi, reml=T)
tab_model(mm_offset, digits=4,   string.p = "p", show.df=T,  p.adjust='none', show.aic = T)


mylist <- list(age_z=seq(-3,3,by=0.5),gender=c('f','m'))
emmip(mm_offset, gender~ age_z | ROI, mult.name = "ROI", at = mylist, 
      CIs=T, ylab='Est. Aperiodic Offset (a.u.)', xlab = 'Age (years)',
      tlab = 'Gender')

tab_model(mm_delta, mm_theta, mm_alpha, mm_beta, 
          mm_extalpha_cf, mm_extalphas_bw, mm_exponent, mm_offset,
          df.method = "satterthwaite",
          dv.labels = c("Delta Pw", "Theta Pw", "Alpha Pw", "Beta Pw",
                        "Dom. Freq", "DFV", "Exponent", "Offset"),
          string.p = "FDR", title = "OPN ComBat-GMM",  digits=3,  show.df=T,  p.adjust='fdr', collapse.ci =T,
          show.aic = T)
# HARMONIZED DATA PLOTS

ggscatter(delta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Delta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(theta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized theta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(alpha_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Alpha Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(beta_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',ylab = "Not-harmonized Beta Power (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")

ggscatter(extalpha_cf_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1,
          ylab = "Not-harmonized Dom. Freq (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")

ggscatter(extalphas_bw_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized DF Bandwidth (Hz)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")



ggscatter(exponent_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Exponent (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "loess",  color = "black")


ggscatter(offset_roi, x = "age", y = "value_roi",
          conf.int = T,
          color = 'center',size=1, alpha = 1, facet.by='ROI',
          ylab = "Not-harmonized Offset (a.u.)",
          xlab = "Age (years)", ggtheme = theme_pubr(border = T, base_family = "mono"))+
  geom_smooth(method = "lm",  color = "black")



