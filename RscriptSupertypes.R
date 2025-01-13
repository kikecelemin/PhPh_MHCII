# Supertypes
#install.packages("ggrepel")
#install.packages("Rtsne")

library(Rtsne)
library(ggplot2)
library(ggrepel)
library(cluster)

setwd("C:/Users/kikec/Desktop/KIKE/PhD_Potsdam_Universität/Harbour_Porpoise/Whole_genome_resequencing/WGS_Data_Analysis/3rd_paper/04_Supertypes")

####______________ DQB___________ ####
#Cleaning: removing duplicates
#zDQB <- read.delim2("02_DQB_z_matrix_noDupl")
#zDQB$Allele=NULL
#duplicate_rows <- zDQB[duplicated(zDQB), ] #these are the duplicates
#duplicate_rows
# Rows 11=6, 15=12, 21=5 are duplicates! Rename, Remove. 
#zDQB <- read.delim2("03_DQB_z_matrix_noDupl")
#zDQB$Allele=NULL
#duplicate_rows <- zDQB[duplicated(zDQB), ] #these are the duplicates
#duplicate_rows # no duplicates left
zDQB <- read.delim2("03_DQB_z_matrix_noDupl")

#Normalising df
# 1. Extract numeric columns (excluding the first column if it contains identifiers)
numeric_colDQB <- zDQB[, sapply(zDQB, is.numeric)]

# 2. Center and Normalize the Data
scaled_DQBraw <- scale(numeric_colDQB)

# 3. Convert the Scaled Data Back to a Dataframe
scaled_DQB <- as.data.frame(scaled_DQBraw)
scaled_DQB$Allele=zDQB$Allele
head(scaled_DQB)
zDQB=scaled_DQB

###ANALYSIS DQB
#t-SNE

DQB_tsne <- Rtsne(zDQB, dims = 2, perplexity=3) #PERPLEXITY
tsne_df_DQB <- as.data.frame(DQB_tsne$Y)
tsne_df_DQB$Allele <- zDQB$Allele  # Add labels from the first column

pdf("DQB_tSNE.pdf")
ggplot(tsne_df_DQB, aes(x = V1, y = V2, label = Allele)) +
  geom_point() +
  geom_text_repel(size = 3) +  # Add text labels with automatic repulsion
  labs(title = "t-SNE Plot", x = "Dimension 1", y = "Dimension 2") +
  theme_classic()
dev.off()

DQB_tsne <- Rtsne(zDQB, dims = 2, perplexity=3)
plot(DQB_tsne$Y, main = "t-SNE Plot")
text(DQB_tsne$Y, labels = zDQB$Allele, pos = 3, col = "red", cex = 0.6)

####PCA###
pca_result <- prcomp(zDQB[, -6], scale. = TRUE)  # Exclude the first column (if it contains non-numeric identifiers)
pca_result$rotation
pca_result$sdev
pca_scores <- as.data.frame(pca_result$x)  # Convert to data frame for further analysis or visualization
pca_scores$ID <- zDQB[, 6]  # Assuming the first column contains identifiers

pdf("DQB_PCA.pdf")
ggplot(pca_scores, aes(x = PC1, y = PC2, label = ID)) +
  geom_point() +
  geom_text_repel(size=2.5) +  # Add text labels with automatic repulsion
  labs(title = "PCA Plot", x = "PC1 (37.5%)", y = "PC2 (32%)") +
  theme_classic()
dev.off()

standard_deviations <- pca_result$sdev
variance_explained <- (standard_deviations^2) / sum(standard_deviations^2)
print(variance_explained)


arrow_data <- data.frame(
  Variable = colnames(pca_result$rotation),
  PC1_end = pca_result$rotation[, 1],
  PC2_end = pca_result$rotation[, 2])
arrow_data$ID=c("Z1","Z2","Z3","Z4","Z5")

ggplot(pca_scores, aes(x = PC1, y = PC2, label = ID)) +
  geom_point() +
  geom_text_repel() +  # Add text labels with automatic repulsion
  geom_segment(data = arrow_data, aes(x = 0, y = 0, xend = PC1_end, yend = PC2_end),
               arrow = arrow(length = unit(0.2, "inches")), color = "blue") +
  labs(title = "PCA Plot with Variable Loadings", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()


#### Hierarchical clusters

dis<-dist(zDQB[,-6])
clusters <- hclust(dis,method = "ward.D2")
clusters$labels <- zDQB$Allele

pdf("DQB_Hierarchichal_clustering.pdf")
plot(clusters,main = "Dendrogram of Hierarchical Clustering")
dev.off()

?hclust



#number of clusters elbow method
wcss <- rep(NA, length = nrow(zDQB))
for (i in 1:length(wcss)) {
  # Cut the dendrogram at height i
  clusters <- cutree(hc, h = i)
  # Calculate the within-cluster sum of squares
  wcss[i] <- sum((dist_matrix^2) * as.numeric(clusters == clusters[i]))
}

# 5. Plot WCSS against Number of Clusters
plot(1:length(wcss), wcss, type = "b", pch = 19, 
     xlab = "Number of Clusters", ylab = "Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Method for Hierarchical Clustering")
#Elbow increasing -> bad solutions:
#Data Scaling: Ensure that your data is scaled appropriately before clustering. Hierarchical clustering is sensitive to the scale of the variables, so standardizing or normalizing your data may be necessary.
#Distance Metric: Check the distance metric used to compute the distance matrix. Depending on the nature of your data, different distance metrics (e.g., Euclidean distance, Manhattan distance, correlation distance) may be more suitable.

#different mehtod "silhuette"
hc <- hclust(dis, method = "complete")  # Use complete linkage method
silhouette_scores <- c()
for (i in 2:(nrow(data) - 1)) {
  # Cut the dendrogram at height i
  clusters <- cutree(hc, h = i)
  # Calculate silhouette score
  silhouette_scores[i] <- mean(silhouette(clusters, dist_matrix)$widths)
}

# 5. Plot Silhouette Scores
plot(2:(nrow(data) - 1), silhouette_scores[2:(nrow(data) - 1)], type = "b",
     xlab = "Number of Clusters", ylab = "Average Silhouette Width",
     main = "Silhouette Method for Hierarchical Clustering")




####______________ DRB2 _____________####
#Rm duplicates
#zDRB2 <- read.delim2("~/Schreibtisch/MHC_Phph_Paper/Supertypes/02_DRB2_z_matrix")
#zDRB2$Allele=NULL
#duplicate_rows <- zDRB2[duplicated(zDRB2), ] #these are the duplicates
#duplicate_rows # 3=4; 11=22; 16=23 -> remove fronm file

#zDRB2 <- read.delim2("~/Schreibtisch/MHC_Phph_Paper/Supertypes/03_DRB2_z_matrix_noDupl")
#zDRB2$Allele=NULL
#duplicate_rows <- zDRB2[duplicated(zDRB2), ] #these are the duplicates
#duplicate_rows # no left
zDRB2 <- read.delim2("03_DRB2_z_matrix_noDupl")

#Normalize
# 1. Extract numeric columns (excluding the first column if it contains identifiers)
numeric_colDRB2 <- zDRB2[, sapply(zDRB2, is.numeric)]

# 2. Center and Normalize the Data
scaled_DRB2raw <- scale(numeric_colDRB2)

# 3. Convert the Scaled Data Back to a Dataframe
scaled_DRB2 <- as.data.frame(scaled_DRB2raw)
scaled_DRB2$Allele=zDRB2$Allele
head(scaled_DRB2)
zDRB2=scaled_DRB2


###ANALYSIS DRB2
#t-SNE

DRB2_tsne <- Rtsne(zDRB2, dims = 2, perplexity=3)
tsne_df_DRB2 <- as.data.frame(DRB2_tsne$Y)
tsne_df_DRB2$Allele <- zDRB2$Allele  # Add labels from the first column

pdf("DRB1_tSNE.pdf")
ggplot(tsne_df_DRB2, aes(x = V1, y = V2, label = Allele)) +
  geom_point() +
  geom_text_repel(size = 3) +  # Add text labels with automatic repulsion
  labs(title = "t-SNE Plot", x = "Dimension 1", y = "Dimension 2") +
  theme_classic()
dev.off()


DRB2_tsne <- Rtsne(zDRB2, dims = 2, perplexity=8)
plot(DRB2_tsne$Y, main = "t-SNE Plot")
text(DRB2_tsne$Y, labels = zDRB2$Allele, pos = 3, col = "red", cex = 0.6)


####PCA###
pca_result <- prcomp(zDRB2[, -6], scale. = TRUE)  # Exclude the first column (if it contains non-numeric identifiers)
pca_result$rotation
pca_result$sdev
pca_scores <- as.data.frame(pca_result$x)  # Convert to data frame for further analysis or visualization
pca_scores$ID <- zDRB2[, 6]  # Assuming the first column contains identifiers

pdf("DRB1_PCA.pdf")
ggplot(pca_scores, aes(x = PC1, y = PC2, label = ID)) +
  geom_point() +
  geom_text_repel(size=2.5) +  # Add text labels with automatic repulsion
  labs(title = "PCA Plot", x = "PC1 68%", y = "PC2 20%") +
  theme_classic()
dev.off()

standard_deviations <- pca_result$sdev
variance_explained <- (standard_deviations^2) / sum(standard_deviations^2)
print(variance_explained)

#hier clust
dis<-dist(zDRB2[,-6])
clusters <- hclust(dis,method = "ward.D2")
clusters$labels <- zDRB2$Allele

pdf("DRB2_Hierarchichal_clustering.pdf")
plot(clusters,main = "Dendrogram of Hierarchical Clustering")
dev.off()

