## ---------------------------
## 0. Packages & Palette
## ---------------------------

# Make sure these are installed once (uncomment if needed)
# install.packages(c("ggplot2", "reshape2", "GGally", "khroma", "dendextend"))

library(ggplot2)
library(reshape2)
library(GGally)
library(khroma)
library(dendextend)

# Create a nice colourblind-friendly palette
make_palette <- function() {
  pal <- colour("okabeito")(8)
  pal <- pal[c(2:8, 1)]    # rotate so greys come last
  names(pal) <- NULL
  pal
}


## ---------------------------
## 1. Data Loading & Indexing
## ---------------------------

load_restaurants <- function(path = "MichelinNY.csv") {
  Restaurants <- read.csv(path)
  
  # Assume column 1 is Michelin indicator: 1 = in Michelin, 0 = not
  mich_idx    <- which(Restaurants[, 1] == 1)
  no_mich_idx <- which(Restaurants[, 1] == 0)
  
  list(
    data      = Restaurants,
    mich      = mich_idx,
    no_mich   = no_mich_idx
  )
}


## ---------------------------
## 2. Correlation & Covariance
## ---------------------------

compute_cor_cov <- function(Restaurants, mich, no_mich) {
  # Numeric columns of interest: Food, Decor, Service, Price (assumed 3:6)
  num_cols <- 3:6
  
  cat("Correlation (Michelin only):\n")
  print(cor(Restaurants[mich,   num_cols]))
  
  cat("\nCorrelation (Non-Michelin only):\n")
  print(cor(Restaurants[no_mich, num_cols]))
  
  cat("\nCovariance (all restaurants):\n")
  print(cov(Restaurants[, num_cols]))
}


## ---------------------------
## 3. PCA Helpers
## ---------------------------

run_pca <- function(Restaurants, rows = NULL, cols = 3:6, scale. = TRUE) {
  # rows = NULL means use all rows
  if (is.null(rows)) {
    dat <- Restaurants[, cols]
  } else {
    dat <- Restaurants[rows, cols]
  }
  prcomp(dat, scale. = scale.)
}

plot_scree <- function(pca, title = "Scree Plot") {
  # Proportion of variance explained
  var_prop <- (pca$sdev^2) / sum(pca$sdev^2)
  
  plot(
    var_prop,
    type = "b",
    main = title,
    xlab = "No. of components",
    ylab = "Proportion of variance explained",
    xaxt = "n"
  )
  axis(1, at = seq_along(var_prop))
}

plot_pca_loadings <- function(pca, var_names, main_prefix = "") {
  # Barplots of loadings for each PC
  par(mfrow = c(2, 2))
  for (i in 1:ncol(pca$rotation)) {
    barplot(
      pca$rotation[, i],
      names.arg = var_names,
      ylim = c(-1, 1),
      col = as.factor(var_names),
      main = paste(main_prefix, "PC", i)
    )
  }
  par(mfrow = c(1, 1))
}

ggplot_pca_loadings <- function(pca, var_names, palette) {
  rot_df      <- as.data.frame(pca$rotation)
  rot_df$Var  <- var_names
  rot_melt    <- melt(rot_df, id.vars = "Var", variable.name = "PC", value.name = "Coefficient")
  
  ggplot(rot_melt, aes(x = Var, y = Coefficient, fill = Var)) +
    geom_col() +
    facet_wrap(~PC, scales = "free") +
    theme(legend.position = "none") +
    ylim(-1, 1) +
    labs(title = "Variable Coefficients for each Principal Component",
         x = NULL, y = "Loading") +
    scale_fill_manual(values = palette)
}


## ---------------------------
## 4. Simple Plots (Food vs Price, Pairs)
## ---------------------------

plot_food_price <- function(Restaurants, palette) {
  Food  <- Restaurants[, 3]
  Price <- Restaurants[, 6]
  
  # Assuming column 1 is Michelin indicator
  michelin_factor <- factor(Restaurants[, 1], labels = c("No Michelin", "Michelin"))
  
  plot(
    Food, Price,
    col  = palette[as.integer(michelin_factor)],
    pch  = as.integer(michelin_factor),
    main = "Scatter Plot of Food vs Price",
    xlab = "Food",
    ylab = "Price"
  )
  legend("topleft",
         legend = levels(michelin_factor),
         col    = palette[1:2],
         pch    = 1:2,
         bty    = "n")
}

base_pairs_plot <- function(Restaurants) {
  # Colours by Michelin indicator (assume col 1)
  michelin_factor <- factor(Restaurants[, 1], labels = c("No Michelin", "Michelin"))
  pairs(
    Restaurants[, 3:6],
    col  = as.integer(michelin_factor),
    pch  = as.integer(michelin_factor),
    main = "Pairs Plot for Restaurants (Base R)"
  )
}

ggpairs_plot <- function(Restaurants) {
  michelin_factor <- factor(Restaurants[, 1], labels = c("No Michelin", "Michelin"))
  ggpairs(
    Restaurants,
    columns = 3:6,
    mapping = aes(colour = michelin_factor),
    title   = "Pairs Plot for Restaurants (GGally)"
  )
}


## ---------------------------
## 5. Hierarchical Clustering
## ---------------------------

run_hclustering <- function(Restaurants, palette) {
  # Standardise each numeric variable (columns 3:6) across restaurants
  num_data <- Restaurants[, 3:6]
  num_sd   <- apply(num_data, 2, sd)
  num_mean <- apply(num_data, 2, mean)
  
  standard_data <- sweep(num_data, 2, num_mean, "-")
  standard_data <- sweep(standard_data, 2, num_sd,   "/")
  
  # Distances
  dis_euclidean <- dist(standard_data, method = "euclidean")
  dis_manhattan <- dist(standard_data, method = "manhattan")
  
  # Hierarchical clustering
  clust_complete <- hclust(dis_euclidean, method = "complete")
  clust_average  <- hclust(dis_euclidean, method = "average")
  
  # Basic dendrograms
  par(mfrow = c(1, 2))
  plot(clust_complete, labels = FALSE, main = "Complete-Linkage Dendrogram")
  plot(clust_average,  labels = FALSE, main = "Average-Linkage Dendrogram")
  par(mfrow = c(1, 1))
  
  # Cut into clusters
  k_complete <- 6
  cut_complete <- cutree(clust_complete, k = k_complete)
  
  # Colour branches
  dend_complete <- as.dendrogram(clust_complete)
  dend_complete <- color_branches(
    dend_complete,
    clusters = cut_complete[clust_complete$order],
    col      = rep(palette, length.out = length(unique(cut_complete)))
  )
  
  plot(dend_complete, leaflab = "none", main = "Coloured Dendrogram (Complete)")
  
  invisible(list(
    clust_complete = clust_complete,
    clust_average  = clust_average,
    cut_complete   = cut_complete,
    dend_complete  = dend_complete
  ))
}


## ---------------------------
## 6. Main Driver
## ---------------------------

run_analysis <- function(csv_path = "MichelinNY.csv") {
  pal <- make_palette()
  
  # Load data and indices
  loaded <- load_restaurants(csv_path)
  Restaurants <- loaded$data
  mich        <- loaded$mich
  no_mich     <- loaded$no_mich
  
  # Basic correlations & covariance
  compute_cor_cov(Restaurants, mich, no_mich)
  
  # PCA: all restaurants
  pca_all <- run_pca(Restaurants, rows = NULL, cols = 3:6)
  cat("\nPCA (all restaurants):\n")
  print(summary(pca_all))
  
  # PCA: Michelin vs non-Michelin
  pca_mich    <- run_pca(Restaurants, rows = mich,    cols = 3:6)
  pca_no_mich <- run_pca(Restaurants, rows = no_mich, cols = 3:6)
  
  cat("\nPCA (Michelin only):\n")
  print(summary(pca_mich))
  
  cat("\nPCA (Non-Michelin only):\n")
  print(summary(pca_no_mich))
  
  # Scree plot for all restaurants
  plot_scree(pca_all, title = "Scree Plot - All Restaurants")
  
  # PCA loadings (barplots)
  var_names <- c("Food", "Decor", "Service", "Price")
  plot_pca_loadings(pca_all, var_names, main_prefix = "All Restaurants")
  
  # PCA loadings (ggplot)
  print(ggplot_pca_loadings(pca_all, var_names, pal))
  
  # Scatter and pairs plots
  plot_food_price(Restaurants, pal)
  base_pairs_plot(Restaurants)
  print(ggpairs_plot(Restaurants))
  
  # Hierarchical clustering
  run_hclustering(Restaurants, pal)
}

## ---------------------------
## Run
## ---------------------------
run_analysis("MichelinNY.csv")
