# ====================================================
# Author:         Lina Möller
# Date:           25.10.2024
# Description:    processing of the data, build datasets, calculate diversity indices, performing of MLE (via cma-es), analysis of parameters
# Ongoing Todos:  result check
# ====================================================


# ====================================================
# 1. Load packages
# ====================================================

#Packages in the order in which they are used
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
install.packages("janitor")
library(janitor)
install.packages("readxl")
library(readxl)
install.packages("mice")
library(mice)
install.packages("ggplot2")
library(ggplot2)
install.packages("factoextra")
library(factoextra)
install.packages("cmaes")
library(cmaes)


# ====================================================
# 2. Load data
# ====================================================

# coordinates
coordinates <- read.csv("/Users/linamoller/Documents/Uni/Master/Masterarbeit/R Skripte/Flexible Neighborhood Analysis/InputDaten/ident_Plot_coordinates.csv")
head(coordinates)
load("/Users/linamoller/Documents/Uni/Master/Masterarbeit/R Skripte/Flexible Neighborhood Analysis/InputDaten/IDENT-FR_trees_array.Rdata")

# productivity
inventory <- read_excel("/Users/linamoller/Documents/Uni/Master/Masterarbeit/R Skripte/Flexible Neighborhood Analysis/InputDaten/IDENT_inventory 2019-20_JF.xlsx", 
                        col_types = c("numeric", "numeric", "numeric", 
                                      "text", "text", "numeric", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric", "text")) %>% 
  separate(pos., c("nrow", "ncol")) %>%       # Split column with position into two columns with coordinates for x and y axes
  mutate(ncol = as.numeric(ncol)) %>%         # Changes the two newly created columns to data type=numeric
  mutate(nrow = as.numeric(nrow)) %>% 
  clean_names() %>%                           # cleans the column names (uniform style: lowercase letters, no special characters or spaces)
  drop_na(d1_mm) %>%                          # Removes all lines where d1_mm is empty
  select(1:7, 10)                             # selects rows (block, plot, pos. (2x), species, mixture, d1_mm)



# bark beetle infestation
infest <- as_tibble(read.csv("/Users/linamoller/Documents/Uni/Master/Masterarbeit/R Skripte/Flexible Neighborhood Analysis/InputDaten/InfestationInventory_Winter2020_FS_corrected.csv",
                             na.strings=c("","na"))) %>%
  # format data set to tibble
  # na.strings determines which character strings are to be evaluated as NA
  select(-Date.) %>%                          # remove the date column
  rename("nrow"=X.Pos) %>%                    # rename: X.Pos to nrow
  rename("ncol"=Y.Pos) %>% 
  clean_names()                               # Cleaning up column names




# ====================================================
# 3. Hilfsfunktionen 
# ====================================================

# ======================== 1. ======================== 
# function converts 3D array to dataframe, 
# each dimension is broken down into columns
# allows us to put coordinates into tabular form
array_to_dataframe <- function(data){
  ##create output dataframe
  result <- data.frame("dim" = rep(1:dim(data)[3],
                                   each = dim(data)[1] * dim(data)[2]),
                       "ncol" = rep(1:dim(data)[2],
                                    dim(data)[1] * dim(data)[3]),
                       "nrow" = rep(1:dim(data)[1], each = dim(data)[2],
                                    dim(data)[3]),
                       "value" = numeric(prod(dim(data))))
  datalength <- dim(data)[1] * dim(data)[2]
  ##fill the dataframe with values from the array
  for(i in seq_len(dim(data)[3])){
    result$value[(i - 1) * datalength + 1:datalength] <-
      pivot_longer(as.data.frame(data[, ,i]),
                   cols = seq_len(NCOL(data[,,i])))$value
    
  }
  return(result)
}

# ======================== 2. ======================== 
# function to calculate Shannon Index per block
calculate_shannon <- function(data) {
  data %>%
    group_by(block) %>%                              # grouping by block
    summarise(
      Shannon_Index = -sum(prop.table(table(species)) * 
                             log(prop.table(table(species))), na.rm = TRUE)
    )
}

# ======================== 3. ======================== 
# Function for calculating the Shannon index across all data
calculate_shannon_overall <- function(data, species_col) {
  # Calculation of the absolute frequencies of each species
  species_counts <- data %>%
    group_by(across(all_of(species_col))) %>% 
    summarise(count = n(), .groups = "drop")
  
  # Calculate the total number of individuals
  total_count <- sum(species_counts$count)
  
  # Calculate relative frequencies (p_i) for each type
  species_counts <- species_counts %>%
    mutate(p = count / total_count)
  
  # calculate Shannon-Index
  shannon_index <- -sum(species_counts$p * log(species_counts$p))
  
  return(shannon_index)
}

# ======================== 4. ======================== 
# weighted Shannon-Index (WS)
# calculates local diversity for each tree with weighting of distance

weighted_shannon_index <- function(data, 
                                   x_col, 
                                   y_col, 
                                   species_col, 
                                   incline = 1, 
                                   maxdist = Inf) {
  # 1) Naming the columns 
  # Does NOT change ‘data’ itself, but only creates copies
  # Simplifies further processing because you can always access the same column names
  df <- data
  colnames(df)[colnames(df) == x_col] <- "x"
  colnames(df)[colnames(df) == y_col] <- "y"
  colnames(df)[colnames(df) == species_col] <- "species"
  
  # 2) Factor for species (for faster grouping)
  # preparation for grouping by tree species
  df$species <- as.factor(df$species)
  sp_levels <- levels(df$species)        # all possible types as a vector of strings
  sp_codes  <- as.integer(df$species)    # creates a numeric vector that assigns a factor code to each row
  
  n <- nrow(df)
  
  # 3) Coordinate matrix for distance calculation
  coords <- cbind(df$x, df$y)            # matrix with 2 columns
  
  # 4) Calculate distance matrix (NxN)
  #    Warning: potentially very memory intensive
  dist_mat <- as.matrix(dist(coords, method = "euclidean"))   # calculates all pairwise Euclidean distances between rows
  
  # 5) If maxdist < Inf: set all distances > maxdist to NA (so that they are ignored)
  #    so trees that are too far away no longer matter
  if (!is.infinite(maxdist)) {
    dist_mat[dist_mat > maxdist] <- NA
  }
  
  # 6) calculate weights
  weights_mat <- exp(-dist_mat * incline)
  weights_mat[is.na(weights_mat)] <- 0
  # Ensure that the custom tree always has a weight of 1
  diag(weights_mat) <- 1
  
  # 7) For every tree i: 
  #    - Summing up the weights by type
  #    - Relative frequencies = proportions per species
  #    - Shannon-Index = exp( - sum( p * log(p) ) )
  #    - t_density = sum of all weights
  # shannon and density as numerical vectors of length n
  Shannon_vec <- numeric(n)
  t_density_vec <- numeric(n)
  
  for (i in seq_len(n)) {
    w_i <- weights_mat[i, ]     # extracts the weights that tree i assigns to all other trees
    # summing up weights per type
    # sp_codes[j] contains the type code of tree j
    # sum w_i[j] for every type
    sp_sums <- tapply(w_i, sp_codes, sum)         # totals weights by type
    # tapply calculates the sum of the weights for type 1, type 2, ...
    
    w_total <- sum(sp_sums)       # Sum of all weights of tree i
    # Prevent division by zero
    if (w_total == 0) {
      # If w_total is 0, set a default value
      Shannon_vec[i] <- 1  
      t_density_vec[i] <- 0
    } else {
      t_density_vec[i] <- w_total
      p <- sp_sums / w_total
      # 0 * log(0) is treated as 0, thus preventing NAs later on
      Shannon_vec[i] <- exp(-sum(ifelse(p == 0, 0, p * log(p))))
    }
  }
  
  # 8) Write result back to original dataset
  df$Shannon   <- Shannon_vec
  df$t_density <- t_density_vec
  
  return(df)
}





# ======================== 4. ======================== 
# weighted functional Shanno-Index (WFS)
# Based on the Shannon index, spatial and functional levels are taken into account
# 

weighted_fd_index <- function(data,
                              x_col,
                              y_col,
                              species_col,
                              incline = 1,
                              maxdist = Inf,
                              lookup = NULL) {
  # 0) Lookup species -> plant_type (A/G)
  if (is.null(lookup)) {
    lookup <- data.frame(
      species    = c("ACPL","ACSA","BEPA","BEPE","LADE","LALA",
                     "PIAB","PIGL","QURO","QURU","PISY","PIST"),
      plant_type = c("A","A","A","A","G","G",
                     "G","G","A","A","G","G"),
      stringsAsFactors = FALSE
    )
  }
  
  # 1) Rename colums
  df <- data
  colnames(df)[colnames(df) == x_col]       <- "x"
  colnames(df)[colnames(df) == y_col]       <- "y"
  colnames(df)[colnames(df) == species_col] <- "species"
  
  # 2) Adding plant type
  df <- merge(df, lookup, by = "species", all.x = TRUE, sort = FALSE)
  
  if (any(is.na(df$plant_type))) {
    warn_sp <- unique(df$species[is.na(df$plant_type)])
    warning("Kein plant_type-Mapping für: ", paste(warn_sp, collapse = ", "))
    df <- df[!is.na(df$plant_type), ]
  }
  
  n <- nrow(df)
  if (n == 0L) return(df)
  
  # 3) Distance matrix
  coords   <- cbind(df$x, df$y)
  dist_mat <- as.matrix(dist(coords, method = "euclidean"))
  
  # 4) apply maxdist
  if (!is.infinite(maxdist)) {
    dist_mat[dist_mat > maxdist] <- NA
  }
  
  # 5) Spatial weights
  weights_mat <- exp(-dist_mat * incline)
  weights_mat[is.na(weights_mat)] <- 0
  diag(weights_mat) <- 1
  
  # 6) calculation FD_i
  pt_fac <- factor(df$plant_type, levels = c("A", "G"))
  
  FD         <- numeric(n)
  t_density  <- numeric(n)
  
  for (i in seq_len(n)) {
    w_i <- weights_mat[i, ]
    
    # W'_{i,p} = Total weight by plant type
    pt_sums <- tapply(w_i, pt_fac, sum, na.rm = TRUE)
    
    W_tot <- sum(pt_sums, na.rm = TRUE) # W'_i
    t_density[i] <- W_tot
    
    if (W_tot > 0) {
      p_prime <- pt_sums / W_tot       # p'_{i,p}
      p_prime[is.na(p_prime)] <- 0
      H_prime <- -sum(ifelse(p_prime == 0, 0, p_prime * log(p_prime)))
      FD[i] <- exp(H_prime)
    } else {
      FD[i] <- 1
    }
  }
  
  # 7) return
  df$FD        <- FD
  df$t_density <- t_density
  return(df)
}


# ====================================================
# 4. Spatial data preparation
# ====================================================
head(coordinates)
# calculate coordinates in meters
# Grid positions of tree locations within a plot
# Explanation of the exact values in continuous text of master thesis
# multiplicators: step size between trees in meters 
#     (6 * 0.45) describes the width (or length) of the 7 trees per row or column within a plot 
#     (6 spacings between 7 trees at 45 cm), 
#     Plot spacings of 1.8 m
coordinates$x_m <- (6 * 0.45 + 1.8) * (coordinates$Plot_coord_X - 1)   
coordinates$y_m <- (6 * 0.45 + 1.8) * (coordinates$Plot_coord_Y - 1)

species <- array_to_dataframe(trees_array)

species_pos  <- cbind(species,
                      x_m = rep(coordinates$x_m,each = 49) +
                        0.45 * 0:6,  
                      y_m = rep(coordinates$y_m, each = 49)+
                        rep(0.45 * 0:6, each = 7)) %>% 
  # species is supplemented by 2 additional columns x_m and y_m
  # x_m contains spatial X coordinates in meters
  # Repeats the x&y coordinates of a plot for each of the 49 trees within the plot (each = 49)
  # Positioning takes into account the arrangement of trees in a 7x7 grid per plot
  rename("species"= value) %>% 
  rename("plot" = dim) 

# add the column "block" to species_pos
species_pos <- species_pos %>%
  mutate(block = rep(coordinates$Block, each = 49))  # Block für jeden der 49 Bäume wiederholen

head(species_pos)     #fast check
print(species_pos)    #long check :D

# ====================================================
# 5a. join data (prod und infest)
# ====================================================

# List of tree species by abbreviations
species_list <- c("ACPL","LADE","ACSA","LALA","BEPE","PIAB",
                  "BEPA","PIGL","QURO","PISY","QURU","PIST")

# infestation
df_infest <- infest %>% 
  mutate(infestation = case_when(
    infestation_age %in% c("f", "a") ~ "1",   # if infestation_age "f" or "a", then 1
    infestation_age == "n" ~ "0",             # if infestation_age "n", then 0
    infestation_age == "na" ~ "na")) %>% 
  filter(species %in% species_list)

# productivity 
df_prod <- inventory %>% 
  filter(species %in% species_list) %>% 
  filter(d1_mm < 200) # prevent incorrect entries

# Insight into both data sets
head(df_infest) 
nrow(df_infest)
head(df_prod)
nrow(df_prod)

# Merging the two data sets
# Merge all rows here, even if there is no match in the other table.
df_all <- full_join(df_prod, df_infest, by = c("block", "plot", "nrow", "ncol"))
# If values are missing in one of the tables, NA is inserted.
head(df_all)
print(df_all)
nrow(df_all)



# ===============================================================================
# 5b. Merging coordinates and df_all
# ===============================================================================

df_all$block <- as.numeric(df_all$block)
# block 3 is divided in "3" and "3b". I consider both to be block "3"
# Remove “b” from the values in the block column
species_pos$block <- gsub("b", "", species_pos$block)
species_pos$block <- as.numeric(species_pos$block)
unique(species_pos$block)
# Merge df_all and species_pos based on plot, ncol, nrow 
df_all <- inner_join(df_all, species_pos, by = c("plot", "ncol", "nrow", "block"))
head(df_all)


# ===============================================================================
# 6. Complete entries from Mixture 
# ===============================================================================

# Create column "Zusammensetzung"
df_all <- df_all %>%
  group_by(plot) %>%  # Gruppieren nach Plot
  mutate(Zusammensetzung = paste(sort(unique(species.x)), collapse = "+")) %>%
  ungroup()  # Gruppierung wieder aufheben

# Show result
df_all <- as.data.frame(df_all)
# Column "composition"Zusammensetzung" does not match the original "mixture" in all rows. 
# However, my calculation should be accurate.

df_all <- subset(df_all, select = -mixture)           # outdated one
df_all <- subset(df_all, select = -species.y)         # we get "species" later through dataset "species_pos"
df_all <- subset(df_all, select = -species.x)
df_all <- subset(df_all, select = -comments)          # bwe don't need it
df_all <- subset(df_all, select = -infestation_age)   # already covered by infestation
head(df_all)
colnames(df_all) <- c("block", "plot", "nrow", "ncol", "d1_mm", "h1_mm", "dead", "infestation", "species", "x", "y", "Zusammensetzung")





# ===============================================================================
# 7. Use MICE to try to import the missing data from df_all
# ===============================================================================

# Overview: how many NA do I have per column
missing_summary <- colSums(is.na(df_all))
print(missing_summary)


# =======================
# adjust data types for MICE
# =======================

str(df_all)

# convert chr variables to factor
df_all$dead <- as.factor(df_all$dead)
df_all$species <- as.factor(df_all$species)
df_all$Zusammensetzung <- as.factor(df_all$Zusammensetzung)
df_all$infestation <- as.numeric(df_all$infestation)
# Check, if converting was successfull
str(df_all)

# Levels of the factors - fits so far, don't need to remove anything unnecessary
levels(df_all$dead)
levels(df_all$species)
levels(df_all$Zusammensetzung)
unique(df_all$infestation)

# ===================
# Start imputation
# ===================
# Methods: pmm for num, polyreg for factor
imputed_data <- mice(df_all, 
                     method = c("pmm", "pmm", "pmm", "pmm", "pmm", "pmm", "logreg", "pmm", "polyreg", "pmm", "pmm", "polyreg"), 
                     m = 5,       # Number of iterations
                     maxit = 5,  # 5 iterations max
                     seed = 456)  # For reproducibility

# Summary of the result
summary(imputed_data)
print(imputed_data)
completed_data <- complete(imputed_data)
str(completed_data)
head(completed_data)

# Check: How many NA do I have in each column?
missing_summary_2 <- colSums(is.na(completed_data))
print(missing_summary_2)

# ==============
# 8. Check: How good was imputation? Visualization
# ==============

# Comparison of the distribution (original and imputed values), histogram of d1_mm
ggplot(completed_data, aes(x = d1_mm)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5) +
  geom_histogram(data = df_all, aes(x = d1_mm), binwidth = 1, fill = "red", alpha = 0.5) +
  labs(title = "Distributions of original and imputed values", x = "d1_mm", y = "Frequency")
# we see that the data is distributed very similarly -> “good” choice of imputed values

# Comparison of the distribution of original and imputed values in the histogram of h1_mm
ggplot(completed_data, aes(x = h1_mm)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5) +
  geom_histogram(data = df_all, aes(x = h1_mm), binwidth = 1, fill = "red", alpha = 0.5) +
  labs(title = "Distributions of original and imputed values", x = "h1_mm", y = "Frequency") +
  coord_cartesian(xlim = c(0, 1000))
# 

# Comparison of the frequencies of categories in original and imputed data for the parameter “dead”
# not very meaningful
ggplot() +
  geom_bar(data = df_all, aes(x = dead, fill = "Original"), position = "dodge", alpha = 0.5) +
  geom_bar(data = completed_data, aes(x = dead, fill = "Imputiert"), position = "dodge", alpha = 0.5) +
  labs(title = "Frequencies between original and imputed data",
       x = "Dead", y = "Häufigkeit") +
  scale_fill_manual(name = "Datensatz", values = c("Original" = "red", "Imputiert" = "blue"))

# Comparison of the distribution of original and imputed values in the histogram of infestation
ggplot(completed_data, aes(x = infestation)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5) +
  geom_histogram(data = df_all, aes(x = infestation), binwidth = 1, fill = "red", alpha = 0.5) +
  labs(title = "Distributions of original and imputed values", x = "infestation", y = "Häufigkeit")
# joa... We see that the data is distributed similarly -> “good” choice of imputed values


# create productivity column based on tree trunk volume formula
completed_data$productivity <- (pi / 4) * (completed_data$d1_mm^2) * completed_data$h1_mm * 0.5
# unit converted from mm to m
# for productivity, then m^3 as a unit
completed_data$productivity <- completed_data$productivity / 1e9
head(completed_data)

max(completed_data$productivity)


# ====================================================================
# 9. Create different Datasets
# ====================================================================

# Check that there are no outliers! No tree is taller than 2 m, which is realistic.
max(completed_data$h1_mm)

# =======================
# dataset CONIFERS
# for looking at infestation
# =======================
list_conifers <- c("PIST", "PISY", "PIGL", "PIAB", "LALA", "LADE")
# only conifers
df_conifers <- df_all[df_all$species %in% list_conifers, ]
head(df_conifers)
# species has 10 levels ... ?
# wondered first. 
# However, this appears to be because levels are retained in filters and entries are filtered but still displayed in levels
str(df_conifers)
unique(df_conifers$species)
# makes sense, that 3906 rows remain 
nrow(df_conifers)

# lösche irrelevante spalten 
#df_conifers <- df_conifers %>% select(-d1_mm, -h1_mm)

missing_summary_conifers <- colSums(is.na(df_conifers))
print(missing_summary_conifers)
# see that 1/9 of infestation and dead data is missing
# IMPUTATION
imputed_data_conifers <- mice(df_conifers, 
                              method = c("pmm", "pmm", "pmm", "pmm", "pmm", "pmm", "logreg", "pmm", "polyreg", "pmm", "pmm", "polyreg"), 
                              m = 5,       # Number of the iterations
                              maxit = 5,   # 5 iterations
                              seed = 456)  # For reproducibility
# Summary of the result
summary(imputed_data_conifers)
print(imputed_data_conifers)
completed_data_conifers <- complete(imputed_data_conifers)
str(completed_data_conifers)
head(completed_data_conifers)
# Check: How many NA do I have in each coolumn?
missing_summary_conifers_2 <- colSums(is.na(completed_data_conifers))
print(missing_summary_conifers_2)
# Renaming
df_conifers <- completed_data_conifers
head(df_conifers)
# create productivity column
df_conifers$productivity <- (pi / 4) * (df_conifers$d1_mm^2) * df_conifers$h1_mm * 0.5
# set unit to m
# productivity has uni m^3
df_conifers$productivity <- df_conifers$productivity / 1e9
head(df_conifers)
# =======================
# dataset CORES
# for looking at prod. parameters 
# =======================
# consider only cores of the plots
df_cores <- df_all %>% filter(nrow > 1, nrow < 7, ncol > 1, ncol < 7)
head(df_cores)
nrow(df_cores)
missing_summary_cores <- colSums(is.na(df_cores))
print(missing_summary_cores)
# delete irrelevant columns
str(df_cores)
# IMPUTATION
imputed_data_cores <- mice(df_cores, 
                           method = c("pmm", "pmm", "pmm", "pmm", "pmm", "pmm", "logreg", "pmm", "polyreg", "pmm", "pmm", "polyreg"), 
                           m = 5,       # Number of the iterations
                           maxit = 5,   # 5 iterations
                           seed = 456)  # For reproducibility
# Summary of the result
summary(imputed_data_cores)
print(imputed_data_cores)
completed_data_cores <- complete(imputed_data_cores)
str(completed_data_cores)
head(completed_data_cores)
# Check: How many NA do I have in each column?
missing_summary_cores_2 <- colSums(is.na(completed_data_cores))
print(missing_summary_cores_2)
# Renaming
df_cores <- completed_data_cores
head(df_cores)
# Productivity
df_cores$productivity <- (pi / 4) * (df_cores$d1_mm^2) * df_cores$h1_mm * 0.5
# Set unit to m
# Unit for productivity is m^3
df_cores$productivity <- df_cores$productivity / 1e9
# log und standardized prod. 
df_cores$productivity <- as.numeric(scale(log(df_cores$productivity)))

head(df_cores)

# Create historgramm
# hier noch mal die achsen bearbeiten, dass man alles sehen kann
x <- na.omit(df_cores$productivity)
# histogramm of the data
hist(x, breaks = 30, col = "gray85", border = "white",
     main = "Histogram: productivity (log + scaled)",
     xlab = "standardized log(productivity)",
     freq = FALSE)  # Important: Histogram as density, not absolute frequency
# Calculate mean and standard deviation
mu <- mean(x)
sigma <- sd(x)
# Sequence of values across the range of data
xx <- seq(min(x), max(x), length.out = 200)
# Calculate normal density
yy <- dnorm(xx, mean = mu, sd = sigma)
# Plot a normal distribution curve
lines(xx, yy, col = "blue", lwd = 2)
# Mark mean
abline(v = mu, col = "red", lwd = 2, lty = 2)

# =======================
# dataset only cores + only conifers
# =======================
#df_cores_conifers <- inner_join(df_prod, df_infest, by = c("block", "plot", "nrow", "ncol"))
df_cores_conifers <- df_all %>% filter(nrow > 1, nrow < 7, ncol > 1, ncol < 7) %>% filter(species %in% list_conifers)
# delete all entries with NA
missing_summary_cores_conifers <- colSums(is.na(df_cores_conifers))
print(missing_summary_cores_conifers)
df_cores_conifers <- df_cores_conifers %>% drop_na()
nrow(df_cores_conifers)
head(df_cores_conifers)
df_cores_conifers$productivity <- (pi / 4) * (df_cores_conifers$d1_mm^2) * df_cores_conifers$h1_mm * 0.5
# Set unit to m
# Productivity unit m^3
df_cores_conifers$productivity <- df_cores_conifers$productivity / 1e9

# log und standardized prod. 
df_cores_conifers$productivity <- as.numeric(scale(log(df_cores_conifers$productivity)))

head(df_cores_conifers)

max(df_cores_conifers$productivity)
mean(df_cores_conifers$productivity, na.rm = TRUE)
sd(df_cores_conifers$productivity, na.rm = TRUE)

# Create historgramm
x <- na.omit(df_cores_conifers$productivity)
# histogramm of the data
hist(x, breaks = 30, col = "gray85", border = "white",
     main = "Histogram: productivity (log + scaled)",
     xlab = "standardized log(productivity)",
     freq = FALSE)  # Important: Histogram as density, not absolute frequency
# Calculate mean and standard deviation
mu <- mean(x)
sigma <- sd(x)
# Sequence of values across the range of data
xx <- seq(min(x), max(x), length.out = 200)
# Calculate normal density
yy <- dnorm(xx, mean = mu, sd = sigma)
# Plot a normal distribution curve
lines(xx, yy, col = "blue", lwd = 2)
# Mark mean
abline(v = mu, col = "red", lwd = 2, lty = 2)




##########################################
# 10. MLE function definition with weighted Shannon-Index (WS) via CMA-ES
##########################################


####### optional : gewichte für likelihoods berechnen 
suggest_eta_shannon_min <- function(p, data) {
  beta0 <- p[1]; beta1 <- p[2]
  alpha0 <- p[3]; alpha1 <- p[4]
  sigma  <- max(exp(p[5]), 1e-6)
  incline <- p[6]; maxdist <- p[7]
  
  d <- weighted_shannon_index(
    data, "x", "y", "species",
    incline = incline, maxdist = maxdist
  )
  
  if (!all(c("Shannon","infestation","productivity") %in% names(d))) return(0.5)
  
  yi <- d$infestation
  if (is.logical(yi)) yi <- as.integer(yi)
  if (is.factor(yi) || is.character(yi)) suppressWarnings(yi <- as.numeric(as.character(yi)))
  
  prod <- d$productivity
  Shan <- d$Shannon
  
  keep <- is.finite(yi) & is.finite(prod) & is.finite(Shan)
  if (!any(keep)) return(0.5)
  
  yi <- yi[keep]; prod <- prod[keep]; Shan <- Shan[keep]
  
  ll_inf  <- dbinom(yi, size = 1, prob = plogis(beta0 + beta1 * Shan), log = TRUE)
  ll_prod <- dnorm(prod, mean = alpha0 + alpha1 * Shan, sd = sigma, log = TRUE)
  
  s_inf <- sd(ll_inf); s_prod <- sd(ll_prod)
  if (!is.finite(s_inf) || !is.finite(s_prod) || (s_inf + s_prod) == 0) return(0.5)
  
  as.numeric(s_inf / (s_inf + s_prod))
}
eta_conifers <- suggest_eta_shannon_min(init_pars, df_cores_conifers)
print(eta_conifers)

########gewichte berechnen ende


## ---- Negative Composite Log-Likelihood
neg_loglik <- function(p, data) {
  # p = (beta0, beta1, alpha0, alpha1, log_sigma_p, incline, maxdist)
  beta0    <- p[1]
  beta1    <- p[2]
  alpha0   <- p[3]
  alpha1   <- p[4]
  sigma_p  <- exp(p[5])
  incline  <- p[6]
  maxdist  <- p[7]
  
  # calculate diversity index (dependent on incline, maxdist)
  sh <- weighted_shannon_index(
    data,
    x_col = "x",
    y_col = "y",
    species_col = "species",
    incline = incline,
    maxdist = maxdist
  )
  S <- sh$Shannon
  
  # Target variables
  y_infest     <- as.numeric(data$infestation)  # 0/1
  productivity <- data$productivity
  
  # --- Part 1: Infestation (Logit) ---
  eta   <- beta0 + beta1 * S
  # Numerically stable: plogis + dbinom(log=TRUE)
  q_inf <- plogis(eta)
  loglik_infest <- sum(dbinom(y_infest, size = 1, prob = q_inf, log = TRUE))
  
  # --- Part 2: Productivity (Normal) ---
  mu_p  <- alpha0 + alpha1 * S
  # Numerically stable: dnorm(log=TRUE)
  loglik_p <- sum(dnorm(productivity, mean = mu_p, sd = sigma_p, log = TRUE))
  
  # Sum of the partial Loglikelihoods (weight = 1 for both parts)
  -(0.1007358*loglik_infest + (1-0.1007358)*loglik_p)  # negative, because of minimization
}

## ---- CMA-ES Setup ----------------------------------------------------------


# (2) Wrapper
neg_loglik_cma <- function(p) neg_loglik(p, data = df_cores_conifers)

# (3) Starting values + Box-Constraints
init_pars <- c(
  beta0       = 0,
  beta1       = 0.1,
  alpha0      = mean(df_cores_conifers$productivity, na.rm = TRUE),
  alpha1      = 0.1,
  log_sigma_p = log(sd(df_cores_conifers$productivity, na.rm = TRUE)),
  incline     = 1,     # in [0, 7]
  maxdist     = 10     # in (0, 30], 0 leads to only-own tree
)


lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, -Inf,   0,  0.1)  # maxdist > 0 makes sense
upper_bounds <- c( Inf,  Inf,  Inf,  Inf,  Inf,   7, 30)

# (4) CMA-ES Call
set.seed(1)
res_cma <- cma_es(
  fn    = neg_loglik_cma,
  par   = init_pars,
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(
    sigma  = 1.0,
    lambda = 20,
    maxit  = 150,
    # optional: Trace-Output
    verbose = TRUE                    # Controls how much progress output CMA-ES writes to the console.
  )
)

# (5) Result
cat("Found parameter set:\n")
print(res_cma$par)
cat("neg. Log-Likelihood:\n")
print(res_cma$value)






##########################################
# 11. MLE function definition with weighted functional Shannon-Index (WFSI) via CMA-ES
# Current status: läuft nicht gut, da incline parameter nach oben crasht
##########################################




# ===================================================================================
# GEWICHT FÜR LIKELIHOODS BERECHNEN - START
# ===================================================================================
suggest_eta <- function(p, data) {
  beta0 <- p[1]; beta1 <- p[2]
  alpha0 <- p[3]; alpha1 <- p[4]
  sigma  <- exp(p[5])
  incline <- p[6]; maxdist <- p[7]
  
  fd <- weighted_fd_index(data, "x","y","species", incline=incline, maxdist=maxdist)
  
  FD  <- fd$FD
  yi  <- fd$infestation
  if (is.logical(yi)) yi <- as.integer(yi)
  if (is.factor(yi) || is.character(yi)) yi <- as.numeric(as.character(yi))
  prod <- fd$productivity
  
  keep <- is.finite(FD) & is.finite(yi) & is.finite(prod)
  if (!any(keep)) return(0.5)  # Fallback
  
  FD <- FD[keep]; yi <- yi[keep]; prod <- prod[keep]
  
  # per-Observation Log-Likelihoods:
  ll_inf_i  <- dbinom(yi, size=1, prob=plogis(beta0 + beta1*FD), log=TRUE)
  ll_prod_i <- dnorm(prod, mean = alpha0 + alpha1*FD, sd = pmax(sigma, 1e-6), log=TRUE)
  
  s_inf  <- sd(ll_inf_i)
  s_prod <- sd(ll_prod_i)
  
  if (!is.finite(s_inf) || !is.finite(s_prod) || (s_inf + s_prod) == 0) return(0.5)
  # Heuristik: beide Teile gleich stark "schwingen" lassen
  eta <- s_inf / (s_inf + s_prod)
  as.numeric(eta)
}


eta0 <- suggest_eta(init_pars, df_cores)
print(eta0)
# ===================================================================================
# GEWICHT FÜR LIKELIHOODS BERECHNEN - ENDE
# ===================================================================================



## ---- Negative Composite Log-Likelihood
neg_loglik <- function(p, data) {
  # p = (beta0, beta1, alpha0, alpha1, log_sigma_p, incline, maxdist)
  beta0   <- p[1]
  beta1   <- p[2]
  alpha0  <- p[3]
  alpha1  <- p[4]
  sigma_p <- max(exp(p[5]), 0.25)   # forces sigma >= 0.25
  incline <- p[6]
  maxdist <- p[7]
  
  # calculate WFS
  fd <- weighted_fd_index(
    data,
    x_col = "x", y_col = "y", species_col = "species",
    incline = incline, maxdist = maxdist
  )
  
  #
  FD <- fd$FD
  
  # infestation safely in 0/1
  yi <- fd$infestation
  if (is.factor(yi) || is.character(yi)) yi <- as.numeric(as.character(yi))
  if (is.logical(yi)) yi <- as.integer(yi)  # TRUE/FALSE -> 1/0
  
  prod <- fd$productivity
  
  # complete cases
  keep <- is.finite(FD) & is.finite(yi) & is.finite(prod)
  if (!any(keep)) return(Inf)
  
  # target variables
  FD   <- FD[keep]
  yi   <- yi[keep]
  prod <- prod[keep]
  
  # --- Part 1: Infestation logit ---
  eta  <- beta0 + beta1 * FD
  # Numerically stable: plogis + dbinom(log=TRUE)
  loglik_infest <- sum(dbinom(yi, size = 1, prob = plogis(eta), log = TRUE))
  
  # --- Part 2: Productivity ---
  mu   <- alpha0 + alpha1 * FD
  # Numerically stable: dnorm(log=TRUE)
  loglik_prod <- sum(dnorm(prod, mean = mu, sd = sigma_p, log = TRUE))
  
  # Sum of the partial-Loglikelihoods (weight = 1 for both parts)
  -(0.07251564*loglik_infest + (1-0.07251564)*loglik_prod) # negative, because of minimization
}


## ---- CMA-ES Setup ----------------------------------------------------------


# (2) Wrapper
neg_loglik_cma <- function(p) neg_loglik(p, data = df_cores)

# (3) Starting Values + Box-Constraints
init_pars <- c(
  beta0       = 0,
  beta1       = 0.1,
  alpha0      = mean(df_cores$productivity, na.rm = TRUE),
  alpha1      = 0.1,
  log_sigma_p = log(sd(df_cores$productivity, na.rm = TRUE)),
  incline     = 1,     # in [0, 7]
  maxdist     = 10     # in (0, 30], 0 leads to only-own tree
)


lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, -Inf,   0,  0.1)  # maxdist > 0 makes sense
upper_bounds <- c( Inf,  Inf,  Inf,  Inf,  Inf,   7, 30)   # checkt upper bounds for incline

# (4) CMA-ES Call
set.seed(1)
res_cma <- cma_es(
  fn    = neg_loglik_cma,
  par   = init_pars,
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(
    sigma  = 1.0,
    lambda = 20,
    maxit  = 150,
    # optional: Trace-Output
    verbose = TRUE                    # Controls how much progress output CMA-ES writes to the console
  )
)

# (5) Ergebnis
cat("Found parameterset:\n")
print(res_cma$par)
cat("neg. Log-Likelihood:\n")
print(res_cma$value)



# ================================
# TEST ANFANG
# dieser test soll mir die plausibilität der ermittelten Parameterwerte aus dem MLE FD prüfen
# dabei ist: optimal wäre eine trichterförmige fläche, die darauf hinweisen würde, dass es ein klar identifizierbares minimum gibt.
# hier ist mein minimum abaer am rande des such-grids - was darauf schließen lässt, dass der suchbereich zu eng ist. 
# roter punkt sind meine schätzwerte, die liegen hier außerhalb der fläche - schlecht
# ================================


# --- Eingaben: dein bestes Ergebnis ---
p_hat <- res_cma$par           # c(beta0,beta1,alpha0,alpha1,log_sigma, incline, maxdist)
inc0  <- p_hat[6]; dmax0 <- p_hat[7]

# --- Grobes Grid um das Optimum (anpassen, falls nötig) ---
incs  <- seq(inc0 - 2, inc0 + 2, length.out = 15)
dmaxs  <- seq(dmax0 - 6, dmax0 + 6, length.out = 15)

grid   <- expand.grid(incline = incs, maxdist = dmaxs)
f_eval <- function(inc, dmax){
  p <- p_hat
  p[6] <- inc
  p[7] <- dmax
  neg_loglik(p, data = df_cores)   # deine Funktion
}

# --- Likelihood-Landschaft berechnen (schnell, keine Optimierung) ---
grid$nll <- mapply(f_eval, grid$incline, grid$maxdist)

# --- Interpretation: Delta zur besten Stelle im Grid ---
best_nll <- min(grid$nll, na.rm = TRUE)
grid$delta <- grid$nll - best_nll

# 1) Schneller Eindruck (kleine Tabelle):
subset(grid[order(grid$delta), ], delta <= 2)[1:10, ]

# 2) Konturplot (zeigt Trichter vs. Ridge):

ggplot(grid, aes(incline, maxdist, z = delta)) +
  geom_contour_filled(breaks = c(0, 0.5, 1, 2, 5, 10)) +
  geom_point(aes(x = inc0, y = dmax0), color = "red", size = 3) +
  labs(title = "2D-Profil der neg. Log-Likelihood (Δ zum Grid-Minimum)",
       x = "incline", y = "maxdist") +
  theme_minimal()


# ================================
# TEST ENDE
# ================================













# ==============================================================================
# 12. Calculate Diversity indices
# ==============================================================================

# ================================
# 12a. normal Shannon-Index
# ================================

shannon_result_completed_data <- calculate_shannon(completed_data)
print(shannon_result_completed_data)

shannon_result_df_conifers <- calculate_shannon(df_conifers)
print(shannon_result_df_conifers)

shannon_result_df_cores <- calculate_shannon(df_cores)
print(shannon_result_df_cores)

shannon_result_df_cores_conifers <- calculate_shannon(df_cores_conifers)
print(shannon_result_df_cores_conifers)


# ===============
# 12b. Shannon-Index across all data
# ===============
# Calculate Shannon-Index for completed_data
shannon_overall_completed_data <- calculate_shannon_overall(
  data = completed_data,   
  species_col = "species"  # Column with species
)
cat("Global Shannon-Index for the entire data set:", shannon_overall_completed_data, "\n")



# ===============
# 12c. weighted shannon index (WSI)
# ===============

# weighted Shannon-Index for df_cores_conifers
# with values for maxdist and incline of MLE
df_cores_conifers_with_shannon <- weighted_shannon_index(
  data = df_cores_conifers,
  x_col = "x",
  y_col = "y",
  species_col = "species",
  incline = 0.200464894,
  maxdist = 11.030479963
)
head(df_cores_conifers_with_shannon)



# Plot of weighted Shannon-Index
x <- na.omit(df_cores_conifers_with_shannon$Shannon)
# Histogramm der Daten
hist(x, breaks = 30, col = "gray85", border = "white",
     main = "weighted Shannon-Index",
     xlab = "Shannon-Index",
     freq = FALSE)  # Important: Histogram as density, not absolute frequency
# calculate mean and standard deviation
mu <- mean(x)
sigma <- sd(x)
# Sequence of values across the range of data
xx <- seq(min(x), max(x), length.out = 200)
# Calculate normal density
yy <- dnorm(xx, mean = mu, sd = sigma)
# Plot a normal distribution curve
lines(xx, yy, col = "blue", lwd = 2)
# Mark mean value
abline(v = mu, col = "red", lwd = 2, lty = 2)











# ===============
# 12d. weighted functional Shannon-Index (WFSI)
# ===============



# weighted functional Shannon-Index for df_cores
# With values for incline and maxdist from MLE 
df_cores_with_fd <- weighted_fd_index(
  data = df_cores,
  x_col = "x",                      # Column with X-coordinates
  y_col = "y",                      # Column with Y-coordinates
  species_col = "species",          # Column with tree species
  incline = 1.250407789,             # Incline
  maxdist = 10.034484631             # Maximal distance in meters
)
head(df_cores_with_fd)











# Plot of weighted functional Shannon-Index
x <- na.omit(df_cores_with_fd$FD)
# Histogramm der Daten
hist(x, breaks = 30, col = "gray85", border = "white",
     main = "weighted functional Shannon-Index",
     xlab = "weighted functional Shannon-Index",
     freq = FALSE)  # Important: Histogram as density, not absolute frequency
# calculate mean and standard deviation
mu <- mean(x)
sigma <- sd(x)
# Sequence of values across the range of data
xx <- seq(min(x), max(x), length.out = 200)
# Calculate normal density
yy <- dnorm(xx, mean = mu, sd = sigma)
# Plot a normal distribution curve
lines(xx, yy, col = "blue", lwd = 2)
# Mark mean value
abline(v = mu, col = "red", lwd = 2, lty = 2)




# =========================================
# Plot incline parameter WSI
# =========================================
incline <- 0.200464894
maxdist <- 11.030479963
r50 <- log(2) / incline      #
print(r50)
df <- data.frame(
  r = seq(0, maxdist + 2, by = 0.05)
)
df$w <- ifelse(df$r <= maxdist, exp(-incline * df$r), 0)

ggplot(df, aes(r, w)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = maxdist, linetype = "dashed", color = "red") +
  geom_vline(xintercept = r50,     linetype = "dotted", color = "blue") +
  labs(
    title = "Decrease in neighborhood weights",
    subtitle = expression(w(r) == e^{-I %.% r} ~ ",  incline = 0.280416243,  maxdist = 10.471021439m"),
    x = "distance r (m)",
    y = "weight w(r)"
  )




# =========================================
# Plot incline parameter WFSI
# =========================================
incline <- 1.250407789
maxdist <-  10.034484631
r50 <- log(2) / incline      # ~3.47 m
print(r50)
df <- data.frame(
  r = seq(0, maxdist + 2, by = 0.05)
)
df$w <- ifelse(df$r <= maxdist, exp(-incline * df$r), 0)

ggplot(df, aes(r, w)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = maxdist, linetype = "dashed", color = "red") +
  geom_vline(xintercept = r50,     linetype = "dotted", color = "blue") +
  labs(
    title = "Decrease in neighborhood weights",
    subtitle = expression(w(r) == e^{-I %.% r} ~ ",  incline = 1.250407789,  maxdist = 10.034484631m"),
    x = "distance r (m)",
    y = "weight w(r)"
  )







# ===================================================================
# 13. Connection: weighted Shannon-Index and parameters (infest, prod)
# ===================================================================

# log regression for infestation
log_model <- glm(infestation ~ Shannon, data = df_cores_conifers_with_shannon, family = binomial)
summary(log_model)
# Shannon coefficient is negative -> The higher the Shannon index, the lower the risk of infestation.
# !!!!!!!!! I want to draw further analyses from here !!!!!


# =====================
# Visualization
# =====================
# =====================
# 13.1 Connection Shannon and infestation
# =====================

# log. reg. infestation ~ Shannon
# with function stat_smooth()
# first scatter plot of the original data points, then overlay with regression curve (using statt_smooth())

ggplot(df_cores_conifers_with_shannon, aes(x = Shannon, y = infestation)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.25, color="blue") +  # 0/1-Punkte sichtbar machen
  stat_smooth(method = "glm",
              method.args = list(family = binomial(link = "logit")),
              se = TRUE, size = 1.1, color = "red") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(
    title = "Probability of infestation vs. Shannon index",
    x = "Shannon-Index",
    y = "Probability (Infestation = 1)"
  ) +
  theme_minimal()





# =====================
# 13.2 Connection productivity and Shannon-Index
# =====================


#check der productivity

hist(df_cores_conifers_with_shannon$productivity,
     main = "Histogramm of the Productivity",
     xlab = "Productivity (log-transformed, standardized)",
     ylab = "Frequency",
     col  = "lightgray",
     border = "white"
     )
plot(df_cores_conifers_with_shannon$Shannon, df_cores_conifers_with_shannon$productivity)
plot(df_cores_conifers_with_shannon$Shannon,
     df_cores_conifers_with_shannon$productivity,
     xlab = "Shannon-Index",
     ylab = "Productivity (log-transformed, standardized)",
     main = "Zusammenhang zwischen Diversität und Produktivität",
     pch  = 19,        # dots, filled
     col  = rgb(0,0,0,0.2))  # transparent dots


# hier noch mal ran... noch nicht fertig
# Linear Mixed Model for Prod
library(lme4)
m1 <- lmer(productivity ~ Shannon + (1|species) + (1|plot), 
           data = df_cores_conifers_with_shannon)
summary(m1)










# I don't see any connection between productivity and Shannon-Index
ggplot(df_cores_conifers_with_shannon, aes(x = Shannon, y = productivity)) +
  geom_point(alpha = 0.5, color = "blue") +  # plot dots
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # regression line with confidence interval
  labs(
    title = "Relationship: Shannon-Index and Productivity",
    x = "Shannon-Index",
    y = "Productivity"
  ) +
  theme_minimal()


# Correlation between productivity and Shannon-Index for each species
# In some species, the index has a greater influence on productivity than in others
ggplot(df_cores_conifers_with_shannon, aes(x = Shannon, y = productivity)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  facet_wrap(~species, scales = "free") +  # Erstellt separate Plots für jede Art
  labs(
    title = "Relationship: Shannon-Index and Productivity by tree species",
    x = "Shannon-Index",
    y = "Productivity"
  ) +
  theme_minimal()


#läuft noch nicht so doll, hier könnte man noch mal dran arbeiten
# Zusammenhang produktivität und shannon in nichtlinearer regression
ggplot(df_cores_conifers_with_shannon, aes(x = Shannon, y = productivity)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_smooth(
    method = "nls",
    formula = y ~ SSlogis(x, Asym, xmid, scal),
    method.args = list(
      algorithm = "port",                          
      lower = c(Asym = 0, xmid = -Inf, scal = 1e-6), # lower bounds
      control = nls.control(
        maxiter = 500,
        minFactor = 1e-12,    # allows little step sizes
        warnOnly = TRUE       # gives Warning instead of error
      )
    ),
    se = FALSE, color = "red"
  ) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Sigmoidal fit of productivity vs. Shannon by species",
       x = "Shannon index", y = "Productivity") +
  theme_minimal()


# =============
# 13.3 Spatial visualization of the weighted Shannon-Index in the coordinate system -> heatmap
# =============

ggplot(df_cores_conifers_with_shannon, aes(x = x, y = y, color = Shannon)) +
  geom_point(size = 4, alpha = 0.8) +  # Plotting points with color scale
  scale_color_viridis_c(option = "plasma") +  # Color scale for Shannon-Index
  labs(title = "Shannon-Index heatmap of the tree positions",
       x = "X-coordinate", y = "Y-coordinate", color = "Shannon-Index") +
  theme_minimal()




#============================
# spatial visualization of infestation
#============================

# vielleicht hier noch mal machen mit completed data für bessere übersicht
ggplot(df_cores_conifers_with_shannon, aes(x = x, y = y, color = infestation)) +
  geom_point(size = 4, alpha = 0.8) +                     # Plot points with color scale
  scale_color_viridis_c(option = "plasma") +              # color scale for Shannon-Index
  labs(title = "Infestation Heatmap of the tree positions",
       x = "X-coordinate", y = "Y-coordinate", color = "Infestation") +
  theme_minimal()












#============================
# spatial visualization of produktivity
#============================

# auch hier noch mal machen mit completed data für bessere übersicht
ggplot(df_cores_conifers_with_shannon, aes(x = x, y = y, color = productivity)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Productivity Heatmap (log-skaliert)",
       x = "X-Koordinate (m)", y = "Y-Koordinate (m)", color = "log&stand. Prod") +
  theme_minimal()












# ===================================================================
# 14. Connection: weighted functional Shannon-Index and parameters (infest, prod)
# ===================================================================

# log regression for infestation
log_model <- glm(infestation ~ FD, data = df_cores_with_fd, family = binomial)
summary(log_model)

#infestation and WFS-Index

ggplot(df_cores_with_fd, aes(x = FD, y = infestation)) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.25, color="blue") +  # make 0/1-dots visible
  stat_smooth(method = "glm",
              method.args = list(family = binomial(link = "logit")),
              se = TRUE, size = 1.1, color = "red") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(
    title = "Probability of infestation vs. FD index",
    x = "FD-Index",
    y = "Probability (Infestation = 1)"
  ) +
  theme_minimal()





# WFS x prod

# statistische untersuchung
cor.test(df_cores_with_fd$FD, df_cores_with_fd$productivity)
cor.test(df_cores_with_fd$FD, df_cores_with_fd$productivity, method = "spearman")
lm_prod <- lm(productivity ~ FD, data = df_cores_with_fd)
summary(lm_prod)
summary(lm(df_cores_with_fd$productivity ~ df_cores_with_fd$FD))$r.squared


ggplot(df_cores_with_fd, aes(x = FD, y = productivity)) +
  geom_point(alpha = 0.5, color = "blue") +               # Plot points
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Regression line with confidence interval
  labs(
    title = "Relationship: WFS-Index and Productivity",
    x = "WFS-Index",
    y = "Productivity"
  ) +
  theme_minimal()


# Connection Productivity and WFS-Index for every species
ggplot(df_cores_with_fd, aes(x = FD, y = productivity)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  facet_wrap(~species, scales = "free") +  # Create different plots for each species
  labs(
    title = "Relationship: WFS-Index and Productivity by tree species",
    x = "WFS-Index",
    y = "Productivity"
  ) +
  theme_minimal()
# ich sehe hier ähnliche trends bei gleichen arten, also dass sich europ. und NA species ähnlich verhalten






# =============
# Spatial visualization of the WFS-Index in the coordinate system -> heatmap
# =============

ggplot(df_cores_with_fd, aes(x = x, y = y, color = FD)) +
  geom_point(size = 3, alpha = 0.8) +                            # Plot points with color scale
  scale_color_viridis_c(option = "plasma") +                     # Color scale for Shannon-Index
  labs(title = "WFS-Index heatmap of the tree positions",
       x = "X-coordinate", y = "Y-coordinate", color = "WFS-Index") +
  theme_minimal()




#============================
# spatial visualization of infestation
#============================

ggplot(df_cores_with_fd, aes(x = x, y = y, color = infestation)) +
  geom_point(size = 2, alpha = 0.8) +                   # Plot points with color scale
  scale_color_viridis_c(option = "plasma") +            # Color Scale for Shannon-Index
  labs(title = "Infestation Heatmap of the tree positions",
       x = "X-coordinate", y = "Y-coordinate", color = "Infestation") +
  theme_minimal()












#============================
# räumliche Visualisierung von Produktivität
#============================

ggplot(df_cores_with_fd, aes(x = x, y = y, color = productivity)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Productivity Heatmap (log-skaliert)",
       x = "X-Koordinate (m)", y = "Y-Koordinate (m)", color = "log&stand. Prod") +
  theme_minimal()





















# TEST CLIC

# ergbnisse von cma es noch mla neu eingeben

res_cma <- list()
res_cma$par <- c(
  1.355252924,  # beta0
  -0.874058698,  # beta1
  0.470231540,  # alpha0
  -0.122897464,  # alpha1
  -0.003324612,  # log_sigma_p
  0.200464894,  # incline
  11.030479963   # maxdist
)


# Pakete
if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
if (!requireNamespace("stats", quietly = TRUE)) install.packages("stats")

library(numDeriv)

# ---- 0) Hooks auf DEINE Objekte ---------------------------------------------
# Annahmen: Folgende Objekte/Funktionen sind bereits in deinem Workspace vorhanden:
# - df_cores: data.frame mit Spalten x,y,species,infestation,productivity
# - weighted_fd_index(): wie oben von dir definiert
# - neg_loglik(p, data): gibt NEGATIVE composite log-likelihood zurück
# - res_cma$par: MCL-Schätzer (Vektor-Reihenfolge: beta0,beta1,alpha0,alpha1,log_sigma_p,incline,maxdist)

theta_hat <- as.numeric(res_cma$par)

# Positiv formulierte CL:
loglik_total <- function(par) {
  -neg_loglik(par, data = df_cores)
}

# ---- 1) l_C(hat) -------------------------------------------------------------
lC_hat <- loglik_total(theta_hat)

# ---- 2) H-Schätzer (Sensitivität) via numerischer Hesse ----------------------
# H = - E[∂^2 l_C / ∂θ∂θ^T], wir approximieren bei hat-θ:
H_hat <- -hessian(func = loglik_total, x = theta_hat)

# Sanfte Regularisierung, falls H numerisch knapp singulär ist:
H_hat <- (H_hat + t(H_hat))/2
eigH  <- eigen(H_hat, symmetric = TRUE)
if (min(eigH$values) < 1e-8) {
  # kleine Ridge, damit invertierbar
  H_hat <- H_hat + diag(1e-6, nrow(H_hat))
}

# ---- 3) J-Schätzer (Variabilität) via Window-Subsampling --------------------
# Idee nach Heagerty & Lumley / Varin & Vidoni:
# Wir bilden K räumliche Fenster (Cluster) und berechnen je Fenster den Score (Gradient) der CL,
# aber nur über Daten dieses Fensters -> ungefähre Unkorreliertheit der Score-Beiträge.
# J ~ Var(Score_window) mit geeigneter Skalierung.

# Hilfsfunktion: Composite-Loglik NUR für einen Index-Vektor (Fenster)
# -> Wir filtern df_cores auf Indizes I und rechnen die CL nur mit diesen Beobachtungen
loglik_window <- function(par, I) {
  # Subset der Daten (Achtung: deine FD-Funktion nutzt globale Nachbarschaften.
  # Für ein sauberes Window-CL ignorieren wir Cross-Fenster-Interaktionen (Standard in Subsampling-Ansätzen))
  d_sub <- df_cores[I, , drop = FALSE]
  -neg_loglik(par, data = d_sub)
}

# Score (Gradient) für ein Fenster
score_window <- function(par, I) {
  grad(func = function(p) loglik_window(p, I), x = par)
}

# Fensterbildung: K-means auf Koordinaten oder Gitter. Hier: K-means (robust und einfach).
set.seed(42)
K <- 25  # Anzahl Fenster (anpassen; 20–50 ist oft ok, je nach N)
coords <- as.matrix(df_cores[, c("x","y")])
km <- kmeans(coords, centers = K, iter.max = 50)
groups <- split(seq_len(nrow(df_cores)), km$cluster)

# Score je Fenster
S_list <- lapply(groups, function(I) score_window(theta_hat, I))
S_mat  <- do.call(cbind, S_list)   # d x K (d = length(theta))

# Schätze J als Kovarianzmatrix der Fenster-Scores (mit nützlicher Skalierung)
# Varianz der Summen ~ Kovarianz der Fensterbeiträge, skaliert mit K
# Wir nehmen die emp. Kovarianz * K (entspricht Sum-of-squares um Mittelwert)
S_centered <- S_mat - rowMeans(S_mat)
J_hat <- (S_centered %*% t(S_centered)) / (K - 1)

# Numerische Symmetrisierung
J_hat <- (J_hat + t(J_hat)) / 2

# ---- 4) CLIC: l_C + tr(J H^{-1}) --------------------------------------------
H_inv <- solve(H_hat)
trace_term <- sum(diag(J_hat %*% H_inv))
CLIC <- lC_hat + trace_term

cat("\n=== Composite Likelihood Information Criterion (CLIC) ===\n")
cat("l_C(hat):        ", sprintf("%.6f", lC_hat), "\n")
cat("trace(J H^{-1}): ", sprintf("%.6f", trace_term), "\n")
cat("CLIC:            ", sprintf("%.6f", CLIC), "\n")

# ---- 5) (Optional) Parametrisches Bootstrap für J ---------------------------
# Gerüst: Nur falls du eine datengenerierende Prozedur hast, um aus deinem
# zusammengesetzten Modell zu simulieren (Infestation via Bernoulli, Productivity via Normal,
# beide konditioniert auf FD(par)). Dann:
#
# simulate_data <- function(par, data_template) {
#   # 1) FD mit current par berechnen (auf Basis der template-Koordinaten/Species)
#   fd <- weighted_fd_index(data_template, "x","y","species",
#                           incline = par[6], maxdist = par[7])
#   FD   <- fd$FD
#   n    <- nrow(fd)
#   # 2) Infestation ~ Bernoulli(plogis(beta0 + beta1*FD))
#   eta  <- par[1] + par[2]*FD
#   p    <- 1/(1+exp(-eta))
#   yi   <- rbinom(n, size=1, prob=p)
#   # 3) Productivity ~ Normal(alpha0 + alpha1*FD, sigma_p)
#   mu   <- par[3] + par[4]*FD
#   sd_p <- max(exp(par[5]), 0.1)
#   prod <- rnorm(n, mean=mu, sd=sd_p)
#   # 4) Zurückschreiben in Kopie
#   d_sim <- data_template
#   d_sim$infestation <- yi
#   d_sim$productivity <- prod
#   d_sim
# }
#
# B <- 200
# S_boot <- matrix(0, nrow=length(theta_hat), ncol=B)
# set.seed(1)
# for (b in 1:B) {
#   d_b <- simulate_data(theta_hat, df_cores)
#   # Score auf GANZE Daten (oder auf Fenster wie oben)
#   S_boot[, b] <- grad(function(p) -neg_loglik(p, data=d_b), x=theta_hat)
# }
# J_hat_boot <- cov(t(S_boot))
# J_hat <- (J_hat_boot + t(J_hat_boot))/2
# CLIC_boot <- lC_hat + sum(diag(J_hat %*% H_inv))
# cat("CLIC (bootstrap-J):", sprintf("%.6f", CLIC_boot), "\n")
