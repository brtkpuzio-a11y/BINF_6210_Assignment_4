#### Introduction  ####


# Assignment 4 - Bartek Puzio 1146247
# Dr. Karl Cottenie BINF*6210
# Project idea #3
# Bombus COI phylogeography: QC, ML tree, and spatial patterns
#
# Purpose:
#   Combine Bombus COI-5P sequences (NCBI) with geographic data (BOLD) to explore
#   how phylogenetic relationships (sister taxa, clades) relate to geography.
#
# Methods:
#   - Clean and summarise BOLD coordinates to check spatial coverage and outliers.
#   - Download, length-filter, and summarise NCBI COI sequences for QC.
#   - Build a GTR ML tree (from aligned sequences), identify sister taxa and
#     patristic-distance clade groups, and link them to BOLD records.
#   - Produce maps (points + 2D density) and an NMDS ordination coloured by clade.
#
# Main outputs:
#   - QC’d Bombus COI-5P sequence + geographic distribution dataset,
#   - GTR ML tree with sister groups and clades,
#   - Figures showing clade/sister-taxon geographic distributions and NMDS structure.

#### Github commands ####

#usethis::git_remotes()
#usethis::git_sitrep()
#usethis::pr_init(branch = "barteks_branch")
#usethis::pr_push()
#usethis::pr_pull()
#usethis::pr_finish()
#usethis::git_sitrep()

#gitcreds::gitcreds_get()
#GITHUB_PAT=

styler::style_file("R.R") # for consistency in styling!
# install and/or load packages

# install.packages("tidyverse")
library(tidyverse)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("maps")
library(maps)
# install.packages("rentrez")
library(rentrez)
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)
# install.packages("dplyr")
library(dplyr)
# install.packages("styler")
library(styler)
# BiocManager::install("DECIPHER")
library(DECIPHER)
# install.packages("ape")
library(ape)
# install.packages("phangorn")
library(phangorn)
# install.packages("sf")
library(sf)
# install.packages("rnaturalearthdata")
library(rnaturalearthdata)
# install.packages("rnaturalearth")
library(rnaturalearth)
# install.packages("cowplot")
library(cowplot)
# install.packages(viridis)
library(viridis)
# install.packages(vegan)
library(vegan)


# ---- BOLD: load raw occurrence data (coords + taxonomy) from TSV export ----


data_api <- read_tsv("../Data/result (1).tsv")

# Split the combined "coord" field into latitude/longitude while keeping the original
# This makes later debugging easier if we need to trace back weird coordinates.
data_api_clean <- data_api %>%
  separate(coord, into = c("lat", "long"), sep = ",", remove = FALSE)

# BOLD sometimes wraps coordinates in square brackets; strip those before parsing
data_api_clean$lat <- gsub("\\[", "", data_api_clean$lat)
data_api_clean$long <- gsub("\\]", "", data_api_clean$long)

# Parse coordinates as numeric so we can summarise and plot them
data_api_clean$lat <- as.numeric(data_api_clean$lat)
data_api_clean$long <- as.numeric(data_api_clean$long)

# Keep only records that are actually usable for spatial and biological analyses:
# require non-missing coordinates, country/ocean, and species ID.
data_api_cleaner <- data_api_clean %>%
  filter(
    !is.na(`lat`),
    !is.na(`long`),
    !is.na(`country/ocean`),
    !is.na(`species`)
  )
rm(data_api_clean) # drop intermediate object

# Summarise coordinate ranges to check that the data cover plausible lat/long values and to report basic spatial coverage in the methods/results.
data_api_cleaner %>%
  summarise(
    n_records = n(),
    lat_min   = min(lat, na.rm = TRUE),
    lat_q1    = quantile(lat, 0.25, na.rm = TRUE),
    lat_med   = median(lat, na.rm = TRUE),
    lat_q3    = quantile(lat, 0.75, na.rm = TRUE),
    lat_max   = max(lat, na.rm = TRUE),
    long_min  = min(long, na.rm = TRUE),
    long_q1   = quantile(long, 0.25, na.rm = TRUE),
    long_med  = median(long, na.rm = TRUE),
    long_q3   = quantile(long, 0.75, na.rm = TRUE),
    long_max  = max(long, na.rm = TRUE)
  )

# Compact summary() calls for quick sanity checks
summary(data_api_cleaner$lat)
summary(data_api_cleaner$long)

# Reshape lat/long into long format so we can plot both variables in one figure
coords_long <- data_api_cleaner %>%
  select(lat, long) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  )

# Histograms of latitude and longitude to visualise spatial coverage and detect any obvious artefacts
ggplot(coords_long, aes(x = value)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  facet_wrap(~variable, scales = "free") +
  labs(
    title = "Distribution of Bombus occurrence coordinates (BOLD)",
    x     = "Value",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Boxplot-based outliers
boxplot.stats(data_api_cleaner$lat)$out
boxplot.stats(data_api_cleaner$long)$out

# Explicitly list clearly implausible coordinates (outside valid lat/long bounds). These are candidates for data entry errors.
data_api_cleaner %>%
  filter(lat < -60 | lat > 80 | long < -180 | long > 180)


# ---- NCBI: download and QC COI-5P sequences for Bombus --------------------


# Query combines Bombus taxon, COI/CO1/cox1 gene names, mitochondrial filter, and COI-5P/barcode terms. Pattern follows examples from the rentrez vignette:
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
query <- '"Bombus"[Organism]
AND (COI[Gene] OR CO1[Gene] OR cox1[Gene])
AND mitochondrion[Filter]
AND ("COI-5P"[All Fields] OR barcode[All Fields])'

# Limit retmax to 300 to keep the dataset tractable for alignment and ML tree building
search_res <- entrez_search(
  db = "nucleotide",
  term = query,
  retmax = 300
)

# Fetch the sequences in FASTA format
seqs_fasta <- entrez_fetch(
  db      = "nucleotide",
  id      = search_res$ids,
  rettype = "fasta",
  retmode = "text"
)

# Save a copy of the raw FASTA
writeLines(seqs_fasta, "../Data/Bombus_COI.fasta")
dna_all <- readDNAStringSet("../Data/Bombus_COI.fasta")

# Length filter: keep sequences in a typical COI-5P barcode range.
# This step removes very short fragments and long records.
keep_len <- width(dna_all) >= 500 & width(dna_all) <= 750
dna_filt <- dna_all[keep_len]
length(dna_all)
length(dna_filt)
writeXStringSet(dna_filt, "../Data/Bombus_COI_500-750bp.fasta")

# All NCBI Bombus COI-5P sequences
len_all <- Biostrings::width(dna_all)

# Filtered set (500–750 bp)
len_filt <- Biostrings::width(dna_filt)

# Summaries to show how the length filter affects the distribution of sequence lengths
summary(len_all)
summary(len_filt)

# Combine lengths into a single tidy data frame for plotting and comparison
len_df <- data.frame(
  length = c(len_all, len_filt),
  dataset = c(
    rep("All NCBI sequences", length(len_all)),
    rep("Filtered 500–750 bp", length(len_filt))
  )
)

# Tidy summary: to report min/median/max lengths for each set
len_df %>%
  group_by(dataset) %>%
  summarise(
    n        = n(),
    min_len  = min(length),
    q1       = quantile(length, 0.25),
    median   = median(length),
    q3       = quantile(length, 0.75),
    max_len  = max(length)
  )

# Boxplots emphasize outliers and show how the filter narrows the length distribution
# This figure doesn't generate sometimes, run it again if a boxplot does not generate. I could not tell you why sorry! It works when i run it again hehe.
ggplot(len_df, aes(x = dataset, y = length)) +
  geom_boxplot(outlier.colour = "red") +
  coord_flip() +
  labs(
    title = "COI-5P sequence lengths: all vs filtered",
    x     = "",
    y     = "Sequence length (bp)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Histograms for all vs filtered sequences
# These panels help confirm that most retained sequences cluster around the expected barcode length.
ggplot(len_df, aes(x = length)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  facet_wrap(~dataset, ncol = 1, scales = "free_x") +
  labs(
    title = "Distribution of COI-5P sequence lengths (NCBI)",
    x     = "Sequence length (bp)",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


box_all <- boxplot.stats(len_all)
box_all$out # vector of suspicious lengths

# Inspect a few extreme sequences directly to decide if any should be excluded
which(len_all %in% head(box_all$out, 10))
dna_all[which(len_all %in% head(box_all$out, 10))]


#### Exploratory graphs from BOLD ####


# Build a simple world basemap to contextualize specimen coordinates
world_map <- map_data("world")

# Global map of all BOLD occurrences to check broad spatial coverage
ggplot() +
  geom_polygon(
    data = world_map, aes(x = long, y = lat, group = group),
    fill = "lightblue", color = "gray70", alpha = 0.5
  ) +
  geom_point(
    data = data_api_cleaner, aes(x = long, y = lat),
    color = "red", alpha = 0.5
  ) +
  labs(title = "Specimen Locations on Globe From BOLD", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


#### NCBI data cleaning ####


# Extract FASTA headers and split into tokens so we can pull out genus + species
hdr <- names(dna_filt)
parts <- strsplit(hdr, " ")

# Derive a "Genus_species" label from tokens 2 and 3 when possible.
# Headers that don't follow this convention are marked as NA and removed later.
species_vec <- sapply(parts, function(x) {
  if (length(x) >= 3) {
    paste(x[2], x[3], sep = "_")
  } else {
    NA_character_
  }
})

# Drop sequences where we couldn't reliably parse a species name
keep_species <- !is.na(species_vec)
dna_filt <- dna_filt[keep_species]
species_vec <- species_vec[keep_species]

# Fix the random seed so the species subsampling step is reproducible
set.seed(123)

# Identify unique species across the filtered sequences
unique_species <- unique(species_vec)
length(unique_species) # quick check on species richness in NCBI data

# Limit the tree to at most 300 species to keep alignment/ML optimisation tractable
n_target <- min(300, length(unique_species))
chosen_species <- sample(unique_species, size = n_target)
chosen_species

# Map each sequence to its species and index, so we can select representatives by row index
seq_df <- data.frame(
  idx = seq_along(species_vec),
  species = species_vec,
  stringsAsFactors = FALSE
)

# For each chosen species, keep the first encountered sequence as a single representative.
# This ensures "one sequence per species" for the tree without biasing towards species that have many records.
one_per_species <- seq_df %>%
  filter(species %in% chosen_species) %>%
  group_by(species) %>%
  summarise(idx = idx[1], .groups = "drop")

# Indices of sequences to retain for phylogenetic analysis
keep_idx <- one_per_species$idx

# Subset to the representative sequences and use species names as labels
dna_tree <- dna_filt[keep_idx]
names(dna_tree) <- one_per_species$species

# Ensure all sequence names are unique so downstream phylo functions do not fail
names(dna_tree) <- make.unique(names(dna_tree))


#### Sequence Alignment ####


# Align representative COI-5P sequences before tree building.
dna_aln <- AlignSeqs(dna_tree, processors = NULL)

dna_aln # useful to print briefly in an interactive session to confirm alignment succeeded
writeXStringSet(dna_aln, "../Data/Bombus_COI_aligned.fasta") # save aligned sequences for reuse/QC


# Convert alignment to phyDat, the format expected by phangorn for distance and ML methods.
aln_mat <- as.matrix(dna_aln)
class(aln_mat) <- "matrix"

phydat <- phyDat(aln_mat, type = "DNA")

# Compute ML-based pairwise distances and build a neighbour-joining tree to use as a starting tree for later maximum likelihood optimisation.
dm <- dist.ml(phydat)
tree_nj <- NJ(dm)

# Initialise an ML fit object using the NJ tree; this provides starting branch lengths and topology for subsequent model optimisation (e.g., GTR).
set.seed(125)
fit <- pml(tree_nj, phydat)


#### GTR ####

# Re-use the same seed so that stochastic tree rearrangements in optim.pml() are reproducible across runs (important for debugging and reporting).
set.seed(125)

# Optimise a maximum likelihood tree under a GTR substitution model.
# - 'fit' provides the starting tree and branch lengths (from NJ).
# - 'rearrangement = "stochastic"' allows tree topology changes beyond simple NNI,  improving the chance of escaping local optima but increasing runtime.
fit_GTRg <- optim.pml(
  fit,
  model = "GTR",
  rearrangement = "stochastic",
  control = pml.control(trace = 0)
) # This will take time


#### GTR Visual ####


# Plot the final ML tree under the GTR model to visualise overall topology and branch lengths.
# This is the main phylogeny used in downstream biogeographic and clade/sister-group analyses.
par(mar = c(2, 2, 5, 1))
plot(
  fit_GTRg$tree,
  type            = "phylogram",
  use.edge.length = TRUE,
  cex             = 0.7,
  font            = 3,
  edge.width      = 2,
  no.margin       = FALSE
)
title("ML tree for Bombus (GTR)", line = 2)

# Use the NJ tree as a separate QC dendrogram to screen for unusual sequences without relying on the full ML optimisation.
qc_tree <- tree_nj

par(mar = c(2, 2, 1, 1))
plot(
  qc_tree,
  type = "phylogram",
  use.edge.length = TRUE,
  cex = 0.8,
  main = "QC dendrogram of Bombus COI-5P sequences"
)

# Root-to-tip distances approximate how divergent each tip is from the rest of the tree. Outliers here can flag sequences that are mis-annotated, non-homologous, or low quality.
tip_depths <- node.depth.edgelength(qc_tree)[1:length(qc_tree$tip.label)]

# Summary of root-to-tip distances to understand the typical range.
summary(tip_depths)

# Identify statistically extreme branch-length values using the boxplot outlier rule.
depth_stats <- boxplot.stats(tip_depths)
depth_stats$out

# Map those extreme distances back to their corresponding tip labels (species). These are candidates for closer manual inspection rather than automatic removal.
outlier_tips <- qc_tree$tip.label[tip_depths %in% depth_stats$out]
outlier_tips

# Colour code tips so that long-branch outliers stand out in red.
tip_cols <- ifelse(qc_tree$tip.label %in% outlier_tips, "red", "black")

par(mar = c(2, 2, 2, 1))
plot(
  qc_tree,
  type = "phylogram",
  use.edge.length = TRUE,
  tip.color = tip_cols,
  cex = 0.8,
  main = "QC dendrogram – long-branch outliers in red"
)


#### Merging NCBI data with BOLD ####


# Create a version of the species name that matches your tree labels
data_api_cleaner <- data_api_cleaner %>%
  mutate(
    species_tree = gsub(" ", "_", species)
  )

head(data_api_cleaner$species) # sanity check


# Species present in the tree
# These are the only taxa for which we have both sequence data and a position in the ML tree.
tree_species <- fit_GTRg$tree$tip.label
tree_species[1:10] # inspect a subset of tip labels

# Keep only BOLD records for species that are actually in the tree. This restricts spatial analyses to taxa represented in the phylogeny.
bold_tree <- data_api_cleaner %>%
  filter(species_tree %in% tree_species)

length(unique(bold_tree$species_tree)) # number of tree-matched species in BOLD
table(bold_tree$species_tree %in% tree_species) # confirm that all records match tree tips

# Convert to sf object for mapping and spatial operations, keeping original lat/long columns.
occ_sf <- st_as_sf(
  bold_tree,
  coords = c("long", "lat"),
  crs    = 4326,
  remove = FALSE
)


#### Sister Taxa  ####


# Work on the final ML tree to identify sister relationships at the tip level
tree <- fit_GTRg$tree
n_tips <- length(tree$tip.label)

# allocate a group ID for each tip; NA = not assigned to a simple sister pair
sister_group <- rep(NA_integer_, n_tips)
group_id <- 1L

# Loop over internal nodes and find nodes whose two species are both tips.
# These represent simple sister-species pairs in the tree.
for (node in (n_tips + 1):(n_tips + tree$Nnode)) {
  children <- tree$edge[tree$edge[, 1] == node, 2]

  # if this node has exactly 2 children and both are tips, they are a sister pair
  if (length(children) == 2 && all(children <= n_tips)) {
    idx <- children

    # only assign if they don't already have a group
    if (all(is.na(sister_group[idx]))) {
      sister_group[idx] <- group_id
      group_id <- group_id + 1L
    }
  }
}

# Lookup table linking each tree tip to its sister group (or NA if unpaired)
sister_df <- data.frame(
  species_tree = tree$tip.label,
  sister_group = sister_group
)

# Drop sf geometry so we can use dplyr joins and summaries on a plain data frame
occ_df <- sf::st_drop_geometry(occ_sf)

# sanity checks on the subset of BOLD records that intersect the tree
summary(occ_df$lat)
summary(occ_df$long)

# Attach sister-group IDs to each occurrence record based on matching species_tree labels
occ_df <- occ_df %>%
  left_join(sister_df, by = "species_tree")

# Keep only occurrences that belong to a named sister group (i.e., part of a tip–tip pair)
occ_sis <- occ_df %>%
  dplyr::filter(!is.na(sister_group))

# Build a readable summary of sister groups: one row per sister_group, listing how many tips it contains and which taxa are included.
sister_list <- sister_df %>%
  as.data.frame() %>%
  filter(!is.na(sister_group)) %>%
  group_by(sister_group) %>%
  summarise(
    n_tips = n(),
    taxa = paste(sort(species_tree), collapse = ", "),
    .groups = "drop"
  )

sister_list # table of sister groups and their member species


#### Clade Groups ####


# Work with the final ML tree to define broader "clade groups" based on overall distances among tips.
tree <- fit_GTRg$tree
n_tips <- length(tree$tip.label)

# distance matrix between tips
# cophenetic() returns pairwise distances along the tree, capturing total evolutionary divergence between each pair of taxa.
dist_mat <- cophenetic(tree)

# Hierarchical clustering on those distances. This clusters tips into groups of closely related taxa using "average" linkage.
hc <- hclust(as.dist(dist_mat), method = "average")

# Cut into k relatedness/clade groups
# Here k = 5 is an analysis choice that balances resolution and interpretability.
k <- 5
clade_group <- cutree(hc, k = k)

# Make lookup table: species -> clade
# This table is the bridge between tree tip labels and group-level mapping/plots.
clade_df <- data.frame(
  species_tree = names(clade_group),
  clade_group  = factor(clade_group)
)

head(clade_df)
table(clade_df$clade_group)

# Vector of species for each group; useful for reporting membership and debugging.
clade_members <- split(clade_df$species_tree, clade_df$clade_group)

clade_members

# Assign a distinct colour to each clade group using a qualitative palette.
# These colours are used to colour tree tips.
clade_cols <- setNames(
  RColorBrewer::brewer.pal(5, "Set1"),
  levels(clade_df$clade_group)
)

# Map clade colours onto tree tip order so that plot() can colour tips correctly.
tip_cols <- clade_cols[clade_df$clade_group][match(
  tree$tip.label,
  clade_df$species_tree
)]


#### Clade Visual ####


# Visualise the ML tree with tips coloured by clade_group to show how the distance-based clusters correspond to major Bombus lineages.
plot(tree,
  tip.color = tip_cols, cex = 0.7,
  main = "ML tree for Bombus (GTR) - Labelled by Clade"
)

for (g in levels(clade_df$clade_group)) {
  cat("\nClade", g, ":\n")
  spp <- sort(clade_df$species_tree[clade_df$clade_group == g])
  print(spp)
}


#### Sister taxa ML tree ####


# Helper function to plot the ML tree with:
# - each sister-group pair in its own solid colour
# - unpaired tips shown in semi-transparent grey to de-emphasize them
plot_sister_tree <- function(tree,
                             sister_group,
                             alpha_unpaired = 0.2,
                             alpha_paired = 1,
                             palette = "Dark2") {
  n_tips <- length(tree$tip.label)

  # Encode sister groups as a factor so we can map each group to a specific colour
  # (NA correspond to unpaired tips).
  sister_factor <- factor(sister_group)
  n_groups <- nlevels(sister_factor)

  # Draw a base colour palette with one colour per sister group.
  base_cols <- hcl.colors(n_groups, palette)

  # Initialise a colour vector for all tips
  tip_cols <- rep(NA_character_, n_tips)

  # Logical indices for tips that are part of a sister pair vs unpaired
  paired <- !is.na(sister_group)
  unpaired <- is.na(sister_group)

  # Colours for paired tips: one colour per sister_group, with full opacity.
  tip_cols[paired] <- adjustcolor(
    base_cols[as.integer(sister_factor[paired])],
    alpha.f = alpha_paired
  )

  # Colours for unpaired tips: neutral grey with transparency
  tip_cols[unpaired] <- adjustcolor("grey40", alpha.f = alpha_unpaired)

  # Plot the phylogram with tip colours encoding sister-groups.
  par(mar = c(2, 2, 4, 2))
  plot(
    tree,
    type = "phylogram",
    use.edge.length = TRUE,
    tip.color = tip_cols,
    cex = 0.7
  )
  title("ML tree for Bombus (GTR) – Sister taxa coloured", line = 2)
}

# Use the optimised GTR ML tree for the sister-taxa visualisation
tree <- fit_GTRg$tree

# Plot ML tree with coloured sister taxa; unpaired tips are semi-transparent grey.
plot_sister_tree(tree, sister_group,
  alpha_unpaired = 0.2,
  alpha_paired   = 1
)

# Drop geometry and attach clade IDs to each occurrence record for downstream spatial analyses.
occ_df <- sf::st_drop_geometry(occ_sf)

occ_df <- occ_df %>%
  left_join(clade_df, by = "species_tree")


#### Sister taxa plotting ####


# plot 2D density heatmaps of sister-group occurrences.
plot_sister_density_region <- function(occ_sis, world_map, xlim, ylim, title_suffix = "") {
  ggplot(occ_sis, aes(x = long, y = lat)) +
    # base map
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = "lightblue",
      color = "gray70",
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    # 2D binned density of occurrence counts, facetted by sister group
    stat_bin_2d(
      aes(fill = ..count..),
      bins = 100,
      alpha = 0.85
    ) +
    facet_wrap(~sister_group) +
    labs(
      title = paste0(
        "Density of BOLD occurrences for Bombus sister taxa (COI ML tree)",
        title_suffix
      ),
      x = "Longitude",
      y = "Latitude",
      fill = "N records"
    ) +
    scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim)
}

# Split sister groups into two sets so that each figure contains a manageable number of facets.
g1 <- c(1, 2, 3, 4, 5, 6)
g2 <- setdiff(sort(unique(occ_sis$sister_group)), g1)

# Europe – sister groups 1–6
plot_sister_density_region(
  occ_sis = dplyr::filter(occ_sis, sister_group %in% g1),
  world_map = world_map,
  xlim = c(-10, 30),
  ylim = c(30, 70),
  title_suffix = " in Europe (groups 1–6)"
)

# Europe – the rest
plot_sister_density_region(
  occ_sis = dplyr::filter(occ_sis, sister_group %in% g2),
  world_map = world_map,
  xlim = c(-10, 30),
  ylim = c(30, 70),
  title_suffix = " in Europe (groups 7–13)"
)

# North America – sister groups 1–6
plot_sister_density_region(
  occ_sis = dplyr::filter(occ_sis, sister_group %in% g1),
  world_map = world_map,
  xlim = c(-170, -50),
  ylim = c(8, 90),
  title_suffix = " in North America (groups 1–6)"
)

# North America – the rest
plot_sister_density_region(
  occ_sis = dplyr::filter(occ_sis, sister_group %in% g2),
  world_map = world_map,
  xlim = c(-170, -50),
  ylim = c(8, 90),
  title_suffix = " in North America (groups 7–13)"
)


#### Density of clade group occurence ####


# plot clade-level occurrence density as 2D binned heatmaps.
plot_clade_density_region <- function(occ_df, world_map, xlim, ylim, title_suffix = "") {
  ggplot(occ_df, aes(x = long, y = lat)) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = "lightblue",
      color = "gray70",
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
    # 2D binned density of occurrence counts, faceted by clade group
    stat_bin_2d(
      aes(fill = ..count..),
      bins = 150,
      alpha = 0.8
    ) +
    facet_wrap(~clade_group) +
    labs(
      title = paste0(
        "BOLD occurrences for Bombus clade groups in COI ML tree",
        title_suffix
      ),
      x = "Longitude",
      y = "Latitude",
      fill = "N records"
    ) +
    scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim)
}

# North America-specific clade density; same function, different geographic window
plot_clade_density_region(
  occ_df = occ_df,
  world_map = world_map,
  xlim = c(-170, -50),
  ylim = c(8, 90),
  title_suffix = " in North America"
)

# Europe-specific clade density for comparison with the North American patterns
plot_clade_density_region(
  occ_df = occ_df,
  world_map = world_map,
  xlim = c(-10, 30),
  ylim = c(30, 70),
  title_suffix = " in Europe"
)


#### Clade Groups Visualization Global ####


# plot raw occurrence points coloured by clade group
# This complements the binned density plots by showing the actual point locations.
plot_clade_map <- function(occ_df, world_map, xlim, ylim, title_suffix = "") {
  ggplot() +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = "lightblue",
      color = "gray70",
      alpha = 0.7
    ) +
    geom_point(
      data = occ_df,
      aes(x = long, y = lat, colour = clade_group),
      alpha = 0.8,
      size = 2
    ) +
    labs(
      title = paste0(
        "BOLD occurrences for Bombus clade groups in COI ML tree",
        title_suffix
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim)
}

# World map of all Bombus occurrences coloured by clade group
plot_clade_map(
  occ_df = occ_df,
  world_map = world_map,
  xlim = c(-180, 180),
  ylim = c(-90, 90),
  title_suffix = ""
)


#### NMDS ####


# Run a 2D non-metric multidimensional scaling (NMDS) on the patristic distance matrix.
nmds <- metaMDS(dist_mat, k = 2, trymax = 100)

# Extract NMDS scores for each tip into a plain data frame
scores_df <- as.data.frame(scores(nmds))
scores_df$species_tree <- rownames(scores_df)

# Join NMDS coordinates with clade membership so we can colour/hull by clade in ordination space. Using base merge() here avoids adding a dplyr dependency to this specific step.
scores_df <- merge(scores_df, clade_df, by = "species_tree", all.x = TRUE)
scores_df$clade_group <- factor(scores_df$clade_group)

# Split the NMDS scores into one data frame per clade_group for hull calculation.
split_list <- split(scores_df, scores_df$clade_group)

# For each clade, calculate the hull in ordination space (NMDS1 x NMDS2).
# Hulls approximate the outer envelope of each clade's spread in the NMDS plot.
hulls_list <- lapply(split_list, function(df) {
  if (nrow(df) >= 3) {
    df[chull(df$NMDS1, df$NMDS2), ]
  } else {
    df
  }
})

# combine hull coordinates into a single data frame for plotting polygons per clade.
hulls <- do.call(rbind, hulls_list)


#### NMDS Plotting ####


# NMDS ordination of Bombus tips, coloured by clade_group.
# Hulls and ellipses show how each clade occupies ordination space based on patristic distances.
ggplot(scores_df, aes(x = NMDS1, y = NMDS2, colour = clade_group)) +
  # hull polygons outline the outer boundary of each clade in NMDS space, highlighting group separation and overlap.
  geom_polygon(
    data = hulls,
    aes(x = NMDS1, y = NMDS2, fill = clade_group, group = clade_group),
    alpha = 0.2,
    colour = NA
  ) +
  # Ellipses capture the core spread of each clade (80% confidence),
  stat_ellipse(
    aes(group = clade_group),
    level = 0.8,
    linewidth = 0.5,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  # Each point represents a tip (species) positioned according to its NMDS coordinates.
  geom_point(size = 5, alpha = 0.9) +
  coord_equal() +
  labs(
    title = "NMDS of Bombus ML tree (patristic distances)",
    x     = "NMDS1",
    y     = "NMDS2"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )




