library(tidyverse)
library(ComplexHeatmap)

# Table with hailstorm presence in each sample
dat <- read_tsv("data/All_samples_hailstorms_inspected_table.tsv")

# Metadata ----
metadata <- read_tsv("data/sample_ids.tsv")

# Create a rename vector for PD_ID -> BR_ID
rename_vector <- setNames(metadata$BR_ID, metadata$PD_ID)

# Select patients that have PDXs
ids <- metadata |>
  group_by(BR_Sample_ID) |>
  filter(n() > 1) |>
  ungroup() |>
  mutate(passage = ifelse(Sample_passage == "Patient", "P", Sample_passage)) |>
  select(BR_Sample_ID, BR_ID, passage, PD_ID)


# Create matrix using only the samples with pdx
mat <- dat |>
  column_to_rownames("Arm") |>
  select(one_of(ids$PD_ID)) |>
  filter_all(any_vars(. != "No"))
colnames(mat) <- rename_vector[colnames(mat)]
mat[mat == "No"] <- ""

chr_order <- c("1p", "1q", "3q", "4p", "4q", "5p", "7p", "10p", "11q", "12q", "22q")

# Complex heatmaps ----
col <- c(
  "Yes" = "#00A087FF",
  "Manual" = "#91D1C2FF")

# Alteration functions
alter_fun <- list(
  background = function(x, y, width, height) {
    grid.rect(x, y, width, height, gp = gpar(fill = "white", col = "grey90", lwd = 1))
  },
  Yes = function(x, y, width, height) {
    grid.rect(x, y, width * 0.95, height * 0.95, gp = gpar(fill = col["Yes"], col = NA))
  },
  Manual = function(x, y, width, height) {
    grid.rect(x, y, width * 0.95, height * 0.95, gp = gpar(fill = col["Manual"], col = NA))
  }
)

# Ensure the metadata is in the same order as the heatmap columns
anno_ids <- ids |> 
  filter(BR_ID %in% colnames(mat)) |>
  arrange(match(BR_ID, colnames(mat)))

# Create an alternating color scheme based on BR_Sample_ID
unique_ids <- unique(anno_ids$BR_Sample_ID)  # Get unique sample groups
group_colors <- setNames(rep(c("black", "white"), length.out = length(unique_ids)), unique_ids)  # Assign alternating colors
bar_colors <- group_colors[anno_ids$BR_Sample_ID]  # Assign the colors to each column

# Create a vector to store BR_Sample_ID labels, with only the first occurrence per group
unique_labels <- rep("", length(anno_ids$BR_Sample_ID))

# Find the first index of each unique BR_Sample_ID group
unique_labels[!duplicated(anno_ids$BR_Sample_ID)] <- anno_ids$BR_Sample_ID[!duplicated(anno_ids$BR_Sample_ID)]


# Define column annotation
top_annotation <- HeatmapAnnotation(
  BR_Sample_ID = anno_text(
    unique_labels,
    which = "column",  # Position the label at the top of each column
    just = "left",
    location = 0,
    gp = gpar(fontsize = 12, fontface = "bold", col = "black")
  ),
  GroupBar = anno_simple(factor(bar_colors, levels = c("black", "white")), 
                         col = c("black" = "black", "white" = "white"),
                         border = TRUE,
                         gp = gpar(col = NA)),
  Passage = anno_text(anno_ids$passage, gp = gpar(fontsize = 10, fontface = "bold"),
                      rot = 0, just = "center"),
  show_annotation_name = FALSE,
  gap = unit(2, "mm")
)

sample_order <- anno_ids$BR_ID[order(anno_ids$BR_ID)]
# Oncoprint ----
onco <- oncoPrint(mat,
                  alter_fun = alter_fun,
                  col = col,
                  show_column_names = TRUE,
                  row_names_side = "left",
                  show_pct = FALSE,
                  right_annotation = NULL,
                  column_order = sample_order,
                  row_order = chr_order,
                  top_annotation = top_annotation,
                  heatmap_legend_param = list(title = "Hailstorm presence", at = c("Yes", "Manual"),
                                              labels = c("Detected", "Manual inspection"))
)

pdf("figures/Figure_4d2.pdf", height = 5, width = 10)
draw(onco)
dev.off()