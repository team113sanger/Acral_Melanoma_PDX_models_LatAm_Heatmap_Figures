library(tidyverse)
library(ComplexHeatmap)

# Metadata ----
metadata <- read_tsv("data/sample_ids.tsv")

# Sample order
sample_order <- metadata$BR_ID[order(metadata$BR_ID)]
# Vector to change PD_ID -> BR_ID
rename_vector <- setNames(metadata$BR_ID, metadata$PD_ID)

# Driver alteration data ----
dat <- read_tsv("data/all_samples_driver_alterations.tsv")

# Plotting Samples + multiple passages ----
# Select Patients that have more than 1 pdx
multi_ids <- metadata |>
  group_by(BR_Sample_ID) |>
  filter(n() > 2)

# Filter alteration list by sample 
multi_compare <- dat |>
  filter(Sample %in% multi_ids$PD_ID)

comp_mat <- multi_compare |>
  pivot_wider(names_from = Sample, values_from = Alteration, values_fill = "Neutral") |>
  column_to_rownames("Gene") |>
  mutate(across(everything(), ~gsub("Neutral;|Neutral", "", .)))

# Rename to BR_ID
colnames(comp_mat) <- rename_vector[colnames(comp_mat)]

# Complex heatmap ----
# Subset sample_order with plotted samples only
col_order <- intersect(sample_order, colnames(comp_mat))

# Split alterations by ;
get_type_fun = function(x) strsplit(x, ";")[[1]]

col <- c(
  "Missense_Mutation" = "#00A087B2",
  "Frame_Shift_Ins" = "#4DBBD5B2", 
  "In_Frame_Ins" = "#7E6148B2", 
  "Gain" = "#E64B35B2",
  "Loss" = "#8491B4B2",
  "Amp" = "#DC0000B2",
  "Del" = "#3C5488B2"
)

# Alteration functions
alter_fun <- list(
  background = function(x, y, width, height) {
    grid.rect(x, y, width, height, gp = gpar(fill = "grey90", col = "white", lwd = 1))
  },
  Missense_Mutation = function(x, y, width, height) {
    grid.rect(x, y, width * 0.9, height * 0.9, gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, width, height) {
    grid.rect(x, y, width * 0.9, height * 0.9, gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Ins = function(x, y, width, height) {
    grid.rect(x, y, width * 0.9, height * 0.9, gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  Amp = function(x, y, width, height) {
    grid.polygon(
      x = unit(c(x - width * 0.3, x + width * 0.3, x), "npc"),
      y = unit(c(y - height * 0.3, y - height * 0.3, y + height * 0.3), "npc"),
      gp = gpar(fill = col["Amp"], col = NA)
    )
  },
  Del = function(x, y, width, height) {
    grid.polygon(
      x = unit(c(x - width * 0.3, x + width * 0.3, x), "npc"),
      y = unit(c(y + height * 0.3, y + height * 0.3, y - height * 0.3), "npc"),
      gp = gpar(fill = col["Del"], col = NA)
    )
  },
  Gain = function(x, y, width, height) {
    grid.polygon(
      x = unit(c(x - width * 0.2, x + width * 0.2, x), "npc"),
      y = unit(c(y, y, y + height * 0.4), "npc"),
      gp = gpar(fill = col["Gain"], col = NA)
    )
  },
  Loss = function(x, y, width, height) {
    grid.polygon(
      x = unit(c(x - width * 0.2, x + width * 0.2, x), "npc"),
      y = unit(c(y, y, y - height * 0.4), "npc"),
      gp = gpar(fill = col["Loss"], col = NA)
    )
  }
)

# Annotations ----
anno_dat <- multi_ids |>
  mutate(passage = ifelse(Sample_passage == "Patient", "P", Sample_passage)) |>
  dplyr::select(BR_Sample_ID, BR_ID, passage) |>
  arrange(match(BR_ID, colnames(comp_mat)))

# Create an alternating color scheme based on BR_Sample_ID
unique_ids <- unique(anno_dat$BR_Sample_ID)  # Get unique sample groups
unique_ids <- unique_ids[order(unique_ids)]
group_colors <- setNames(rep(c("black", "white"), length.out = length(unique_ids)), unique_ids)  # Assign alternating colors
bar_colors <- group_colors[anno_dat$BR_Sample_ID]  # Assign the colors to each column

# Create a vector to store BR_Sample_ID labels, with only the first occurrence per group
unique_labels <- rep("", length(anno_dat$BR_Sample_ID))

# Find the first index of each unique BR_Sample_ID group
unique_labels[!duplicated(anno_dat$BR_Sample_ID)] <- anno_dat$BR_Sample_ID[!duplicated(anno_dat$BR_Sample_ID)]

# Define column annotation
top_annotation <- HeatmapAnnotation(
  BR_Sample_ID = anno_text(
    unique_labels,
    which = "column",  # Position the label at the top of each column
    just = "left",
    location = 0,
    gp = gpar(fontsize = 10, fontface = "bold", col = "black")
  ),
  GroupBar = anno_simple(factor(bar_colors, levels = c("black", "white")), 
                         col = c("black" = "black", "white" = "white"),
                         border = TRUE,
                         gp = gpar(col = NA)),
  Passage = anno_text(anno_dat$passage, gp = gpar(fontsize = 10, fontface = "bold"),
                      rot = 0, just = "center"),
  show_annotation_name = FALSE,
  gap = unit(2, "mm")
)

# Oncoprint ----
onco <- oncoPrint(comp_mat,
          alter_fun = alter_fun,
          col = col,
          show_column_names = TRUE,
          row_names_side = "left",
          show_pct = FALSE,
          right_annotation = NULL,
          column_order = col_order,
          top_annotation = top_annotation
          )
pdf("figures/Figure_S8b.pdf", height = 8, width = 10)
draw(onco)
dev.off()
