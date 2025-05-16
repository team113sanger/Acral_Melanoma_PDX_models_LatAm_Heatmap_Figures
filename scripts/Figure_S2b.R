
library(tidyverse)
library(ComplexHeatmap)

select <- dplyr::select

# Get matrix from Gistic data
dat <- read.csv("data/all_lesions.conf_75.txt", header = TRUE, sep = "\t")

samples <- colnames(dat)[grep("PD", colnames(dat))]

mat <- dat |>
  filter(Amplitude.Threshold != "Actual Copy Change Given") |>
  filter(q.values <= 0.1) |>
  select(Unique.Name, Descriptor, all_of(samples))

# Identify the sample columns
sample_cols <- grep("^PD", names(mat), value = TRUE)

# Modify the sample columns based on 'Unique.Name'
new_mat <- mat %>%
  mutate(across(all_of(sample_cols), 
                ~ ifelse(. != 0, ifelse(grepl("Amplification", Unique.Name), "Amp", "Del"), 0))) |>
  select(-1) |>
  column_to_rownames("Descriptor")

# Convert new_mat into a matrix (ensure it's character format)
mat <- as.matrix(new_mat)
mat[mat == "0"] <- ""

# Mutation data
mut_data <- read.csv("data/sample_mutation_types.csv", header = T, sep = ";", stringsAsFactors = F)
mut <- filter(mut_data, Tumor_Sample_Barcode %in% samples)

# Master table for renaming
metadata <- read_tsv("data/sample_ids.tsv")
metadata <- metadata |>
  filter(PD_ID %in% samples)
rename_vector <- setNames(metadata$BR_ID, metadata$PD_ID)

# CNV alteration percetages table
alter_dat <- read.delim("data/cnv_percentages.tsv")

alter_bar <- alter_dat |>
  select(Sample, total_amp_plus_gain, total_loss_plus_del) |>
  filter(Sample %in% samples)

# Plotting on complex heatmap

alter_fun <- list(
  background = function(x, y, width, height) {
    grid.rect(x, y, width, height, gp = gpar(fill = "grey90", col = "white", lwd = 1))
  },
  Amp = function(x, y, width, height) {
    grid.rect(x, y, width * 0.9, height * 0.9, gp = gpar(fill = "#DC0000FF", col = NA))
  },
  Del = function(x, y, width, height) {
    grid.rect(x, y, width * 0.9, height * 0.9, gp = gpar(fill = "#3C5488FF", col = NA))
  }
)

# Define color mapping for legend
col_list <- c("Amp" = "#DC0000FF", "Del" = "#3C5488FF")

# Define annotation colors for mutations
mutation_colors <- c(
  "NRAS" = "#00A087FF",
  "WT"   = "#E64B35FF",
  "KIT"  = "#4DBBD5FF",
  "HRAS" = "#91D1C2FF",
  "KRAS" = "#B09C85FF",
  "NF1"  = "#F39B7FFF"
)

# Create column annotation
column_ha <- columnAnnotation(
  Gained = anno_barplot(alter_bar$total_amp_plus_gain, gp = gpar(fill = "#DC0000FF")),
  Lost = anno_barplot(alter_bar$total_loss_plus_del, gp = gpar(fill = "#3C5488FF")),
  Mutation = mut$Mutation,  # Values from mut dataframe
  col = list(Mutation = mutation_colors),  # Define colors
  annotation_legend_param = list(title = "Mutation Type"),
  gap = unit(2, "mm")
  
)

colnames(mat) <- rename_vector[colnames(mat)]
mut$Tumor_Sample_Barcode <- rename_vector[mut$Tumor_Sample_Barcode]
alter_bar$Sample <- rename_vector[alter_bar$Sample]

# Create the oncoplot
pdf("figures/Figure_S2b.pdf" ,width = 8, height = 6)
oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col_list,
  row_names_gp = gpar(fontsize = 10),  # Adjust row font size
  column_names_gp = gpar(fontsize = 10),  # Adjust column font size
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",
  show_pct = TRUE,
  pct_side = "right",
  heatmap_legend_param = list(title = "Alteration"),
  top_annotation = column_ha,
  right_annotation = NULL
)
dev.off()

