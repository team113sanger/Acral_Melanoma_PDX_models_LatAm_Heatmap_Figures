# Molecular and functional profiling of Brazilian acral melanoma reveals aberrant genetic landscape and therapeutic vulnerabilities
This repository is a sub-repository and is included with other repositories in the [Acral Melanoma PDX models from Latin America](https://github.com/team113sanger/Acral_Melanoma_PDX_models_from_Latin_America) parent repository.

## Overview
This repository contains R scripts used to plot figures in the manuscript. Some figures were further edited using the software **Affinity Designer 2** [here](https://affinity.serif.com/en-gb/designer/) in order to get to their final versions used in the manuscript.

## Directory structure
```bash
tree -L 2
.
├── README.md
├── data
│   ├── All_samples_hailstorms_inspected_table.tsv
│   ├── all_lesions.conf_75.txt
│   ├── all_samples_driver_alterations.tsv
│   ├── cnv_percentages.tsv
│   ├── sample_ids.tsv
│   └── sample_mutation_types.csv
├── figures
│   ├── Figure_4d2.pdf
│   ├── Figure_S2b.pdf
│   ├── Figure_S8a.pdf
│   └── Figure_S8b.pdf
└── scripts
    ├── Figure_4d2.R
    ├── Figure_S2b.R
    ├── Figure_S8a.R
    └── Figure_S8b.R

3 directories, 15 files

```

## Required Software

The following software is required to recreate the figures:

- **R** (version 4.2.2 or later) [Download here](https://cran.r-project.org/)
- **R Packages**:
  - tidyverse
  - ComplexHeatmap

If you're interested in [reproducing the R environment using `renv`](https://rstudio.github.io/renv/reference/index.html) please follow the official documentation.

## Contact
  If you have any questions or comments about this repository, please contact:
  - Antonio Facciolo : (<af31@sanger.ac.uk>)