# Seurat DoubletFinder Workflow

This repository provides a full workflow for filtering out doublets in single-cell RNA-seq data using Seurat and DoubletFinder in R.

## Features

- End-to-end R script for doublet detection and filtering in scRNA-seq datasets
- PBMC 10k 3' NextGem example (Cell Ranger output)
- Quality control, filtering, clustering, and DoubletFinder pipeline
- Inline comments to guide customization for your own data

## Requirements

- R (>= 4.0)
- Seurat
- DoubletFinder
- ggplot2
- tidyverse

Install required R packages:
```r
install.packages(c("Seurat", "ggplot2", "tidyverse"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

## Usage

1. Place your 10x Genomics output files in the expected path (`data/raw_feature_bc_matrix/`).
2. Run the provided R script:
   ```
   Rscript doubletfinder_workflow.R
   ```
3. Follow inline comments to customize thresholds or adapt to your dataset.

## References

- [Seurat](https://satijalab.org/seurat/)
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)

## License

MIT License

---

For questions or contributions, please open an issue or pull request.