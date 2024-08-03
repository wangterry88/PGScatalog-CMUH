# Polygenic Score (PGS) Analysis for a hospital cohort

This project involves a series of scripts for Polygenic Score (PGS) analysis. The following three steps are included:

## Step 1: PGS-Catalog-Grep.sh

### Description

This script is designed to extract relevant information from a PGS catalog using grep. It helps filter and organize data based on specific criteria.

### Usage

```bash
bash PGS-Catalog-Grep.sh [options] input_file output_file
```

Replace `[options]`, `input_file`, and `output_file` with your specific parameters.

## Step 2: PGS-6-DistPlot.R

### Description

The R script `PGS-6-DistPlot.R` generates a distribution plot for a given PGS dataset. It provides insights into the distribution of polygenic scores across the analyzed samples.

### Dependencies

- R (>=3.5.0)

### Usage

```R
Rscript PGS-6-DistPlot.R input_file output_plot
```

Replace `input_file` and `output_plot` with your specific file paths.

## Step 3: PGS-AUC-Calculation.R

### Description

This R script calculates the Area Under the Curve (AUC) for a given PGS dataset. AUC is a metric commonly used to assess the performance of polygenic scores.

### Dependencies

- R (>=3.5.0)
- Required R packages (specified in the script)

### Usage

```R
Rscript PGS-AUC-Calculation.R input_file output_auc
```

Replace `input_file` and `output_auc` with your specific file paths.

## Contributing

If you'd like to contribute to this project, please follow the standard GitHub Fork and Pull Request workflow. Your contributions are welcome!

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
