# phytometer 0.0.0.9000

## Initial Release

* Initial development version of phytometer package
* Package setup with MIT license and GitHub repository

## Functions

* Added `calculate_gdd()` - Calculate growing degree days from temperature data
* Added `showtidydata()` - Exploratory data analysis with distribution plots and correlation heatmaps
  - Supports both Pearson and Spearman correlation methods
  - Converts data to long format for visualization
  - Returns list with long data, distribution plot, and correlation matrix

## Data

* Added `envirodata` - Environmental monitoring dataset with 120 observations
  - Includes temperature, precipitation, soil moisture, VPD, PAR
  - Drought and herbivory indices
* Added `phenotydata` - Phenotypic measurements from 50 oak saplings
  - Morphological traits (height, DBH, leaf area, biomass)
  - Physiological measurements (stomatal conductance, photosynthesis, WUE)
* Added `genotydata` - Genotypic and molecular data from 50 individuals
  - SNP genotypes, epigenetic marks
  - Gene expression for drought and herbivory responses
  - Microbiome diversity metrics

## Dependencies

* Added dependencies: dplyr, tidyr, ggplot2, corrplot
* Added pipe operator support (`%>%`)

## Documentation

* Created README with installation instructions and examples
* Added function documentation with roxygen2
* Added dataset documentation