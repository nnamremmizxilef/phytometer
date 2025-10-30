# phytometer

## Overview

phytometer provides tools for analyzing multi-omics data integrated with environmental measurements in ecological and evolutionary studies.

## Installation

Install from GitHub:
```r
devtools::install_github("nnamremmizxilef/phytometer")
```

## Quick Start
```r
library(phytometer)

# Load example datasets
data(envirodata)
data(phenotydata)
data(genotydata)

# Explore data
result <- showtidydata(envirodata, cor_method = "spearman")
print(result$distribution_plot)
```

## Features

- Multi-omics integration
- Environmental data analysis
- Distribution plots and correlation heatmaps
- Example datasets included

## License

MIT © Felix Zimmermann