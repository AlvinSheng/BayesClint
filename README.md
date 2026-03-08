# BayesClint

This is the BayesClint R package, which implements the method proposed in the manuscript Sheng, A., Chekouo, T., Safo, S. E. (2026). BayesClint: Bayesian multi-scale clustering and multi-sample integration with feature selection for spatial transcriptomics data.
It uses high-performance C++ via `Rcpp`, `RcppArmadillo`, and the GNU Scientific Library (GSL).

License: GPL (>= 3).

For more information, please contact tchekouo@umn.edu.

# Installation

This package contains compiled C++ code and requires the **GNU Scientific Library (GSL)** to be installed on your system.

## 1. Prerequisites

Before installing the R package, please follow the instructions for your Operating System:

### macOS (Sequoia and later)
The easiest way to install GSL is via [Homebrew](https://brew.sh/):
```bash
brew install gsl
```

*Note: If you are on an Apple Silicon Mac (M1/M2/M3), ensure your R installation is configured to look in /opt/homebrew/include.*

### Windows

1. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) (match the version to your R version).

2. Rtools includes the GSL headers by default, so no further manual GSL installation is typically required.

### Linux (Ubuntu/Debian)
```bash
sudo apt-get install libgsl-dev
```

## 2. Install from GitHub

Once the system dependencies are met, you can install the development version of BayesClint from GitHub using devtools:

```r
# Install devtools if you don't have it
if (!require("devtools")) install.packages("devtools")

# Install the package
devtools::install_github("AlvinSheng/BayesClint")
```

# Quick Start

See the example located at the bottom of the documentation for the main function, `BayesClint_run`, which applies BayesClint to a small example dataset, taking less than 10 minutes to analyze it.

```r
help(BayesClint_run)
```

# Troubleshooting

## Speed Issues

This package is optimized for performance using `-O3`. If the code feels slow, ensure you are not running R in a "Debug" mode and that your compiler supports the flags defined in `src/Makevars` or `src/Makevars.win`.

## Compilation Errors

If you see an error like `gsl/gsl_matrix.h: No such file or directory`, it means the GSL library was not found.

* **macOS:** Try running `brew link gsl` in the terminal.
* **Windows:** Ensure Rtools is in your system PATH.


