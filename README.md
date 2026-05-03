# dd-cfDNA Meta-Analysis

Bayesian HSROC threshold-regression meta-analysis of donor-derived cell-free DNA (dd-cfDNA) for the diagnosis of acute heart-allograft rejection.

The analysis fits a hierarchical SROC model with a threshold-regression term across 14 included studies, propagates uncertainty from raw 2×2 counts through to clinical decision thresholds, and reports posterior summaries, leave-one-out diagnostics, and decision-curve analysis under a Gaussian-copula construction for repeat-test strategies.

## What's in here

- **`dd_cfdna_analysis_only.R`** — the analysis pipeline.
- **`DNA.xlsx`** — per-study 2×2 counts, decision thresholds (% donor fraction), and patient counts. Place this in the working directory before sourcing the script.

## Methods

A companion methods document describing the rationale behind every modelling choice — cluster-correction of effective sample size, prior specification, PSIS-LOO with exact-refit fallback, Gaussian-copula decision-curve construction, and the threshold-marginalised summary point — will be published alongside the manuscript.

## Requirements

R (≥ 4.1) with the following packages:

```r
install.packages(c("rjags", "loo", "pbivnorm", "readxl", "coda", "mada",
                   "MASS", "ggrepel", "dplyr", "ggplot2", "patchwork",
                   "scales", "tidyr"))
```

`rjags` requires the JAGS system binary. Install JAGS first, then the R package:

- macOS: `brew install jags`
- Ubuntu: `sudo apt-get install jags`
- Windows: <https://sourceforge.net/projects/mcmc-jags/>

## Running the analysis

```r
setwd("path/to/this/repo")   # so DNA.xlsx is in the working directory
source("dd_cfdna_analysis_only.R")
```

The script fits the HSROC model at four cluster-correlation values (ρ ∈ {0.00, 0.30, 0.50, 0.70}) with ρ = 0.50 as the primary analysis, then prints the comparison table, the per-study Pareto-k diagnostics, the final results summary, and the six figures (HSROC, sens/spec vs. threshold, two forest plots, sensitivity-anchor posterior, strategy-based DCA, ρ_repeat sensitivity, clinical decision map).

Expected runtime: roughly 45–90 minutes for the full ρ sweep on a modern laptop with four cores.

## Author

Jabez David John

## License

Released under the MIT License — see `LICENSE`. (Note: CC-BY 4.0 is the standard for open-access *manuscripts and figures*, but the Creative Commons organisation recommends against CC licenses for code; MIT is the conventional choice for research code repositories on GitHub.)
