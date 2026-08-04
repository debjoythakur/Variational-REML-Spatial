# Variational-REML-Spatial

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Official implementation of **Variational Restricted Maximum Likelihood Estimation (VRMLE)** for Gaussian Intrinsic Conditional Autoregressive (ICAR) spatial models.

This repository accompanies the paper:

> **Variational Approximated Restricted Maximum Likelihood Estimation for Spatial Data**
>
> Debjoy Thakur

---

# Overview

Classical Restricted Maximum Likelihood (REML) estimation for Gaussian ICAR spatial models becomes computationally demanding for large spatial datasets because it repeatedly requires sparse matrix factorization and inversion.

This repository implements a **Variational Restricted Maximum Likelihood Estimation (VRMLE)** framework that replaces the computationally expensive restricted likelihood optimization with a scalable **Evidence Lower Bound (ELBO)** optimization using coordinate-ascent variational inference.

The repository also contains a low-rank covariance approximation for large spatial datasets together with theoretical guarantees and reproducible simulation studies.

---

# Features

- Variational Restricted Maximum Likelihood (VRMLE)
- Gaussian ICAR spatial models
- Coordinate-ascent ELBO optimization
- Closed-form parameter updates
- Low-rank covariance approximation
- Convergence guarantees
- Simulation framework
- Real-data analysis
- Reproducible research code

---

# Repository Structure

```
Variational-REML-Spatial/
│
├── R/
│   ├── utilities.R
│   ├── graph_generation.R
│   ├── data_generation.R
│   ├── exact_mle.R
│   ├── vrmle.R
│   ├── lowrank_vrmle.R
│   ├── inla.R
│   ├── prediction.R
│   ├── posterior.R
│   └── plotting.R
│
├── simulations/
│   ├── run_simulation.R
│   └── simulation_parameters.R
│
├── results/
│   ├── figures/
│
├── manuscript/
│   └── Variational_REML.pdf
│
├── examples/
│   └── VREML_Simulation.Rmd
│
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```

---
# Installation

Clone the repository

```bash
git clone https://github.com/debjoythakur/Variational-REML-Spatial.git
cd Variational-REML-Spatial
```

Open R or RStudio and install the required packages.

Install the required R packages:
```r
install.packages(
  c(
    "MASS",
    "Matrix",
    "ggplot2",
    "dplyr",
    "RANN"
  )
)

#### Install INLA separately:

install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  ),
  dependencies = TRUE
)

```
---

# Example

```r
source("R/utilities.R")
source("R/graph_generation.R")
source("R/data_generation.R")
source("R/vrmle.R")

set.seed(123)

# Construct an 8 x 8 lattice, giving n = 64 spatial units
scenario <- build_scenario(m = 8)

# Generate Gaussian ICAR spatial data
simulated_data <- simulate_icar_data(
  X = scenario$X,
  H = scenario$H,
  Kinv = scenario$Kinv,
  beta = scenario$beta_true,
  sigma2_e = scenario$sigma2_e_true,
  sigma2_u = scenario$sigma2_u_true
)

# Fit the VRMLE model
fit <- fit_variational_reml_train(
  y_tr = simulated_data$y,
  X_tr = scenario$X,
  H_tr = scenario$H,
  K = scenario$Kmat
)

# Estimated observation-error variance
fit$sigma2_e

# Estimated spatial-effect variance
fit$sigma2_u

# Estimated regression coefficients
fit$beta

# Algorithm diagnostics
fit$converged
```

---

# Contents

The repository includes

- VRMLE algorithm
- Low-rank covariance VRMLE
- Simulation studies
- Real data example
- Figures used in the manuscript
- Reproducible scripts

---

# Citation

If you use this software, please cite

```
Debjoy Thakur.

Variational Approximated Restricted Maximum Likelihood Estimation
for Spatial Data.
```

(BibTeX will be added after publication.)

---

# License

This project is distributed under the MIT License.

See the LICENSE file for details.

---

# Contact

**Debjoy Thakur**

School of Arts and Sciences  
Ahmedabad University, India

Email:
- debjoy.thakur@ahduni.edu.in
- debjoythakur@outlook.com

GitHub:
https://github.com/debjoythakur
