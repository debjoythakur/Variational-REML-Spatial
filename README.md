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

├── R/
│   ├── vrmle.R
│   ├── lowrank_vrmle.R
│   ├── simulation.R
│   ├── utilities.R
│   └── plotting.R
│
├── data/
│
├── examples/
│
├── results/
│
├── manuscript/
│
├── README.md
├── LICENSE
└── CITATION.cff
```

---

# Installation

Clone the repository

```bash
git clone https://github.com/debjoythakur/Variational-REML-Spatial.git
```

Open R or RStudio and install the required packages.

---

# Example

```r
source("R/vrmle.R")

result <- vrmle(
    Y = Y,
    X = X,
    W = W
)
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
