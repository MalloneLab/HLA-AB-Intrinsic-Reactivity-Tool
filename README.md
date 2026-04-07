<p align="center">
  <img src="logo.png" width="400">
</p>

<p align="center">
  <b>HLA-antibodies Intrinsic Reactivity Tool</b><br>
  MalloneLab
</p>

---

## Description

Shiny tool and analysis code for quantifying assay-intrinsic HLA reactivity in single-antigen bead (SAB) assays.

This repository accompanies the manuscript:

**Reclassifying Assay-Intrinsic HLA Reactivity Expands Donor Availability Without Short-Term Allograft Risk**

---

## Overview

This tool implements a bead-level probabilistic model to estimate assay-intrinsic (non-alloimmune) HLA reactivity and supports interpretation of SAB results.

---

## Access the tool

👉 Web application:  
https://mallonelabimmunotools.shinyapps.io/hla_unspecificity_app/

---

## Repository contents

- `app.R` — Shiny application
- `.Rmd` files — analysis and model code
- `.csv` files — bead-level reference data
- `example_data/` — example input/output for reproducibility

---

## Requirements

- R (version 4.5.2)
- RStudio (version 2025.09.2+418)

Main packages:
- tidyverse
- dplyr
- ggplot2
- ggrepel
- caret
- pROC
- shiny

---

## Run locally

1. Open `app.R` in RStudio  
2. Install required packages  
3. Click **Run App**

---

## Example data

Example input and expected outputs are provided in the `example_data/` folder.

---

## Expected runtime

- App launch: < 10 seconds  
- Analysis: seconds to minutes depending on input  

---

## License

MIT License
