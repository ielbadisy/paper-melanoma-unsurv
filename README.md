# Survival Trajectory Phenotypes in Advanced Melanoma Under Systemic Therapy

This repository contains the R analysis code for the manuscript:

**Survival Trajectory Phenotypes in Advanced Melanoma Under Systemic Therapy**

The analysis combines deep survival learning, unsupervised clustering of predicted survival trajectories, and dynamic restricted mean survival time (RMST) summaries to identify and describe progression-free survival phenotypes in an advanced melanoma cohort.

## Repository contents

- `code.R`: complete R analysis pipeline used for model fitting, survival-trajectory clustering, RMST summaries, tables, and figures.

The patient-level dataset is not included in this public repository.

## Data availability

The analysis script expects a local file named `melanoma2.csv` in the repository root. This file is intentionally excluded because it contains patient-level study data.

To run the code with an authorized local dataset, provide a CSV file with the variables used by `code.R`, including:

- `PFS`
- `Event_of_PFS`
- `Sex`
- `Age`
- `Histology`
- `Comorbidity`
- `PS`
- `Lung_Metastasis`
- `GG_metastasis`
- `Bone_metastasis`
- `Other`
- `therapy_type`
- `Therapeutic_line`
- `Total_doses`
- `Response_type`

Expected treatment labels in the raw input are `dacarbazine` and `pembrolizumab`. Expected therapy-line labels are `first_line` and `second_line`.

## R packages

The script uses the following R packages:

- `survival`
- `survdnn`
- `unsurv`
- `tvrmst`
- `ggplot2`
- `dplyr`
- `tidyr`
- `table1`
- `patchwork`
- `pacman`

## Running the analysis

Install the required packages, place the authorized `melanoma2.csv` file in the repository root, then run:

```r
source("code.R")
```

The script writes manuscript tables and figures to the working directory.

## Notes

The analysis is descriptive and intended to reproduce the statistical workflow reported in the manuscript. The dataset is not distributed with this repository.
