# modelTools

This repository contains a collection of MATLAB functions for performing contrast analysis on treatment-coded statistical models. The toolbox enables users to construct design matrices, compute contrast and condition weight matrices, run hypothesis tests (t-tests and F-tests), and generate detailed HTML reports summarizing model results.

## Overview

The functions in this toolbox are designed to help researchers and data analysts work with generalized linear models (GLMs) and ANOVA frameworks. Key capabilities include:

- **Design Matrix Construction:** Build design matrices that encode experimental conditions and interactions.
- **Contrast Computation:** Compute contrast matrices for main effects and interactions from treatment-coded models.
- **Hypothesis Testing:** Execute t-tests and F-tests based on computed contrasts to assess model effects.
- **Condition Estimates:** Generate predicted responses and 95% confidence intervals for each experimental condition.
- **Reporting:** Produce an HTML report summarizing model outputs, contrasts, and figures.

## Repository Contents

- **getTreatWeights.m**  
  Computes contrast matrices and condition weight vectors based on a treatment-coded model.  
  *(See documentation within the file for detailed usage and examples.)*

- **addXterms.m**  
  Augments a design matrix with interaction terms by computing element-wise products of main effect columns.

- **runFCon.m**  
  Computes F contrast test statistics (and corresponding p-values, degrees of freedom, and effect direction) for a specified contrast matrix and fitted model.

- **getXtermStruct.m**  
  Parses coefficient names to determine the structure of interactions, identifying which terms are main effects versus interactions.

- **getXEUL.m**  
  Constructs a design matrix from provided predictors, computes estimated responses, and derives 95% confidence intervals (CIs) from a fitted model.

- **runHCons.m**  
  Applies hypothesis contrasts to a fitted model and augments a table of contrast matrices with test types, statistics, and p-values.

- **getCondEsts.m**  
  Calculates condition estimates and corresponding 95% CIs for each experimental condition based on a fitted model.

- **runExample.m**  
  Provides a complete, integrated example pipeline that demonstrates how to use the toolbox. This script:
  - Generates a synthetic model and factor structure.
  - Computes ANOVA contrast vectors and condition weights.
  - Runs hypothesis tests and computes condition estimates.
  - Plots the effect of a continuous predictor.
  - Generates an HTML report summarizing the analysis.

- **reportResults.m**  
  Generates an HTML report that compiles model summaries, contrast tables, condition estimates, and figures into a single, formatted document.

## Getting Started

### Prerequisites

- **MATLAB:** Ensure that MATLAB is installed. While the toolbox is self-contained, using the Statistics and Machine Learning Toolbox may be beneficial for model fitting.
- **Operating System:** Compatible with Windows, macOS, and Linux.

### Installation

1. Clone or download this repository to your local machine.
2. Add the repository folder to your MATLAB path:
   ```matlab
   addpath('path/to/repository');
   ```
3. TEST
