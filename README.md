# Masked Bounded Nuclear Norm Regularization (MBNNR)

## Introduction
This repository provides the implementation of Masked Bounded Nuclear Norm Regularization (maskBNNR) for predicting post-translational modification (PTM)-disease associations.  
The framework integrates heterogeneous biological similarities and applies an external masking mechanism to distinguish known from unknown associations, ensuring unbiased evaluation.

## Code execution
You can run `MBNNR_5_fold.m` to perform 100 rounds of five-fold cross-validation.  
You can run `plot_models_evaluation.m` to get evaluation metrics for MBNNR, including AUC, PR_AUC, F1, and ACC.  
In `predict.m`, you can predict PTM (protein) related to the disease you are interested in.

## Requirements
- MATLAB version 2022b

