# Causal Framework for Quantiles of Residual Lifetime

This repository provides the R implementation for the paper: 
**"A Causal Framework for Estimating Quantiles of Residual Lifetime conditional on Surviving a Landmark Time"**.

The repository contains functions to estimate the **Observed Survivor Quantile Contrast (OSQC)** using a doubly robust (DR) approach and the **Principal Survivor Quantile Contrast (PSQC)** via reweighting.

---

## 🚀 Key Features
* **Doubly Robust Estimation:** Consistent quantile estimation if either the treatment or outcome model is correct.
* **Informative Censoring:** Integrated IPCW (Inverse Probability of Censoring Weighting).
* **Principal Stratification:** Supplementary tools to disentangle compositional selection from causal efficacy.
* **Bootstrap Inference:** Nonparametric bootstrap for variance estimation of non-smooth quantile functions.
  
---

## 📂 Folder Overview
* `/simulations`: Scripts to replicate the bias, SE, and CP results from the paper.
* `/applications`: Analysis scripts for the SUPPORT and NSABP B-14 datasets. Note that B-14 dataset is not publicaly available.

---

## Reference
If you use this code, please cite the following article:  

Hong, T., Bae, W., LEE, SK., Choi, D. & Jeong, JH. (2026). A causal framework for quantile residual lifetime. Arxiv
