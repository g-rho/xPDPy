eXplAInability
=============================================================================

Numerical and visual tools to analyse the explainability of partial dependence functions.

Partial dependence plots are a popular tool to analyze black box machine learning models. 
The repo provides an R package that computes the *explainability* of a PDP, which is a measure to quantify how far a PDP is able to explain a model.

Supported functionalities:
- ...computation of explainability,
- ...matchplot of PDP vs. model predictions,
- ...computation of a forward variable selection based on explainability,
- ...visualization of 2D PDP vs. unexplained residual predictions,
- ...scatterplot matrix of 2D partial dependence plots.  

Details are described in this [paper](https://arxiv.org/pdf/1910.13376.pdf). The examples are taken from the paper.  