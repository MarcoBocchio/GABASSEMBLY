# GABASSEMBLY
Pipeline to analyse 2p calcium imaging data from mouse hippocampal cells.

## Structure
### 1) PreProcessing.m
a)  Movement correction using NoRMCorre code (Pnevmatikakis & Giovannucci, 2017; https://doi.org/10.1016/j.jneumeth.2017.07.031)
b)  Semi-automatic ROI detection, denoising, demixing and spike inference using the CaImAn pipeline (Giovannucci et al., 2019; https://doi.org/10.7554/eLife.38173)

### 2) mvmSCE.m
Analysis of inferred spike activity in relation to movement and Synchronous Calcium Events (SCEs)

### 3) Assemblies.m
Assembly detection with the following methods
a) seqAssembly: set of codes for k-means clustering of SCEs and cell assignment to significant clusters (Malvache et al., 2016)
b) SGCassembly: similarity graph clustering (Avitan et al., 2017, Moelner et al., 2018; https://doi.org/10.1186/s12915-018-0606-4)
c) ICAssembly: PCA/ICA method (Lopes Dos Santos et al., 2013; https://doi.org/10.1186/s12915-018-0606-4)
d) Cross-correlation of single cell activities to assembly activities detected with PCA/ICA method

## Required third-party codes (all included)
- NoRMCorre (https://github.com/flatironinstitute/NoRMCorre
- CaImAn (https://github.com/flatironinstitute/CaImAn-MATLAB)
- plot_staggered (Marco Bilucaglia)
- sepblockfun (Matt Jacobson)
- abfload (F. Collman) https://github.com/fcollman/abfload
- toolbox (PCA/ICA method, Vitor Lopes-Dos Santos, vtlsantos@gmail.com)
- SCEkmAssembly (k-means method, Arnaud Malvache/Yannick Bollmann)
- SGCAssembly (similarity graph clustering method, Jan Moelter, part of the neural assembly detection package: https://github.com/GoodhillLab/neural-assembly-detection)

Marco Bocchio, updated 25/11/19

