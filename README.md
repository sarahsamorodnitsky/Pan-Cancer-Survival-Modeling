# Pan-Cancer-Survival-Modeling

This is the repository with code used in the Pan-Cancer Survival Modeling project, titled "A Pan-Cancer and Polygenic Bayesian Hierarchical Model for the Effect of Somatic Mutations on Survival." In this repository, one will find code for (1) data importation and gene set selection, (2) the four models we tried and the calculation of the cross-validated out-of-sample log-posterior likelihoods to compare them, (3) the forward selection procedure we employed, (4) the survival plots we created, and (5) the simulation we used to validate our in-house Gibbs sampler.

# Description of each file 
(1) ImportData.R: In this file, one can find code on importing data through the TCGA2STAT pipeline and how we selected genes for our study. One will need to save the data from this code in order to run the rest of the analysis. At the bottom of this script there is a save() statement and at the top of each additional analysis script there is a load() statement to facilitate this. 

(2) The four models we tried were the log-normal, normal, exponential, and Weibull. The corresponding scripts for these models are as follows: PL5FCVLogNorm50Genes.R, PL5FCVNorm50Genes.R, JagsModelExponential.R, JagsModelWeibull.R. 



