# Pan-Cancer-Survival-Modeling

This is the repository for code used in the Pan-Cancer Survival Modeling project. You can find our article [here](https://arxiv.org/abs/1910.03447). This repository contains code for (1) data importation and gene set selection, (2) the four models we tried and the calculation of the cross-validated out-of-sample log-posterior likelihoods to compare them, (3) the forward selection procedure we employed, (4) the survival plots we created, and (5) the simulation we used to validate our in-house Gibbs sampler.

# Description of each file 
(1) ImportData.R: This file contains code for importing data through the TCGA2STAT pipeline and how we selected genes for our study. One will need to save the data from this code in order to run the rest of the analysis. At the bottom of this script there is a save() statement and at the top of each additional analysis script there is a load() statement to facilitate this. 

(2) GeneCorrelations.R: This file contains code for computing the correlation between mutation burdens of each gene across all cancer types. Code for plotting the correlations based on significance is also included in this script. 

(3) The four models we tried were the log-normal, normal, exponential, and Weibull. The corresponding scripts for these models are as follows: PL5FCVLogNorm50Genes.R, PL5FCVNorm50Genes.R, JagsModelExponential.R, JagsModelWeibull.R. Based on the results in these, we used the log-normal model to continue our analysis. 

(4) To do forward selection, first we ran the forward selection procedure without any genes. This served as our baseline. The code for this can be found in ForwardSelectionNoGenes.R. Once this has been run, then one can move onto comparing single-gene models, two-gene models, and so on. 

ForwardSelection_5FoldCV.R contains the forward selection procedure. Running parallel loops in R with the foreach package requires the user to import everything in the outside environment into the parallel environment which is a bit tricky. To avoid this, we chose not write the entire forward selection procedure in a function. This also encourages the user to come back and check every once in awhile where the procedure is at and see how the results are progressing. For this reason, in order to run the one-gene model and the two-gene model and so on, copy-and-paste the chunk of code between the ### comments. Every line with a ### between these two comments must be changed at each step of the forward selection procedure. 

PL5CVLogNormEstCoefByCancerSeparately.R contains the code for estimating the posteriors for each parameter by cancer type separately (without borrowing information across cancer types). 

(5) Survival_Curves_ACC.R: This script contains code we used to create the survival curves for a patient with ACC depending on various ages and combinations of mutation statuses. The script imports the posteriors from fitting the final model, the FinalModelFromForwardSelectionFAT4_TP53.rda file of which can be found in this repository. This will load in a variable called "FinalModelPosteriors."

(6) SimulateCoverage.R: This is the simulation we used to assess the coverage rates of our in-house Gibbs sampler. We generate 1000 datasets based on fixed true values of the parameters, calculate the posteriors for the parameters using our model, calculate the credible intervals for each, and then check if the true value of each parameter used to generate the data is contained in the credible interval. 



