# Extensions of Health Economic Evaluation in R for Microsoft Excel Users: A Tutorial for Incorporating Heterogeneity and Conducting Value of Information Analyses

> This is the repo for the journal paper _Extensions of Health Economic Evaluation in R for Microsoft Excel Users: A Tutorial for Incorporating Heterogeneity and Conducting Value of Information Analyses_ by 
Nichola R. Naylor & Jack Williams (joint first author), Nathan Green, Felicity Lamrock, Andrew Briggs.

## Abstract

Advanced health economic analysis techniques currently performed in Excel, such as incorporating heterogeneity, time-dependent transitions and value of information analysis, can be easily transferred to R. Whilst previous R tutorials for health economists have focused on setting up Markov models and probabilistic sensitivity analyses, there is a need for more detailed exploration and explanation of how to develop more advanced Markov models in R. This tutorial provides a step-by-step guide of how to incorporate heterogeneity and value of information with Markov models in R, through side-by-side comparisons of how such analyses are done in Excel. Though there are papers highlighting the efficiency benefits of R in comparison to Excel, and general tutorials of how to perform health economic analyses in R, we hope that this paper can act as a reference point to those switching more complex models from Excel to R. We provide open-access code and data, suitable for future adaptation.

## Contents

File | Description
-----|------------
[THR_Model.R](https://github.com/Excel-R-tutorials/Markov_Extensions/blob/main/THR_Model.R) | R script of the main model functions
[THR_Model_VOI.R](https://github.com/Excel-R-tutorials/Markov_Extensions/blob/main/THR_Model_VOI.R) | R script containing the value of information analysis functions
[hazardfunction.csv](https://github.com/Excel-R-tutorials/Markov_Extensions/blob/main/hazardfunction.csv) | CSV file containing the output of the survival analysis
[cov55.csv](https://github.com/Excel-R-tutorials/Markov_Extensions/blob/main/cov55.csv) | CSV file containing the corresponding covariance matrix of the suvival analysis
[life-table.csv](https://github.com/Excel-R-tutorials/Markov_Extensions/blob/main/life-table.csv) | CSV file containing the life-table data

## How to Cite this Code

Nichola R. Naylor* & Jack Williams*, Nathan Green, Felicity Lamrock, Andrew Briggs (2022) Extensions of Health Economic Evaluation in R for Microsoft Excel Users: A Tutorial for Incorporating Heterogeneity and Conducting Value of Information Analyses. GitHub (https://github.com/Excel-R-tutorials/Markov_Extensions)

*joint first authors that wrote the majority of this code.

## 👂 Feedback

Please feel free to raise an issue on this GitHub, or email any feedback to chilgithub@lshtm.ac.uk, either Nichola or Jack can then respond to any queries.
