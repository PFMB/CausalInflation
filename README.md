### Repository for the Reproduction of *Estimating the Effect of Central Bank Independence on Inflation Using Longitudinal Targeted Maximum Likelihood Estimation*

This repository is provided for the reproduction of the results of the research paper "Estimating the Effect of Central Bank Independence on Inflation Using Longitudinal Targeted Maximum Likelihood Estimation". It contains all relevant data and code for the main results. **Please install the following packages** prior to running the files listed below: ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, magrittr, randomForest, earth, gbm, gam , mgcv, reshape2, dplyr, data.table, plyr, tibble, scales, polspline, BaBooN, simcausal (currently only on CRAN archive). 

The repository cotains the following files:

* `Simulation1.R` and `Simulation2.R` reproduce the simulation studies presented in the paper

* `EstimateData.R` reproduces the data analysis

* `CalcBiasCP.R` takes the output from the run on the cluster and calculates the abs. bias and the coverage probabilites reported in the paper

* `NodesFormInterv.RData` contains g/Q-formula and further arguments needed for the specification of `ltmle::ltmle`

* `WeightsSummary.R` contains functions that extract the weights of the learners which are used during the run of the `SuperLearner`

* `LearnerLibrary.R` and `SelectedLearners.RData` contain customized sets of learners that are used for the estimations made by `EstimateData.R`, `Simulation1.R` and `Simulation2.R`

* `InflData.RData` contains the 5 imputed data sets that were used in `EstimateData.R`

