### Repository for the Reproduction of *Estimating the Effect of Central Bank Independence on Inflation Using Longitudinal Targeted Maximum Likelihood Estimation*

This repository is provided for the reproduction of the results of the research paper ["Estimating the Effect of Central Bank Independence on Inflation Using Longitudinal Targeted Maximum Likelihood Estimation"](https://arxiv.org/abs/2003.02208). It contains all relevant data and code for the main results. **Please install the following packages** prior to running the files listed below: ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, magrittr, randomForest, earth, gbm, gam , mgcv, reshape2, dplyr, data.table, plyr, tibble, scales, polspline, BaBooN, simcausal (currently only on CRAN archive). 

The repository cotains the following files:

* `Simulation1.R` and `Simulation2.R` reproduce the simulation studies presented in the paper.

* `EstimateData.R` reproduces the data analysis.

* `PrepareLTMLE.R` holds the specifiations for every run of `ltmle::ltmle`. In particular, `g-form`, `Q-form` for every interventon and the coresponding `L/A/Y-nodes`.

* `CalcBiasCP.R` takes the output from the run on the cluster and calculates the abs. bias and the coverage probabilites reported in the paper.

* `NodesFormInterv.RData` contains g/Q-formula and further arguments needed for the specification of `ltmle::ltmle`.

* `WeightsSummary.R` contains functions that extract the weights of the learners which are used during the run of the `SuperLearner`.

* `LearnerLibrary.R` and `SelectedLearners.RData` contain customized sets of learners that are used for the estimations made by `EstimateData.R`, `Simulation1.R` and `Simulation2.R`.

* `causalinfl.RData` contains the 5 imputed data sets that were used in `EstimateData.R`. It is composed of the merged and processed data sets `infl_XYI_imputed.RData` and `infl_Z.RData`.

* `/plots` holds the plots presented in the paper and `/results` contains the results from the [HPC cluster](https://scicomp.ethz.ch/wiki/Euler) run.

* `CausalOrder.csv` contains the topological ordering that has to be assigned based on the DAG presented in the manuscript.

* `causal_log` holds the details (i.e. p2s) of the run of `EstimateData.R` on the [HPC cluster](https://scicomp.ethz.ch/wiki/Euler).

* `DescriptivePlots.R` reproduces the plots and tables presented in the main text and the appendix of the paper.

The run for the results in the paper was conducted on the 8th of March, 2021 and is based on an imputation run from the very same day containing 60 countries. The file `EstimateData.R` was run with `bsub -n 5 -W 24:00 -N -B -R "rusage[mem=16384]" 'R --vanilla --slave` `< $HOME/CausalInflation/EstimateData.R > causal_log'`.
