### Baumann, Schomaker, Rossi: __Estimating the Effect of Central Bank Independence on Inflation Using Longitudinal Targeted Maximum Likelihood Estimation__, Journal of Causal Inference, 2021, in press

The code and data in this repository can be used freely given that attribution to the authors is given. Please cite our paper if you use the data or code. Send an email to `baumann@kof.ethz.ch` to obtain the password for the `data.zip` folder. Use the subdirectory `/Users/flipst3r/RStHomeDir/GitHub/CausalInflation/data` from the unzipped `data.zip` and set the paths accordingly to reproduce the results. **Please install the following packages** prior to running the files listed below: ltmle, SuperLearner, vcd, arm, rpart, nnet, glmnet, magrittr, randomForest, ggsci, metR, cowplot, xtable, grid, gridextra, earth, gbm, gam , mgcv, reshape2, dplyr, data.table, plyr, mFilter, RColorBrewer, tibble, scales, polspline, BaBooN, simcausal (currently only on CRAN archive). 

The repository cotains the following files:

* `data/CausInfl.RData` contains the 5 imputed data sets that were used in `code/33EstimateData.R` which contains the main analyses. It is composed of the merged and processed data sets `3_infl_XYI_imp.RData` and `5_infl_Z.RData` which were imputed and processed in a different repository. Merging and preprocessing is conducted in `code/1PrepareData.R`

* `code/2PrepareLTMLE.R` holds the specifiations for every run of `ltmle::ltmle`. In particular, `g-form`, `Q-form` for every interventon and the coresponding `L/A/Y-nodes`.

* `code/34CBITest.R` performs an unpaired Welch Two Sample t-test.

* `code/31WeightsSummary.R` contains functions that extract the weights of the learners which are used during the run of the `SuperLearner`.

* `data/CausalOrder.csv` contains the topological ordering that has to be assigned based on the DAG presented in the manuscript.

* `data/NodesFormInterv.RData` contains g/Q-formula and further arguments needed for the specification of `ltmle::ltmle`.

* `code/4DescriptivePlots.R` reproduces the plots and tables presented in the main text and the appendix of the paper.

* `code/32LearnerLibrary.R` and `data/SelectedLearners.RData` contain customized sets of learners that are used for the estimations made by `code/33EstimateData.R`, `code/Simulation1.R` and `code/Simulation2.R`.

* `code/5Simulation1.R` and `code/6Simulation2.R` reproduce the simulation studies presented in the paper.

* `code/7CalcBiasCP.R` takes the output from the run on the cluster and calculates the abs. bias and the coverage probabilites reported in the paper.

* `/plots` holds the plots presented in the paper and `/results` contains the results from the [HPC](https://scicomp.ethz.ch/wiki/Euler) run.

* `results/causal_log` holds the details (i.e. p2s) of the run of `33EstimateData.R`.

The run for the results in the paper was conducted on the **8/9th of March, 2021** and is based on an imputation run from **8th of March** containing 60 countries. The file `code/33EstimateData.R`, which contains the main analyses, was run with:

`module load new gcc/4.8.2 r/3.6.0`
`bsub -n 5 -W 24:00 -N -B -R "rusage[mem=16384]" 'R --vanilla --slave < $HOME/CausalInflation/code/33EstimateData.R > causal_log'`
