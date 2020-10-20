# Double-robust estimation in longitudinal data with time-varying confounding: a simulation study

This is a repository of code related to the work-in-progress above, nicknamed "drsim."  

Three programs are provided under drsim/code:
  1. drsim_ltmle.R : Generates simulations and analyzes data using Longitudinal Targeted Minimum Loss-Based Estimation (LTMLE).
  2. drsim_ipw.R : Generates simulations and analyzes data using inverse probability weighting (IPW).
  3. drsim_truth.R : Generates 2 large simulations, identical except that in the first exposure is set to and in the second exposure is set to 0, and obtains the true risk differences.
  
We further provide the following files under drsim/results:
  1. drsim_truth.txt : The output of drsim_truth.R, a file containing the true risk difference by maximum number of time points included.
