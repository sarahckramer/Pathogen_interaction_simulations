# Statistical methods for inferring pathogen-pathogen interactions

Code to test 1) Pearson correlation coefficients, 2) generalized additive models (GAMs), 3) Granger causality, 4) transfer entropy, and 5) convergent cross-mapping (CCM) as methods to determine whether or not an interaction exists between two pathogens

Directory Structure
-------------------
* src
    * 01_dependencies: Code required for running all analyses
    * 01_methods: Code files for each of the individual methods tested
    * 02_dependencies: Code required for processing results
    * sensitivity: Code to run and process results from various sensitivity analyses

Dependencies
------------

All code was run in R version 4.4.0. The following packages are used by this repository:

* 

All models are run using the package "pomp," which has additional dependencies, depending on your operating system. Detailed installation instructions can be found [here](https://kingaa.github.io/pomp/install.html).

Additionally, in order to calculate transfer entropy, users will need to install both Java, and the Java Information Dynamics Toolkit. Installation instructions can be found [here](https://github.com/jlizier/jidt/tree/master).

Running Analysis
----------------

The full analysis is conducted by running 01_run_analyses.R, which generates synthetic data and performs all five statistical methods in turn. I ran the correlation coefficients (because of their quick run time), granger causality (similar), and transfer entropy (because of the multiple dependencies) analyses locally, but ran the other methods on a high performance computing cluster. To run remotely, the code should be run once each for jobid 1-18, so that all combinations of interaction strength and duration are covered; run_local should be set as FALSE. For local running, there is code at the top of 01_run_analyses.R that can be uncommented and pasted into the console.

Various sensitivity analyses can be run by changing the variable sens in 01_run_analyses.R. Options for the sensitivity analyses can be seen in generate_data.R, lines 117-154. To run sensitivity/check_ccm_alpha.R, additional sensitivity analyses were conducted by changing the value of alpha in 01_methods/CCM.R lines 191-192.

After all analyses have been run, results can be formatted and plotted using 02_process_results.R (or sensitivity/process_sens_results.R, for sensitivity analyses).
