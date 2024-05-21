# Uncertainty Quantification in Instrumental Variables Estimation
This is the R code for numerical analysis in the paper "Uncertainty Quantification in Instrumental Variables Estimation". To replicate the simulation results, readers need to set their own working direction and input the corresponding $h_{2}^{p}, h_{2}^{e}$ values as stated in the main paper.

In the direction "code", we have:
* "sim.R" is the main code of our simulation.
* "ComPACT_support_V20.R" is some necessary functions used in "sim.R"
* "sub.lsf" is the submission file to server, readers can ignore it.
* "powerplot_gwas.R" is the code to generate power plots.
* "bias_coverage_gwas.R" is the code to generate MAB, MSE and coverage probability plots.
* "Fig2c_qq_plot_type1.R" is the code to generate QQ-plot.

In the direction "data", we have the data of simulation results for each of 10 scenarios described in the main paper. As for the data we used in this work, it is from UK Biobank dataset, and the individual data can be accessed from the UK Biobank study at https://www.ukbiobank.ac.uk.
