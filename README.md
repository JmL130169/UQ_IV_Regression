# Uncertainty Quantification in Instrumental Variables Estimation
This is the R code for numerical analysis in the paper "Uncertainty Quantification in Instrumental Variables Estimation". To replicate the simulation results, readers need to set their own working direction and input the corresponding $h_{2}^{p}, h_{2}^{e}$ values as stated in the main paper.

In the direction "code", we have:
* "sim.R" is the main code of our simulation.
* "ComPACT_support_V20.R" is some necessary functions used in "sim.R"
* "sub.lsf" is the submission file to server, readers can ignore it.

In the direction "data", we have:
* "MAF011_VIM" is the genotyoe data of gene VIM, with MAF > 0.01.
* "1000-Genome" is the reference panel, where we do not use it in this project but only use it to select common SNPs between genotype data and reference panel, in case we need to deal with summary data in the future.
