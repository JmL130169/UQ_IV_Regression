# Uncertainty Quantification in Instrumental Variables Estimation
This is the R code for numerical analysis in the paper "Uncertainty Quantification in Instrumental Variables Estimation", which discusses the influence of uncertainty in Stage 1 estimators that on the estimation of effect size in two-stage least square (2SLS). To replicate the simulation results, readers need to set their own working direction and input the corresponding $h_{2}^{p}, h_{2}^{e}$ values as stated in the main paper.

## Dataset
### Description
In the direction "data", we have the data of simulation results for each of 10 scenarios described in the main paper. As for the data we used in this work, it is from UK Biobank dataset. The individual data can be accessed from the UK Biobank study at https://www.ukbiobank.ac.uk.

### License
The UK Biobank dataset is provided under a specific set of terms and conditions which all users must adhere to. By using this dataset, you agree to the following terms:

- The data can only be used for health-related research that is in the public interest.
- Data access is granted under the terms of the Material Transfer Agreement (MTA) and the specific Data Access Application (DAA) that was approved.
- The dataset cannot be shared with third parties or used for commercial purposes without explicit permission from UK Biobank.

For more detailed information on the terms of use and licensing, please refer to the UK Biobank's official documentation:
- [UK Biobank Terms and Conditions](https://www.ukbiobank.ac.uk/enable-your-research/approved-research/legal-agreements)
- [Material Transfer Agreement (MTA)](https://www.ukbiobank.ac.uk/media/g0jjlycd/uk-biobank-material-transfer-agreement.pdf)
- [Data Access Application (DAA)](https://www.ukbiobank.ac.uk/enable-your-research/approved-research/register)

### Acknowledgement
Please ensure to acknowledge UK Biobank in any publications or presentations that use this dataset as follows:
"This research has been conducted using the UK Biobank Resource under Application Number [your application number]."

For more detailed information on how to acknowledge UK Biobank, please refer to their [publication guidelines](https://www.ukbiobank.ac.uk/media/uz5f3zpp/publications-policy.pdf).

## Code
### Description
In the direction "code", we have:
* "sim.R" is the main code of our simulation.
* "ComPACT_support_V20.R" is some necessary functions used in "sim.R"
* "sub.lsf" is the submission file to server, readers can ignore it.
* "powerplot_gwas.R" is the code to generate power plots.
* "bias_coverage_gwas.R" is the code to generate MAB, MSE and coverage probability plots.
* "Fig2c_qq_plot_type1.R" is the code to generate QQ-plot.

### Dependencies
- R/4.3.1

### Usage
All analyses were performed on our server, utilizing parallel computing to efficiently process the simulations. Each simulation was allocated a single CPU core and 4 GB of memory, with jobs distributed across 500 cores. The submission file is "sub.lsf" in direction "code".

