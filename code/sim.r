# Working directory
setwd("") # Set your own working direction

###########
# Package #
###########

# Package
library(BEDMatrix)
library(compiler)
suppressMessages(library(data.table))
suppressMessages(library(ddpcr))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
library(boot)
library(caret)
library(lassosum)
library(Matrix)
library(optparse)
library(psych)
library(parallel)
library(Rcpp)
library(stringr)

# Option
option_list <- list(
make_option("--h2_e", type = "numeric", default = NA, action = "store"
),
make_option("--h2_p", type = "numeric", default = NA, action = "store"
),
make_option("--causal", type = "numeric", default = NA, action = "store"
),
make_option("--name_batch", type = "character", default = NA, action = "store"
),
make_option("--gene", type = "character", default = NA, action = "store"
),
make_option("--seed", type = "numeric", default = NA, action = "store"
),
make_option("--runs_association", type = "numeric", default = 20, action = "store"
)
)

opt <- parse_args(OptionParser(option_list = option_list))

h2.e       <- opt$h2_e
h2.p       <- opt$h2_p
causal     <- opt$causal
name.batch <- opt$name_batch
gene       <- opt$gene
seed       <- opt$seed
runs       <- opt$runs_association

out.id = as.numeric(gsub("sim","",name.batch))

# Make folders
paste0("mem/", name.batch) %>% dir.create(., showWarnings = FALSE, recursive = TRUE)

# My R functions
source("ComPACT_support_V20.r")

# Elastic Net
weights.enet <- function(genos, pheno, alpha = 0.5) {
  
    eff.wgt = matrix(0 , ncol = 1 , nrow = ncol(genos))
    # remove monomorphics
    sds = apply( genos , 2 , sd )
    keep = sds != 0 & !is.na(sds)
    enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
    
    # Added: exclude models with n.snp < 10
    enet$cvm[enet$nzero <= 10] <- 99999999 
    enet$lambda.min <- enet$lambda[which.min(enet$cvm)]
    
    eff.wgt[keep] = coef(enet , s = "lambda.min")[2 : (sum(keep)+1)]
    # return(eff.wgt)
    return(list(weight = eff.wgt, lam = enet$lambda.min))
  
}

weights.lasso <- function(genos, pheno, alpha = 1) {
  
  eff.wgt = matrix(0 , ncol = 1 , nrow = ncol(genos))
  # remove monomorphics
  sds = apply( genos , 2 , sd )
  keep = sds != 0 & !is.na(sds)
  enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
  
  # Added: exclude models with n.snp < 10
  enet$cvm[enet$nzero <= 10] <- 99999999 
  enet$lambda.min <- enet$lambda[which.min(enet$cvm)]
  
  eff.wgt[keep] = coef(enet , s = "lambda.min")[2 : (sum(keep)+1)]
  # return(eff.wgt)
  return(list(weight = eff.wgt, lam = enet$lambda.min))
  
}

#################################
# Exterior executables and data #
#################################

path.Plink <- "executables/plink"

################
# Prepare data #
################

index.train <- 1:10000
index.test  <- (10000 + 1):(10000 + 10000)
index.add   <- (20000 + 1):(20000 + 10000)
index.ref   <- (30000 + 1):(30000 + 500)

N <- 10000 + 10000 + 10000 + 500

index.train.0 <- index.train

# Load the processed genotype data
# Keep only the SNPs with MAF > 0.01

bim <- paste0("MAF001_", gene, ".bim") %>% fread(., data.table = FALSE)

genotype = BEDMatrix(paste0("MAF001_", gene), simple_names = TRUE)

genotype = genotype[1:N, ]

chromosome = bim[1, 1]

tmpbim = fread(paste0("1000-Genome/1000G.EUR.ALLSNP.QC.CHR", chromosome,".bim")) %>% as.data.frame()

# only used SNPs in 1000G
genotype = genotype[, colnames(genotype) %in% tmpbim$V2]
bim = bim[bim$V2 %in% tmpbim$V2, ]

rownames(tmpbim) = tmpbim$V2
tmpbim = tmpbim[bim$V2, ]

# we further remove ambignous and unmapped SNPs to avoid any potential problems
qcinfo = allele.qc(bim$V5, bim$V6, tmpbim$V5, tmpbim$V6)
indx = !qcinfo$flip & qcinfo$keep
bim = bim[indx, ]
genotype = genotype[, colnames(genotype) %in% bim$V2]

genotype <- genotype %>% PatchUp() %>% scale()
size <- ncol(genotype) # number of SNPs

##################
# Main iteration #
##################

for (i in 1:20) {
  
    #################
    # Simulate data #
    #################
    
    # Set seed
    seed.used = seed * 40 - 40 + i + 10000 * out.id
    set.seed(seed.used)
    
    ###########################
    # Added: new causal rules #
    ###########################
    
    size <- size  # Size of the vector
    min_distance <- 3  # Minimum distance between non-zero weights
    
    # initialization of generated w
    w = numeric(size)

    # Function to check if the proposed index is at least min_distance away from all selected indices
    is_far_enough <- function(selected, candidate, min_distance) {
        all(abs(selected - candidate) >= min_distance)
    }
    
    # Pick SNPs with enforced minimum distance
    picks <- c()
    while(length(picks) < causal) {
        candidate <- sample(1:length(w), 1)
        if (is_far_enough(picks, candidate, min_distance)) {
            picks <- c(picks, candidate)
        }
    }

    # generate IV-expression data, while control the heritability
    w[picks] <- rnorm(length(picks), 0, 1) # Stage-1 true parameters
    var_e <- var(genotype %*% w) * ((1 - h2.e)/(h2.e))
    noise <- rnorm(nrow(genotype), 0, sqrt(var_e))
    exp_true <- genotype %*% w
    exp_obs <- exp_true + noise

    # estimation
    set.seed(seed.used)
    
    genotypetrain = genotype[index.train, ] # %>% scale()
    genotest = genotype[index.add, ]
    exp_train <- exp_obs[index.train]
    
    # 2SLS estimator, Stage-I
    ls.result <- weights.enet(
        genos = genotypetrain,
        pheno = exp_train,
        alpha = 0.5
    )
    
    #########################
    # Correction parameters #
    #########################
    
    # Ground truth
    weight_est <- as.numeric(ls.result$weight)
    pred <- genotest %*% weight_est
    sigma_z <- var(pred)
    
    sigma_zt <- var(genotest %*% w)
    sigma_ze <- cov(genotest %*% w, pred - genotest %*% w)
    sigma_err <- var(pred - genotest %*% w)
    rr <- (sigma_zt + sigma_ze)/sigma_z
    
    # 2SLS, Stage-II
    runT1e = runs * 10
    runPower = runs * 5
    p.t1e <- matrix(0, nrow = 3, ncol = runT1e)
    p.power <- array(0, dim = c(5, runPower, 3))
    bias <- matrix(0, nrow = 3, ncol = 5)
    mse <- matrix(0, nrow = 3, ncol = 5)
    beta_est <- array(0, dim = c(5, runPower, 3))
    lb <- array(0, dim = c(5, runPower, 3))
    ub <- array(0, dim = c(5, runPower, 3))
    cover <- array(0, dim = c(5, runPower, 3))
    
    #######
    # T1E #
    #######
    
    cat("Start T1E simulation\n")
    set.seed(seed.used)
    for (rep in 1:runT1e) {
        
        # estimation of stage-I
        w_est = as.numeric(ls.result$weight)
        
        # Generate beta
        beta = 0
        
        # Generate alpha
        genotype_train = genotype[index.add, ] # %>% scale()
        n_train = length(index.add)
        
        expression <- exp_obs[index.add]
        expression_true <- exp_true[index.add]
        
        # Generate phenotype levels
        # When H0 is true, no need to control the heti, just make sure we have condounfers
        corr_level <- 0.3
        noise_p0 <- rnorm(n_train, 0, sqrt(var(noise)))
        noise_p <- noise[index.add] * corr_level + 
          noise_p0 * sqrt(1 - corr_level^2)
        
        phenotype_train <- expression * beta + noise_p
        
        # Stage-II estimation, OLS
        Y.temp <- phenotype_train
        fit0 <- lm(Y.temp ~ (genotype_train %*% w_est))
        summary_fit0 <- summary(fit0)
        p.t1e[1, rep] <- summary_fit0$coefficients[2, 4]
        
        # Stage-II estimation, truth
        Y.temp <- phenotype_train
        fit1 <- lm(Y.temp ~ expression_true)
        summary_fit1 <- summary(fit1)
        p.t1e[2, rep] <- summary_fit1$coefficients[2, 4]
        
        # Stage-II estimation, bias corrected using truth
        reli_ratio <- rr
        coef_beta <- coef(fit0)
        beta_corr <- coef_beta[2] / reli_ratio # I need this
        res_var <- var(noise_p + beta * noise[index.add])
        V <- ((sigma_zt^2 * sigma_err) - sigma_ze^2) / ((sigma_z)^2)
        beta_var <- (((beta)^2) * V +
          res_var / (sigma_z)) / n_train
        beta_var_corr <- beta_var / (reli_ratio ^ 2) # I need this
        p.t1e[3, rep] <- 2 * (pnorm(abs(beta_corr / sqrt(beta_var_corr)), lower.tail = FALSE))
        
    }
    
    #########
    # Power #
    #########
    
    cat("Start Power simulation\n")
    set.seed(seed.used)
    beta = c(-0.1, -0.05, 0, 0.05, 0.1)
    size_heti <- length(beta)
    
    for (rep in 1:runPower) {
        
        for (pw in 1:size_heti) {
          
            w_est = as.numeric(ls.result$weight)
            
            # Generate alpha
            genotype_train = genotype[index.add, ] # %>% scale()
            n_train = length(index.add)
            
            expression <- exp_obs[index.add]
            expression_true <- exp_true[index.add]
            
            # Generate phenotype levels
            # When H0 is true, no need to control the heti, just make sure we have condounfers
            corr_level <- 0.3
            var_level <- var(expression * beta[1]) * ((1 - h2.p) / h2.p)
            noise_p0 <- rnorm(n_train, 0, var(noise))
            noise_p1 <- noise[index.add] * corr_level + 
              noise_p0 * sqrt(1 - corr_level^2)
            noise_p <- (noise_p1 / sqrt(var(noise_p1))) * sqrt(var_level)
            
            phenotype_train <- expression * beta[pw] + noise_p
            
            # Stage-II estimation, OLS
            Y.temp <- phenotype_train
            fit0 <- lm(Y.temp ~ (genotype_train %*% w_est))
            summary_fit0 <- summary(fit0)
            p.power[pw, rep, 1] <- summary_fit0$coefficients[2, 4]
            beta_est[pw, rep, 1] <- summary_fit0$coefficients[2, 1]
            CI <- confint(fit0, level = 0.95)
            lb[pw, rep, 1] <- CI[2, 1]
            ub[pw, rep, 1] <- CI[2, 2]
            
            if (beta[pw] < CI[2, 2] && beta[pw] > CI[2, 1]) {
                cover[pw, rep, 1] <- 1
            } else {
                cover[pw, rep, 1] <- 0
            }
                    
            # Stage-II estimation, truth
            Y.temp <- phenotype_train
            fit1 <- lm(Y.temp ~ expression_true)
            summary_fit1 <- summary(fit1)
            p.power[pw, rep, 2] <- summary_fit1$coefficients[2, 4]
            beta_est[pw, rep, 2] <- summary_fit1$coefficients[2, 1]
            CI <- confint(fit1, level = 0.95)
            lb[pw, rep, 2] <- CI[2, 1]
            ub[pw, rep, 2] <- CI[2, 2]
            
            if (beta[pw] < CI[2, 2] && beta[pw] > CI[2, 1]) {
                cover[pw, rep, 2] <- 1
            } else {
                cover[pw, rep, 2] <- 0
            }
            
            # Stage-II estimation, bias corrected using truth
            reli_ratio <- rr
            coef_beta <- coef(fit0)
            beta_corr <- coef_beta[2] / reli_ratio # I need this
            res_var <- var(noise_p + beta[pw] * noise[index.add])
            V <- ((sigma_zt * sigma_err) - sigma_ze^2) / ((sigma_z)^2)
            beta_var <- (((beta[pw])^2) * V +
                           res_var / (sigma_z)) / n_train
            beta_var_corr <- beta_var / (reli_ratio ^ 2)
            
            p.power[pw, rep, 3] <- 2 * (pnorm(abs(beta_corr / sqrt(beta_var_corr)), lower.tail = FALSE))
            beta_est[pw, rep, 3] <- beta_corr
            lb[pw, rep, 3] <- beta_corr - qnorm(0.975) * sqrt(beta_var_corr)
            ub[pw, rep, 3] <- beta_corr + qnorm(0.975) * sqrt(beta_var_corr)
            
            if (beta[pw] < beta_corr + qnorm(0.975) * sqrt(beta_var_corr) && 
                beta[pw] > beta_corr - qnorm(0.975) * sqrt(beta_var_corr)) {
                cover[pw, rep, 3] <- 1
            } else {
                cover[pw, rep, 3] <- 0
            }

        }
              
    }
    
    for (pw in 1:size_heti) {
        
        bias[1, pw] <- mean(abs(beta_est[pw, , 1] - beta[pw]))
        bias[2, pw] <- mean(abs(beta_est[pw, , 2] - beta[pw]))
        bias[3, pw] <- mean(abs(beta_est[pw, , 3] - beta[pw]))
        
        mse[1, pw] <- mean((beta_est[pw, , 1] - beta[pw])^2)
        mse[2, pw] <- mean((beta_est[pw, , 2] - beta[pw])^2)
        mse[3, pw] <- mean((beta_est[pw, , 3] - beta[pw])^2)
        
    }
    
    # Power calculation; those also quick; save results in a similiar format such that the previous codes can be used
    cat("Finish !\n")
    
    # save the results
    save(
        p.t1e,
        p.power,
        cover,
        bias,
        lb,
        ub,
        beta_est,
        mse, 
        file = paste0("mem/", name.batch, "/", seed.used, "-", length(index.train), ".RData")
        
    )
    
}

