# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })
    
    return(M)
}

PatchUp <- cmpfun(PatchUp)

# Standardize
Standardize <- function(M) {
    # Centralize
    M <- M - matrix(rep(colMeans(M), times = nrow(M)), nrow = nrow(M) , ncol = ncol(M), byrow = T)
    
    # Standardize
    M <- sweep(M, 2, sqrt(apply(M, 2, crossprod) / nrow(M)), "/")
    
    return(M)
}

Standardize <- cmpfun(Standardize)



# allele.qc
allele.qc <- function(a1, a2, ref1, ref2) {
    ref <- ref1
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip1 <- flip
    ref <- ref2
    flip <- ref
    flip[ref == "A"] <- "T"
    flip[ref == "T"] <- "A"
    flip[ref == "G"] <- "C"
    flip[ref == "C"] <- "G"
    flip2 <- flip
    snp <- list()
    snp[["keep"]] <- !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
    snp[["flip"]] <- (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
    return(snp)
}


SetUpJobs <- function(n.jobs) {
    null.object <- NULL
    if (!dir.exists(paste0(name.batch, "_quick-submit"))) {
        dir.create(paste0(name.batch, "_quick-submit"))
    }
    
    if (!dir.exists(paste0(name.batch, "_quick-submit/to-do"))) {
        dir.create(paste0(name.batch, "_quick-submit/to-do"))
        dir.create(paste0(name.batch, "_quick-submit/working"))
        dir.create(paste0(name.batch, "_quick-submit/done"))
        dir.create(paste0(name.batch, "_quick-submit/error"))
        for (i in 1:n.jobs) {
            save(null.object, file = paste0(name.batch, "_quick-submit/to-do/", i))
        }
    }
    
    invisible(NULL)
}

SetUpJobs <- cmpfun(SetUpJobs)

GetJob <- function() {
    null.object <- NULL
    
    if (length(dir(paste0(name.batch, "_quick-submit/to-do/"))) == 0) {
        print("All jobs are done!")
        
        return(NULL)
    } else {
        to_do <- as.numeric(dir(paste0(name.batch, "_quick-submit/to-do/")))
        
        if (length(to_do) == 1) {
            i <- to_do[1]
        } else {
            i <- sample(to_do, 1)
        }
        
        save(null.object, file = paste0(name.batch, "_quick-submit/working/", i))
        file.remove(paste0(name.batch, "_quick-submit/to-do/", i))
    }
    
    print(paste0("Start working on job ", i, "!"))
    
    return(i)
}

GetJob <- cmpfun(GetJob)

FinishJob <- function(i) {
    null.object <- NULL
    
    file.remove(paste0(name.batch, "_quick-submit/working/", i))
    save(null.object, file = paste0(name.batch, "_quick-submit/done/", i))
    
    invisible(NULL)
}

FinishJob <- cmpfun(FinishJob)


RecordError <- function(i) {
    null.object <- NULL
    
    file.remove(paste0(name.batch, "_quick-submit/working/", i))
    save(null.object, file = paste0(name.batch, "_quick-submit/error/", i))
    
    invisible(NULL)
}

RecordError <- cmpfun(RecordError)


PUMAS.II.FUN.I <- function(GWAS, Nt, LD, LD.block, trait_name, k, chr){
    
    #GWAS=matched_data$gwas_matched
    #Nt=floor(partitions[2]*min(matched_data$gwas_matched$N))
    #LD = matched_data$new_ld_blocks
    #LD.block=matched_data$ld_block_size
    #trait_name=trait_name
    #k=k
    #chr=chr
    
    LD = LD[[1]]
    # parameter initialization
    N.samplesize <- GWAS$N
    Ntr <- N.samplesize - Nt
    
    #MAF <- GWAS$MAF
    X.VAR <- rep(1,dim(GWAS)[1]) #2*MAF*(1-MAF)
    
    SE <- GWAS$SE
    
    beta <- GWAS$BETA # this is standardized beta (for standardized X and Y)
    
    # statistics initialization
    Var.Y <- 1 #quantile(SE^2*N.samplesize*X.VAR,probs=seq(0,1,0.1))[10]
    XtY <- beta * (N.samplesize - 1) # Note that X'X is N-1 not N; both X and Y are standardized
    XtY_mu <- (Nt/ N.samplesize) * XtY
    
    # sampling covariance matrix for XtY_t
    #project_mat <- mclapply(1:length(LD.block), function(j){
    j = 1
    SE.tmp <- SE
    X.VAR.tmp <- X.VAR
    X.VAR.mat.tmp <- sqrt(X.VAR.tmp %*% t(X.VAR.tmp))
    N.samplesize.tmp <- N.samplesize
    NminusNv.tmp <- N.samplesize.tmp - Nt
        
    SE_mat_1 <- matrix(rep(SE.tmp * sqrt(N.samplesize.tmp),length(SE.tmp)),length(SE.tmp),length(SE.tmp))
    
    SE_mat_2 <- matrix(rep(SE.tmp * sqrt(N.samplesize.tmp),each=length(SE.tmp)),length(SE.tmp),length(SE.tmp))
    
    epsilon_mat <- SE_mat_1
    epsilon_mat[SE_mat_1>=SE_mat_2] <- SE_mat_2[SE_mat_1>=SE_mat_2]
        
    min_matrix <- outer(N.samplesize, N.samplesize, FUN = "pmin")
    XtY_sigma <- floor(Nt - Nt^2 / min_matrix) * (epsilon_mat^2) * LD
    eigen_decom <- eigen(XtY_sigma)
    eigen_decom_value <- pmax(eigen_decom$value,0)
    project_block <- eigen_decom$vector %*% diag(sqrt(eigen_decom_value))
    project_mat = project_block
    
    # sample summary statistics
    Cor.squared <- c()
    for (i in 1:k){
        set.seed(chr + i)
        XtY.test.temp <- c()

        XtY_mu_tmp <- XtY_mu
        if (length(XtY_mu)==1) {
                XtY.test.temp <- c(XtY.test.temp,unlist(XtY_mu_tmp + project_mat * rnorm(n=1,0,1)))
        } else {
                XtY.test.temp <- c(XtY.test.temp,unlist(XtY_mu_tmp + project_mat %*% rnorm(ncol(LD),0,1)))
        }
        
        
        ## calculate training summary statistics
        XtY.train.temp <- XtY - XtY.test.temp
        beta.train.temp <- XtY.train.temp/(Ntr-1)
        beta.v.temp     <- XtY.test.temp/(Nt-1)
        SE.train.temp <- sqrt(N.samplesize/Ntr) * GWAS$SE
        SE.v.temp     <- sqrt(N.samplesize/Nt) * GWAS$SE
        Zscore.temp <- beta.train.temp/ SE.train.temp
        pvalue.temp <- 2*pnorm(abs(Zscore.temp),lower.tail=F)
        
        #Z = beta.v.temp / SE.v.temp
        #beta.v.temp - Z / sqrt(Z^2 + Nt - 2)
        
        ## write new beta, SE, and p-value in the GWAS sumstats
        GWAS.tmp <- GWAS
        GWAS.tmp$`BETA` <- beta.train.temp
        GWAS.tmp$`SE` <- SE.train.temp
        GWAS.tmp$`P` <- pvalue.temp
        GWAS.tmp$`N` <- Ntr
        GWAS.tmp$`Z` <- Zscore.temp

        # this is training data set
        fwrite(GWAS.tmp,paste0(output_path,trait_name,".gwas.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
        
        ## write XtY
        ##fwrite(data.frame(CHR=GWAS$CHR,SNP=GWAS$SNP,A1=GWAS$A1,A2=GWAS$A2,test=XtY.test.temp),paste0(output_path,trait_name,".xty.ite",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)

        ## write sumstat for tuning 
        fwrite(data.frame(CHR=GWAS$CHR,SNP=GWAS$SNP,A1=GWAS$A1,A2=GWAS$A2,BETA=beta.v.temp,SE=SE.v.temp, N=Nt),paste0(output_path,trait_name,".valid",i,chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
    }
    
    return(Var.Y)
}

main.subsampling <- function(SS, matrix.LD){
    matched_data <- list(
    gwas_matched = as.data.table(SS),
    new_rs_blocks = list(SS$SNP),
    new_ld_blocks = list(matrix.LD),
    ld_block_size = list(c(1, nrow(SS)))
    )

    # PUMAS R2, k times
    #GWAS=matched_data$gwas_matched
    #Nt=floor(partitions[2]*min(matched_data$gwas_matched$N))
    #LD = matched_data$new_ld_blocks
    #LD.block=matched_data$ld_block_size
    #trait_name=trait_name
    #k=k
    #chr=chr
    
    pumas_tmp <- PUMAS.II.FUN.I(GWAS=matched_data$gwas_matched, Nt=floor(partitions[2]*min(matched_data$gwas_matched$N)),LD = matched_data$new_ld_blocks, LD.block=matched_data$ld_block_size, trait_name=trait_name, k=k, chr=chr)
    
    #fwrite(data.frame(var.Y=pumas_tmp,N.t=floor(partitions[2]*min(matched_data$gwas_matched$N))),paste0(output_path,trait_name,".forEVAL",chr,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
}


# TODO: to make it more faster; read X directly...
quasicors2 <- function(testr, N, penalizedBetas, refPanel, ldadjust =NULL, useld=FALSE){
    if(length(N) == 1) N = rep(N, length(betas))
    
    penalizedBetas = as.matrix(penalizedBetas)
    
    if(useld) {
        corX = ldadjust
        corX = corX[,SS.original$SNP]
        corX = corX[SS.original$SNP,]
    } else {
        X = BEDMatrix(refPanel, simple_names = TRUE) %>% PatchUp()
        X = X[, SS.original$SNP]
        
        # calculate the S and S inverse first
        X = scale(X)            # standardise
        X[is.na(X)] = 0            # mean imputation for missing data
        N.ref = nrow(X)
        
        corX = cor(X)
        #s = 0.1
        #corX = s * diag(rep(1,dim(corX)[1])) + (1-s) * corX
    }
    
    quasicorrelations = rep(0, ncol(penalizedBetas))
    num = rep(0, ncol(penalizedBetas))
    denom = rep(0, ncol(penalizedBetas))
    
    N[is.na(N)] = mean(N,na.rm=TRUE)
    
    #STAT = testBetas/ testSes
    r = testr #STAT / sqrt(STAT^2 + N - 2)    # for continuous phenos, convert Z to r;
    
    for(i in 1:length(quasicorrelations)){
        
        penalizedBetasTemp = penalizedBetas[,i]
        num[i] = sum( r * penalizedBetasTemp)
        
        denom[i] = sqrt(t(penalizedBetasTemp) %*% corX %*% penalizedBetasTemp)
        quasicorrelations[i] = num[i] / denom[i]
    }
    
    toReturn = list(quasicorrelations = quasicorrelations)
    return(toReturn)
}



main.evaluation <- function(){
    
    pumas.cor2 <- c()

    for(j in 1:k) {
        snp.w <- as.data.frame(fread(paste0(weight_path,trait_name,".",prs_method,".ite",j,".txt"),header=T))

        validdat = as.data.frame(fread(paste0(weight_path,trait_name,".valid",j,chr,".txt")))

        matrix.LDused <- matrix.LD.0[SS$SNP, SS$SNP]
        result.temp <- quasicors2(
            testr = validdat$BETA,
            N = validdat$N,
            penalizedBetas = snp.w[, 5:ncol(result)],
            refPanel = paste0(temp.outdir, gene),
            ldadjust = matrix.LDused,
            useld = TRUE
        )

        result.temp <- result.temp[[1]]
        result.temp[is.na(result.temp)] <- 0
        result.QuasiCorr <- result.temp
        
        pumas.cor2 <- rbind(pumas.cor2,result.QuasiCorr)

    }


    suffix <- colnames(snp.w)[-c(1:4)]
    colnames(pumas.cor2) <- suffix
    rownames(pumas.cor2) <- 1:k
    # store R2
    write.table(pumas.cor2,paste0(output_path,trait_name,".",prs_method,".txt"),col.names = T,row.names=F,quote=F,sep="\t")
}


normalize <- function(genotypes) {
    .Call('_lassosum_normalize', PACKAGE = 'lassosum', genotypes)
}

train_SUMMIT <- function(SS, matrix.LD,r=NULL) {
    
    if(is.null(r)) {
        SS["Z"] <- SS$BETA / SS$SE
        r <- SS$Z / sqrt(SS$N - 1 + SS$Z ^2)
    }
    
    matrix.LD <- matrix.LD[SS$SNP, SS$SNP]
    
    ############
    # r vector #
    ############
    
    # Compute r vector
    #
    
    # Get size
    size <- nrow(SS)
    
    ##########################
    # Iterate By parameter s #
    ##########################
    
    result <- SS[, c("CHR", "SNP", "A1", "A2")]
    
    for (method in c("ElNet")) {#, "LASSO", "MCP", "MNet", "SCAD"
        for (j in c(1,5,9)) {
            s <- 0.1 * j
            
            #########################################
            # Constructing the initial lambda array #
            #########################################
            
            # Get the maximum for lambda
            z.temp <- numeric(size)
            for (m in 1:size) {
                z.temp[m] <- abs(r[m])
            }
            
            lambda.max   <- max(z.temp)
            lambda.min   <- lambda.max * 1E-3
            lambda.array <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 10000)
            rm(z.temp)
            
            ###########################################
            # Explore the possible minimum for lambda #
            ###########################################
            
            ######################################
            # What do 0, 1, and 2 stand for?     #
            # 0: Not tested                      #
            # 1: tested and no problem detected  #
            # 2: tested and problem is detected  #
            ######################################
            
            problem <- numeric(10000)
            start   <- 1
            end     <- 10000
            OVERFITTING <- FALSE
            
            max.iteration <- 300
            
            threshold <- 5
            alpha     <- 0.5
            
            if (method == "SCAD") {
                gamma <- 3.7
            } else {
                gamma <- 3
            }
            
            dummy.detect <- TRUE
            
            while (dummy.detect) {
                # Stop if all lambdas are OK
                if (sum(problem == rep(1, 10000)) == 10000) {
                    #print("All lambdas check out!")
                    break
                }
                
                start <- min(which(problem == 0))
                end   <- max(which(problem == 0))
                
                try    <- ceiling((start + end) /2)
                lambda <- lambda.array[try]
                
                if (method == "MCP") {
                    problem[try] <- MCPDetect(r, as.matrix(matrix.LD), lambda.array, s, gamma, max.iteration, threshold, try)
                }
                if (method == "LASSO") {
                    problem[try] <- ElNetDetect(r, as.matrix(matrix.LD), lambda.array, s, 1, max.iteration, threshold, try)
                }
                if (method == "ElNet") {
                    problem[try] <- ElNetDetect(r, as.matrix(matrix.LD), lambda.array, s, 0.5, max.iteration, threshold, try)
                }
                if (method == "MNet") {
                    problem[try] <- MNetDetect(r, as.matrix(matrix.LD), lambda.array, s, alpha, gamma, max.iteration, threshold, try)
                }
                if (method == "SCAD") {
                    problem[try] <- SCADDetect(r, as.matrix(matrix.LD), lambda.array, s, gamma, max.iteration, threshold, try)
                }
                
                if (problem[try] == 1) {
                    problem[1:try] <- 1
                } else {
                    problem[try:10000] <- 2
                }
                
                if (sum(problem == 2) >= 9999) {
                    print("Unresolved overfitting!")
                    
                    OVERFITTING <- TRUE
                    dummy.detect <- FALSE
                } else {
                    if (problem[try] == 1 & problem[min(10000, try + 1)] == 2) {
                        print(paste0("Problem detected for lambda=", lambda.array[try + 1], ". Finer interval has been saved!"))
                        lambda.min   <- lambda
                        dummy.detect <- FALSE
                    }
                    
                    if (problem[try] == 2 & problem[try - 1] == 1) {
                        print(paste0("Problem detected for lambda=", lambda.array[try], ". Finer interval has been saved!"))
                        lambda.min   <- lambda.array[try - 1]
                        dummy.detect <- FALSE
                    }
                }
            }
            
            # Update the lambda array
            lambda.array <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 100)
            
            #######################################################################
            # Optimizing with the chosen penalty (using the updated lambda array) #
            #######################################################################
            
            #####################
            # Iterate by lambda #
            #####################
            
            beta <- t(matrix(0, nrow = 1, ncol = size))
            
            max.iteration <- 300
            
            threshold <- 1e-6
            alpha     <- 0.5
            
            if (method == "SCAD") {
                gamma <- 3.7
            } else {
                gamma <- 3
            }
            
            if (!OVERFITTING) {
                if (method == "MCP") {
                    res <- MCP(r, as.matrix(matrix.LD), lambda.array, s, gamma, max.iteration, threshold)
                }
                if (method == "LASSO") {
                    res <- ElNet(r, as.matrix(matrix.LD), lambda.array, s, 1, max.iteration, threshold)
                }
                if (method == "ElNet") {
                    res <- ElNet(r, as.matrix(matrix.LD), lambda.array, s, 0.5, max.iteration, threshold)
                }
                if (method == "MNet") {
                    res <- MNet(r, as.matrix(matrix.LD), lambda.array, s, alpha, gamma, max.iteration, threshold)
                }
                if (method == "SCAD") {
                    res <- SCAD(r, as.matrix(matrix.LD), lambda.array, s, gamma, max.iteration, threshold)
                }
                
                result <- cbind(result, res)
            } else {
                result <- cbind(result, matrix(0, nrow = size, ncol = 100))
            }
        }
    }
    
    return(result)
}


pseudoGic <- function(penalizedBetas, betas, ses, N, refPanel,ldadjust =NULL, useld=FALSE){
    
    if(length(N) == 1) N = rep(N, length(betas))
    
    penalizedBetas = as.matrix(penalizedBetas)
    
    fam = read.table(paste0(refPanel,".fam"))
    n = nrow(fam)
    
    if(useld) {
        corX = ldadjust
        corX = corX[,SS.original$SNP]
        corX = corX[SS.original$SNP,]
    } else {
        
        X = BEDMatrix(refPanel, simple_names = TRUE) %>% PatchUp()
        X = X[, SS.original$SNP]
        
        # calculate the S and S inverse first
        X = scale(X)            # standardise
        X[is.na(X)] = 0            # mean imputation for missing data
        N.ref = nrow(X)
        
        corX = cor(X)
    }
    
    eig = eigen(corX)
    lambda = eig$values
    Q = eig$vectors
    
    
    cum.perc = cumsum(lambda / sum(lambda) * 100)
    prune.thresh = 99 #default and do not want to change this
    
    keep = 1:min(which(cum.perc >= prune.thresh))
    
    locN =  mean(N,na.rm=TRUE)
    locK = length(keep)
    
    N[is.na(N)] = mean(N,na.rm=TRUE)
    
    # cap the number of PCs at a proportion of the (lowest) GWAS input sample size (cdl 16/3)
    max.prop.K = 0.75
    if (!is.null(max.prop.K) && !is.na(max.prop.K)) {
        max.K = floor(max.prop.K * min(N))
        if (length(keep) > max.K) {
            keep = keep[1:max.K]
            Q = Q[,keep]
            lambda = lambda[keep]
        }
    }
    
    
    STAT = betas/ ses
    r = STAT / sqrt(STAT^2 + N - 2)    # for continuous phenos, convert Z to r;
    
    alpha = Q[,keep] %*% diag(1/lambda[keep]) %*% t(Q[,keep]) %*% r
    
    delta = diag(c(sqrt(lambda[keep]))) %*% t(Q[,keep]) %*% alpha    # using keep to filter out dropped PCs (since this one is updated every time a PC is dropped)
    eta = t(r) %*% alpha; eta = (locN-1) / (locN-locK-1) * (1-eta)
    sigma = eta/(locN-1) # eta is the variance of the noise~N(0,eta)
    
    # get h2
    h2.obs = 1 - (1 - t(delta) %*% delta) * ((locN-1) / (locN-locK-1))
    
    SSEvec = NULL
    qVec = NULL
    aicVec = NULL
    bicVec = NULL
    bxxb = NULL
    bxy = NULL
    gicVec = NULL
    gicVec2 = NULL
    
    for(k in 1:ncol(penalizedBetas)){
        
        penalizedBetasTemp = penalizedBetas[,k]
        
        qElast = sum(penalizedBetasTemp!=0)
        qVec = c(qVec,qElast)
        
        
        bxyTemp = t(penalizedBetasTemp)%*% r
        bxy = c(bxy,bxyTemp[1,1])
        
        
        bxxbWeight = t(penalizedBetasTemp) %*% corX %*% penalizedBetasTemp
        
        bxxb = c(bxxb, bxxbWeight[1,1])
        
        SSEest = 1 - (2 * bxyTemp) + bxxbWeight #ytyEst = 1
        SSEvec = c(SSEvec,SSEest[1,1])
        
        logLik = -SSEest*median(N) / (2 *  eta)
        
        aicTemp = (2 * qElast) - (2 * logLik)
        bicTemp = (log(median(N)) * qElast) - (2 * logLik)
        gicTemp = log(median(N)) * log(nrow(penalizedBetas)) * qElast - (2 * logLik)
        p = nrow(penalizedBetas)
        
        gicTemp2 = 2 * (log(p) + log(log(p))) * qElast - (2 * logLik)
        
        aicVec = c(aicVec, aicTemp[1,1])
        bicVec = c(bicVec, bicTemp[1,1])
        gicVec = c(gicVec, gicTemp[1,1])
        gicVec2 = c(gicVec2, gicTemp2[1,1])
    }
    
    toReturn = list(aic = aicVec, bic=bicVec, gic=gicVec, gic2 = gicVec2, eta = eta, h2.obs = h2.obs, nPC = length(keep))
    return(toReturn)
}

