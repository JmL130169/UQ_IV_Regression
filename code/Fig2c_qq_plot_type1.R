setwd("")

library(data.table)
library(lattice)
library(ggplot2)
library(ggsci)
library(dplyr)

qqunif.plot <- function(pvalues,
                        should.thin = T, thin.obs.places = 4, thin.exp.places = 4,
                        xlab = expression(paste("Expected (", -log[10], " p-value)")),
                        ylab = expression(paste("Observed (", -log[10], " p-value)")),
                        draw.conf = TRUE, conf.points = 1000, conf.col = "lightgray", conf.alpha = .05,
                        already.transformed = FALSE, pch = 20, aspect = "fill", prepanel = prepanel.qqunif,
                        par.settings = list(superpose.symbol = list(pch = pch)), ...) {


    # error checking
    if (length(pvalues) == 0) stop("pvalue vector is empty, can't draw plot")
    if (!(class(pvalues) == "numeric" ||
        (class(pvalues) == "list" && all(sapply(pvalues, class) == "numeric")))) {
        stop("pvalue vector is not numeric, can't draw plot")
    }
    if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed == FALSE) {
        if (any(unlist(pvalues) == 0)) stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues) < 0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
    }

    grp <- NULL
    n <- 1
    exp.x <- c()
    if (is.list(pvalues)) {
        nn <- sapply(pvalues, length)
        rs <- cumsum(nn)
        re <- rs - nn + 1
        n <- min(nn)
        if (!is.null(names(pvalues))) {
            grp <- factor(rep(names(pvalues), nn), levels = names(pvalues))
            names(pvalues) <- NULL
        } else {
            grp <- factor(rep(1:length(pvalues), nn))
        }
        pvo <- pvalues
        pvalues <- numeric(sum(nn))
        exp.x <- numeric(sum(nn))
        for (i in 1:length(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method = "first") - .5) / nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvo[[i]], ties.method = "first") - .5) / (nn[i] + 1))
            }
        }
    } else {
        n <- length(pvalues) + 1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method = "first") - .5) / n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n - rank(pvalues, ties.method = "first") - .5) / n)
        }
    }


    # this is a helper function to draw the confidence interval
    panel.qqconf <- function(n, conf.points = 1000, conf.col = "gray", conf.alpha = .05, ...) {
        require(grid)
        conf.points <- min(conf.points, n - 1)
        mpts <- matrix(nrow = conf.points * 2, ncol = 2)
        for (i in seq(from = 1, to = conf.points)) {
            mpts[i, 1] <- -log10((i - .5) / n)
            mpts[i, 2] <- -log10(qbeta(1 - conf.alpha / 2, i, n - i))
            mpts[conf.points * 2 + 1 - i, 1] <- -log10((i - .5) / n)
            mpts[conf.points * 2 + 1 - i, 2] <- -log10(qbeta(conf.alpha / 2, i, n - i))
        }
        grid.polygon(x = mpts[, 1], y = mpts[, 2], gp = gpar(fill = conf.col, lty = 0), default.units = "native")
    }

    # reduce number of points to plot
    if (should.thin == T) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(
                pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places),
                grp = grp
            ))
            grp <- thin$grp
        } else {
            thin <- unique(data.frame(
                pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places)
            ))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }
    gc()

    prepanel.qqunif <- function(x, y, ...) {
        A <- list()
        A$xlim <- range(x, y) * 1.02
        A$xlim[1] <- 0
        A$ylim <- A$xlim
        return(A)
    }
    
    # Define the colors
    cols <- c('#33A02C', '#FF7F00', 'hotpink2', 'mediumpurple1', 'lightskyblue', 'yellow3', 'orange')
        
    # Define the par.settings for the legend
    legend.settings <- list(
            superpose.symbol = list(
                col = cols,  # Set the color
                pch = pch    # Keep the same point shape
            )
    )
        
        larger_point_size <- 1.6 # You can adjust this value

        # Define larger font size for axis labels
        larger_font_size <- 1.2  # You can adjust this value
        
        larger_title_size = 1.4

        # draw the plot
          xyplot(pvalues ~ exp.x,
              groups = grp,  xlab = list(label = expression(paste("Expected (", -log[10], " p-value)")), cex = larger_title_size),
              ylab = list(label = expression(paste("Observed (", -log[10], " p-value)")), cex = larger_title_size), aspect = aspect,
              prepanel = prepanel,  scales = list(axs = "i", x = list(cex = larger_font_size), y = list(cex = larger_font_size)), pch = pch, xlim = c(0, 5), ylim = c(0, 5),
              panel = function(x, y, subscripts, groups, ...) {
                  if (draw.conf) {
                      panel.qqconf(n,
                          conf.points = conf.points,
                          conf.col = conf.col, conf.alpha = conf.alpha
                      )
                  }
                  col.index <- as.integer(groups[subscripts])
                  panel.xyplot(x, y, col=cols[col.index],alpha = 0.9, cex=larger_point_size,...)
                  panel.abline(0, 1)
              }, par.settings = legend.settings,  # Use the updated settings
              ...
          )
}


# wd = "/Users/jliu28/Dropbox/MDACC_research/Undergoing/scTWAS/shared/simulations2/"
# temp <- fread(paste0(wd, "map_file.csv"), data.table = FALSE)
# 
# p1 = NULL
# p2 = NULL
# p3 = NULL
# p4 = NULL
# p5 = NULL
# 
# for(i in 1:nrow(temp)) {
#     tryCatch({
#         setwd(paste0(wd,temp$name_batch[i]))
#     
#     files <- dir()
#     files <- files[grepl("-",files)]
#     
#     p1.tmp = NULL
#     p2.tmp = NULL
#     p3.tmp = NULL
#     p4.tmp = NULL
#     p5.tmp = NULL
#     
#     for (j in 1:length(files)) {
#         n <- files[j] %>% gsub(".RData", "", .) %>% str_split(., "-") %>% unlist()
#         n <- n[2]
#         
#         load(files[j])
#         
#         if(n == "50") {
#             p1.tmp = c(p1.tmp,p.t1e[1,])
#         } else if (n == "100") {
#             p2.tmp = c(p2.tmp,p.t1e[1,])
# 
#         } else if (n == "200") {
#             p3.tmp = c(p3.tmp,p.t1e[1,])
# 
#         } else if (n == "500") {
#             p4.tmp = c(p4.tmp,p.t1e[1,])
# 
#         } else {
#             p5.tmp = c(p5.tmp,p.t1e[1,])
#             if(sum(p.t1e[1,]<1e-3)>2) {
#                 cat("out\n")
#                 #stop()
#             }
#         }
#     }
#     
#     p1 = c(p1, p1.tmp)
#     p2 = c(p2, p2.tmp)
#     p3 = c(p3, p3.tmp)
#     p4 = c(p4, p4.tmp)
#     p5 = c(p5, p5.tmp)
#    
#    cat("finisih ",i, "  ", sum(p1.tmp<0.05),"\n")
#    cat("finisih ",i, "  ", sum(p5.tmp<0.05),"\n")
# 
#     }, error=function(e){})
#    
# }
# 
# p1 = p1[!is.na(p1)]
# p2 = p2[!is.na(p2)]
# p3 = p3[!is.na(p3)]
# p4 = p4[!is.na(p4)]
# p5 = p5[!is.na(p5)]

# Set the working directory to the folder containing the .RData files
setwd("")
file_names <- dir(pattern = "\\.RData$")
load(file_names[1])
combined_data = p.t1e

# Loop through the list of .RData files
for(file in 2:length(file_names)) {
  
    load(file_names[file])
    combined_data <- cbind(combined_data, p.t1e)
    
}

combined_data <- combined_data[, !apply(is.na(combined_data), 2, any)]

# Your combined data is now in 'combined_data'

# my.pvalue.list <- list("old" = combined_data[1, ], "data-driven" = combined_data[2, ], "true h2" = combined_data[3, ],
#                        "true var" = combined_data[4, ], "non-zero data-driven" = combined_data[5, ],
#                        "non-zero true h2" = combined_data[6, ], "non-zero true var" = combined_data[7, ])
# my.pvalue.list <- list("old" = combined_data[1, ], "data-driven" = combined_data[2, ], "true h2" = combined_data[3, ],
#                        "non-zero data-driven" = combined_data[4, ],
#                        "non-zero true h2" = combined_data[5, ])
# my.pvalue.list <- list("old" = combined_data[1, ], "data-driven" = combined_data[2, ], "true h2" = combined_data[3, ],
#                        "true var" = combined_data[4, ], "individual all SNPs" = combined_data[5, ],
#                        "individual valid SNPs" = combined_data[6, ])
# my.pvalue.list <- list("data-driven" = combined_data[2, ], "true h2" = combined_data[3, ],
#                        "individual all SNPs" = combined_data[5, ],
#                        "individual valid SNPs" = combined_data[6, ])
# my.pvalue.list <- list("Summary" = combined_data[1, ], "Individual" = combined_data[2, ])
my.pvalue.list <- list("2SLS" = combined_data[1, ], "Oracle" = combined_data[2, ],
                       "2SLS-Corrected" = combined_data[3, ])

# my.pvalue.list <- lapply(my.pvalue.list, na.omit)
remove_zero <- function(x) {
  x[x != 0]
}
my.pvalue.list <- lapply(my.pvalue.list, remove_zero)


setwd("")
pdf(".pdf")
qqunif.plot(my.pvalue.list, auto.key = list(corner = c(.95, .05)), thin.obs.places = 2, thin.exp.places =2)
dev.off()
