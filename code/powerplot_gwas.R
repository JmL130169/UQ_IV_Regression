# Set the working directory to the folder containing the .RData files
setwd("")
file_names <- dir(pattern = "\\.RData$")
load(file_names[1])
combined_data = vector("list", length = 6)
null_data = p.t1e


# Loop through the list of .RData files
for(file in 2:length(file_names)) {
  
  load(file_names[file])
  
  for (i in 1:6) {
    
    combined_beta[[i]] = cbind(combined_beta[[i]], beta_estimation[ , , i])
    
  }
  
}



# Loop through the list of .RData files

# type 1 error
for(file in 2:length(file_names)) {
  
    load(file_names[file])
    null_data <- cbind(null_data, p.t1e)
  
}
null_data <- null_data[, !apply(is.na(null_data), 2, any)]


# power
for (i in 1:6) {
  
    combined_data[[i]] = p.power[ , , i]
    
}


# Loop through the list of .RData files
for(file in 2:length(file_names)) {
  
    load(file_names[file])
  
    for (i in 1:6) {
      
        combined_data[[i]] = cbind(combined_data[[i]], p.power[ , , i])
        
    }
  
}



# plot function
thr = 0.05/10000
power_thr = matrix(0, ncol = 5, nrow = 6)

for (i in 1:nrow(power_thr)) {
    
    for (j in 1:ncol(power_thr)) {
      
        pw = combined_data[[i]][j, ]
        power_thr[i, j] = length(pw[pw < thr])/length(pw)
        
    }
  
}

# Preprocess the data for ggplot
x = c(-0.1, -0.05, 0, 0.05, 0.1)
power_melt = data.frame(eff_size = rep(x, 3), p = as.vector(t(power_thr[1:3, ])))
power_melt$method = c(rep("m1", 5), rep("m2", 5), rep("m3", 5))

colors <- c('#33A02C', '#FF7F00', 'hotpink2')
shapes <- c(16, 3, 2)
line_types <- c("solid", "dotdash", "dotted")

pdf(".pdf")

ggplot(power_melt, aes(x = eff_size, y = p, group = method)) +
  geom_line(aes(color = method, linetype = method), size = 1) +
  geom_point(aes(color = method, shape = method), size = 3.5) +
  ylim(0, 1) + 
  theme(aspect.ratio = 1) +
  theme(
    axis.title.x = element_text(size = 18),  # Customize x-axis title
    axis.text.x = element_text(size = 15), # Customize x-axis labels (ticks)
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 15)
  ) +
  theme(legend.position = c(0.9, 0.15)) +
  #theme_minimal(base_size = 14) +
  scale_shape_manual(values = shapes,
                     labels = c("2SLS", "Oracle", "2SLS-Corrected"),
                     name = "Method") +
  scale_color_manual(values = colors,
                     labels = c("2SLS", "Oracle", "2SLS-Corrected"),
                     name = "Method") +
  scale_linetype_manual(values = line_types,
                        labels = c("2SLS", "Oracle", "2SLS-Corrected"),
                        name = "Method") +
  scale_x_continuous(breaks = x,
                     limits = c(min(power_melt$eff_size), max(power_melt$eff_size))) +
  # labs(x = "Phenotypic heritability", y = "Statistical power")
  labs(x = "True effect size", y = "Statistical power")

dev.off()


