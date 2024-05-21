# set dir
setwd("")
file_names <- dir(pattern = "\\.RData$")
load(file_names[1])
combined_data = vector("list", length = 6)
null_data = bias
null_data_mse = mse

# Loop through the list of .RData files

# bias
for(file in 2:length(file_names)) {
  
    load(file_names[file])
    null_data <- null_data + bias
    null_data_mse <- null_data_mse + mse
  
}

bias_1 <- abs(null_data / length(file_names))
mse_1 <- abs(null_data_mse / length(file_names))


# coverage
for (i in 1:6) {
    combined_data[[i]] = cover[ , , i]
}

# Loop through the list of .RData files
for(file in 2:length(file_names)) {
    load(file_names[file])
    for (i in 1:6) {
        combined_data[[i]] = cbind(combined_data[[i]], cover[ , , i])
    }
}


# plot function
cover_prob = matrix(0, ncol = 5, nrow = 6)

for (i in 1:nrow(cover_prob)) {
    
    for (j in 1:ncol(cover_prob)) {
      
        cover_per = combined_data[[i]][j, ]
        cover_prob[i, j] = sum(cover_per)/length(cover_per)
        
    }
  
}

# Preprocess the data for ggplot

x = c(-0.1, -0.05, 0, 0.05, 0.1)
cover_melt = data.frame(eff_size = rep(x, 3), p = as.vector(t(cover_prob[1:3, ])))
cover_melt$method = c(rep("m1", 5), rep("m2", 5), rep("m3", 5))

bias_melt = data.frame(eff_size = rep(x, 3), b = as.vector(t(bias_1[1:3, ])))
bias_melt$method = c(rep("m1", 5), rep("m2", 5), rep("m3", 5))
                     
mse_melt = data.frame(eff_size = rep(x, 3), b = as.vector(t(mse_1[1:3, ])))
mse_melt$method = c(rep("m1", 5), rep("m2", 5), rep("m3", 5))

colors <- c('#33A02C', '#FF7F00', 'hotpink2')
shapes <- c(16, 3, 2)
line_types <- c("solid", "dotdash", "dotted")


pdf(".pdf")

ggplot(cover_melt, aes(x = eff_size, y = p, group = method)) +
  geom_line(aes(color = method, linetype = method), size = 1) +
  geom_point(aes(color = method, shape = method), size = 3.5) +
  ylim(0.5, 1) + 
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
                     limits = c(min(cover_melt$eff_size), max(cover_melt$eff_size))) +
  # labs(x = "Phenotypic heritability", y = "Statistical power")
  labs(x = "True effect size", y = "Coverage probability")

dev.off()

pdf(".pdf")

ggplot(bias_melt, aes(x = eff_size, y = b, group = method)) +
  geom_line(aes(color = method, linetype = method), size = 1) +
  geom_point(aes(color = method, shape = method), size = 3.5) +
  # ylim(0, 0.05) + 
  theme(aspect.ratio = 1) +
  theme(
    axis.title.x = element_text(size = 18),  # Customize x-axis title
    axis.text.x = element_text(size = 15), # Customize x-axis labels (ticks)
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 15)
  ) +
  theme(legend.position = c(0.9, 0.25)) +
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
                     limits = c(min(bias_melt$eff_size), max(bias_melt$eff_size))) +
  # labs(x = "Phenotypic heritability", y = "Statistical power")
  labs(x = "True effect size", y = "Mean Absolute Error")

dev.off()

pdf(".pdf")

ggplot(mse_melt, aes(x = eff_size, y = b, group = method)) +
  geom_line(aes(color = method, linetype = method), size = 1) +
  geom_point(aes(color = method, shape = method), size = 3.5) +
  # ylim(0, 0.05) + 
  theme(aspect.ratio = 1) +
  theme(
    axis.title.x = element_text(size = 18),  # Customize x-axis title
    axis.text.x = element_text(size = 15), # Customize x-axis labels (ticks)
    axis.title.y = element_text(size = 18),  
    axis.text.y = element_text(size = 15)
  ) +
  theme(legend.position = c(0.9, 0.25)) +
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
                     limits = c(min(mse_melt$eff_size), max(mse_melt$eff_size))) +
  # labs(x = "Phenotypic heritability", y = "Statistical power")
  labs(x = "True effect size", y = "Mean Square Error")

dev.off()


