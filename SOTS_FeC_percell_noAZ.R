setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment") #setwd
#loading packages
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(ggtext)
library(patchwork)
library(car)
library(emmeans)
library(visreg)
library(writexl)
library(rstan)
library(rstanarm)

#-----------import the data --------------------------------------------------------------------------------------
FeC <- read_csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/sots_FeC.csv") %>%
  mutate(light = as.factor(light), 
         inhibitor = as.factor(inhibitor) %>% fct_recode("+AZ" = "AZ"),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         size = as.factor(size) %>% fct_recode("0.2-2" = "0.2", "2-20" = "2", ">20" = "20"),
         bottle = as.factor(bottle))

#remove straight uptake data and fraction data
FeC <- FeC %>% select(-Fe_up, -C_up, -FeC)

# Reorder  levels
FeC$treatment <- factor(FeC$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
FeC$inhibitor <- factor(FeC$inhibitor, levels = c("Control", "+AZ"), ordered = TRUE)

# log transform Fe and C variables
FeC$Fe_cell_ln <- log(FeC$Fe_cell)
FeC$C_cell_ln <- log(FeC$C_cell)
# square root transform Fe and C variables
FeC$Fe_cell_sqrt <- sqrt(FeC$Fe_cell)
FeC$C_cell_sqrt <- sqrt(FeC$C_cell)

#check the structure
str(FeC)

#FeC_noAZ 
FeC_noAZ <- FeC %>% filter(inhibitor == "Control")

# ---Compute the mean and variance for each group defined by light, treatment, and size-----------------------
FeC_noAZ_summary <- FeC_noAZ %>%
  group_by(treatment, size) %>%
  summarise(mean_Fe_cell = mean(Fe_cell, na.rm = TRUE),
            var_Fe_cell = var(Fe_cell, na.rm = TRUE),
            mean_C_cell = mean(C_cell, na.rm = TRUE),
            var_C_cell = var(C_cell, na.rm = TRUE))

# Function to create scatterplots of variance vs. mean
create_scatterplot <- function(data, x_var, y_var, title) {
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point() +
    theme_bw() +
    xlab("Mean") +
    ylab("Variance") +
    ggtitle(title)
}

# Create scatterplots for Fe_up, C_up, and FeC variables
scatterplot_Fe_cell <- create_scatterplot(FeC_noAZ_summary, "mean_Fe_cell", "var_Fe_cell", "Fe Variance vs. Mean")
scatterplot_C_cell <- create_scatterplot(FeC_noAZ_summary, "mean_C_cell", "var_C_cell", "C Variance vs. Mean")

# Display scatterplots
scatterplot_Fe_cell | scatterplot_C_cell 

#------FITTING MODELS :  -----------------------------------------------------------------------------------
Fe_model  <- glm((Fe_cell) ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model <- glm((C_cell) ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))

Fe_model_reduced <- glm(Fe_cell ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model_reduced <- glm(C_cell ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"))

#---- Fitting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_cell ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced <- stan_glm(Fe_cell ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_cell ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced <- stan_glm(C_cell ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for each model
Fe_loo_full <- loo(Fe_stan_full, k_threshold = 0.7)
Fe_loo_reduced <- loo(Fe_stan_reduced, k_threshold = 0.7)
C_loo_full <- loo(C_stan_full, k_threshold = 0.7)
C_loo_reduced <- loo(C_stan_reduced, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(Fe_loo_full, Fe_loo_reduced)
loo_compare(C_loo_full, C_loo_reduced)

# -----Create a list of models with their corresponding names--------------------------
models_list <- list(
  Fe = Fe_model,
  C = C_model,
  Fe_2 = Fe_model_reduced,
  C_2 = C_model_reduced)

# ---- QQ plots for the gamma models ----------------------------------------------------------
#function to calculate standardised residuals
calculate_standardized_residuals <- function(model) {
  pearson_residuals <- residuals(model, type = "pearson")
  return(pearson_residuals)
}

# Function to calculate shape and rate for a given model
calculate_shape_rate <- function(model) {
  fitted_values <- fitted(model)
  dispersion <- summary(model)$dispersion
  shape <- 1 / dispersion
  rate <- shape / mean(fitted_values)
  return(list(shape = shape, rate = rate))
}

# Calculate shape and rate for each model in models_list
shape_rate_list <- lapply(models_list, calculate_shape_rate)

# Function to generate QQ plots with shape and rate
create_qq_plot <- function(residuals, model_name, shape, rate) {
  ggplot(data = data.frame(residuals), aes(sample = residuals)) +
    geom_qq(distribution = qgamma, dparams = list(shape = shape, rate = rate)) +
    geom_qq_line(distribution = qgamma, dparams = list(shape = shape, rate = rate)) +
    theme_bw() +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    ggtitle(paste(model_name, "QQ Plot (Gamma)"))
}

# Generate standardized residuals for each model in models_list
std_residuals_list <- lapply(models_list, calculate_standardized_residuals)

# Create a list of QQ plots for each model in models_list, incorporating shape and rate
qq_plot_list <- mapply(function(x, y) {
  create_qq_plot(std_residuals_list[[x]], x, y$shape, y$rate)
}, x = names(std_residuals_list), y = shape_rate_list, SIMPLIFY = FALSE)

# Convert the list of plots to a patchwork object
qq_plot_patchwork <- wrap_plots(qq_plot_list, ncol = 2)

# Print the plot
qq_plot_patchwork

# quick diagnostics
par(mfrow = c(4, 2))
plot(Fe_model, which = 1:4)
plot(Fe_model_reduced, which = 1:4)
plot(C_model, which = 1:4)
plot(C_model_reduced, which = 1:4)

# ----Function to create residual vs fitted plot for a model------------------------------------------------
res_fit_plot <- function(model, season, response) {
  res_fit <- broom::augment(model)
  
  plot <- ggplot(data = res_fit, aes(x = .fitted, y = .resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(paste(season, response, "Res-vs-Fit")) +
    theme_bw()
  
  return(plot)
}

# Create a list of Residuals vs Fitted plots
res_fit_plot_list <- lapply(names(models_list), function(x) {
  res_fit_plot(models_list[[x]], x, "Residuals")
})

# Convert the list of plots to a patchwork object
res_plot_patchwork <- wrap_plots(res_fit_plot_list, ncol = 2)
res_plot_patchwork

# --- interrogating the model results using EMMEANS ---------------------------------------------------
# emmeans for Fe uptake
emm_all_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ treatment : size))
emm_trt_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ treatment | size), type = "unlink")
emm_size_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ size | treatment), type = "unlink")
# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model, specs = pairwise ~ treatment : size))
emm_trt_C <- summary(emmeans(C_model, specs = pairwise ~ treatment | size), type = "unlink")
emm_size_C <- summary(emmeans(C_model, specs = pairwise ~ size | treatment), type = "unlink")

# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, trt = emm_trt_Fe$contrasts, size = emm_size_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, trt = emm_trt_C$contrasts, size = emm_size_C$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_noAZ_PERCELL.xlsx")

#--plotting them up ----------------------------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_noAZ_summary <- FeC_noAZ %>%
  group_by(treatment, size) %>%
  summarise(
    Fe_mean = mean(Fe_cell),
    Fe_se = sd(Fe_cell) / sqrt(n()),
    C_mean = mean(C_cell),
    C_se = sd(C_cell) / sqrt(n()))

format_y_labels <- function(x) {
  ifelse(x > 1, round(x), x)
}

# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, show_xlab = TRUE, y_min = NULL, y_max = NULL) {
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), fill = size)) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    geom_point(size = 3, shape = 21, colour = "black") +
    scale_fill_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
    scale_x_discrete(breaks = unique(data$treatment)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Size (μm)" else NULL,
         y = paste(variable_name)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"), 
          panel.grid.minor = element_blank()) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    ) 
  return(p)
}

# Generate plots for each variable
Fe_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 60) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "top") +
  labs(y = expression(paste("Fe uptake per cell")))
Fe_plot_errorbars_low <- generate_plots_with_errorbars(FeC_noAZ_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 3) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none") 
C_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "C_mean", "C_se", "", y_min = 0, y_max = 6) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none") +
  labs(y = expression(paste("C uptake per cell")))
C_plot_errorbars_low <- generate_plots_with_errorbars(FeC_noAZ_summary, "C_mean", "C_se", "", y_min = 0, y_max = 0.3) +
  theme(axis.title.x = element_blank(), legend.position = "none") 

# ----- plotting model results in vireg -------------------------------------------------------------------
Fe_asterisks <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(-1.2, -1.4, -0.5, 0.6, 0.7, 1.7, 3.5, 4, 4.5), 
  label = c("a", "a", "b", "a", "a", "b", "a", "b", "c"), colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D","#007B3D"))

C_asterisks <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(-1, -1.5, -1, -1.7, -2, -1.4, 2.1, 2.2, 2.1), 
  label = c("a", "b", "c", "a", "b", "a", "a", "a", "a"), colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D","#007B3D"))

Fe_visreg <- visreg(Fe_model, "treatment", by = "size", data = FeC_noAZ, 
                        scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33")) +
  scale_y_continuous(limits = c(-4, 5), breaks = seq(-4, 5, by = 3), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe uptake)")), breaks = seq(-4, 5, by = 3))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = Fe_asterisks, aes(x = x_pos, y = y_pos, label = label), 
            colour = Fe_asterisks$colours, size = 4, inherit.aes = FALSE)  +
  xlab("") +
  ylab("")

C_visreg <- visreg(C_model, "treatment", by = "size", data = FeC_noAZ, 
                       scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-5, 3, by = 2), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")), breaks = seq(-5, 3, by = 2))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = C_asterisks, aes(x = x_pos, y = y_pos, label = label), 
            colour = C_asterisks$colours, size = 4, inherit.aes = FALSE) +
  xlab("Treatment") +
  ylab("")

#----- Combining all in one plot ----------------------------------------------------------------------------------------------
Fe_data <- (Fe_plot_errorbars / Fe_plot_errorbars_low) +
  plot_layout(nrow = (2), heights = c(2,1)) & 
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "none") 

C_data <- (C_plot_errorbars / C_plot_errorbars_low) + 
  plot_layout(nrow = (2), heights = c(2,1)) & 
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "none") 
        
# Create combined plots without legends for each data type
Fe_comb_plot <- Fe_data | Fe_visreg + theme(legend.position = "none", axis.text.x = element_blank())
C_comb_plot <- C_data | C_visreg + theme(legend.position = "none")

FeC_percell_plot <-  guide_area() / Fe_comb_plot / C_comb_plot + 
  plot_annotation(tag_levels = list(c('A','','B','C','','D'))) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) & 
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.02, 0.96)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plot
ggsave("FeC_percell_plot.svg", FeC_percell_plot, width = 7.2, height = 8)
