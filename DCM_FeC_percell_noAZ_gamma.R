setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment") #setwd
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
FeC <- read_csv("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment/dcm_FeC.csv") %>%
  mutate(light = as.factor(light), 
         inhibitor = as.factor(inhibitor) %>% fct_recode("+AZ" = "AZ"),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         size = as.factor(size) %>% fct_recode("0.2-2" = "0.2", "2-20" = "2", ">20" = "20"),
         bottle = as.factor(bottle))

#remove straight uptake data and fraction data
FeC <- FeC %>% select(-Fe_up, -C_up, -Fe_fract, -C_fract, -FeC)

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
  group_by(light, treatment, size) %>%
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
Fe_model  <- glm((Fe_cell) ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model <- glm((C_cell) ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))

Fe_model_reduced <- glm(Fe_cell ~ light + treatment + size + light:treatment + light:size + 
                          treatment:size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model_reduced <- glm(C_cell ~ light + treatment + size + light:treatment + light:size + 
                         treatment:size, data = FeC_noAZ, family = Gamma(link = "log"))

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_cell ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced <- stan_glm(Fe_cell ~ light + treatment + size + light:treatment + light:size + treatment:size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_cell ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced <- stan_glm(C_cell ~ light + treatment + size + light:treatment + light:size + treatment:size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)

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
emm_all_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ treatment : light : size))
emm_trt_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ treatment | light | size), type = "unlink")
emm_light_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ light | treatment | size), type = "unlink")
emm_size_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ size | treatment | light), type = "unlink")
# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model, specs = pairwise ~ treatment : light : size))
emm_trt_C <- summary(emmeans(C_model, specs = pairwise ~ treatment | light | size), type = "unlink")
emm_light_C <- summary(emmeans(C_model, specs = pairwise ~ light | treatment | size), type = "unlink")
emm_size_C <- summary(emmeans(C_model, specs = pairwise ~ size |  treatment | light ), type = "unlink")

# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, trt = emm_trt_Fe$contrasts, light = emm_light_Fe$contrasts,
                        size = emm_size_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, trt = emm_trt_C$contrasts, light = emm_light_C$contrasts,
                        size = emm_size_C$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_noAZ_PERCELL.xlsx")

# ----- plotting model results BY SIZE ------------------------------------------------------------------------------------------------------
# ---------------Plots of the data without inhibitor ----------------------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_noAZ_summary <- FeC_noAZ %>%
  group_by(treatment, size, light) %>%
  summarise(
    Fe_mean = mean(Fe_cell),
    Fe_se = sd(Fe_cell) / sqrt(n()),
    C_mean = mean(C_cell),
    C_se = sd(C_cell) / sqrt(n()))

format_y_labels <- function(x) {
  ifelse(x > 1, round(x), x)
}

# Generate function for facet grid with point and error bars
generate_plots <- function(data, size_value, variable_mean, variable_se, variable_name,  
                           show_legend = FALSE, show_xlab = TRUE, y_min = NULL, y_max = NULL,  
                           y_lab = NULL, legend_position = "top") {
  # Subset the data for the specified size
  data_subset <- data[data$size == size_value, ]
  
  p <- ggplot(data_subset, aes(x = treatment, y = !!sym(variable_mean), color = light)) +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se),
                      ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
    scale_fill_manual(values = c("#D5131733", "#0094CD33")) +
    scale_x_discrete(breaks = unique(data_subset$treatment)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
         y = y_lab) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.position = legend_position,
          plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", 
                                              halign = 0.5, lineheight = 0.9, linewidth = 0.25, 
                                              linetype = 1, padding = margin(4.5, 3, 4.5, 3), 
                                              box.color = "black")) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE))
  
  # If no y-axis label is provided, remove the axis title entirely
  if (is.null(y_lab) || y_lab == "") {
    p <- p + theme(axis.title.y = element_blank())
  }
  
  return(p)
}


Fe_0.2_plot <- generate_plots(FeC_noAZ_summary, "0.2-2", "Fe_mean", "Fe_se","Fe uptake", y_min = 0, y_max = 0.08, 
                              y_lab = expression(atop("Fe uptake per cell", "(amol cell"^-1~" day"^-1~")")), show_legend = TRUE, show_xlab = FALSE) +
                            ggtitle("0.2-2 μm") + theme(axis.title.x = element_blank())
Fe_2_plot <- generate_plots(FeC_noAZ_summary, "2-20", "Fe_mean", "Fe_se","Fe uptake", y_min = 0, y_max = 1.5, 
                            y_lab = "", show_legend = TRUE, show_xlab = TRUE) +
                           ggtitle("2-20 μm") + theme(axis.title.x = element_blank())
Fe_20_plot <- generate_plots(FeC_noAZ_summary, ">20", "Fe_mean", "Fe_se","Fe uptake", y_min = 0, y_max = 60, 
                             y_lab = "", show_legend = TRUE, show_xlab = TRUE) +
                              ggtitle("&gt;20 μm")
C_0.2_plot <- generate_plots(FeC_noAZ_summary, "0.2-2", "C_mean", "C_se","C uptake", y_min = 0, y_max = 0.002, 
                             y_lab = expression(atop("C uptake per cell", "(pmol cell"^-1~" day"^-1~")")), show_legend = TRUE, show_xlab = FALSE) +
                              ggtitle("0.2-2 μm")+ theme(axis.title.x = element_blank())
C_2_plot <- generate_plots(FeC_noAZ_summary, "2-20", "C_mean", "C_se","C uptake", y_min = 0, y_max = 0.1, 
                           y_lab = "", show_legend = TRUE, show_xlab = TRUE) +
                              ggtitle("2-20 μm")+ theme(axis.title.x = element_blank())
C_20_plot <- generate_plots(FeC_noAZ_summary, ">20", "C_mean", "C_se","C uptake", y_min = 0, y_max = 12, 
                            y_lab = "", show_legend = TRUE, show_xlab = TRUE) +
                              ggtitle("&gt;20 μm")

#---adding asterisks for significance ----------------------------------------------------------------
Feasterisks_0.2 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(3.37, 4.45, 4.98, 0.08, -0.24, 0.32),
  label = c("a", "b", "b", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_2 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(4.2, 6.84, 7.89, 6.56, 4.1, 5.18),
  label = c("a", "b", "c", "a", "b", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_20 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(10.6, 11.3, 11.8, 6.48, 7.52, 8.69),
  label = c("a,a", "a,b", "b,b", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by SIZE -----------------------------------------------
Fe_0.2_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                        cond = list(size="0.2-2"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-1, 5), breaks = seq(-1, 5, by = 2),
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe uptake)")), breaks = seq(-1, 5, by = 2))) +
  ggtitle("0.2-2 μm") +
  geom_text(data = Feasterisks_0.2, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_0.2$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") 


Fe_2_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                         cond=list(size="2-20"), scale = "linear", partial = TRUE, 
                         overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(4, 8), breaks = seq(4, 8, by = 1),
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe uptake)")))) +
  ggtitle("2-20 μm") +
  geom_text(data = Feasterisks_2, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_2$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") 

Fe_20_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                       cond=list(size=">20"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(6, 12), breaks = seq(6, 12, by = 2),
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe uptake)")))) +
  ggtitle("&gt;20 μm") +
  geom_text(data = Feasterisks_20, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_20$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("Treatment")

#---adding asterisks for significance for plots by size C uptake ----------------------------------------------------------------
Casterisks_0.2 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(0.1, -0.3, 1, -3.9, -3.45, -3.4),
  label = c("a", "a", "b", "a", "b", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_2 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(4.03, 3.55, 4.6, 1.5, 1.7, 1.29), 
  label = c("a", "b", "a", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_20 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(9.66, 9.2, 9.78, 6.5, 6.7, 6.2),
  label = c("a", "b", "a", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
C_0.2_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                        cond = list(size="0.2-2"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-4, 1), breaks = seq(-4, 1, by = 1),
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")))) +
  ggtitle("0.2-2 μm") +
  geom_text(data = Casterisks_0.2, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_0.2$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") 


C_2_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                      cond=list(size="2-20"), scale = "linear", partial = TRUE, 
                      overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(), axis.title.x = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, by = 1),
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")))) +
  ggtitle("2-20 μm") +
  geom_text(data = Casterisks_2, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_2$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") 

C_20_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                       cond=list(size=">20"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(6, 10), breaks = seq(6, 10, by = 1),
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.text.y.left = element_blank(),
        legend.position = "none", axis.ticks.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  ggtitle("&gt;20 μm") +
  geom_text(data = Casterisks_20, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_20$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("Treatment")

#----- Combining all in one plot ----------------------------------------------------------------------------------------------
FeC_percell_Feplot <- guide_area() / (Fe_0.2_plot| Fe_0.2_visreg) / (Fe_2_plot| Fe_2_visreg) / (Fe_20_plot| Fe_20_visreg) + 
                      plot_annotation(tag_levels = list(c('A','B','C','D','E','F'))) +
                      plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) & 
                      theme(axis.text = element_text(size = 10, colour = "black"),
                            legend.title = element_blank(), legend.text = element_text(size = 9), 
                            legend.direction = "horizontal", legend.position = "top",
                            plot.tag.position = c(0.02, 0.97)) & 
                      guides(fill = "none")

FeC_percell_Cplot <- guide_area() / (C_0.2_plot| C_0.2_visreg) / (C_2_plot| C_2_visreg) / (C_20_plot| C_20_visreg) +
  plot_annotation(tag_levels = list(c('A','B','C','D','E','F'))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) & 
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.02, 0.97)) & 
  guides(fill = "none")

# printing both the by-treatment and by-size plots
FeC_percell_Feplot
FeC_percell_Cplot
# Save the plots 
ggsave("FeC_percell_Fe_noAZ.svg", FeC_percell_Feplot, width = 7.2, height = 8)
ggsave("FeC_percell_C_noAZ.svg", FeC_percell_Cplot, width = 7.2, height = 8)

FeC_percell__allplot <- guide_area()  / (Fe_0.2_plot|Fe_2_plot|Fe_20_plot) / (C_0.2_plot| C_2_plot|C_20_plot) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) & 
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.02, 0.97)) & 
  guides(fill = "none")

ggsave("FeC_percell_noAZ.svg", FeC_percell__allplot, width = 7.2, height = 4.9)
