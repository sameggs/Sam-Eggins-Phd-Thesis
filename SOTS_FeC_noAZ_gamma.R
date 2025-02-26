setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment") #setwd#loading packages
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
FeC <- read_csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/SOTS_FeC.csv") %>%
  mutate(inhibitor = as.factor(inhibitor) %>% fct_recode("+AZ" = "AZ"),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         size = as.factor(size) %>% fct_recode("0.2-2" = "0.2", "2-20" = "2", ">20" = "20"),
         bottle = as.factor(bottle))

#remove per cell data and fraction data
FeC <- FeC %>% select(-Fe_cell, -C_cell)

# Reorder  levels
FeC$treatment <- factor(FeC$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
FeC$inhibitor <- factor(FeC$inhibitor, levels = c("Control", "+AZ"), ordered = TRUE)

# log transform Fe and C variables
FeC$Fe_up_ln <- log(FeC$Fe_up)
FeC$C_up_ln <- log(FeC$C_up)
FeC$FeC_ln <- log(FeC$FeC)
# square root transform Fe and C variables
FeC$Fe_up_sqrt <- sqrt(FeC$Fe_up)
FeC$C_up_sqrt <- sqrt(FeC$C_up)
FeC$FeC_sqrt <- sqrt(FeC$FeC)

#check the structure
str(FeC)

#FeC_noAZ 
FeC_noAZ <- FeC %>% filter(inhibitor == "Control")

# ---------------Plots of the data without inhibitor ----------------------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_noAZ_summary <- FeC_noAZ %>%
  group_by(treatment, size, light) %>%
  summarise(
    Fe_mean = mean(Fe_up),
    Fe_se = sd(Fe_up) / sqrt(n()),
    C_mean = mean(C_up),
    C_se = sd(C_up) / sqrt(n()),
    FeC_mean = mean(FeC),
    FeC_se = sd(FeC) / sqrt(n())
  )

# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, show_xlab = TRUE, y_min = NULL, y_max = NULL) {
  p <- ggplot(data, aes(x = size, y = !!sym(variable_mean))) +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    facet_grid(~ treatment) +
    scale_x_discrete(breaks = unique(data$size)) +
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
Fe_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 100) +
  theme(axis.title.x = element_blank(), legend.position = "top") +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " day"^-1, ")")))
C_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "C_mean", "C_se", "", y_min = 0, y_max = 5) +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " day"^-1, ")")))
FeC_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0, y_max = 40) +
  theme(legend.position = "none")

# Arrange plots in a grid
grid_plot_errorbars <- (guide_area() / Fe_plot_errorbars / C_plot_errorbars / FeC_plot_errorbars) + 
  plot_annotation(tag_levels = list(c('A','B','C'))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) &
  theme(panel.spacing.x = unit(8, "pt"), plot.margin = margin(t=5.5, b=5.5, r=5.5, l=2.5),
        legend.position = 'top', legend.text = element_text(size = 9), 
        legend.title = element_blank(), legend.direction = "horizontal",
        plot.tag.position = c(0.02, 0.96))
grid_plot_errorbars
# Save the plot to a file
ggsave("FeC_dataplots_noAZ.svg", grid_plot_errorbars, width = 7.2, height = 8)

#------ An alternative with size all folded together--------------------------------------------------------------------------------
# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, show_xlab = TRUE, y_min = NULL, y_max = NULL) {
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), color = size)) +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    scale_colour_manual(values = c("#EFD500FF", "#95C11FFF", "#007B3DFF")) +
    scale_fill_manual(values = c("#EFD50033", "#95C11F33", "#007B3D33"))  +
    scale_x_discrete(breaks = unique(data$treatment)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
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
Fe_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 100) +
  theme(axis.title.x = element_blank(), legend.position = "top") +
  geom_hline(yintercept = 43.3, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 15.8, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 9.0, linetype = "dashed", color = "#007B3DFF") +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " d"^-1, ")")))
C_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "C_mean", "C_se", "", y_min = 0, y_max = 5) +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 0.86, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 1.83, linetype = "dashed", color = "#007B3DFF") +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " d"^-1, ")")))
FeC_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0, y_max = 60) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 15)) +
  geom_hline(yintercept = 56.5, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 18.3, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 5.8, linetype = "dashed", color = "#007B3DFF") +
  theme(legend.position = "none")

#------FITTING MODELS :  -----------------------------------------------------------------------------------
Fe_model  <- glm((Fe_up) ~  treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model <- glm((C_up) ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
FeC_model  <- glm((FeC) ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))

Fe_model_reduced  <- glm((Fe_up) ~  treatment + size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model_reduced <- glm((C_up) ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"))
FeC_model_reduced  <- glm((FeC) ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"))

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_up ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced <- stan_glm(Fe_up ~ treatment + size , data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_up ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced <- stan_glm(C_up ~ treatment + treatment, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
FeC_stan_full <- stan_glm(FeC ~ treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
FeC_stan_reduced <- stan_glm(FeC ~ treatment + size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for each model
Fe_loo_full <- loo(Fe_stan_full, k_threshold = 0.7)
Fe_loo_reduced <- loo(Fe_stan_reduced, k_threshold = 0.7)
C_loo_full <- loo(C_stan_full, k_threshold = 0.7)
C_loo_reduced <- loo(C_stan_reduced, k_threshold = 0.7)
FeC_loo_full <- loo(FeC_stan_full, k_threshold = 0.7)
FeC_loo_reduced <- loo(FeC_stan_reduced, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(Fe_loo_full, Fe_loo_reduced)
loo_compare(C_loo_full, C_loo_reduced)
loo_compare(FeC_loo_full, FeC_loo_reduced)

# -----Create a list of models with their corresponding names--------------------------
models_list <- list(
  Fe = Fe_model,
  C = C_model,
  FeC = FeC_model,
  Fe_noint = Fe_model_reduced,
  C_noint = C_model_reduced,
  FeC_noint = FeC_model_reduced
)

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
qq_plot_patchwork <- wrap_plots(qq_plot_list, ncol = 3)

# Print the plot
qq_plot_patchwork


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

# Arrange the plots in a 3x3 grid
plot_grid(
  plotlist = res_fit_plot_list, align = "hv",
  ncol = 3, nrow = 2
)


# --- interrogating the model results using EMMEANS ---------------------------------------------------
# emmeans for Fe uptake
emm_all_Fe <- summary(emmeans(Fe_model_reduced, specs = pairwise ~ treatment : size), type = "unlink")
emm_trt_Fe <- summary(emmeans(Fe_model_reduced, specs = pairwise ~ treatment | size), type = "unlink")
emm_size_Fe <- summary(emmeans(Fe_model_reduced, specs = pairwise ~ size | treatment), type = "unlink")
# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model_reduced, specs = pairwise ~ treatment : size), type = "unlink")
emm_trt_C <- summary(emmeans(C_model_reduced, specs = pairwise ~ treatment | size), type = "unlink")
emm_size_C <- summary(emmeans(C_model_reduced, specs = pairwise ~ size | treatment), type = "unlink")
# emmeans for Fe:C ratio
emm_all_FeC <- summary(emmeans(FeC_model_reduced, specs = pairwise ~ treatment : size), type = "unlink")
emm_trt_FeC <- summary(emmeans(FeC_model_reduced, specs = pairwise ~ treatment | size), type = "unlink")
emm_size_FeC <- summary(emmeans(FeC_model_reduced, specs = pairwise ~ size | treatment), type = "unlink")

# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, trt = emm_trt_Fe$contrasts, size = emm_size_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, trt = emm_trt_C$contrasts, size = emm_size_C$contrasts)
emmeans_FeC <- bind_rows(emmeans = emm_all_FeC$emmeans, trt = emm_trt_FeC$contrasts, size = emm_size_FeC$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C,
  emmeans_FeC = emmeans_FeC
)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_noAZ_gamstan_red.xlsx")


# ----- plotting model results in vireg -------------------------------------------------------------------
#---adding asterisks for significance ----------------------------------------------------------------
Fe_asterisks_size <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(4.9, 4.5, 4, 3.7, 2.9, 2.5, 2.3, 1.8, 1.4), 
  label = c("a", "b", "c", "a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#000000", "#000000", "#000000", "#164194", "#164194","#164194"))

Fe_asterisks_trt <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(3.4, 3.5, 4.9, 2.7, 2.8, 4.2, 1.5, 1.6, 2.9), 
  label = c("a", "a", "b", "a", "a", "b", "a", "a", "b"), colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D","#007B3D"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
Fe_visreg_size <- visreg(Fe_model_reduced, "size", by = "treatment", data = FeC_noAZ, 
                        scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_color_manual(values = c("#164194FF", "#000000FF", "#D51317FF")) +
  scale_fill_manual(values = c("#16419433", "#00000033", "#D5131733"))  +
  scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, by = 1), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe uptake)")), breaks = seq(1, 5, by = 1))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = Fe_asterisks_size, aes(x = x_pos, y = y_pos, label = label), 
            colour = Fe_asterisks_size$colours, size = 4, inherit.aes = FALSE)  +
  xlab("")
  

Fe_visreg_trt <- visreg(Fe_model_reduced, "treatment", by = "size", data = FeC_noAZ, 
                         scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33")) +
  scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, by = 1), 
                   sec.axis = sec_axis(~., name = expression(paste("ln(",Delta,"Chl:C)")), breaks = seq(1, 5, by = 1))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  
  geom_text(data = Fe_asterisks_trt, aes(x = x_pos, y = y_pos, label = label), 
            colour = Fe_asterisks_trt$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the plots into a single row
Fe_row <- (Fe_plot_errorbars | Fe_visreg_trt | Fe_visreg_size) & theme(axis.title.x = element_blank())

#---adding asterisks for significance ----------------------------------------------------------------
C_asterisks_size <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(1.4, 1.6, 1.7, -0.7, -0.7, 0, 0.7, 0.5, 1.1), 
  label = c("a", "a", "a", "a", "a", "b", "a", "a", "b"), colours = c("#D51317", "#D51317", "#D51317", "#000000", "#000000", "#000000", "#164194", "#164194","#164194"))

C_asterisks_trt <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(-0.2, -0.7, 0.5, -0.4, -0.9, 0.7, 1.2, 0.7, 1.7), 
  label = c("a", "b", "c", "a", "b", "c", "a", "b", "c"), colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D","#007B3D"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
C_visreg_size <- visreg(C_model, "size", by = "treatment", data = FeC_noAZ, 
                         scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_color_manual(values = c("#164194FF", "#000000FF", "#D51317FF")) +
  scale_fill_manual(values = c("#16419433", "#00000033", "#D5131733"))  +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")), breaks = seq(-2, 2, by = 1))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = C_asterisks_size, aes(x = x_pos, y = y_pos, label = label), 
            colour = C_asterisks_size$colours, size = 4, inherit.aes = FALSE)  +
  xlab("")

C_visreg_trt <- visreg(C_model, "treatment", by = "size", data = FeC_noAZ, 
                        scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33")) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
      axis.title.y.right = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(C uptake)")), breaks = seq(-2, 2, by = 1))) +
  geom_text(data = C_asterisks_trt, aes(x = x_pos, y = y_pos, label = label), 
            colour = C_asterisks_trt$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the plots into a single row
C_row <- (C_plot_errorbars | C_visreg_trt | C_visreg_size) & theme(axis.title.x = element_blank())

#---adding asterisks for significance ----------------------------------------------------------------
FeC_asterisks_size <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(3.1, 3.4, 2.4, 3.9, 2.4, 1.8, 2.3, 1.7, 0.8), 
  label = c("a", "b", "c", "a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#000000", "#000000", "#000000", "#164194", "#164194","#164194"))

FeC_asterisks_trt <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), y_pos = c(3, 3.9, 3.8, 1.8, 2.4, 2.7, 0.8, 1.2, 1.85), 
  label = c("a", "a", "b", "a", "a", "b", "a", "a", "b"), colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D","#007B3D"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
FeC_visreg_size <- visreg(FeC_model, "size", by = "treatment", data = FeC_noAZ, 
                         scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_color_manual(values = c("#164194FF", "#000000FF", "#D51317FF")) +
  scale_fill_manual(values = c("#16419433", "#00000033", "#D5131733"))  +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(Fe:C)")), breaks = seq(0, 4, by = 1))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = FeC_asterisks_size, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeC_asterisks_size$colours, size = 4, inherit.aes = FALSE)  +
  xlab("Size (μm)")

FeC_visreg_trt <- visreg(FeC_model, "treatment", by = "size", data = FeC_noAZ, 
                        scale = "linear", overlay = TRUE, partial = TRUE, gg = TRUE, points.par = list(shape = 1, size = 2)) +
  theme_bw() +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33")) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1), 
                   sec.axis = sec_axis(~., name = expression(paste("ln(Fe:C)")), breaks = seq(0, 4, by = 1))) +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y.right = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank()) +
  geom_text(data = FeC_asterisks_trt, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeC_asterisks_trt$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("Treatment")

# Combine the plots into a single row
FeC_row <- (FeC_plot_errorbars | FeC_visreg_trt | FeC_visreg_size) 

#----- Combining all in one plot ----------------------------------------------------------------------------------------------
FeC_models_plot <- guide_area() / Fe_row / C_row / FeC_row + 
  plot_annotation(tag_levels = list(c('A','B','','C','D','','E','F',''))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(2,10,10,10)) & 
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.04, 0.96)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("FeC_models_noAZ_SOTS.svg", FeC_models_plot, width = 7.2, height = 8.26)


#----- NOT SIZE FRACTIONATED-----------
# Step 1: Sum the values for each replicate (each bottle) over all sizes
FeC_noAZ_reps <- FeC_noAZ %>%
  group_by(treatment, light, bottle) %>%
  summarise(
    Fe_total  = sum(Fe_up, na.rm = TRUE),
    C_total   = sum(C_up, na.rm = TRUE),
    FeC_total = sum(FeC, na.rm = TRUE),
    .groups   = "drop"
  )

# Step 2: For each treatment and light combination, calculate the mean and SE across replicates
FeC_noAZ_summary <- FeC_noAZ_reps %>%
  group_by(treatment, light) %>%
  summarise(
    Fe_mean   = mean(Fe_total, na.rm = TRUE),
    Fe_se     = sd(Fe_total, na.rm = TRUE) / sqrt(n()),
    C_mean    = mean(C_total, na.rm = TRUE),
    C_se      = sd(C_total, na.rm = TRUE) / sqrt(n()),
    FeC_mean  = mean(FeC_total, na.rm = TRUE),
    FeC_se    = sd(FeC_total, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )
