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

#remove per cell data and fraction data
FeC <- FeC %>% select(-Fe_cell, -C_cell, -Fe_fract, -C_fract)

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
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), color = light)) +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
    scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
    facet_grid(~ size) +
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
Fe_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 60) +
  theme(axis.title.x = element_blank(), legend.position = "top") +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " d"^-1, ")")))
C_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "C_mean", "C_se", "", y_min = 0, y_max = 2) +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " d"^-1, ")")))
FeC_plot_errorbars <- generate_plots_with_errorbars(FeC_noAZ_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0, y_max = 125) +
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
ggsave("FeC_dataplots_noAZ2.svg", grid_plot_errorbars, width = 7.2, height = 8)


# ---Compute the mean and variance for each group defined by light, treatment, and size-----------------------
FeC_noAZ_summary <- FeC_noAZ %>%
  group_by(light, treatment, size) %>%
  summarise(mean_Fe_up = mean(Fe_up, na.rm = TRUE),
            var_Fe_up = var(Fe_up, na.rm = TRUE),
            mean_C_up = mean(C_up, na.rm = TRUE),
            var_C_up = var(C_up, na.rm = TRUE),
            mean_FeC = mean(FeC, na.rm = TRUE),
            var_FeC = var(FeC, na.rm = TRUE))

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
scatterplot_Fe_up <- create_scatterplot(FeC_noAZ_summary, "mean_Fe_up", "var_Fe_up", "Fe Variance vs. Mean")
scatterplot_C_up <- create_scatterplot(FeC_noAZ_summary, "mean_C_up", "var_C_up", "C Variance vs. Mean")
scatterplot_FeC <- create_scatterplot(FeC_noAZ_summary, "mean_FeC", "var_FeC", "Fe:C Variance vs. Mean")

# Display scatterplots
scatterplot_Fe_up | scatterplot_C_up | scatterplot_FeC 

#------FITTING MODELS :  -----------------------------------------------------------------------------------
Fe_model  <- glm((Fe_up) ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model <- glm((C_up) ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))
FeC_model  <- glm((FeC) ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"))

Fe_model_reduced <- glm(Fe_up ~ light + treatment + size + light:treatment + light:size + 
                          treatment:size, data = FeC_noAZ, family = Gamma(link = "log"))
C_model_reduced <- glm(C_up ~ light + treatment + size + light:treatment + light:size + 
                         treatment:size, data = FeC_noAZ, family = Gamma(link = "log"))
FeC_model_reduced <- glm(FeC ~ light + treatment + size + light:treatment + light:size + 
                           treatment:size, data = FeC_noAZ, family = Gamma(link = "log"))

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_up ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced <- stan_glm(Fe_up ~ light + treatment + size + light:treatment + light:size + treatment:size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_up ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced <- stan_glm(C_up ~ light + treatment + size + light:treatment + light:size + treatment:size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
FeC_stan_full <- stan_glm(FeC ~ light * treatment * size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)
FeC_stan_reduced <- stan_glm(FeC ~ light + treatment + size + light:treatment + light:size + treatment:size, data = FeC_noAZ, family = Gamma(link = "log"), iter = 5000)

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

#posterior predictive checks
Fe_pp_full <- pp_check(Fe_stan_full, plotfun = "stat_grouped", stat = "mean", group = "treatment") +
  theme(legend.position = c(0.05, 0.95), legend.justification = c(0, 1)) +
  theme_bw() 
Fe_pp_reduced <- pp_check(Fe_stan_reduced, plotfun = "stat_grouped", stat = "mean", group = "treatment") +
  theme(legend.position = "none") +
  theme_bw() 
plot_grid(Fe_pp_full, Fe_pp_reduced, nrow = 2, align = "hv", axis = "v", labels = c("A", "B"))

#mean vs stdev
pp_check(Fe_stan_full, plotfun = "stat_2d", stat = c("mean", "sd")) +
  theme_bw()

#compare Bayes to frequentist model
coefficients(Fe_model)
coefficients(Fe_stan_full)

# -----Create a list of models with their corresponding names--------------------------
models_list <- list(
  Fe = Fe_model,
  C = C_model,
  FeC = FeC_model,
  Fe_2 = Fe_model_reduced,
  C_2 = C_model_reduced,
  FeC_2 = FeC_model_reduced
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
emm_all_Fe <- summary(emmeans(Fe_stan_full, specs = pairwise ~ treatment : light : size), type = "unlink")
emm_trt_Fe <- summary(emmeans(Fe_stan_full, specs = pairwise ~ treatment | light | size), type = "unlink")
emm_light_Fe <- summary(emmeans(Fe_stan_full, specs = pairwise ~ light | treatment | size), type = "unlink")
emm_size_Fe <- summary(emmeans(Fe_stan_full, specs = pairwise ~ size | treatment | light), type = "unlink")
# emmeans for C uptake
emm_all_C <- summary(emmeans(C_stan_reduced, specs = pairwise ~ treatment : light : size), type = "unlink")
emm_trt_C <- summary(emmeans(C_stan_reduced, specs = pairwise ~ treatment | light | size), type = "unlink")
emm_light_C <- summary(emmeans(C_stan_reduced, specs = pairwise ~ light | treatment | size), type = "unlink")
emm_size_C <- summary(emmeans(C_stan_reduced, specs = pairwise ~ size |  treatment | light ), type = "unlink")
# emmeans for Fe:C ratio
emm_all_FeC <- summary(emmeans(FeC_stan_full, specs = pairwise ~ treatment : light : size), type = "unlink")
emm_trt_FeC <- summary(emmeans(FeC_stan_full, specs = pairwise ~ treatment | light | size), type = "unlink")
emm_light_FeC <- summary(emmeans(FeC_stan_full, specs = pairwise ~ light | treatment | size), type = "unlink")
emm_size_FeC <- summary(emmeans(FeC_stan_full, specs = pairwise ~ size | treatment | light), type = "unlink")

# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, trt = emm_trt_Fe$contrasts, light = emm_light_Fe$contrasts,
                        size = emm_size_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, trt = emm_trt_C$contrasts, light = emm_light_C$contrasts,
                        size = emm_size_C$contrasts)
emmeans_FeC <- bind_rows(emmeans = emm_all_FeC$emmeans, trt = emm_trt_FeC$contrasts, light = emm_light_FeC$contrasts,
                        size = emm_size_FeC$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C,
  emmeans_FeC = emmeans_FeC
)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_noAZ_gamstan.xlsx")


# ----- plotting model results in vireg -------------------------------------------------------------------

#---adding asterisks for significance ----------------------------------------------------------------
Feasterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(3, 1.4, 1.6, 0.2, -1.1, -2.8), 
  label = c("a", "b", "b", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(4, 2.2, 2.5, -0.9, -2.4, -2), 
  label = c("a", "b", "b", "a", "b", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(4.7, 3.9, 4.9, -0.2, -0.7, -0.3), 
  label = c("a", "b", "a", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
Fe_DFB_visreg <- visreg(Fe_model, "size", by = "light", data = FeC_noAZ, 
                        cond = list(treatment="+DFB"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("+DFB") +
  geom_text(data = Feasterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(Fe uptake)"))) +
  xlab("")

Fe_Ctrl_visreg <- visreg(Fe_model, "size", by = "light", data = FeC_noAZ, 
                      cond=list(treatment="Control"), scale = "linear", partial = TRUE, 
                      overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("Control") +
  geom_text(data = Feasterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

Fe_Fe_visreg <- visreg(Fe_model, "size", by = "light", data = FeC_noAZ, 
                       cond=list(treatment="+Fe"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=-1),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("+Fe") +
  geom_text(data = Feasterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the two rows into a single grid
Fe_row <- (Fe_DFB_visreg | Fe_Ctrl_visreg | Fe_Fe_visreg) & theme(axis.title.x = element_blank())

#---adding asterisks for significance for plots by size C uptake ----------------------------------------------------------------
Casterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(-0, 0, 1.15, -4.35, -4.4, -3.7), 
  label = c("a", "a", "b", "a", "a", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(-0.6, -0.8, 0.8, -4.4, -4.65, -3.7), 
  label = c("a", "a", "b", "a", "a", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(1.2, 0.85, 2.9, -4.4, -4.8, -3.2), 
  label = c("a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
C_DFB_visreg <- visreg(C_model_reduced, "size", by = "light", data = FeC_noAZ, 
                        cond = list(treatment="+DFB"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("+DFB") +
  geom_text(data = Casterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(C uptake)"))) +
  xlab("")

C_Ctrl_visreg <- visreg(C_model_reduced, "size", by = "light", data = FeC_noAZ, 
                         cond=list(treatment="Control"), scale = "linear", partial = TRUE, 
                         overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.ticks.y = element_blank(),
        legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("Control") +
  geom_text(data = Casterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

C_Fe_visreg <- visreg(C_model_reduced, "size", by = "light", data = FeC_noAZ, 
                       cond=list(treatment="+Fe"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=-1),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("+Fe") +
  geom_text(data = Casterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the two rows into a single grid
C_row <- (C_DFB_visreg | C_Ctrl_visreg | C_Fe_visreg) & theme(axis.title.x = element_blank())

#---adding asterisks for significance Fe:C -----------------------------------------------------------------------
FeCasterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(2.2, 0.7, -0.3, 5.5, 4.3, 1.9), 
  label = c("a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

FeCasterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(5.5, 4, 2.7, 2.7, 1.7, 0.8), 
  label = c("a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

FeCasterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(2.7, 2.2 ,1.2, 5, 5.2, 4), 
  label = c("a", "a", "b", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for FeC ratio by treatment -----------------------------------------------
FeC_DFB_visreg <- visreg(FeC_model, "size", by = "light", data = FeC_noAZ, 
                        cond = list(treatment="+DFB"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("+DFB") +
  geom_text(data = FeCasterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(Fe:C ratio)"))) +
  xlab("")

FeC_Ctrl_visreg <- visreg(FeC_model, "size", by = "light", data = FeC_noAZ, 
                         cond=list(treatment="Control"), scale = "linear", partial = TRUE, 
                         overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0), 
        legend.position = "none", axis.text.y = element_blank(),  
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("Control") +
  geom_text(data = FeCasterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("Size (μm)")

FeC_Fe_visreg <- visreg(FeC_model, "size", by = "light", data = FeC_noAZ, 
                       cond=list(treatment="+Fe"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=0),
        legend.position = "none", axis.text.y = element_blank(),  
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("+Fe") +
  geom_text(data = FeCasterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine into a single grid
FeC_row <- (FeC_DFB_visreg | FeC_Ctrl_visreg | FeC_Fe_visreg)

#----- Combining all in one plot ----------------------------------------------------------------------------------------------
FeC_treatment_plot <- guide_area() / Fe_row / C_row / FeC_row + 
                      plot_annotation(tag_levels = list(c('A','','','B','','','C','',''))) +
                      plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) & 
                      theme(axis.text = element_text(size = 10, colour = "black"),
                      legend.title = element_blank(), legend.text = element_text(size = 9), 
                      legend.direction = "horizontal", legend.position = "top",
                      plot.tag.position = c(0.02, 0.96)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


# ----- plotting model results NOW BY SIZE INSTEAD------------------------------------------------------------------------------------------------------

#---adding asterisks for significance ----------------------------------------------------------------
Feasterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(3, 4, 4.7, 0.2, -0.9, -0.2),
  label = c("a", "b", "c", "a,a", "b,b", "a,b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(1.4, 2.2, 3.9, -1.1, -2.4, -0.7),
  label = c("a", "b", "c", "a", "b", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Feasterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(1.6, 2.5, 4.9, -2.8, -2, -0.3),
  label = c("a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by SIZE -----------------------------------------------
Fe_0.2_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                        cond = list(size="0.2-2"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("0.2-2 μm") +
  geom_text(data = Feasterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(Fe uptake)"))) +
  xlab("")

Fe_2_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                         cond=list(size="2-20"), scale = "linear", partial = TRUE, 
                         overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("2-20 μm") +
  geom_text(data = Feasterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

Fe_20_visreg <- visreg(Fe_model, "treatment", by = "light", data = FeC_noAZ, 
                       cond=list(size=">20"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=-1),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-3, 5), breaks = seq(-3, 5, by = 2)) +
  ggtitle("&gt;20 μm") +
  geom_text(data = Feasterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = Feasterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the two rows into a single grid
Fe_row <- (Fe_0.2_visreg | Fe_2_visreg | Fe_20_visreg) & theme(axis.title.x = element_blank())

#---adding asterisks for significance for plots by size C uptake ----------------------------------------------------------------
Casterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(0, -0.6, 1.2, -4.35, -4.4, -4.4),
  label = c("a", "b", "c", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(0, -0.8, 0.85, -4.4, -4.65, -4.8), 
  label = c("a", "b", "c", "a,a", "a,b", "b,b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

Casterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(1.15, 0.8, 2.9, -3.7, -3.7, -3.2),
  label = c("a", "b", "c", "a", "a", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for Fe uptake by treatment -----------------------------------------------
C_0.2_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                       cond = list(size="0.2-2"), scale = "linear", partial = TRUE, 
                       overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("0.2-2 μm") +
  geom_text(data = Casterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(C uptake)"))) +
  xlab("")

C_2_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                        cond=list(size="2-20"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.ticks.y = element_blank(),
        legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("2-20 μm") +
  geom_text(data = Casterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

C_20_visreg <- visreg(C_model_reduced, "treatment", by = "light", data = FeC_noAZ, 
                      cond=list(size=">20"), scale = "linear", partial = TRUE, 
                      overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=-1),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-5, 3), breaks = seq(-4, 2, by = 2)) +
  ggtitle("&gt;20 μm") +
  geom_text(data = Casterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = Casterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine the two rows into a single grid
C_row <- (C_0.2_visreg | C_2_visreg | C_20_visreg) & theme(axis.title.x = element_blank())

#---adding asterisks for significance Fe:C -----------------------------------------------------------------------
FeCasterisks_DFB <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(2.2, 5.5, 2.7, 5.5, 2.7, 5),
  label = c("a", "b", "c", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

FeCasterisks_Ctrl <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(0.7, 4, 2.2, 4.3, 1.7, 5.2), 
  label = c("a", "b", "b", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

FeCasterisks_Fe <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(-0.3, 2.7 , 1.2, 1.9, 0.8, 4),
  label = c("a", "b", "b", "a", "b", "c"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

# ----- plotting model means and residuals for FeC ratio by treatment -----------------------------------------------
FeC_0.2_visreg <- visreg(FeC_model, "treatment", by = "light", data = FeC_noAZ, 
                         cond = list(size="0.2-2"), scale = "linear", partial = TRUE, 
                         overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(margin = margin(t = 0, r = 7.5, b = 0, l = 0)),
        legend.position = "none",
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=5.5),
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("0.2-2 μm") +
  geom_text(data = FeCasterisks_DFB, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_DFB$colours, size = 4, inherit.aes = FALSE)  +
  ylab(expression(paste("ln(Fe:C ratio)"))) +
  xlab("")

FeC_2_visreg <- visreg(FeC_model, "treatment", by = "light", data = FeC_noAZ, 
                          cond=list(size="2-20"), scale = "linear", partial = TRUE, 
                          overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        plot.margin = margin(t=5.5, r=5.5, b=5.5, l=0), 
        legend.position = "none", axis.text.y = element_blank(),  
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("2-20 μm") +
  geom_text(data = FeCasterisks_Ctrl, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_Ctrl$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("Treatment")

FeC_20_visreg <- visreg(FeC_model, "treatment", by = "light", data = FeC_noAZ, 
                        cond=list(size=">20"), scale = "linear", partial = TRUE, 
                        overlay = TRUE, gg = TRUE, points = list(shape = 19)) +
  theme_bw() +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = -15, b = 0, l = 0)),
        plot.margin = margin(t=5.5, b=5.5, r=5.5, l=0),
        legend.position = "none", axis.text.y = element_blank(),  
        panel.grid.minor = element_blank(),
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, linewidth = 0.25,
                                            linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  scale_y_continuous(limits = c(-0.5, 5.5), breaks = seq(0, 5, by = 1)) +
  ggtitle("&gt;20 μm") +
  geom_text(data = FeCasterisks_Fe, aes(x = x_pos, y = y_pos, label = label), 
            colour = FeCasterisks_Fe$colours, size = 4, inherit.aes = FALSE)  +
  ylab("") +
  xlab("")

# Combine into a single grid
FeC_row <- (FeC_0.2_visreg | FeC_2_visreg | FeC_20_visreg) 

#----- Combining all in one plot ----------------------------------------------------------------------------------------------
FeC_size_plot <- guide_area() / Fe_row / C_row / FeC_row + 
                      plot_annotation(tag_levels = list(c('A','','','B','','','C','',''))) +
                      plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) & 
                      theme(axis.text = element_text(size = 10, colour = "black"),
                            legend.title = element_blank(), legend.text = element_text(size = 9), 
                            legend.direction = "horizontal", legend.position = "top",
                            plot.tag.position = c(0.04, 0.96)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# printing both the by-treatment and by-size plots
FeC_treatment_plot
FeC_size_plot
# Save the plots 
ggsave("FeC_visreg_treatment_noAZ.svg", FeC_treatment_plot, width = 7.2, height = 8)
ggsave("FeC_visreg_size_noAZ.svg", FeC_size_plot, width = 7.2, height = 8)
