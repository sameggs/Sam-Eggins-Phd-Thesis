setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment") #setwd
#loading packages
library(tidyverse)
library(broom)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(cowplot)
library(svglite)
library(car)
library(emmeans)
library(visreg)
library(writexl)
library(rstan)
library(rstanarm)

#-----------import the data --------------------------------------------------------------------------------------
FeC <- read_csv("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment/FeC data/dcm_FeC.csv") %>%
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

# ---------------Plots of the data with inhibitor ----------------------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_summary <- FeC %>%
  group_by(treatment, size, light, inhibitor) %>%
  summarise(
    Fe_mean = mean(Fe_up),
    Fe_se = sd(Fe_up) / sqrt(n()),
    C_mean = mean(C_up),
    C_se = sd(C_up) / sqrt(n()),
    FeC_mean = mean(FeC),
    FeC_se = sd(FeC) / sqrt(n())
  )

# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, legend_to_show = NULL, show_xlab = TRUE, show_x_axis_labels = TRUE, y_min = NULL, y_max = NULL) {
  # Add a new column for x-axis position
  data$x_pos <- as.numeric(as.factor(data$size)) + ifelse(data$size == "2", -0.1, 0.1)
  
  p <- ggplot(data, aes(x = x_pos, y = data[[variable_mean]], color = light, shape = inhibitor)) +
    facet_grid(~ treatment) +
    geom_point(size = 3) + 
    geom_line(aes(group = interaction(size, light, treatment)), linewidth = 0.5, show.legend = FALSE) +
    geom_errorbar(aes(ymin = data[[variable_mean]] - data[[variable_se]], ymax = data[[variable_mean]] + data[[variable_se]]), width = 0.1) +
    scale_x_continuous(breaks = 1:length(unique(data$size)), labels = c("0.2" = "0.2-2", "2" = "2-20", "20" = ">20"), limits = c(0.5, length(unique(data$size)) + 0.5)) +
    scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
    scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
    scale_shape_manual(values = c("Control" = 1, "+AZ" = 19)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Size (μm)" else NULL,
         y = paste(variable_name),
         color = NULL, shape = NULL) + 
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"), 
          legend.text = element_text(size = 9),
          legend.direction = "horizontal",
          legend.justification = c(0, 1),
          legend.title = element_blank()) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    ) 
}

#---adding asterisks for significance ----------------------------------------------------------------
asterisks_Fe <- data.frame(
  x_pos = c(1.25, 2.25, 2.25, 3.25, 1.25, 2.25, 1.25, 2.25), y_pos = c(0.95, 0.964, 0.1, 0.115, 0.68, 0.18, 0.94, 33.94), 
  label = c("\u2193", "\u2193", "\u2193", "\u2193", "\u2193", "\u2193", "\u2193", "\u2191"), treatment = c("+DFB", "+DFB", "+DFB", "+DFB", "Control", "Control", "+Fe", "+Fe"))

asterisks_C <- data.frame(
  x_pos = c(1.25, 2.25, 3.25, 2.25, 2.25, 3.25), y_pos = c(0.23, 0.165, 0.69, 0.015, 1.31, 0.13), 
  label = c("\u2193", "\u2193", "\u2193", "\u2193", "\u2191", "\u2191"), treatment = c("+DFB", "+DFB", "+DFB", "+DFB", "+Fe", "+Fe"))

asterisks_FeC <- data.frame(
  x_pos = c(1.25, 1.25, 2.25, 3.25, 2.25, 1.25, 2.25, 3.25), y_pos = c(55.53, 30.19, 7.18, 2.6, 0.436, 28.48, 36.37, 12.1), 
  label = c("\u2191", "\u2193", "\u2193", "\u2191", "\u2193", "\u2193", "\u2193", "\u2193"), treatment = c("+DFB", "+DFB", "+DFB", "+DFB", "Control", "+Fe", "+Fe", "+Fe"))

lightsig_FeC <- data.frame(
  x_pos = c(1.42, 1.4, 2.4, 3.42, 2.4, 1.4, 2.4, 3.4), y_pos = c(50.53, 30.19, 7.18, 2.35, 0.436, 28.48, 36.37, 12.1), 
  label = c("HL", "LL", "LL", "HL", "LL", "LL", "LL", "LL"), treatment = c("+DFB", "+DFB", "+DFB", "+DFB", "Control", "+Fe", "+Fe", "+Fe"))

asterisks_Fe$treatment <- factor(asterisks_Fe$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
asterisks_C$treatment <- factor(asterisks_C$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
asterisks_FeC$treatment <- factor(asterisks_FeC$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
lightsig_FeC$treatment <- factor(asterisks_FeC$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)

#-----------------------log scale plots-----------------------------------------------------
# Define label formatting function
format_y_labels <- function(x) {
  ifelse(x > 1, round(x), x)
}
# plots
Fe_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "Fe_mean", "Fe_se", "", y_min = 0.06, y_max = 130) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_Fe, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  annotation_logticks(sides = "l") +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " day"^-1, ")")))
C_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "C_mean", "C_se", "", y_min = 0.006, y_max = 13) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_C, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  annotation_logticks(sides = "l") +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " day"^-1, ")")))
FeC_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0.6, y_max = 130) +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_FeC, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  geom_text(data = lightsig_FeC, aes(x = x_pos, y = y_pos, label = label), size = 3, inherit.aes = FALSE) +
  annotation_logticks(sides = "l")

# Arrange plots in a grid
grid_plot_errorbars_log <- (guide_area() / Fe_plot_errorbars_log / C_plot_errorbars_log / FeC_plot_errorbars_log) + 
  plot_annotation(tag_levels = list(c('A','B','C'))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) &
  theme(panel.spacing.x = unit(8, "pt"), plot.margin = margin(t=5.5, b=5.5, r=5.5, l=2),
        legend.position = 'top', legend.text = element_text(size = 9), 
        legend.title = element_blank(), legend.direction = "horizontal",
        plot.tag.position = c(0.02, 0.96))
grid_plot_errorbars_log

# save plot
ggsave("FeC_AZ_logscale_plot.svg", plot = grid_plot_errorbars_log, width = 7.2, height = 8)

#FITTING MODELS :  -----------------------------------------------------------------------------------
Fe_model  <- glm((Fe_up) ~ light * treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))
C_model <- glm((C_up) ~ light * treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))
FeC_model  <- glm(log(FeC) ~ light * treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))

Fe_model_reduced_3 <- glm(Fe_up ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment + treatment:inhibitor
                            light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                            data = FeC, family = Gamma(link = "log"))
Fe_model_reduced_2 <- glm(Fe_up ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                            data = FeC, family = Gamma(link = "log"))
C_model_reduced_3 <- glm(C_up ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment + treatment:inhibitor
                            light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                          data = FeC, family = Gamma(link = "log"))
C_model_reduced_2 <- glm(C_up ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                          data = FeC, family = Gamma(link = "log"))
FeC_model_reduced_3 <- glm(FeC ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment + treatment:inhibitor
                            light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                          data = FeC, family = Gamma(link = "log"))
FeC_model_reduced_2 <- glm(FeC ~ light + treatment + size + inhibitor +
                            light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                          data = FeC, family = Gamma(link = "log"))

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_up ~ light * treatment * size * inhibitor, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced_3 <- stan_glm(Fe_up ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment +
                                light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced_2 <- stan_glm(Fe_up ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_up ~ light * treatment * size * inhibitor, 
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced_3 <- stan_glm(C_up ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment +
                                light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_reduced_2 <- stan_glm(C_up ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_full <- stan_glm(FeC ~ light * treatment * size * inhibitor, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_reduced_3 <- stan_glm(FeC ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment +
                                light:treatment:size + light:inhibitor:size + treatment:inhibitor:size + light:treatment:inhibitor,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_reduced_2 <- stan_glm(FeC ~ light + treatment + size + inhibitor +
                                light:treatment + light:size + treatment:size + inhibitor:size + inhibitor:light + inhibitor:treatment,
                                data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without inhibitor
Fe_stan_noIn <- stan_glm(Fe_up ~ light * treatment * size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noIn <- stan_glm(C_up ~ light * treatment * size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noIn <- stan_glm(FeC ~ light * treatment * size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
# models with inhibitor as Main effect only
Fe_stan_mainIn <- stan_glm(Fe_up ~ light * treatment * size + inhibitor, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_mainIn <- stan_glm(C_up ~ light * treatment * size + inhibitor, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_mainIn <- stan_glm(FeC ~ light * treatment * size + inhibitor, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
# models with inhibitor*size 
Fe_stan_sizeIn <- stan_glm(Fe_up ~ light * treatment * size + inhibitor*size, 
                           data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_sizeIn <- stan_glm(C_up ~ light * treatment * size + inhibitor*size, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_sizeIn <- stan_glm(FeC ~ light * treatment * size + inhibitor*size, 
                            data = FeC, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for each model
Fe_loo_full <- loo(Fe_stan_full, k_threshold = 0.7)
Fe_loo_reduced_3 <- loo(Fe_stan_reduced_3, k_threshold = 0.7)
Fe_loo_reduced_2 <- loo(Fe_stan_reduced_2, k_threshold = 0.7)
Fe_loo_noIn <- loo(Fe_stan_noIn, k_threshold = 0.7)
Fe_loo_mainIn <- loo(Fe_stan_mainIn, k_threshold = 0.7)
Fe_loo_sizeIn <- loo(Fe_stan_sizeIn, k_threshold = 0.7)

C_loo_full <- loo(C_stan_full, k_threshold = 0.7)
C_loo_reduced_3 <- loo(C_stan_reduced_3, k_threshold = 0.7)
C_loo_reduced_2 <- loo(C_stan_reduced_2, k_threshold = 0.7)
C_loo_noIn <- loo(C_stan_noIn, k_threshold = 0.7)
C_loo_mainIn <- loo(C_stan_mainIn, k_threshold = 0.7)
C_loo_sizeIn <- loo(C_stan_sizeIn, k_threshold = 0.7)

FeC_loo_full <- loo(FeC_stan_full, k_threshold = 0.7)
FeC_loo_reduced_3 <- loo(FeC_stan_reduced_3, k_threshold = 0.7)
FeC_loo_reduced_2 <- loo(FeC_stan_reduced_2, k_threshold = 0.7)
FeC_loo_noIn <- loo(FeC_stan_noIn, k_threshold = 0.7)
FeC_loo_mainIn <- loo(FeC_stan_mainIn, k_threshold = 0.7)
FeC_loo_sizeIn <- loo(FeC_stan_sizeIn, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(Fe_loo_full, Fe_loo_reduced_3, Fe_loo_reduced_2, Fe_loo_noIn, Fe_loo_mainIn, Fe_loo_sizeIn)
loo_compare(C_loo_full, C_loo_reduced_3, C_loo_reduced_2, C_loo_noIn, C_loo_mainIn, C_loo_sizeIn)
loo_compare(FeC_loo_full, FeC_loo_reduced_3, FeC_loo_reduced_2, FeC_loo_noIn, FeC_loo_mainIn, FeC_loo_sizeIn)

# -----Create a list of models with their corresponding names--------------------------
models_list <- list(
  Fe_full = Fe_model,
  C_full = C_model,
  FeC_full = FeC_model,
  Fe_3 = Fe_model_reduced_3,
  C_3 = C_model_reduced_3,
  FeC_3 = FeC_model_reduced_3,
  Fe_2 = Fe_model_reduced_2,
  C_2 = C_model_reduced_2,
  FeC_2 = FeC_model_reduced_2
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

# Arrange the plots in a grid with 2 rows
plot_grid(
  plotlist = qq_plot_list, align = "hv",
  ncol = 3, nrow = 3
)

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
  ncol = 3, nrow = 3
)

# --- interrogating the model results using EMMEANS ---------------------------------------------------
# emmeans for Fe uptake
emm_all_Fe <- summary(emmeans(Fe_model_reduced_2, specs = pairwise ~ inhibitor : treatment : light : size), type = "unlink")
emm_IinTS_Fe <- summary(emmeans(Fe_model_reduced_2, specs = pairwise ~ inhibitor | treatment : light : size), type = "unlink")

# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model_reduced_3, specs = pairwise ~ inhibitor : treatment : light : size), type = "unlink")
emm_IinTS_C <- summary(emmeans(C_model_reduced_3, specs = pairwise ~ inhibitor | treatment : light : size), type = "unlink")

# emmeans for Fe:C ratio
emm_all_FeC <- summary(emmeans(FeC_model_full, specs = pairwise ~ inhibitor : treatment : light : size), type = "unlink")
emm_IinTS_FeC <- summary(emmeans(FeC_model_full, specs = pairwise ~ inhibitor | treatment : light : size), type = "unlink")


# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, IinTS = emm_IinTS_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, IinTS = emm_IinTS_C$contrasts)
emmeans_FeC <- bind_rows(emmeans = emm_all_FeC$emmeans, IinTS = emm_IinTS_FeC$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C,
  emmeans_FeC = emmeans_FeC
)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_+AZ_gam.xlsx")





# ---------------Plots of the data without size ----------------------------------------------------------
FeC_bottle_totals <- FeC %>%
  group_by(treatment, light, inhibitor, bottle) %>%
  summarise(
    Fe_total = sum(Fe_up, na.rm = TRUE),
    C_total  = sum(C_up, na.rm = TRUE),
    FeC_total = sum(FeC, na.rm = TRUE),
    .groups = "drop"
  )

FeC_summary_nosize <- FeC_bottle_totals %>%
  group_by(treatment, light, inhibitor) %>%
  summarise(
    Fe_mean = mean(Fe_total, na.rm = TRUE),
    Fe_se   = sd(Fe_total, na.rm = TRUE) / sqrt(n()),
    C_mean  = mean(C_total, na.rm = TRUE),
    C_se    = sd(C_total, na.rm = TRUE) / sqrt(n()),
    FeC_mean = mean(FeC_total, na.rm = TRUE),
    FeC_se   = sd(FeC_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


generate_plots_no_facets <- function(data, variable_mean, variable_se, variable_name, 
                                     show_xlab = TRUE, y_min = NULL, y_max = NULL, 
                                     log_scale = FALSE) {
  # Create an x position for dodged points: 
  # add a small offset based on inhibitor (you can adjust the offset as needed)
  data$x_pos <- as.numeric(as.factor(data$treatment)) + ifelse(data$inhibitor == "Control", -0.3, 0.3)
  # Create a new x variable for the line that ignores inhibitor differences:
  data$x_line <- as.numeric(as.factor(data$treatment))
  
  p <- ggplot(data, aes(x = x_pos, y = .data[[variable_mean]], color = light, shape = inhibitor)) +
    # Points and error bars use the dodged x positions:
    geom_point(size = 3, position = position_dodge(width = 0)) +
    geom_errorbar(aes(ymin = .data[[variable_mean]] - .data[[variable_se]],
                      ymax = .data[[variable_mean]] + .data[[variable_se]]), 
                  width = 0.1, position = position_dodge(width = 0)) +
    # Lines use the non-dodged x value to join the points across inhibitor groups (grouped by treatment and light)
    geom_line(aes(x = x_line, group = interaction(treatment, light)), 
              linewidth = 0.5, show.legend = FALSE) +
    scale_x_continuous(breaks = unique(data$x_line), labels = unique(data$treatment)) +
    scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
    scale_shape_manual(values = c("Control" = 1, "+AZ" = 19)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
         y = variable_name, color = NULL, shape = NULL) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 9),
          legend.direction = "horizontal",
          legend.justification = c(0, 1),
          legend.title = element_blank()) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
           shape = guide_legend(nrow = 1, byrow = TRUE))
  
  # If log_scale is TRUE, add a log10 scale with log ticks:
  if(log_scale) {
    p <- p + scale_y_log10(labels = scales::label_number_auto()) +
      annotation_logticks(sides = "l")
  }
  
  return(p)
}


# For Fe uptake
Fe_plot <- generate_plots_no_facets(FeC_summary_nosize, 
                                    variable_mean = "Fe_mean", 
                                    variable_se = "Fe_se", 
                                    variable_name = expression(paste("Fe uptake (pmol L"^-1, " day"^-1, ")")), 
                                    y_min = 0.4, y_max = 400, 
                                    log_scale = TRUE)

# For C uptake
C_plot <- generate_plots_no_facets(FeC_summary_nosize, 
                                   variable_mean = "C_mean", 
                                   variable_se = "C_se", 
                                   variable_name = expression(paste("C uptake (nmol L"^-1, " day"^-1, ")")), 
                                   y_min = 0, y_max = 12, 
                                    log_scale = FALSE)

# For the Fe:C ratio
FeC_plot <- generate_plots_no_facets(FeC_summary_nosize, 
                                     variable_mean = "FeC_mean", 
                                     variable_se = "FeC_se", 
                                     variable_name = "Fe:C ratio (μmol:mol)", 
                                     y_min = 0, y_max = 200, 
                                      log_scale = FALSE)

community_AZ_plot <- (guide_area() / (Fe_plot | C_plot) / (FeC_plot | FeC_plot)) + 
  plot_annotation(tag_levels = list(c('a','b','c','d'))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) &
  theme(panel.spacing.x = unit(8, "pt"), plot.margin = margin(t=5.5, b=5.5, r=5.5, l=2),
        legend.position = 'top', legend.text = element_text(size = 9), 
        legend.title = element_blank(), legend.direction = "horizontal",
        plot.tag.position = c(0.02, 0.96))
community_AZ_plot

# save plot
ggsave("FeC_AZ_community_plot.svg", plot = community_AZ_plot, width = 7.2, height = 7.6)

Fe_model  <- glm((Fe_total) ~ light * treatment * inhibitor, data = FeC_bottle_totals, family = Gamma(link = "log"))
C_model <- glm((C_total) ~ light * treatment * inhibitor, data = FeC_bottle_totals, family = Gamma(link = "log"))
FeC_model  <- glm(log(FeC_total) ~ light * treatment * inhibitor, data = FeC_bottle_totals, family = Gamma(link = "log"))


# emmeans for Fe uptake
emm_total_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ inhibitor | treatment : light ), type = "unlink")
emm_total_C <- summary(emmeans(C_model, specs = pairwise ~ inhibitor | treatment : light), type = "unlink")
emm_total_FeC <- summary(emmeans(FeC_model, specs = pairwise ~ inhibitor | treatment : light), type = "unlink")

# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_total_Fe$emmeans, IinTS = emm_total_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_total_C$emmeans, emm_total_C$contrasts)
emmeans_FeC <- bind_rows(emmeans = emm_total_FeC$emmeans, emm_total_FeC$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C,
  emmeans_FeC = emmeans_FeC
)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_+AZ_sizeaggregated.xlsx")
