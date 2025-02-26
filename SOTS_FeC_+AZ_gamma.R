setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment") #setwd
#loading packages
library(tidyverse)
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
FeC <- read_csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/SOTS_FeC.csv") %>%
  mutate(inhibitor = as.factor(inhibitor) %>% fct_recode("+AZ" = "AZ"),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         size = as.factor(size) %>% fct_recode("0.2-2 μm" = "0.2", "2-20 μm" = "2", ">20 μm" = "20"),
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

# ---------------Plots of the data with inhibitor ----------------------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_summary <- FeC %>%
  group_by(treatment, size, inhibitor) %>%
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
  data$x_pos <- as.numeric(as.factor(data$treatment)) + ifelse(data$inhibitor == "Control", -0.1, 0.1)
  
  p <- ggplot(data, aes(x = x_pos, y = data[[variable_mean]], shape = inhibitor)) +
    facet_grid(~ size) +
    geom_point(size = 3) + 
    geom_line(aes(group = interaction(size, treatment)), linewidth = 0.5, show.legend = FALSE) +
    geom_errorbar(aes(ymin = data[[variable_mean]] - data[[variable_se]], ymax = data[[variable_mean]] + data[[variable_se]]), width = 0.2) +
    scale_x_continuous(breaks = 1:length(unique(data$treatment)), labels = c("+DFB" = "+DFB", "Control" = "Control", "+Fe" = "+Fe"), limits = c(0.5, length(unique(data$treatment)) + 0.5)) +
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
  x_pos = c(2.3, 1.3, 3.3, 2.3), y_pos = c(25.56, 6.46, 40.63, 17.69), 
  label = c("\u2193", "\u2191", "\u2191", "\u2193"), size = c("0.2-2 μm", "2-20 μm", "2-20 μm", ">20 μm"))

asterisks_C <- data.frame(
  x_pos = c(1.3, 2.3, 3.3), y_pos = c(1.054, 0.165, 2.295), 
  label = c("\u2193", "\u2193", "\u2193"), size = c("2-20 μm", "2-20 μm", "2-20 μm"))

asterisks_FeC <- data.frame(
  x_pos = c(2.3, 2.3), y_pos = c(340.7, 7.26), 
  label = c("\u2191", "\u2191"), size = c("2-20 μm", ">20 μm"))

# Reorder the levels of the 'size' variable
asterisks_Fe$size <- factor(asterisks_Fe$size, levels = c("0.2-2 μm", "2-20 μm", ">20 μm"))
asterisks_C$size <- factor(asterisks_C$size, levels = c("0.2-2 μm", "2-20 μm", ">20 μm"))
asterisks_FeC$size <- factor(asterisks_FeC$size, levels = c("0.2-2 μm", "2-20 μm", ">20 μm"))
#-----------------------log scale plots-----------------------------------------------------
# Define label formatting function
format_y_labels <- function(x) {
  ifelse(x > 1, round(x), x)
}
# plots
Fe_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "Fe_mean", "Fe_se", "", y_min = 0.3, y_max = 130) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "top") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_Fe, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  annotation_logticks(sides = "l") +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " d"^-1, ")")))
C_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "C_mean", "C_se", "", y_min = 0.03, y_max = 13) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_C, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  annotation_logticks(sides = "l") +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " d"^-1, ")")))
FeC_plot_errorbars_log <- generate_plots_with_errorbars(FeC_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0.3, y_max = 1300) +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  scale_y_log10(labels = format_y_labels) +
  geom_text(data = asterisks_FeC, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  annotation_logticks(sides = "l")

# Arrange plots in a grid
grid_plot_errorbars_log <- (guide_area() / Fe_plot_errorbars_log / C_plot_errorbars_log / FeC_plot_errorbars_log) + 
  plot_annotation(tag_levels = list(c('A','B','C'))) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) &
  theme(panel.spacing.x = unit(8, "pt"), plot.margin = margin(t=5.5, b=5.5, r=5.5, l=2),
        legend.position = 'top', legend.text = element_text(size = 9), 
        legend.title = element_blank(), legend.direction = "horizontal",
        plot.tag.position = c(0.02, 0.96))

# save plot
ggsave("FeC_AZ_logscale_plot.svg", plot = grid_plot_errorbars_log, width = 7.2, height = 8)

#FITTING MODELS :  -----------------------------------------------------------------------------------
Fe_model  <- glm((Fe_up) ~ treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))
C_model <- glm((C_up) ~ treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))
FeC_model  <- glm(log(FeC) ~ treatment * size * inhibitor, data = FeC, family = Gamma(link = "log"))

Fe_model_reduced_2 <- glm(Fe_up ~ treatment + size + inhibitor +
                            treatment:size + inhibitor:size + inhibitor:treatment,
                          data = FeC, family = Gamma(link = "log"))
C_model_reduced_2 <- glm(C_up ~ treatment + size + inhibitor +
                            treatment:size + inhibitor:size + inhibitor:treatment,
                          data = FeC, family = Gamma(link = "log"))
FeC_model_reduced_2 <- glm(FeC ~ treatment + size + inhibitor +
                            treatment:size + inhibitor:size + inhibitor:treatment,
                          data = FeC, family = Gamma(link = "log"))

C_model_noTI <- glm(C_up ~ treatment + size + inhibitor +
                           treatment:size + inhibitor:size,
                         data = FeC, family = Gamma(link = "log"))

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
Fe_stan_full <- stan_glm(Fe_up ~ treatment * size * inhibitor, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
Fe_stan_reduced_2 <- stan_glm(Fe_up ~ treatment + size + inhibitor +
                                treatment:size + inhibitor:size + inhibitor:treatment,
                              data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_full <- stan_glm(C_up ~ treatment * size * inhibitor, 
                        data = FeC, family = Gamma(link = "log"), iter = 50000)
C_stan_reduced_2 <- stan_glm(C_up ~ treatment + size + inhibitor +
                               treatment:size + inhibitor:size + inhibitor:treatment,
                             data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_full <- stan_glm(FeC ~ treatment * size * inhibitor, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_reduced_2 <- stan_glm(FeC ~ treatment + size + inhibitor +
                                 treatment:size + inhibitor:size + inhibitor:treatment,
                               data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without treatment inhibitor interaction 
Fe_stan_noTI <- stan_glm(Fe_up ~ treatment * size + inhibitor*size, 
                           data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noTI <- stan_glm(C_up ~ treatment * size + inhibitor*size, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noTI <- stan_glm(FeC ~ treatment * size + inhibitor*size, 
                            data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without treatment inhibitor interaction 
Fe_stan_noST <- stan_glm(Fe_up ~ inhibitor*size + inhibitor*size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noST <- stan_glm(C_up ~ inhibitor*size + inhibitor*size, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noST <- stan_glm(FeC ~ inhibitor*size + inhibitor*size, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without size inhibitor interaction
Fe_stan_noSI <- stan_glm(Fe_up ~ treatment*size + inhibitor*treatment, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noSI <- stan_glm(C_up ~ treatment*size + inhibitor*treatment, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noSI <- stan_glm(FeC ~ treatment*size + inhibitor*treatment, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without inhibitor
Fe_stan_noI <- stan_glm(Fe_up ~ treatment*size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noI <- stan_glm(C_up ~ treatment*size, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noI <- stan_glm(FeC ~ treatment*size, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without size
Fe_stan_noS <- stan_glm(Fe_up ~ inhibitor*treatment, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noS <- stan_glm(C_up ~ inhibitor*treatment, 
                       data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noS <- stan_glm(FeC ~ inhibitor*treatment, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)
# models without treatment
Fe_stan_noT <- stan_glm(Fe_up ~ inhibitor * size, 
                        data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_noT <- stan_glm(C_up ~ inhibitor * size, 
                       data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_noT <- stan_glm(FeC ~ inhibitor * size, 
                         data = FeC, family = Gamma(link = "log"), iter = 5000)

# models with inhibitor as Main effect only
Fe_stan_mainI <- stan_glm(Fe_up ~ treatment * size + inhibitor, 
                           data = FeC, family = Gamma(link = "log"), iter = 5000)
C_stan_mainI <- stan_glm(C_up ~ treatment * size + inhibitor, 
                          data = FeC, family = Gamma(link = "log"), iter = 5000)
FeC_stan_mainI <- stan_glm(FeC ~ treatment * size + inhibitor, 
                            data = FeC, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for each model
Fe_loo_full <- loo(Fe_stan_full, k_threshold = 0.7)
Fe_loo_reduced_2 <- loo(Fe_stan_reduced_2, k_threshold = 0.7)
Fe_loo_noSI <- loo(Fe_stan_noSI, k_threshold = 0.7)
Fe_loo_noST <- loo(Fe_stan_noST, k_threshold = 0.7)
Fe_loo_noTI <- loo(Fe_stan_noTI, k_threshold = 0.7)
Fe_loo_noI <- loo(Fe_stan_noI, k_threshold = 0.7)
Fe_loo_noT <- loo(Fe_stan_noT, k_threshold = 0.7)
Fe_loo_noS <- loo(Fe_stan_noS, k_threshold = 0.7)
Fe_loo_mainI <- loo(Fe_stan_mainI, k_threshold = 0.7)

C_loo_full <- loo(C_stan_full, k_threshold = 0.7)
C_loo_reduced_2 <- loo(C_stan_reduced_2, k_threshold = 0.7)
C_loo_noSI <- loo(C_stan_noSI, k_threshold = 0.7)
C_loo_noST <- loo(C_stan_noST, k_threshold = 0.7)
C_loo_noTI <- loo(C_stan_noTI, k_threshold = 0.7)
C_loo_noI <- loo(C_stan_noI, k_threshold = 0.7)
C_loo_noT <- loo(C_stan_noT, k_threshold = 0.7)
C_loo_noS <- loo(C_stan_noS, k_threshold = 0.7)
C_loo_mainI <- loo(C_stan_mainI, k_threshold = 0.7)

FeC_loo_full <- loo(FeC_stan_full, k_threshold = 0.7)
FeC_loo_reduced_2 <- loo(FeC_stan_reduced_2, k_threshold = 0.7)
FeC_loo_noSI <- loo(FeC_stan_noSI, k_threshold = 0.7)
FeC_loo_noST <- loo(FeC_stan_noST, k_threshold = 0.7)
FeC_loo_noTI <- loo(FeC_stan_noTI, k_threshold = 0.7)
FeC_loo_noI <- loo(FeC_stan_noI, k_threshold = 0.7)
FeC_loo_noT <- loo(FeC_stan_noT, k_threshold = 0.7)
FeC_loo_noS <- loo(FeC_stan_noS, k_threshold = 0.7)
FeC_loo_mainI <- loo(FeC_stan_mainI, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(Fe_loo_full, Fe_loo_reduced_2, Fe_loo_noSI, Fe_loo_noST, Fe_loo_noTI, Fe_loo_noI, Fe_loo_noT, Fe_loo_noS, Fe_loo_mainI)
loo_compare(C_loo_full, C_loo_reduced_2, C_loo_noSI, C_loo_noST, C_loo_noTI, C_loo_noI, C_loo_noT, C_loo_noS, C_loo_mainI)
loo_compare(FeC_loo_full, FeC_loo_reduced_2, FeC_loo_noSI, FeC_loo_noST, FeC_loo_noTI, FeC_loo_noI, FeC_loo_noT, FeC_loo_noS, FeC_loo_mainI)

# -----Create a list of models with their corresponding names--------------------------
models_list <- list(
  Fe_full = Fe_model,
  C_full = C_model,
  FeC_full = FeC_model,
  Fe_2 = Fe_model_reduced_2,
  C_2 = C_model_reduced_2,
  FeC_2 = FeC_model_reduced_2,
  C_noTI = C_model_noTI)

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
emm_all_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ inhibitor : treatment : size), type = "unlink")
emm_IinTS_Fe <- summary(emmeans(Fe_model, specs = pairwise ~ inhibitor | treatment : size), type = "unlink")

# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model, specs = pairwise ~ inhibitor : treatment : size), type = "unlink")
emm_IinTS_C <- summary(emmeans(C_model, specs = pairwise ~ inhibitor | treatment : size), type = "unlink")

# emmeans for Fe:C ratio
emm_all_FeC <- summary(emmeans(FeC_model, specs = pairwise ~ inhibitor : treatment : size), type = "unlink")
emm_IinTS_FeC <- summary(emmeans(FeC_model, specs = pairwise ~ inhibitor | treatment : size), type = "unlink")


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
write_xlsx(combined_FeC_emm, "FeC_emmeans_+AZ_gama.xlsx")



#--- OKAY also just doing the whole thing without size fractions --------------------------------
FeC_nosize <- FeC %>%
  group_by(treatment, inhibitor, bottle) %>%
  summarize(
    Fe_up = sum(Fe_up, na.rm = TRUE),
    C_up = sum(C_up, na.rm = TRUE)
  ) %>%
  mutate(
    FeC = Fe_up / C_up
  )

#-- models --------------
Fe_model_nosize  <- glm((Fe_up) ~ treatment * inhibitor, data = FeC_nosize, family = Gamma(link = "log"))
C_model_nosize <- glm((C_up) ~ treatment * inhibitor, data = FeC_nosize, family = Gamma(link = "log"))
FeC_model_nosize  <- glm((FeC) ~ treatment * inhibitor, data = FeC_nosize, family = Gamma(link = "log"))

models_list2 <- list(
  Fe_full = Fe_model_nosize,
  C_full = C_model_nosize,
  FeC_full = FeC_model_nosize)

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
shape_rate_list <- lapply(models_list2, calculate_shape_rate)

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
std_residuals_list <- lapply(models_list2, calculate_standardized_residuals)

# Create a list of QQ plots for each model in models_list, incorporating shape and rate
qq_plot_list <- mapply(function(x, y) {
  create_qq_plot(std_residuals_list[[x]], x, y$shape, y$rate)
}, x = names(std_residuals_list), y = shape_rate_list, SIMPLIFY = FALSE)

# Arrange the plots in a grid with 2 rows
plot_grid(
  plotlist = qq_plot_list, align = "hv",
  ncol = 3, nrow = 1
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
res_fit_plot_list <- lapply(names(models_list2), function(x) {
  res_fit_plot(models_list2[[x]], x, "Residuals")
})

# Arrange the plots in a 3x3 grid
plot_grid(
  plotlist = res_fit_plot_list, align = "hv",
  ncol = 3, nrow = 1
)

# --- interrogating the model results using EMMEANS ---------------------------------------------------
# emmeans for Fe uptake
emm_all_Fe <- summary(emmeans(Fe_model_nosize, specs = pairwise ~ inhibitor : treatment), type = "unlink")
emm_IinT_Fe <- summary(emmeans(Fe_model_nosize, specs = pairwise ~ inhibitor | treatment), type = "unlink")

# emmeans for C uptake
emm_all_C <- summary(emmeans(C_model_nosize, specs = pairwise ~ inhibitor : treatment), type = "unlink")
emm_IinT_C <- summary(emmeans(C_model_nosize, specs = pairwise ~ inhibitor | treatment), type = "unlink")

# emmeans for Fe:C ratio
emm_all_FeC <- summary(emmeans(FeC_model_nosize, specs = pairwise ~ inhibitor : treatment), type = "unlink")
emm_IinT_FeC <- summary(emmeans(FeC_model_nosize, specs = pairwise ~ inhibitor | treatment), type = "unlink")


# Combine emmeans and contrasts for each data frame
emmeans_Fe <- bind_rows(emmeans = emm_all_Fe$emmeans, IinT = emm_IinT_Fe$contrasts)
emmeans_C <- bind_rows(emmeans = emm_all_C$emmeans, IinT = emm_IinT_C$contrasts)
emmeans_FeC <- bind_rows(emmeans = emm_all_FeC$emmeans, IinT = emm_IinT_FeC$contrasts)

# Create a list of combined data frames
combined_FeC_emm <- list(
  emmeans_Fe = emmeans_Fe,
  emmeans_C = emmeans_C,
  emmeans_FeC = emmeans_FeC
)

# Export the combined data frames to an Excel file
write_xlsx(combined_FeC_emm, "FeC_emmeans_+AZ_gama_nosize.xlsx")

#----making that plot --------------------------------------------
# Calculate the mean and standard error for each unique combination
FeC_nosize_summary <- FeC_nosize %>%
  group_by(treatment, inhibitor) %>%
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
  data$x_pos <- as.numeric(as.factor(data$treatment)) + ifelse(data$inhibitor == "Control", -0.1, 0.1)
  
  p <- ggplot(data, aes(x = x_pos, y = data[[variable_mean]], shape = inhibitor)) +
    geom_point(size = 3) + 
    geom_line(aes(group = interaction(treatment)), linewidth = 0.5, show.legend = FALSE) +
    geom_errorbar(aes(ymin = data[[variable_mean]] - data[[variable_se]], ymax = data[[variable_mean]] + data[[variable_se]]), width = 0.1) +
    scale_x_continuous(breaks = 1:length(unique(data$treatment)), labels = c("+DFB" = "+DFB", "Control" = "Control", "+Fe" = "+Fe"), limits = c(0.5, length(unique(data$treatment)) + 0.5)) +
    scale_shape_manual(values = c("Control" = 1, "+AZ" = 19)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
         y = paste(variable_name),
         color = NULL, shape = NULL) + 
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"), 
          legend.text = element_text(size = 9),
          legend.direction = "vertical",
          legend.justification = c(0, 1)) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    ) 
}

#---adding asterisks for significance ----------------------------------------------------------------
asterisks_Fe <- data.frame(
  x_pos = c(2.3), y_pos = c(54.63), 
  label = c("\u2191"))

asterisks_FeC <- data.frame(
  x_pos = c(2.3), y_pos = c(17), 
  label = c("\u2191"))

# Define label formatting function
format_y_labels <- function(x) {
  ifelse(x > 1, round(x), x)
}

Fe_plot_errorbars <- generate_plots_with_errorbars(FeC_nosize_summary, "Fe_mean", "Fe_se", "", y_min = 0, y_max = 200) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  geom_text(data = asterisks_Fe, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) +
  labs(y = expression(paste("Fe uptake (pmol L"^-1, " d"^-1, ")")))
C_plot_errorbars <- generate_plots_with_errorbars(FeC_nosize_summary, "C_mean", "C_se", "", y_min = 0, y_max = 12) +
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, 12, by = 3)) +
  labs(y = expression(paste("C uptake (", mu, "mol L"^-1, " d"^-1, ")")))
FeC_plot_errorbars <- generate_plots_with_errorbars(FeC_nosize_summary, "FeC_mean", "FeC_se", "Fe:C ratio (μmol:mol)", y_min = 0, y_max = 20) +
  theme(panel.grid.minor = element_blank(), legend.position = "none") +
  geom_text(data = asterisks_FeC, aes(x = x_pos, y = y_pos, label = label), size = 4, inherit.aes = FALSE) 

# Arrange plots in a grid
grid_plot_errorbars <- (Fe_plot_errorbars / FeC_plot_errorbars) | (C_plot_errorbars / guide_area()) + 
  plot_annotation(tag_levels = list(c('A','B','C'))) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(10,10)) &
  theme(panel.spacing.x = unit(8, "pt"), plot.margin = margin(t=5.5, b=5.5, r=5.5, l=2),
        legend.position = 'bottom', legend.text = element_text(size = 9), 
        legend.direction = "vertical",
        plot.tag.position = c(0.02, 0.96))

# save plot
ggsave("FeC_AZ_nosize_plot.svg", plot = grid_plot_errorbars, width = 7.2, height = 6.4)
