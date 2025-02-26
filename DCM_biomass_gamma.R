setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment") #setwd
#----------------loading packages--------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(patchwork)
library(emmeans)
library(visreg)
library(writexl)
library(rstan)
library(rstanarm)

#------load the POC PON Chla data  into R---------------------------------------------------------------------------------------
chems <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment/DCM_chems.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(light = as.factor(light), 
         inhibitor = as.factor(inhibitor),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         bottle = as.factor(bottle))
# Replace N/A with NA
chems[chems == "N/A"]  <- NA
# make values numeric
chems <- chems %>% mutate_at(vars(POC:SiC_ratio), as.numeric)
# Reorder  levels
chems$treatment <- factor(chems$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)

# grab data from initial conditions
time_0 <- chems %>% filter(time == "0") %>% mutate_at(vars(POC:Chla), as.numeric)
chems <- chems %>% filter(time == 10)

#Check the structure
str(chems)

#----checking density distributions.can edit the exressions to find the best transformations ------------------------------
p1 <- ggplot(chems, aes(x = log(dPOC))) + geom_density(alpha = .2, fill = "#FF6FF6")
p2 <- ggplot(chems, aes(x = log(dPON))) + geom_density(alpha = .2, fill = "#FF6FF6")
p3 <- ggplot(chems, aes(x = log(dChla))) + geom_density(alpha = .2, fill = "#FF6FF6")
p4 <- ggplot(chems, aes(x = (abs(d13C)))) + geom_density(alpha = .2, fill = "#FF6FF6")
p5 <- ggplot(chems, aes(x = sqrt(dCN))) + geom_density(alpha = .2, fill = "#FF6FF6")
p6 <- ggplot(chems, aes(x = (ChlC_ratio))) + geom_density(alpha = .2, fill = "#FF6FF6")
p7 <- ggplot(chems, aes(x = log(dSiO3))) + geom_density(alpha = .2, fill = "#FF6FF6")
p8 <- ggplot(chems, aes(x = (SiP_ratio))) + geom_density(alpha = .2, fill = "#FF6FF6")
p9 <- ggplot(chems, aes(x = (SiC_ratio))) + geom_density(alpha = .2, fill = "#FF6FF6")
p10 <- ggplot(chems, aes(x = log(dNOx))) + geom_density(alpha = .2, fill = "#FF6FF6")

# Combine plots using patchwork syntax
combined_plot <- p1 + p2 + p3 + p4 + p5 + 
  p6 + p7 + p8 + p9 + p10 +
  plot_layout(ncol = 2)  # arrange plots in 2 columns

# Display combined plot
print(combined_plot)

# ---Compute the mean and variance for each group defined by light, treatment, and size-----------------------
# select columns that are required
vars <- c("dPOC", "dPON", "dChla", "d13C", "CN_ratio", "ChlC_ratio", "dSiO3", "SiP_ratio", 
          "SiC_ratio", "dNOx", "dPO4", "dNH4", "NoxN_ratio", "CP_ratio", "NP_ratio", "NOx", "PO4", "SiO3", "NH4")

#create summary data frame with variable names as "{.fn}_{.col}"
chems_summary <- chems %>%
  group_by(light, treatment) %>%
  summarise(across(all_of(vars), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        var = ~var(.x, na.rm = TRUE)), 
                   .names = "{.fn}_{.col}"))

# Redefine the scatterplot function to take a single variable
create_scatterplot <- function(var) {
  ggplot(chems_summary, aes(x = .data[[paste0("mean_", var)]], y = .data[[paste0("var_", var)]])) +
    geom_point() +
    theme_bw() +
    xlab("Mean") +
    ylab("Variance") +
    ggtitle(paste0(var, " Var vs. Mean"))
}

# Now use map() to apply this function to each variable in `vars`
variance_plots <- purrr::map(vars, create_scatterplot)
names(variance_plots) <- vars

# display plots POC/PON/Chla
(variance_plots$dPOC | variance_plots$dPON | variance_plots$CN_ratio) / (variance_plots$dChla | variance_plots$ChlC_ratio | variance_plots$d13C)
# display plots major nuts
(variance_plots$dNOx | variance_plots$dSiO3) / (variance_plots$dPO4 | variance_plots$dNH4)
# display plots ratios
(variance_plots$CP_ratio | variance_plots$CN_ratio | variance_plots$NP_ratio) / 
(variance_plots$NoxN_ratio | variance_plots$SiP_ratio | variance_plots$SiC_ratio)

#------FITTING MODELS :  -----------------------------------------------------------------------------------
# Define a list of the response variables for Gaussian models
resp_vars_gauss <- c("dPOC", "dPON", "dCN", "d13C", "dChla", "dChlC", "dSiO3", "SiP_ratio", 
                     "SiC_ratio", "dNOx", "dPO4", "dNH4", "NoxN_ratio", "CP_ratio", "NP_ratio")

# Define a function to fit the Gaussian model
fit_gauss_model <- function(var) {
  formula <- as.formula(paste(var, "~ light * treatment"))
  if (var == "dCN") {
    formula <- as.formula(paste("(dCN) ~ light + treatment"))}
  glm(formula, data = chems, family = gaussian())
}

# Use lapply to fit the Gaussian models for all variables and store the results in a list
gauss_models_list <- setNames(lapply(resp_vars_gauss, fit_gauss_model), resp_vars_gauss)

# Define a list of the response variables for gamma models
resp_vars <- c("dPOC", "dPON", "dCN", "d13C", "dChla", "dChlC", "dSiO3", "SiP_ratio", 
               "SiC_ratio", "dNOx", "dPO4", "dNH4", "NoxN_ratio", "CP_ratio", "NP_ratio")

# Define a function to fit the model
fit_model <- function(var) {
  formula <- as.formula(paste(var, "~ light * treatment"))
  if (var == "d13C") {
    formula <- as.formula(paste("abs(d13C) ~ light * treatment"))}
  if (var == "dCN") {
    formula <- as.formula(paste("(dCN) ~ light + treatment"))}
  glm(formula, data = chems, family = Gamma(link = "log"))
}

# Use lapply to fit the models for all variables and store the results in a list
gamma_models_list <- setNames(lapply(resp_vars, fit_model), resp_vars)

#--------------- Function to generate QQ plots with normal distribution --------------------------------------------------------------------------------------------
create_qq_plot <- function(residuals, model_name) {
  residuals <- unlist(residuals)
  ggplot(data = data.frame(residuals), aes(sample = residuals)) +
    geom_qq() +
    geom_qq_line() +
    theme_bw() +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    ggtitle(paste(model_name, "QQ Plot (Normal)"))
}

# Generate standardized residuals for each model in models_list and models_list_reduced
std_residuals_list <- gauss_models_list %>% 
  map(~rstandard(.x)) %>%
  set_names(c("dPOC", "dPON", "dCN", "d13C", "dChla", "ChlC_ratio", "dSiO3", "SiP_ratio", 
              "SiC_ratio", "dNOx", "dPO4", "dNH4", "NoxN_ratio", "CP_ratio", "NP_ratio"))

# Create a list of QQ plots for each model in models_list and models_list_reduced
qq_plot_list <- map(names(std_residuals_list), 
                    ~create_qq_plot(std_residuals_list[.x], paste0(.x, " Full")))

# Convert the list of plots to a patchwork object
qq_plot_patchwork <- wrap_plots(qq_plot_list, ncol = 3)

# Print the plot
print(qq_plot_patchwork)

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
shape_rate_list <- lapply(gamma_models_list, calculate_shape_rate)

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
std_residuals_list <- lapply(gamma_models_list, calculate_standardized_residuals)

# Create a list of QQ plots for each model in models_list, incorporating shape and rate
qq_plot_list <- mapply(function(x, y) {
  create_qq_plot(std_residuals_list[[x]], x, y$shape, y$rate)
}, x = names(std_residuals_list), y = shape_rate_list, SIMPLIFY = FALSE)

# Convert the list of plots to a patchwork object
qq_plot_patchwork <- wrap_plots(qq_plot_list, ncol = 3)

# Print the plot
qq_plot_patchwork

# ----Function to create residual vs fitted plot for a model------------------------------------------------
res_fit_plot <- function(model, model_name) {
  res_fit <- broom::augment(model)
  
  plot <- ggplot(data = res_fit, aes(x = .fitted, y = .resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(paste(model_name, "Residuals vs Fitted")) +
    theme_bw()
  
  return(plot)
}

# Create a list of Residuals vs Fitted plots for Gaussian models
res_fit_plot_list_gauss <- lapply(names(gauss_models_list), function(x) {
  res_fit_plot(gauss_models_list[[x]], paste0(x, " gaussian"))
})

# Combine the plots for Gaussian models
combined_res_fit_plots_gauss <- wrap_plots(res_fit_plot_list_gauss, ncol = 3)

# Create a list of Residuals vs Fitted plots for Gamma models
res_fit_plot_list_gamma <- lapply(names(gamma_models_list), function(x) {
  res_fit_plot(gamma_models_list[[x]], paste0(x, " gamma"))
})

# Combine the plots for Gamma models
combined_res_fit_plots_gamma <- wrap_plots(res_fit_plot_list_gamma, ncol = 3)

# Print the plots
print(combined_res_fit_plots_gauss)
print(combined_res_fit_plots_gamma)

#---------- Investigate group comparisons for POC PON Chla using emmeans-------------------------------------------------------------------------------------------
# Create a list where model names are matched with the corresponding model list
model_info <- list(
  "dPOC" = list(model_list = gamma_models_list, model_name = "dPOC"),
  "dPON" = list(model_list = gamma_models_list, model_name = "dPON"),
  "dCN" = list(model_list = gauss_models_list, model_name = "dCN"),
  "dChla" = list(model_list = gamma_models_list, model_name = "dChla"),
  "dChlC" = list(model_list = gauss_models_list, model_name = "dChlC"),
  "d13C" = list(model_list = gauss_models_list, model_name = "d13C")
)

# Calculate emmeans for each model
emm_nuts <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ light:treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

emm_nuts_light <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ light|treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

emm_nuts_trt <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ treatment|light)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

# Combine emmeans and contrasts for each data frame
emmeans_dPOC <- bind_rows(emmeans = emm_nuts$dPOC$emmeans, light = emm_nuts_light$dPOC$contrasts, trt = emm_nuts_trt$dPOC$contrasts)
emmeans_dPON <- bind_rows(emmeans = emm_nuts$dPON$emmeans, light = emm_nuts_light$dPON$contrasts, trt = emm_nuts_trt$dPON$contrasts)
emmeans_dCN <- bind_rows(emmeans = emm_nuts$dCN$emmeans, light = emm_nuts_light$dCN$contrasts, trt = emm_nuts_trt$dCN$contrasts)
emmeans_dChla <- bind_rows(emmeans = emm_nuts$dChla$emmeans, light = emm_nuts_light$dChla$contrasts, trt = emm_nuts_trt$dChla$contrasts)
emmeans_dChlC <- bind_rows(emmeans = emm_nuts$dChlC$emmeans, light = emm_nuts_light$dChlC$contrasts, trt = emm_nuts_trt$dChlC$contrasts)
emmeans_d13C <- bind_rows(emmeans = emm_nuts$d13C$emmeans, light = emm_nuts_light$d13C$contrasts, trt = emm_nuts_trt$d13C$contrasts)

# Create a list of combined data frames
combined_CNChl_emm <- list(
  emmeans_dPOC = emmeans_dPOC,
  emmeans_dPON = emmeans_dPON,
  emmeans_dCN = emmeans_dCN,
  emmeans_dChla = emmeans_dChla,
  emmeans_dChlC = emmeans_dChlC,
  emmeans_d13C = emmeans_d13C
)

# Export the combined data frames to an Excel file
write_xlsx(combined_CNChl_emm, "CNChl_emmeans.xlsx")

# ---Calculate the mean and standard error for each unique combination----------------------------------------------------
# select columns that are required
vars <- c("dPOC", "dPON", "dCN", "D13C", "dChla", "dChlC", "dSiO3", "SiP_ratio", 
          "SiC_ratio", "dNOx", "dPO4", "dNH4", "NoxN_ratio", "CP_ratio", "NP_ratio", "NOx", "PO4", "SiO3", "NH4")
#create summary
chems_summary <- chems %>%
  group_by(light, treatment) %>%
  summarise(across(all_of(vars), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))), 
                   .names = "{.col}_{.fn}"))


# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, 
                                          show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, 
                                          legend_position = "top", hline_y) {
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), color = light)) +
    geom_hline(yintercept = hline_y, linetype = "dashed", color = "#6F286AFF") +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
    scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
    scale_x_discrete(breaks = unique(data$treatment)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
         y = y_lab) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"), 
          panel.grid.minor = element_blank(), legend.title = element_blank(),
          legend.direction = "horizontal", legend.position = legend_position) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)) 
  return(p)
}

#---adding asterisks for significance PONPOC ----------------------------------------------------------------
asterisks_dPOC <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.5, 0.858, 0.38, 0.738), y_pos = c(5.96, 5.61, 6.92, 3.2, 3.37, 5.61, 6.92), 
  label = c("a", "b", "c", "a", "a", "\u263C", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "black", "black"))

asterisks_dPON <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.5, 0.858, 0.38, 0.738), y_pos = c(3.79, 3.59, 4.98, 2.47, 2.66, 3.59, 4.98), 
  label = c("a", "b", "c", "a", "a", "\u263C", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "black", "black"))

asterisks_dCN <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.5, 0.858), y_pos = c(7.6, 7.13, 6.92, 5.2, 4.97), 
  label = c("a", "a", "a", "a", "a"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD"))

asterisks_dChla <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.022, 0.38, 0.738), y_pos = c(0.3, 0.2, 2.99, -1.2, -1.08, 0.05, 0.3, 0.2, 2.99), 
  label = c("a", "a", "b", "a", "a", "b","\u263C", "\u263C","\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black", "black", "black"))

asterisks_dChla2 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(0.3, 0.2, 2.99, -1.2, -1.08, 0.05, 2.99), 
  label = c("a", "a", "b", "a", "a", "b","all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_dChlC <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.5, 0.858, 0.738), y_pos = c(7.54, 8.55, 27.2, 2.55, 11, 27.2), 
  label = c("a", "a", "b", "a", "b", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "black"))

# -------------------------------Create plots for POC PON Chla -----------------------------------------------------------------------
# plots for the data...
plot_dPOC <- generate_plots_with_errorbars(chems_summary, "dPOC_mean", "dPOC_se", "dPOC", show_xlab = FALSE, y_min = 0, y_max = 600,
                                           y_lab = expression(paste(Delta,"POC (",mu,"g L"^-1,")")), legend_position = "none", hline_y = -100) + 
                                           scale_y_continuous(breaks = seq(0, 600, by = 150)) 

plot_dPON <- generate_plots_with_errorbars(chems_summary, "dPON_mean", "dPON_se", "dPON", show_xlab = FALSE, y_min = 0, y_max = 120, 
                                           y_lab = expression(paste(Delta,"PON (",mu,"g L"^-1,")")), legend_position = "none", hline_y = -10) + 
                                           scale_y_continuous(breaks = seq(0, 120, by = 30)) 

plot_dCN <- generate_plots_with_errorbars(chems_summary, "dCN_mean", "dCN_se", "dCN", show_xlab = TRUE, y_min = 4, y_max = 8,
                                           y_lab = expression(paste("C:N ratio")), legend_position = "none", hline_y = -10)

plot_Chla <- generate_plots_with_errorbars(chems_summary, "dChla_mean", "dChla_se", "dChla", show_xlab = FALSE, y_min = 0, y_max = 3,
                                          y_lab = expression(paste(Delta,"Chl-",italic("a"),"s",mu,"g L"^-1,"s")), legend_position = "none", hline_y = -10) + 
                                          scale_y_continuous(breaks = seq(0, 3, by = 0.5)) 

plot_ChlC_ratio <- generate_plots_with_errorbars(chems_summary, "dChlC_mean", "dChlC_se", "dChlC", show_xlab = TRUE, y_min = 0, y_max = 30,
                                          y_lab = expression(paste(Delta,"Chl:C (mg g"^-1,")")), legend_position = "none", hline_y = -10) 

# visreg plots...
visreg_POC <- visreg(gamma_models_list$dPOC, "treatment", by="light", overlay = TRUE,
                      gg = TRUE, points.par = list(shape = 19), scale = "response",  top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 8, by = 2),
                     sec.axis = sec_axis(~., name = expression(paste("ln(",Delta, "POC)")))) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dPOC, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dPOC$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_PON <- visreg(gamma_models_list$dPON, "treatment", by="light", overlay = TRUE, 
                     gg = TRUE, points = list(shape = 19), scale = "response",  top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(1, 120), breaks = seq(1, 5, by = 1),
                     sec.axis = sec_axis(~., name = expression(paste("ln(",Delta, "PON)")))) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dPON, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dPON$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_CN_ratio <- visreg(gauss_models_list$dCN, "treatment", by="light", overlay = TRUE, 
                     gg = TRUE, points = list(shape = 19), scale = "linear",  top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(4, 8), breaks = seq(4, 8, by = 1)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dCN, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dCN$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_Chla <- visreg(gamma_models_list$dChla, "treatment", by="light", overlay = TRUE, 
       gg = TRUE, points = list(shape = 19), scale = "response",  top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(0, 3), breaks = seq(-3, 3, by = 2), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(",Delta,"Chl-", italic("a"),")")), breaks = seq(-3, 3, by = 2))) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.direction = "horizontal", 
        legend.position = "top") +
  geom_text(data = asterisks_dChla2, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dChla2$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_ChlC <- visreg(gauss_models_list$dChlC, "treatment", by="light", overlay = TRUE, 
                      gg = TRUE, points = list(shape = 19),  top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 10)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dChlC, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dChlC$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# combined plots and printing
POCPON_plot <- guide_area() / (plot_dPOC | visreg_POC) / (plot_dPON | visreg_PON) / (plot_Chla | visreg_Chla) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

Chl_plot <- guide_area() / (plot_Chla | visreg_Chla) / (plot_ChlC_ratio | visreg_ChlC)  +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("POCPON_model_plot.svg", POCPON_plot, width = 7.2, height = 7.4)
ggsave("Chl_model_plot.svg", Chl_plot, width = 7.2, height = 7.4)

# -------------------------------Create plots for 13C -----------------------------------------------------------------------
# plot for data... PDB
plot_d13C <- generate_plots_with_errorbars(chems_summary, "D13C_mean", "D13C_se", "D13C", show_xlab = TRUE, y_min = -26.5, y_max = -22.5,
                                           y_lab = expression(paste(delta ^13,"C")), legend_position = "none", hline_y = -22.56) + 
  scale_y_continuous(breaks = seq(-23, -26, by = -1))

asterisks_d13C <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.5, 0.858, 0.738), y_pos = c(-4, -4.2, -0.6, -2.5, -3.4, -0.6), 
  label = c("a", "a", "b", "a", "a", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "black"))

# visreg plot...
visreg_d13C <- visreg(gauss_models_list$d13C, "treatment", by="light", overlay = TRUE, 
                          gg = TRUE, points = list(shape = 19), scale = "response", top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(-4.5, -0.5), breaks = seq(-4, -1, by = 1),
                  sec.axis = sec_axis(~., name = expression(paste(Delta ^13,"C")))) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        #axis.text.y = element_blank(), 
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_d13C, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_d13C$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# combine plots and printing
d13C_plot <- guide_area() / (plot_d13C | visreg_d13C)  +
  plot_layout(guides = 'collect', nrow = (2), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("13C_model_plot.svg", d13C_plot, width = 7.2, height = 3.87)

#---------- Investigate group comparisons for NUTRIENTS using emmeans-------------------------------------------------------------------------------------------
# Create a list where model names are matched with the corresponding model list
model_info <- list(
  "dNOx" = list(model_list = gamma_models_list, model_name = "dNOx"),
  "dNH4" = list(model_list = gamma_models_list, model_name = "dNH4"),
  "dPO4" = list(model_list = gauss_models_list, model_name = "dPO4"),
  "dSiO3" = list(model_list = gauss_models_list, model_name = "dSiO3")
)

# Calculate emmeans for each model
emm_nuts <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ light:treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

emm_nuts_light <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ light|treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

emm_nuts_trt <- map(model_info, ~emmeans(.x$model_list[[.x$model_name]], specs = pairwise ~ treatment|light)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(names(model_info))

# Combine emmeans and contrasts for each data frame
emmeans_dNOx <- bind_rows(emmeans = emm_nuts$dNOx$emmeans, light = emm_nuts_light$dNOx$contrasts, trt = emm_nuts_trt$dNOx$contrasts)
emmeans_dNH4 <- bind_rows(emmeans = emm_nuts$dNH4$emmeans, light = emm_nuts_light$dNH4$contrasts, trt = emm_nuts_trt$dNH4$contrasts)
emmeans_dPO4 <- bind_rows(emmeans = emm_nuts$dPO4$emmeans, light = emm_nuts_light$dPO4$contrasts, trt = emm_nuts_trt$dPO4$contrasts)
emmeans_dSiO3 <- bind_rows(emmeans = emm_nuts$dSiO3$emmeans, light = emm_nuts_light$dSiO3$contrasts, trt = emm_nuts_trt$dSiO3$contrasts)

# Create a list of combined data frames
combined_nuts_emm <- list(
  emmeans_dNOx = emmeans_dNOx,
  emmeans_dNH4 = emmeans_dNH4,
  emmeans_dPO4 = emmeans_dPO4,
  emmeans_dSiO3 = emmeans_dSiO3
)

# Export the combined data frames to an Excel file
write_xlsx(combined_nuts_emm, "nuts_emmeans.xlsx")

#---adding asterisks for significance nutrients --------------------------------------------------------------------------
asterisks_dNOx <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), 
  y_pos = c(1.14, 1.09, 3.18, -2.78, -2.7, -2.61, 3.6), 
  label = c("a", "a", "b", "a", "a", "a","all:\u263C"), 
  colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_dNH4 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), 
  y_pos = c(0.1, 0.24, 0.24, -2.1, -3.9, -2.24, 0.8), 
  label = c("a", "a", "a", "a", "b", "a","all:\u263C"), 
  colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))


asterisks_dPO4 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), 
  y_pos = c(0.315, 0.337, 0.857, 0.123, 0.107, 0.12, 1.00), 
  label = c("a", "a", "b", "a", "a", "a","all:\u263C"), 
  colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_dSiO3 <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), 
  y_pos = c(5.23, 6.03, 10.43, 2.20, 2.06, 1.90, 11.50), 
  label = c("a", "a", "b", "a", "a", "a","all:\u263C"), 
  colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

# -------------------------------Create all plots for nutrients -----------------------------------------------------------------------

plot_dNOx <- generate_plots_with_errorbars(chems_summary, "NOx_mean", "NOx_se", "NOx_ratio", show_xlab = FALSE, y_min = 0, y_max = 30,
                                           y_lab = expression(paste(Delta,"Nitrate (",mu,"M)")), legend_position = "none", hline_y = -100)

plot_dPO4 <- generate_plots_with_errorbars(chems_summary, "PO4_mean", "PO4_se", "PO4_ratio", show_xlab = FALSE, y_min = -0.05, y_max = 10,
                                           y_lab = expression(paste(Delta,"Phosphate (",mu,"M)")), legend_position = "none", hline_y = -100)

plot_dSiO3 <- generate_plots_with_errorbars(chems_summary, "SiO3_mean", "SiO3_se", "SiO3_ratio", show_xlab = FALSE, y_min = -0.2, y_max = 30,
                                            y_lab = expression(paste(Delta,"Silicic acid (",mu,"M)")), legend_position = "none", hline_y = -100) + 
                                            scale_y_continuous(breaks = seq(0, 12, by = 3))

plot_dNH4 <- generate_plots_with_errorbars(chems_summary, "NH4_mean", "NH4_se", "NH4_ratio", show_xlab = TRUE, y_min = 0, y_max = 1,
                                           y_lab = expression(paste(Delta,"Ammonia (",mu,"M)")), legend_position = "none", hline_y = -100) + 
                                           scale_y_continuous(breaks = seq(0, 1, by = 0.2))

#visreg plots
visreg_dNOx <- visreg(gamma_models_list$dNOx, "treatment", by="light", overlay = TRUE, 
                      gg = TRUE, points = list(shape = 19), scale = "response", top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(0, 12), breaks = seq(-4, 4, by = 2), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(",Delta,"Nitrate)")))) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dNOx, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dNOx$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_dNH4 <- visreg(gamma_models_list$dNH4, "treatment", by="light", overlay = TRUE, 
                       gg = TRUE, points = list(shape = 19), scale = "response", top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(0, 1), breaks = seq(-4, 1, by = 1), 
                     sec.axis = sec_axis(~., name = expression(paste("ln(",Delta,"Ammonia)")))) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dNH4, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dNH4$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_dSiO3 <- visreg(gauss_models_list$dSiO3, "treatment", by="light", overlay = TRUE, 
                       gg = TRUE, points = list(shape = 19), top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(-0.2, 12.2), breaks = seq(0, 12, by = 3)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dSiO3, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dSiO3$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_dPO4 <- visreg(gauss_models_list$dPO4, "treatment", by="light", overlay = TRUE, 
                      gg = TRUE, points = list(shape = 19), top = 'points', partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dPO4, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_dPO4$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# ----------putting plots together------------------------------------------------------
nuts_plot <- guide_area() / (plot_dNOx | visreg_dPO4) / (plot_dSiO3 | visreg_dSiO3)/ 
            (plot_dNOx | visreg_dNOx) / (plot_dNH4 | visreg_dNH4) +
  plot_layout(guides = 'collect', nrow = (5), heights = c(1,10,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D','E','F','G','H'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top", 
        legend.box = "horizontal", 
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("nuts_model_plot.svg", nuts_plot, width = 7.2, height = 9.8) 

#------------ looking at plots of some of the nutrient ratios...----------------------------

plot_SiP <- generate_plots_with_errorbars(chems_summary, "SiP_ratio_mean", "SiP_ratio_se", "SiP_ratio", show_xlab = FALSE, y_min = 0, y_max = 75,
                                           y_lab = "Si:P ratio (mol:mol)", legend_position = "none", hline_y = -100)

plot_SiC <- generate_plots_with_errorbars(chems_summary, "SiC_ratio_mean", "SiC_ratio_se", "SiC_ratio", show_xlab = FALSE, y_min = 0, y_max = 1,
                                           y_lab = "Si:C ratio (mol:mol)", legend_position = "none", hline_y = -100)

plot_CP <- generate_plots_with_errorbars(chems_summary, "CP_ratio_mean", "CP_ratio_se", "CP_ratio", show_xlab = FALSE, y_min = 0, y_max = 400,
                                            y_lab = "C:P ratio (mol:mol)", legend_position = "none", hline_y = -100) 

plot_CN <- generate_plots_with_errorbars(chems_summary, "dCN_mean", "dCN_se", "dCN_ratio", show_xlab = TRUE, y_min = 0, y_max = 10,
                                           y_lab = "C:N ratio (mol:mol)", legend_position = "none", hline_y = -100) 

plot_NP <- generate_plots_with_errorbars(chems_summary, "NP_ratio_mean", "NP_ratio_se", "NP_ratio", show_xlab = TRUE, y_min = 0, y_max = 75,
                                         y_lab = "N:P ratio (mol:mol)", legend_position = "none", hline_y = -100) 

plot_NoxN <- generate_plots_with_errorbars(chems_summary, "NoxN_ratio_mean", "NoxN_ratio_se", "NoxN_ratio", show_xlab = TRUE, y_min = 0, y_max = 1.2,
                                         y_lab = "Nitrate:PON (mol:mol)", legend_position = "none", hline_y = -100) 

nuts_ratios_plot <- guide_area() / (plot_SiP | plot_SiC) / (plot_CP | plot_CN) / (plot_NP | plot_NoxN) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D','E','F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = "none")
  
ggsave("nuts_ratios_plot.svg", nuts_ratios_plot, width = 7.2, height = 7.4) 

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
# Gaussian full models
dPOC_gaussian_full <- stan_glm(dPOC ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dPON_gaussian_full <- stan_glm(dPON ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dCN_gaussian_full <- stan_glm(dCN ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
Chla_gaussian_full <- stan_glm(dChla ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dChlC_gaussian_full <- stan_glm(dChlC ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
d13C_gaussian_full <- stan_glm(abs(d13C) ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dNOx_gaussian_full <- stan_glm(dNOx ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dPO4_gaussian_full <- stan_glm(dPO4 ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dNH4_gaussian_full <- stan_glm(dNH4 ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dSiO3_gaussian_full <- stan_glm(dSiO3 ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
SiP_ratio_gaussian_full <- stan_glm(SiP_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
SiC_ratio_gaussian_full <- stan_glm(SiC_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
CP_ratio_gaussian_full <- stan_glm(CP_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
CN_ratio_gaussian_full <- stan_glm(CN_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
NP_ratio_gaussian_full <- stan_glm(NP_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
NoxN_ratio_gaussian_full <- stan_glm(NoxN_ratio ~ light * treatment, data = chems, family = gaussian(), iter = 5000)

# Gaussian reduced models
dPOC_gaussian_reduced <- stan_glm(dPOC ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dPON_gaussian_reduced <- stan_glm(dPON ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dCN_gaussian_reduced <- stan_glm(dCN ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
Chla_gaussian_reduced <- stan_glm(dChla ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dChlC_gaussian_reduced <- stan_glm(dChlC ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
d13C_gaussian_reduced <- stan_glm(abs(d13C) ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dNOx_gaussian_reduced <- stan_glm(dNOx ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dPO4_gaussian_reduced <- stan_glm(dPO4 ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dNH4_gaussian_reduced <- stan_glm(dNH4 ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dSiO3_gaussian_reduced <- stan_glm(dSiO3 ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
SiP_ratio_gaussian_reduced <- stan_glm(SiP_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
SiC_ratio_gaussian_reduced <- stan_glm(SiC_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
CP_ratio_gaussian_reduced <- stan_glm(CP_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
CN_ratio_gaussian_reduced <- stan_glm(CN_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
NP_ratio_gaussian_reduced <- stan_glm(NP_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
NoxN_ratio_gaussian_reduced <- stan_glm(NoxN_ratio ~ light + treatment, data = chems, family = gaussian(), iter = 5000)

#mucking around with other transforms on gaussian models
dPOC_sqrt_full <- stan_glm(sqrt(dPOC) ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dPON_sqrt_full <- stan_glm(sqrt(dPON) ~ light * treatment, data = chems, family = gaussian(), iter = 5000)
dPOC_sqrt_reduced <- stan_glm(sqrt(dPOC) ~ light + treatment, data = chems, family = gaussian(), iter = 5000)
dPON_sqrt_reduced <- stan_glm(sqrt(dPON) ~ light + treatment, data = chems, family = gaussian(), iter = 5000)

# Gamma full models
dPOC_gamma_full <- stan_glm(dPOC ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dPON_gamma_full <- stan_glm(dPON ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dCN_gamma_full <- stan_glm(dCN ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChla_gamma_full <- stan_glm(dChla ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChlC_gamma_full <- stan_glm(dChlC ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
d13C_gamma_full <- stan_glm(abs(d13C) ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dNOx_gamma_full <- stan_glm(dNOx ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dPO4_gamma_full <- stan_glm(dPO4 ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dNH4_gamma_full <- stan_glm(dNH4 ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dSiO3_gamma_full <- stan_glm(dSiO3 ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
SiP_ratio_gamma_full <- stan_glm(SiP_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
SiC_ratio_gamma_full <- stan_glm(SiC_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
CP_ratio_gamma_full <- stan_glm(CP_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
CN_ratio_gamma_full <- stan_glm(CN_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
NP_ratio_gamma_full <- stan_glm(NP_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
NoxN_ratio_gamma_full <- stan_glm(NoxN_ratio ~ light * treatment, data = chems, family = Gamma(link = "log"), iter = 5000)

# Gamma reduced models
dPOC_gamma_reduced <- stan_glm(dPOC ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dPON_gamma_reduced <- stan_glm(dPON ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dCN_gamma_reduced <- stan_glm(dCN ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChla_gamma_reduced <- stan_glm(dChla ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChlC_gamma_reduced <- stan_glm(dChlC ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
d13C_gamma_reduced <- stan_glm(abs(d13C) ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dNOx_gamma_reduced <- stan_glm(dNOx ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dPO4_gamma_reduced <- stan_glm(dPO4 ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dNH4_gamma_reduced <- stan_glm(dNH4 ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dSiO3_gamma_reduced <- stan_glm(dSiO3 ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
SiP_ratio_gamma_reduced <- stan_glm(SiP_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
SiC_ratio_gamma_reduced <- stan_glm(SiC_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
CP_ratio_gamma_reduced <- stan_glm(CP_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
CN_ratio_gamma_reduced <- stan_glm(CN_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
NP_ratio_gamma_reduced <- stan_glm(NP_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
NoxN_ratio_gamma_reduced <- stan_glm(NoxN_ratio ~ light + treatment, data = chems, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for Gaussian full models
dPOC_gaussian_full_loo <- loo(dPOC_gaussian_full, k_threshold = 0.7)
dPON_gaussian_full_loo <- loo(dPON_gaussian_full, k_threshold = 0.7)
dCN_gaussian_full_loo <- loo(dCN_gaussian_full, k_threshold = 0.7)
Chla_gaussian_full_loo <- loo(Chla_gaussian_full, k_threshold = 0.7)
dChlC_gaussian_full_loo <- loo(dChlC_gaussian_full, k_threshold = 0.7)
d13C_gaussian_full_loo <- loo(d13C_gaussian_full, k_threshold = 0.7)
dNOx_gaussian_full_loo <- loo(dNOx_gaussian_full, k_threshold = 0.7)
dPO4_gaussian_full_loo <- loo(dPO4_gaussian_full, k_threshold = 0.7)
dNH4_gaussian_full_loo <- loo(dNH4_gaussian_full, k_threshold = 0.7)
dSiO3_gaussian_full_loo <- loo(dSiO3_gaussian_full, k_threshold = 0.7)
SiP_ratio_gaussian_full_loo <- loo(SiP_ratio_gaussian_full, k_threshold = 0.7)
SiC_ratio_gaussian_full_loo <- loo(SiC_ratio_gaussian_full, k_threshold = 0.7)
CP_ratio_gaussian_full_loo <- loo(CP_ratio_gaussian_full, k_threshold = 0.7)
CN_ratio_gaussian_full_loo <- loo(CN_ratio_gaussian_full, k_threshold = 0.7)
NP_ratio_gaussian_full_loo <- loo(NP_ratio_gaussian_full, k_threshold = 0.7)
NoxN_ratio_gaussian_full_loo <- loo(NoxN_ratio_gaussian_full, k_threshold = 0.7)

# Compute LOO estimates for Gaussian reduced models
dPOC_gaussian_reduced_loo <- loo(dPOC_gaussian_reduced, k_threshold = 0.7)
dPON_gaussian_reduced_loo <- loo(dPON_gaussian_reduced, k_threshold = 0.7)
dCN_gaussian_reduced_loo <- loo(dCN_gaussian_reduced, k_threshold = 0.7)
Chla_gaussian_reduced_loo <- loo(Chla_gaussian_reduced, k_threshold = 0.7)
dChlC_gaussian_reduced_loo <- loo(dChlC_gaussian_reduced, k_threshold = 0.7)
d13C_gaussian_reduced_loo <- loo(d13C_gaussian_reduced, k_threshold = 0.7)
dNOx_gaussian_reduced_loo <- loo(dNOx_gaussian_reduced, k_threshold = 0.7)
dPO4_gaussian_reduced_loo <- loo(dPO4_gaussian_reduced, k_threshold = 0.7)
dNH4_gaussian_reduced_loo <- loo(dNH4_gaussian_reduced, k_threshold = 0.7)
dSiO3_gaussian_reduced_loo <- loo(dSiO3_gaussian_reduced, k_threshold = 0.7)
SiP_ratio_gaussian_reduced_loo <- loo(SiP_ratio_gaussian_reduced, k_threshold = 0.7)
SiC_ratio_gaussian_reduced_loo <- loo(SiC_ratio_gaussian_reduced, k_threshold = 0.7)
CP_ratio_gaussian_reduced_loo <- loo(CP_ratio_gaussian_reduced, k_threshold = 0.7)
CN_ratio_gaussian_reduced_loo <- loo(CN_ratio_gaussian_reduced, k_threshold = 0.7)
NP_ratio_gaussian_reduced_loo <- loo(NP_ratio_gaussian_reduced, k_threshold = 0.7)
NoxN_ratio_gaussian_reduced_loo <- loo(NoxN_ratio_gaussian_reduced, k_threshold = 0.7)

#LOOs with other transforms on gaussian models
dPOC_sqrt_full_loo <- loo(dPOC_sqrt_full, k_threshold = 0.7)
dPON_sqrt_full_loo <- loo(dPON_sqrt_full, k_threshold = 0.7)
dPOC_sqrt_reduced_loo <- loo(dPOC_sqrt_reduced, k_threshold = 0.7)
dPON_sqrt_reduced_loo <- loo(dPON_sqrt_reduced, k_threshold = 0.7)

# Compute LOO estimates for Gamma full models
dPOC_gamma_full_loo <- loo(dPOC_gamma_full, k_threshold = 0.7)
dPON_gamma_full_loo <- loo(dPON_gamma_full, k_threshold = 0.7)
dCN_gamma_full_loo <- loo(dCN_gamma_full, k_threshold = 0.7)
dChla_gamma_full_loo <- loo(dChla_gamma_full, k_threshold = 0.7)
dChlC_gamma_full_loo <- loo(dChlC_gamma_full, k_threshold = 0.7)
d13C_gamma_full_loo <- loo(d13C_gamma_full, k_threshold = 0.7)
dNOx_gamma_full_loo <- loo(dNOx_gamma_full, k_threshold = 0.7)
dPO4_gamma_full_loo <- loo(dPO4_gamma_full, k_threshold = 0.7)
dNH4_gamma_full_loo <- loo(dNH4_gamma_full, k_threshold = 0.7)
dSiO3_gamma_full_loo <- loo(dSiO3_gamma_full, k_threshold = 0.7)
SiP_ratio_gamma_full_loo <- loo(SiP_ratio_gamma_full, k_threshold = 0.7)
SiC_ratio_gamma_full_loo <- loo(SiC_ratio_gamma_full, k_threshold = 0.7)
CP_ratio_gamma_full_loo <- loo(CP_ratio_gamma_full, k_threshold = 0.7)
CN_ratio_gamma_full_loo <- loo(CN_ratio_gamma_full, k_threshold = 0.7)
NP_ratio_gamma_full_loo <- loo(NP_ratio_gamma_full, k_threshold = 0.7)
NoxN_ratio_gamma_full_loo <- loo(NoxN_ratio_gamma_full, k_threshold = 0.7)

# Compute LOO estimates for Gamma reduced models
dPOC_gamma_reduced_loo <- loo(dPOC_gamma_reduced, k_threshold = 0.7)
dPON_gamma_reduced_loo <- loo(dPON_gamma_reduced, k_threshold = 0.7)
dCN_gamma_reduced_loo <- loo(dCN_gamma_reduced, k_threshold = 0.7)
dChla_gamma_reduced_loo <- loo(dChla_gamma_reduced, k_threshold = 0.7)
dChlC_gamma_reduced_loo <- loo(dChlC_gamma_reduced, k_threshold = 0.7)
d13C_gamma_reduced_loo <- loo(d13C_gamma_reduced, k_threshold = 0.7)
dNOx_gamma_reduced_loo <- loo(dNOx_gamma_reduced, k_threshold = 0.7)
dPO4_gamma_reduced_loo <- loo(dPO4_gamma_reduced, k_threshold = 0.7)
dNH4_gamma_reduced_loo <- loo(dNH4_gamma_reduced, k_threshold = 0.7)
dSiO3_gamma_reduced_loo <- loo(dSiO3_gamma_reduced, k_threshold = 0.7)
SiP_ratio_gamma_reduced_loo <- loo(SiP_ratio_gamma_reduced, k_threshold = 0.7)
SiC_ratio_gamma_reduced_loo <- loo(SiC_ratio_gamma_reduced, k_threshold = 0.7)
CP_ratio_gamma_reduced_loo <- loo(CP_ratio_gamma_reduced, k_threshold = 0.7)
CN_ratio_gamma_reduced_loo <- loo(CN_ratio_gamma_reduced, k_threshold = 0.7)
NP_ratio_gamma_reduced_loo <- loo(NP_ratio_gamma_reduced, k_threshold = 0.7)
NoxN_ratio_gamma_reduced_loo <- loo(NoxN_ratio_gamma_reduced, k_threshold = 0.7)

# Compare the models using loo_compare()
dPOC_loo_compare <- loo_compare(dPOC_gaussian_full_loo, dPOC_gaussian_reduced_loo, dPOC_gamma_full_loo, dPOC_gamma_reduced_loo)
dPON_loo_compare <- loo_compare(dPON_gaussian_full_loo, dPON_gaussian_reduced_loo, dPON_gamma_full_loo, dPON_gamma_reduced_loo)
dCN_loo_compare <- loo_compare(dCN_gaussian_full_loo, dCN_gaussian_reduced_loo, dCN_gamma_full_loo, dCN_gamma_reduced_loo)
dChla_loo_compare <- loo_compare(Chla_gaussian_full_loo, Chla_gaussian_reduced_loo, dChla_gamma_full_loo, dChla_gamma_reduced_loo)
dChlC_loo_compare <- loo_compare(dChlC_gaussian_full_loo, dChlC_gaussian_reduced_loo, dChlC_gamma_full_loo, dChlC_gamma_reduced_loo)
d13C_loo_compare <- loo_compare(d13C_gaussian_full_loo, d13C_gaussian_reduced_loo, d13C_gamma_full_loo, d13C_gamma_reduced_loo)
dNOx_loo_compare <- loo_compare(dNOx_gaussian_full_loo, dNOx_gaussian_reduced_loo, dNOx_gamma_full_loo, dNOx_gamma_reduced_loo)
dPO4_loo_compare <- loo_compare(dPO4_gaussian_full_loo, dPO4_gaussian_reduced_loo, dPO4_gamma_full_loo, dPO4_gamma_reduced_loo)
dNH4_loo_compare <- loo_compare(dNH4_gaussian_full_loo, dNH4_gaussian_reduced_loo, dNH4_gamma_full_loo, dNH4_gamma_reduced_loo)
dSiO3_loo_compare <- loo_compare(dSiO3_gaussian_full_loo, dSiO3_gaussian_reduced_loo, dSiO3_gamma_full_loo, dSiO3_gamma_reduced_loo)
SiP_ratio_loo_compare <- loo_compare(SiP_ratio_gaussian_full_loo, SiP_ratio_gaussian_reduced_loo, SiP_ratio_gamma_full_loo, SiP_ratio_gamma_reduced_loo)
SiC_ratio_loo_compare <- loo_compare(SiC_ratio_gaussian_full_loo, SiC_ratio_gaussian_reduced_loo, SiC_ratio_gamma_full_loo, SiC_ratio_gamma_reduced_loo)
CP_ratio_loo_compare <- loo_compare(CP_ratio_gaussian_full_loo, CP_ratio_gaussian_reduced_loo, CP_ratio_gamma_full_loo, CP_ratio_gamma_reduced_loo)
CN_ratio_loo_compare <- loo_compare(CN_ratio_gaussian_full_loo, CN_ratio_gaussian_reduced_loo, CN_ratio_gamma_full_loo, CN_ratio_gamma_reduced_loo)
NP_ratio_loo_compare <- loo_compare(NP_ratio_gaussian_full_loo, NP_ratio_gaussian_reduced_loo, NP_ratio_gamma_full_loo, NP_ratio_gamma_reduced_loo)
NoxN_ratio_loo_compare <- loo_compare(NoxN_ratio_gaussian_full_loo, NoxN_ratio_gaussian_reduced_loo, NoxN_ratio_gamma_full_loo, NoxN_ratio_gamma_reduced_loo)

dPOC_loo_compare
dPON_loo_compare 
dCN_loo_compare
dChla_loo_compare 
dChlC_loo_compare 
d13C_loo_compare
dNOx_loo_compare
dPO4_loo_compare
dNH4_loo_compare
dSiO3_loo_compare
SiP_ratio_loo_compare
SiC_ratio_loo_compare
CP_ratio_loo_compare
NP_ratio_loo_compare
NoxN_ratio_loo_compare
