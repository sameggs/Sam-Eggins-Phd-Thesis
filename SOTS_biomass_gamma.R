setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment") #setwd
#----------------loading packages--------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(patchwork)
library(broom)
library(emmeans)
library(visreg)
library(writexl)
library(rstan)
library(rstanarm)

#------load the POC PON Chla data  into R---------------------------------------------------------------------------------------
chems <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/sots_chems.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(inhibitor = as.factor(inhibitor),
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         bottle = as.factor(bottle))
# Replace N/A with NA
chems[chems == "N/A"]  <- NA
# make values numeric
chems <- chems %>% mutate_at(vars(POC:dChlC), as.numeric)
# Reorder  levels
chems$treatment <- factor(chems$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)

# grab data from initial conditions
time_0 <- chems %>% filter(time == "0") %>% mutate_at(vars(POC:Chla), as.numeric)
chems <- chems %>% filter(time == 6)

#Check the structure
str(chems)

#----checking density distributions.can edit the exressions to find the best transformations ------------------------------
p1 <- ggplot(chems, aes(x = (dPOC))) + geom_density(alpha = .2, fill = "#FF6FF6")
p2 <- ggplot(chems, aes(x = (dPON))) + geom_density(alpha = .2, fill = "#FF6FF6")
p3 <- ggplot(chems, aes(x = log(dChla))) + geom_density(alpha = .2, fill = "#FF6FF6")
p4 <- ggplot(chems, aes(x = (abs(d13C)))) + geom_density(alpha = .2, fill = "#FF6FF6")
p5 <- ggplot(chems, aes(x = (dCN))) + geom_density(alpha = .2, fill = "#FF6FF6")
p6 <- ggplot(chems, aes(x = (ChlC_ratio))) + geom_density(alpha = .2, fill = "#FF6FF6")

# Combine plots using patchwork syntax
combined_plot <- p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 2)  # arrange plots in 2 columns

# Display combined plot
print(combined_plot)

# ---Compute the mean and variance for each group defined by light, treatment, and size-----------------------
# select columns that are required
vars <- c("dPOC", "dPON", "dChla", "d13C", "CN_ratio", "ChlC_ratio")

#create summary data frame with variable names as "{.fn}_{.col}"
chems_summary <- chems %>%
  group_by(treatment) %>%
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

#------FITTING MODELS :  -----------------------------------------------------------------------------------
# Define a list of the response variables for Gaussian models
resp_vars <- c("dPOC", "dPON", "dCN", "d13C", "dChla", "dChlC")
# Define a function to fit the Gaussian model
fit_gauss_model <- function(var) {
  formula <- as.formula(paste(var, "~ treatment"))
  if (var == "dPOC") {
    formula <- as.formula(paste("(dPOC^2)/1000 ~ treatment"))}
  if (var == "dPON") {
    formula <- as.formula(paste("(dPON^2)/10 ~ treatment"))}
  glm(formula, data = chems, family = gaussian())
}
# Use lapply to fit the Gaussian models for all variables and store the results in a list
gauss_models_list <- setNames(lapply(resp_vars, fit_gauss_model), resp_vars)

# Define a gamma function to fit the model
fit_model <- function(var) {
  formula <- as.formula(paste(var, "~ treatment"))
  if (var == "d13C") {
    formula <- as.formula(paste("abs(d13C) ~ treatment"))}
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
  set_names(c("dPOC", "dPON", "dCN", "d13C", "dChla", "ChlC_ratio"))

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
  "dPOC" = list(model_list = gauss_models_list, model_name = "dPOC"),
  "dPON" = list(model_list = gauss_models_list, model_name = "dPON"),
  "dCN" = list(model_list = gauss_models_list, model_name = "dCN"),
  "dChla" = list(model_list = gamma_models_list, model_name = "dChla"),
  "dChlC" = list(model_list = gamma_models_list, model_name = "dChlC"),
  "d13C" = list(model_list = gauss_models_list, model_name = "d13C")
)

# Calculate emmeans for each model
emm_nuts <- map(model_info, 
                ~emmeans(.x$model_list[[.x$model_name]], pairwise ~ treatment)) %>% 
  map(~summary(., infer = TRUE)) %>% 
  set_names(names(model_info))

# Combine emmeans and contrasts for each data frame
emmeans_dPOC <- bind_rows(emmeans = emm_nuts$dPOC$emmeans, light = emm_nuts$dPOC$contrasts)
emmeans_dPON <- bind_rows(emmeans = emm_nuts$dPON$emmeans, light = emm_nuts$dPON$contrasts)
emmeans_dCN <- bind_rows(emmeans = emm_nuts$dCN$emmeans, light = emm_nuts$dCN$contrasts)
emmeans_dChla <- bind_rows(emmeans = emm_nuts$dChla$emmeans, light = emm_nuts$dChla$contrasts)
emmeans_dChlC <- bind_rows(emmeans = emm_nuts$dChlC$emmeans, light = emm_nuts$dChlC$contrasts)
emmeans_d13C <- bind_rows(emmeans = emm_nuts$d13C$emmeans, light = emm_nuts$d13C$contrasts)

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
write_xlsx(combined_CNChl_emm, "SOTS_CNChl_emmeans.xlsx")

# ---Calculate the mean and standard error for each unique combination----------------------------------------------------
# select columns that are required
vars <- c("dPOC", "dPON", "dCN", "D13C", "dChla", "dChlC")
#create summary
chems_summary <- chems %>%
  group_by(treatment) %>%
  summarise(across(all_of(vars), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))), 
                   .names = "{.col}_{.fn}"))

# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, 
                                          show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, 
                                          legend_position = "top", hline_y) {
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean))) +
    geom_hline(yintercept = hline_y, linetype = "dashed", color = "#95C11FFF") +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    geom_point(size = 3, shape = 21, colour = "black", fill = "black") +
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

# -------------------------------Create plots for POC PON Chla -----------------------------------------------------------------------
# plots for the data...
plot_dPOC <- generate_plots_with_errorbars(chems_summary, "dPOC_mean", "dPOC_se", "dPOC", show_xlab = FALSE, y_min = 0, y_max = 400,
                                           y_lab = expression(paste(Delta,"POC (",mu,"g L"^-1,")")), legend_position = "none", hline_y = -100) + 
  scale_y_continuous(breaks = seq(0, 400, by = 50)) 

plot_dPON <- generate_plots_with_errorbars(chems_summary, "dPON_mean", "dPON_se", "dPON", show_xlab = FALSE, y_min = 0, y_max = 90, 
                                           y_lab = expression(paste(Delta,"PON (",mu,"g L"^-1,")")), legend_position = "none", hline_y = -100) + 
  scale_y_continuous(breaks = seq(0, 90, by = 15))

plot_dCN <- generate_plots_with_errorbars(chems_summary, "dCN_mean", "dCN_se", "dCN", show_xlab = TRUE, y_min = 5, y_max = 8,
                                          y_lab = expression(paste("C:N ratio")), legend_position = "none", hline_y = -10) +
  scale_y_continuous(breaks = seq(5, 8, by = 0.5))

plot_Chla <- generate_plots_with_errorbars(chems_summary, "dChla_mean", "dChla_se", "dChla", show_xlab = FALSE, y_min = 0, y_max = 8,
                                           y_lab = expression(paste(Delta,"Chl-",italic("a")," (",mu,"g L"^-1,")")), legend_position = "none", hline_y = -10) + 
  scale_y_continuous(breaks = seq(0, 8, by = 1))  

plot_ChlC_ratio <- generate_plots_with_errorbars(chems_summary, "dChlC_mean", "dChlC_se", "dChlC", show_xlab = TRUE, y_min = 0, y_max = 30,
                                                 y_lab = expression(paste(Delta,"Chl:C (mg g"^-1,")")), legend_position = "none", hline_y = -10) +
  scale_y_continuous(breaks = seq(0, 30, by = 5))

#---adding asterisks for significance PONPOC ----------------------------------------------------------------
asterisks_dPOC <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(90, 105, 170), 
  label = c("a", "a", "b"), colours = c("black", "black", "black"))

asterisks_dPON <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(290, 320, 550), 
  label = c("a", "a", "b"), colours = c("black", "black", "black"))

asterisks_dCN <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(7.1, 7.2, 6.9), 
  label = c("a", "a", "a"), colours = c("black", "black", "black"))

asterisks_dChla <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(1.37, 0.85, 2.25), 
  label = c("a", "b", "c"), colours = c("black", "black", "black"))

asterisks_dChlC <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(2.76, 2.24, 3.24), 
  label = c("a", "b", "c"), colours = c("black", "black", "black"))

# visreg plots...
visreg_POC <- visreg(gauss_models_list$dPOC, top = "points", partial = TRUE, gg = TRUE,
                     trans = function(x) sqrt(x * 1000),
                     line.par = list(col = "#706F6FFF"),
                     fill.par = list(fill = "#706F6F", alpha = 0.33),
                     points.par = list(col = "black", shape = 21, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, 25),
                     sec.axis = sec_axis(~., name = expression(Delta * POC^2 / 1000),
                                         breaks = seq(0, 400, 25))) +
  ylab("") + xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dPOC, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_dPOC$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


visreg_PON <- visreg(gauss_models_list$dPON, top = "points", partial = TRUE, gg = TRUE,
                     trans = function(x) sqrt(x * 10),
                     line.par = list(col = "#706F6FFF"),
                     fill.par = list(fill = "#706F6F", alpha = 0.33),
                     points.par = list(col = "black", shape = 21, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 15),
                     sec.axis = sec_axis(~ ., name = expression(Delta * PON^2 / 10),
                                         breaks = seq(0, 90, 15))) +
  ylab("") + xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dPON, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_dPON$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


visreg_CN_ratio <- visreg(gauss_models_list$dCN, "treatment", top = "points", partial = TRUE, gg = TRUE,
                          line.par = list(col = "#706F6FFF"),
                          fill.par = list(fill = "#706F6F", alpha = 0.33),
                          points.par = list(col = "black", shape = 21, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(5, 8), breaks = seq(5, 8, 0.5)) +
  ylab("") + xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_dCN, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_dCN$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


visreg_Chla <- visreg(gamma_models_list$dChla, "treatment", top = "points", partial = TRUE, gg = TRUE,
                      scale = "response",
                      line.par = list(col = "#706F6FFF"),
                      fill.par = list(fill = "#706F6F", alpha = 0.33),
                      points.par = list(col = "black", shape = 21, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, by = 1)) +
  ylab("") +
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.position = "top") +
  geom_text(data = asterisks_dChla, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_dChla$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


visreg_ChlC <- visreg(gamma_models_list$dChlC, "treatment", top = "points", partial = TRUE, gg = TRUE,
                      scale = "response",
                      line.par = list(col = "#706F6FFF"),
                      fill.par = list(fill = "#706F6F", alpha = 0.2),
                      points.par = list(col = "black", shape = 21, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30),
                     breaks = seq(0, 30, by = 5)) +
  ylab("") +
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        legend.position = "top") +
  geom_text(data = asterisks_dChlC, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_dChlC$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))



# combined plots and printing
POCPON_plot <- (plot_dPOC | visreg_POC) / (plot_dPON | visreg_PON) / (plot_dCN | visreg_CN_ratio) +
  plot_layout(nrow = (3), heights = c(10,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd', 'e', 'f'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        plot.tag.position = c(0.01, 0.97))

Chl_plot <- (plot_Chla | visreg_Chla) / (plot_d13C | visreg_d13C) / (plot_dPOC | visreg_POC) +
  plot_layout(nrow = (3), heights = c(10,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        plot.tag.position = c(0.01, 0.97)) 

ggsave("POCPON_model_plot_2.svg", POCPON_plot, width = 7.2, height = 7.1)
ggsave("Chl_model_plot_2.svg", Chl_plot, width = 7.2, height = 7.1)

# -------------------------------Create plots for 13C -----------------------------------------------------------------------
# plot for data... PDB
plot_d13C <- generate_plots_with_errorbars(chems_summary, "D13C_mean", "D13C_se", "D13C", show_xlab = TRUE, y_min = -25, y_max = -22,
                                           y_lab = expression(paste(delta ^13,"C")), legend_position = "none", hline_y = -23.5668) + 
  scale_y_continuous(breaks = seq(-22, -25, by = -1))

asterisks_d13C <- data.frame(
    x_pos = c(0.212, 0.57, 0.928), y_pos = c(1.04, 0.43, 0.69), 
    label = c("a", "b", "ab"), colours = c("black", "black", "black"))

visreg_d13C <- visreg(gamma_models_list$d13C, "treatment", top = "points", partial = TRUE, gg = TRUE,
                     scale = "response",           # Use the linear predictor (log scale)
                      trans = exp, # Transform: negative of the exponential
                      line.par = list(col = "#706F6FFF"),
                      fill.par = list(fill = "#706F6F", alpha = 0.2),
                      points.par = list(col = "black", shape = 23, size = 2)) +
  theme_bw() +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, by = 1)) +
  ylab("") +
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        #axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.direction = "horizontal",
        legend.position = "top") +
  geom_text(data = asterisks_d13C, aes(x = x_pos, y = y_pos, label = label),
            colour = asterisks_d13C$colours, size = 4, inherit.aes = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


# combine plots and printing
d13C_plot <- plot_d13C | visreg_d13C  +
  plot_layout(guides = 'collect', nrow = (1), heights = c(10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("13C_model_plot.svg", d13C_plot, width = 7.2, height = 3.52)

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
# Gaussian full models
dPOC_gaussian_full <- stan_glm(dPOC ~ treatment, data = chems, family = gaussian(), iter = 5000)
dPON_gaussian_full <- stan_glm(dPON ~ treatment, data = chems, family = gaussian(), iter = 5000)
dCN_gaussian_full <- stan_glm(dCN ~ treatment, data = chems, family = gaussian(), iter = 5000)
Chla_gaussian_full <- stan_glm(dChla ~ treatment, data = chems, family = gaussian(), iter = 5000)
dChlC_gaussian_full <- stan_glm(dChlC ~ treatment, data = chems, family = gaussian(), iter = 5000)
d13C_gaussian_full <- stan_glm(abs(d13C) ~ treatment, data = chems, family = gaussian(), iter = 5000)

# Gamma full models
dPOC_gamma_full <- stan_glm(dPOC ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dPON_gamma_full <- stan_glm(dPON ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dCN_gamma_full <- stan_glm(dCN ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChla_gamma_full <- stan_glm(dChla ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
dChlC_gamma_full <- stan_glm(dChlC ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)
d13C_gamma_full <- stan_glm(abs(d13C) ~ treatment, data = chems, family = Gamma(link = "log"), iter = 5000)

# Compute LOO estimates for Gaussian full models
dPOC_gaussian_full_loo <- loo(dPOC_gaussian_full, k_threshold = 0.7)
dPON_gaussian_full_loo <- loo(dPON_gaussian_full, k_threshold = 0.7)
dCN_gaussian_full_loo <- loo(dCN_gaussian_full, k_threshold = 0.7)
Chla_gaussian_full_loo <- loo(Chla_gaussian_full, k_threshold = 0.7)
dChlC_gaussian_full_loo <- loo(dChlC_gaussian_full, k_threshold = 0.7)
d13C_gaussian_full_loo <- loo(d13C_gaussian_full, k_threshold = 0.7)

# Compute LOO estimates for Gamma full models
dPOC_gamma_full_loo <- loo(dPOC_gamma_full, k_threshold = 0.7)
dPON_gamma_full_loo <- loo(dPON_gamma_full, k_threshold = 0.7)
dCN_gamma_full_loo <- loo(dCN_gamma_full, k_threshold = 0.7)
dChla_gamma_full_loo <- loo(dChla_gamma_full, k_threshold = 0.7)
dChlC_gamma_full_loo <- loo(dChlC_gamma_full, k_threshold = 0.7)
d13C_gamma_full_loo <- loo(d13C_gamma_full, k_threshold = 0.7)

# Compare the models using loo_compare()
dPOC_loo_compare <- loo_compare(dPOC_gaussian_full_loo, dPOC_gamma_full_loo)
dPON_loo_compare <- loo_compare(dPON_gaussian_full_loo, dPON_gamma_full_loo)
dCN_loo_compare <- loo_compare(dCN_gaussian_full_loo, dCN_gamma_full_loo)
dChla_loo_compare <- loo_compare(Chla_gaussian_full_loo, dChla_gamma_full_loo)
dChlC_loo_compare <- loo_compare(dChlC_gaussian_full_loo, dChlC_gamma_full_loo)
d13C_loo_compare <- loo_compare(d13C_gaussian_full_loo, d13C_gamma_full_loo)

dPOC_loo_compare
dPON_loo_compare 
dCN_loo_compare
dChla_loo_compare 
dChlC_loo_compare 
d13C_loo_compare
