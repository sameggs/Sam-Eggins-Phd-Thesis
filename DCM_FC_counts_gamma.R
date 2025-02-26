setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment") #setwd
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
cells <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment/dcm_flowcyt.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(light = as.factor(light), 
         treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         bottle = as.factor(bottle))
# Replace N/A with NA
cells[cells == "N/A"]  <- NA
# make values numeric
cells <- cells %>% mutate_at(vars(counts:FB3), as.numeric)
# Reorder  levels
cells$treatment <- factor(cells$treatment, levels = c("85","+DFB", "Control", "+Fe"), ordered = TRUE)
# summing nanophytoplankton  
nano_gates <- c("nano 1", "nano 2", "nano 3")
cells <- cells %>%
  group_by(light, treatment, bottle, time, gate = ifelse(gate %in% nano_gates, "nanos", gate)) %>%
  summarise(across(c(counts, Fsize, FV12, FB3), ~sum(.x, na.rm = TRUE)), .groups = "drop")

# Add log(counts) and sqrt(counts) columns
cells <- cells %>%
  mutate(log_counts = log(counts), sqrt_counts = sqrt(counts))

# grab data from initial conditions
time_0 <- cells %>% filter(treatment == "85") %>% mutate_at(vars(counts:FB3), as.numeric)
cells <- cells %>% filter(time == 10)

#Check the structure
str(cells)

# ----- Investigating the density distributions of the data ------------------------------------------------------
# List of unique gate variables
gate_vars <- unique(cells$gate)

# List of count variables
count_vars <- c("counts", "log_counts", "sqrt_counts")

# Generate a list of plots
plots_list <- purrr::map(gate_vars, ~ {
  gate_var <- .x
  purrr::map(count_vars, ~ {
    count_var <- .x
    data_subset <- cells[cells$gate == gate_var, ]
    ggplot(data_subset, aes(x = !!sym(count_var))) + 
      geom_density(alpha = .2, fill = "#FF6FF6") +
      ggtitle(paste(gate_var, "-", count_var))
  })
}) %>% purrr::flatten()

# Combine all plots in a grid
density_plots <- wrap_plots(plots_list, ncol = length(count_vars))  # adjust ncol and nrow as needed
print(density_plots)

# FOR THE WHOLE COUNT SET OF DATA .................................
# Generate a list of plots
plots_list <- purrr::map(count_vars, ~ {
  count_var <- .x
  ggplot(cells, aes(x = !!sym(count_var))) + 
    geom_density(alpha = .2, fill = "#FF6FF6") +
    ggtitle(paste("Density of", count_var))
})

# Combine all plots in a grid
density_plot_all <- wrap_plots(plots_list, ncol = length(count_vars))  # adjust ncol as needed
print(density_plot_all)

#--------- Reshape data ----------------------------------------------------------------------------------------
counts <- cells %>% 
  select(light, treatment, bottle, gate, counts) %>%
  pivot_wider(names_from = gate, values_from = counts)
counts <- counts %>% select(-all_phyto)
str(counts)

#--- plots of variance against the mean ----------------------------------------------------------------
# Get the count column names
count_cols <- colnames(counts)[4:length(colnames(counts))]

# Calculate the mean and variance for each unique combination of light and treatment
counts_summary <- counts %>%
  group_by(light, treatment) %>%
  summarise(across(all_of(count_cols), list(mean = mean, var = var), na.rm = TRUE))

# Function to create scatterplot
create_scatterplot <- function(var) {
  ggplot(counts_summary, aes(x = !!sym(paste0(var, "_mean")), y = !!sym(paste0(var, "_var")))) +
    geom_point() +
    theme_bw() +
    xlab("Mean") +
    ylab("Variance") +
    ggtitle(paste0(var, " Var vs. Mean"))
}

# Create a scatter plot for each column
plot_list <- map(count_cols, ~ {
  var <- .x
  create_scatterplot(var)
})

# Patchwork the plots together
variance_plot <- wrap_plots(plot_list)
print(variance_plot)

#------FITTING MODELS fo counts:  -----------------------------------------------------------------------------------
# a model for the entire set of counts, with gate as a variable
counts_model <- glm((counts) ~ light * treatment * gate, data = cells, family = Gamma(link = "log"))
# quick diagnostics
par(mfrow = c(3, 2))
plot(counts_model, which = 1:6)

# ---- QQ plots for the gamma model ----------------------------------------------------------
# Function to calculate shape and rate for a given model
calculate_shape_rate <- function(model) {
  fitted_values <- fitted(model)
  dispersion <- summary(model)$dispersion
  shape <- 1 / dispersion
  rate <- shape / mean(fitted_values)
  return(list(shape = shape, rate = rate))
}

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

# Function to calculate standardized residuals
calculate_standardized_residuals <- function(model) {
  resid(model, type = "pearson")
}

# Calculate shape and rate for the counts_model
shape_rate_counts_model <- calculate_shape_rate(counts_model)

# Calculate standardized residuals for the counts_model
std_residuals_counts_model <- calculate_standardized_residuals(counts_model)

# Create the QQ plot
create_qq_plot(std_residuals_counts_model, "Counts Model", shape_rate_counts_model$shape, shape_rate_counts_model$rate)


# -----Function to create residual vs fitted plot for model----------------------------------------------
res_fit_plot <- function(model, model_name) {
  res_fit <- broom::augment(model)
  
  plot <- ggplot(data = res_fit, aes(x = .fitted, y = .resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(paste(model_name, "Residuals vs Fitted")) +
    theme_bw()
  
  return(plot)
}

# Create a Residuals vs Fitted plot for the counts_model
res_fit_plot_counts_model <- res_fit_plot(counts_model, "counts_model")

# Print the plot
print(res_fit_plot_counts_model)

# -------- EMMEANS for the models ----------------------------------------------------------------------------
# Create a list where model names are matched with the corresponding model list
emm_all <- summary(emmeans(counts_model, specs = pairwise ~ treatment : light : gate))
emm_trt <- summary(emmeans(counts_model, specs = pairwise ~ treatment | light | gate))
emm_light <- summary(emmeans(counts_model, specs = pairwise ~ light | treatment | gate))

# Combine emmeans and contrasts for each data frame
counts_contrasts <- bind_rows(emmeans = emm_all$emmeans, 
                            light = emm_trt$contrasts, 
                            trt = emm_light$contrasts)

# Export the combined data frames to an Excel file
write_xlsx(counts_contrasts, "FC_counts_emmeans.xlsx")

# ---Calculate the mean and standard error for each unique combination----------------------------------------------------
# select columns that are required
vars <- c("bacteria", "peuks", "nanos", "micro", "syn", "Cocco")
#create summary
counts_summary <- counts %>%
  mutate(bacteria = bacteria / 1000, 
         nanos = nanos / 1000) %>%
  group_by(light, treatment) %>%
  summarise(across(all_of(vars), 
                   list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))), 
                   .names = "{.col}_{.fn}"))


# Generate function for facet grid with point and error bars
generate_plots <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, 
                           show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, 
                           legend_position = "top", hline_y) {
  
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), color = light)) +
    geom_hline(yintercept = hline_y, linetype = "dashed", color = "#95C11FFF") +
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

# ----- plot the data --------------------------------------------------------------------------------------------------------
plot_peuks <- generate_plots(counts_summary, "peuks_mean", "peuks_se", "peuks", show_xlab = FALSE, y_min = 0, y_max = 4000,
                              y_lab = expression(paste("cells mL"^-1,)),legend_position = "none", hline_y = 317) +
                              scale_y_continuous(breaks = seq(0, 4000, by = 500)) 
plot_nanos <- generate_plots(counts_summary, "nanos_mean", "nanos_se", "nanos", show_xlab = FALSE, y_min = 0, y_max = 30,
                             y_lab = expression(paste("10"^3, "cells mL"^-1)), legend_position = "none", hline_y = 1.060) + 
                            scale_y_continuous(breaks = seq(0, 30, by = 5))
plot_micro <- generate_plots(counts_summary, "micro_mean", "micro_se", "micro", show_xlab = FALSE, y_min = 0, y_max = 1400,
                             y_lab = expression(paste("cells mL"^-1)), legend_position = "none", hline_y = 16) + 
                              scale_y_continuous(breaks = seq(0, 1400, by = 200))
plot_cyano <- generate_plots(counts_summary, "syn_mean", "syn_se", "syn", show_xlab = FALSE, y_min = 0, y_max = 600,
                             y_lab = expression(paste("cells mL"^-1)), legend_position = "none", hline_y = 124) + 
                              scale_y_continuous(breaks = seq(0, 600, by = 200))
plot_cocco <- generate_plots(counts_summary, "Cocco_mean", "Cocco_se", "Cocco", show_xlab = FALSE, y_min = 0, y_max = 500,
                             y_lab = expression(paste("cells mL"^-1)), legend_position = "none", hline_y = 188) + 
                             scale_y_continuous(breaks = seq(0, 500, by = 100))
plot_bacteria <- generate_plots(counts_summary, "bacteria_mean", "bacteria_se", "syn", show_xlab = FALSE, y_min = 0, y_max = 1200,
                                y_lab = expression(paste("10"^3,"cells mL"^-1)), legend_position = "none", hline_y = 264.95) + 
                              scale_y_continuous(breaks = seq(0, 1200, by = 400))


#---adding asterisks for significance ----------------------------------------------------------------
asterisks_peuks <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(8.46, 8.66, 8.3, 6.36, 5.73, 5.3, 8.9), 
  label = c("a", "a", "a", "a", "b", "c", "all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_nano <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(9.94, 9.71, 10.5, 7.5, 7.1, 7.45, 10.9), 
  label = c("a", "a", "b", "a", "a", "a", "all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_micro <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(5.4, 5.56, 7.44, 3.64, 3.53, 3.96, 7.8), 
  label = c("a", "a", "b", "a", "a", "a", "all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_cocco <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(6.53, 6.1, 6.06, 4.2, 3.83, 3.83, 6.9), 
  label = c("a", "b", "b", "a", "a", "a", "all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_cyano <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.38, 0.738), y_pos = c(6.375, 6.355, 6.379, 5.22, 4.99, 5.17, 6.355, 6.379), 
  label = c("a", "a", "a", "a", "a", "a", "\u263C", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black", "black"))

asterisks_bacteria <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.38), y_pos = c(13, 13.7, 13.85, 14, 12.5, 12.7, 13.7), 
  label = c("a", "a", "a", "a", "b", "b", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

# ------------ VISREG plots --------------------------------------------------------------------------------------------------------
visreg_peuks <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                     cond = list(gate="peuks"), gg = TRUE, points = list(shape = 19), 
                     scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 4000), breaks = seq(5, 9, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"pico"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_peuks, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_peuks$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_nanos <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                            cond = list(gate="nanos"), gg = TRUE, points = list(shape = 19), 
                            scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(7, 11, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"nano"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_nano, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_nano$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_micro <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                       cond = list(gate="micro"), gg = TRUE, points = list(shape = 19), 
                       scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 1400), breaks = seq(3, 8, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"micro"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_micro, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_micro$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_cocco <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                       cond = list(gate="Cocco"), gg = TRUE, points = list(shape = 19), 
                       scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(3, 7, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"cocco"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(), axis.text.y = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_cocco, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_cocco$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_cyano <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                       cond = list(gate="syn"), gg = TRUE, points = list(shape = 19), 
                       scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 600), breaks = seq(4, 7, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"cyano"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(), axis.text.y = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_cyano, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_cyano$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_bacteria <- visreg(counts_model, "treatment", by="light", overlay = TRUE, data = cells,
                       cond = list(gate="bacteria"), gg = TRUE, points = list(shape = 19), 
                       scale = "response", partial = TRUE) + 
  theme_bw() + 
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
  ylab("") +  
  xlab("") +
  scale_y_continuous(limits = c(0, 1200000), breaks = seq(12, 15, by = 1),
                     sec.axis = sec_axis(~., name = expression(ln(cells[~"bacteria"])))) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y.left = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_bacteria, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_bacteria$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# -----------putting plots together------------------------------------------------------------------------------------------------------
FC_counts_model_plot <- guide_area() / (plot_peuks | visreg_peuks) / (plot_nanos | visreg_nanos)/ (plot_bacteria | visreg_bacteria) / (plot_cyano | visreg_cyano) +
  #(plot_micro | visreg_micro) / (plot_cocco | visreg_cocco) +
  plot_layout(guides = 'collect', nrow = (5), heights = c(1,10,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D','E','F','G','H'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

FC_bacteria_model_plot <- guide_area() / (plot_bacteria | visreg_bacteria) / (plot_cyano | visreg_cyano) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("FC_counts_model_plot3.svg", FC_counts_model_plot, width = 7.2, height = 9.8) 
ggsave("FC_bacteria_model_plot.svg", FC_bacteria_model_plot, width = 7.2, height = 6.8) 
