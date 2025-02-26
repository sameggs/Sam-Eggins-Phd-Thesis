setwd("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment") #setwd
#----------------loading packages--------------------------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggtext)
library(cowplot)
library(patchwork)
library(broom)
library(emmeans)
library(visreg)
library(writexl)
library(rstan)
library(rstanarm)

#------load the FC data  into R---------------------------------------------------------------------------------------
cells <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/sots_flowcyt.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         domain = as.factor(domain),
         bottle = as.factor(bottle))
# Replace N/A with NA
cells[cells == "N/A"]  <- NA
# make values numeric
cells <- cells %>% mutate_at(vars(counts:FV12), as.numeric)
# Reorder  levels
cells$treatment <- factor(cells$treatment, levels = c("+DFB", "Control", "+Fe"), ordered = TRUE)
# summing nanophytoplankton  
nano_gates <- c("nano 1", "nano 2", "nano 3")
cells <- cells %>%
  group_by(treatment, bottle, time, domain, gate = ifelse(gate %in% nano_gates, "nanos", gate)) %>%
  summarise(across(c(counts, Fsize, FV12, FB3), ~sum(.x, na.rm = TRUE)), .groups = "drop")

# Add log(counts) and sqrt(counts) columns and gate as factor
cells <- cells %>%
  mutate(log_counts = log(counts), sqrt_counts = sqrt(counts),
         gate = as.factor(gate))

cells$gate <- factor(cells$gate, levels = c("bacteria", "syn", "peuks", "nanos", "micro", "Cocco", "all_phyto"), ordered = TRUE)

# grab data from initial conditions
time_0 <- cells %>% filter(time == "0") %>% mutate_at(vars(counts:FV12), as.numeric)
cells <- cells %>% filter(time == 5)

# Remove all_phyto data
cells <- cells %>% filter(gate != "all_phyto")

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
  select(treatment, bottle, gate, counts) %>%
  pivot_wider(names_from = gate, values_from = counts)
str(counts)

#--- plots of variance against the mean ----------------------------------------------------------------
# Get the count column names
count_cols <- colnames(counts)[4:length(colnames(counts))]

# Calculate the mean and variance for each unique combination of light and treatment
counts_summary <- counts %>%
  group_by(treatment) %>%
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
counts_model <- glm((counts) ~ treatment * gate, data = cells, family = Gamma(link = "log"))
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
emm_all <- summary(emmeans(counts_model, specs = pairwise ~ treatment : gate))
emm_trt <- summary(emmeans(counts_model, specs = pairwise ~ treatment | gate))

# Combine emmeans and contrasts for each data frame
counts_contrasts <- bind_rows(emmeans = emm_all$emmeans, 
                              trt = emm_trt$contrasts)

# Export the combined data frames to an Excel file
write_xlsx(counts_contrasts, "SOTS_FC_counts_emmeans.xlsx")

# ---Calculate the mean and standard error for each unique combination----------------------------------------------------
# summarise the Fsize from cells df
counts_summary <- cells %>%
  dplyr::filter(!(gate %in% c("all_phyto"))) %>%
  dplyr::group_by(treatment, gate) %>%
  dplyr::summarize(
    mean_counts = mean(counts, na.rm = TRUE) / 1000,
    se_counts = (sd(counts, na.rm = TRUE) / sqrt(n())) / 1000) %>%
  mutate(
    domain = case_when(
      gate %in% c("nanos", "peuks", "Cocco", "micro") ~ "Eukaryotes",
      gate %in% c("syn", "bacteria") ~ "Prokaryotes",
      TRUE ~ "Unknown"
    ),
    treatment = as.factor(treatment),
    gate = as.factor(gate) %>% fct_recode("Coccolithophores" = "Cocco", "Microeukaryotes" = "micro",
                                          "Nanoeukaryotes" = "nanos", "Picoeukaryotes" = "peuks", 
                                          "Cyanobacteria" = "syn", "Bacteria" = "bacteria"))

counts_summary$gate <- factor(counts_summary$gate, levels = c("Bacteria", "Cyanobacteria", "Picoeukaryotes", "Nanoeukaryotes", 
                                                            "Microeukaryotes", "Coccolithophores"), ordered = TRUE)
counts_summary$domain <- factor(counts_summary$domain, levels = c("Prokaryotes", "Eukaryotes"), ordered = TRUE)

counts_summary <- na.omit(counts_summary)# Remove rows with NA values

options(scipen = 999) # this stops R printing scientific notation on the graphs

counts_prok <- ggplot(counts_summary[counts_summary$domain == "Prokaryotes", ], aes(x = treatment, y = mean_counts, color = gate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.2) +
  geom_hline(yintercept = 62, linetype = "dashed", color = "#D51317FF") +
  geom_hline(yintercept = 4.15, linetype = "dashed", color = "#F39200FF") +
  scale_color_frontiers() +
  scale_x_discrete(breaks = unique(counts_summary$treatment)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  labs(x = "Treatment", y = expression(paste("10"^3,"cells mL"^-1))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.grid.minor = element_blank(),  axis.title.x = element_blank(),
    legend.title = element_blank(),legend.direction = "horizontal",legend.position = "none",
    plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, 
                                        linewidth = 0.25, linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  ggtitle("Prokaryotes") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

counts_euk <- ggplot(counts_summary[counts_summary$domain == "Eukaryotes", ], aes(x = treatment, y = mean_counts, color = gate)) +
  geom_hline(yintercept = 10.88, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 5.63, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 0.13, linetype = "dashed", color = "#007B3DFF") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_counts - se_counts, ymax = mean_counts + se_counts), width = 0.2) +
  scale_color_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_x_discrete(breaks = unique(counts_summary$treatment)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  labs(x = "Treatment", y = expression(paste("10"^3,"cells mL"^-1))) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  ylab("") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y.left = element_blank(),
    panel.grid.minor = element_blank(),  axis.title.x = element_blank(),
    legend.title = element_blank(),legend.direction = "horizontal",legend.position = "none",
    plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, 
                                        linewidth = 0.25, linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  ggtitle("Eukaryotes") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

#---adding asterisks for significance ----------------------------------------------------------------
asterisks_prok <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928), 
  y_pos = c(10.64, 11.00, 10.71, 8.79, 9.26, 8.93), 
  label = c("a", "b", "a", "a", "b", "a"), 
  colours = c("#D51317", "#D51317", "#D51317", "#F39200", "#F39200", "#F39200"))

asterisks_euk <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), 
  y_pos = c(8.21, 8.57, 8.35, 9.97, 9.85, 10.63, 5.93, 5.45, 6.36, 6.87, 6.81, 7.22), 
  label = c("a", "b", "ab", "a", "a", "b", "a", "b", "c", "a", "a", "b"), 
  colours = c("#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D", "#007B3D","#31B7BC","#31B7BC","#31B7BC"))

# -- doing up the visreg plots ---------------------------------------------------------------------------------
visreg_counts <- visreg(counts_model, "treatment", by="gate", data = cells, plot = FALSE) 
visreg_prok <- subset(visreg_counts, gate %in% c("bacteria", "syn"))
visreg_euk <- subset(visreg_counts, gate %in% c("peuks", "nanos", "micro", "Cocco"))
visreg_prok$res$gate <- factor(visreg_prok$res$gate, levels = c("bacteria", "syn"))
visreg_prok$fit$gate <- factor(visreg_prok$fit$gate, levels = c("bacteria", "syn"))
visreg_euk$res$gate <- factor(visreg_euk$res$gate, levels = c("peuks", "nanos", "micro", "Cocco"))
visreg_euk$fit$gate <- factor(visreg_euk$fit$gate, levels = c("peuks", "nanos", "micro", "Cocco"))
  
visreg_counts_prok <- plot(visreg_prok, overlay = TRUE, gg = TRUE, points = list(shape = 1, size =2),
     scale = "linear", partial = TRUE, top = 'points', data = cells) + 
  theme_bw() + 
  ylab(expression(paste("ln(cells mL"^-1,")"))) +  
  xlab("Treatment") +
  scale_y_continuous(limits = c(4.9, 13.1), breaks = seq(5, 13, by = 2)) +
  scale_colour_manual(values = c("#D51317FF","#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_fill_manual(values = c("#D5131733","#F3920033", "#EFD50033","#95C11F33","#007B3D33","#31B7BC33")) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "none",
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, 
                                            linewidth = 0.25, linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  geom_text(data = asterisks_prok, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_prok$colours, size = 4, inherit.aes = FALSE)  +
  ggtitle("Prokaryotes") +
  guides(element_blank()) 

visreg_counts_euk <- plot(visreg_euk, overlay = TRUE, gg = TRUE, points = list(shape = 1, size = 2),
     scale = "linear", partial = TRUE, top = 'points', date = cells) + 
  theme_bw() + 
  ylab("") +  
  xlab("Treatment") +
  scale_y_continuous(limits = c(4.9, 13.1), breaks = seq(5, 13, by = 2)) +
  scale_colour_manual(values = c("#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_fill_manual(values = c("#EFD50033","#95C11F33","#007B3D33","#31B7BC33")) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "none", 
        plot.title = element_textbox_simple(size = 9, color = "black", fill = "lightgrey", halign = 0.5,lineheight = 0.9, 
                                            linewidth = 0.25, linetype = 1, padding = margin(4.5, 3, 4.5, 3), box.color = "black")) +
  geom_text(data = asterisks_euk, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_euk$colours, size = 4, inherit.aes = FALSE) +
  ggtitle("Eukaryotes") +
  guides(element_blank()) 

#--- putting all of the plots together ----------------------------------------------------------------------------------
combined_counts_plots <- guide_area() / (counts_prok |counts_euk) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top", 
        legend.spacing.x = unit(0.1, "cm"), legend.box = "vertical",
        legend.spacing.y = unit(-0.3, "cm"), 
        plot.tag.position = c(0.01, 0.97)) 

SOTS_FC_counts_plots <- combined_counts_plots / (visreg_counts_prok | visreg_counts_euk) +
  plot_layout(nrow = (3), heights = c(2,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        plot.tag.position = c(0.01, 0.97))

ggsave("SOTS_FC_counts_plot2.svg", SOTS_FC_counts_plots, width = 7.2, height = 7.75) 
