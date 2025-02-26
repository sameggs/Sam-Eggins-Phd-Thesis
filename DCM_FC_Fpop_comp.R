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
library(compositions)

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
cells$treatment <- factor(cells$treatment, levels = c("85", "+DFB", "Control", "+Fe"), ordered = TRUE)
# summing nanophytoplankton  
nano_gates <- c("nano 1", "nano 2", "nano 3")
cells <- cells %>%
  group_by(light, treatment, bottle, time, gate = ifelse(gate %in% nano_gates, "nanos", gate)) %>%
  summarise(across(c(counts, Fsize, Fsize_bv, FV12, FB3), ~sum(.x, na.rm = TRUE)), .groups = "drop")

# grab data from initial conditions
time_0 <- cells %>% filter(treatment == "85") %>% mutate_at(vars(counts:FB3), as.numeric)
cells <- cells %>% filter(time == 10)

#Check the structure
str(cells)

# ----- Investigating the density distributions of the data ------------------------------------------------------
# List of unique gate variables
gate_vars <- unique(cells$gate)

# List of count variables
Fsize_vars <- c("Fsize_bv")

# Generate a list of plots
plots_list <- purrr::map(gate_vars, ~ {
  gate_var <- .x
  purrr::map(Fsize_vars, ~ {
    F_var <- .x
    data_subset <- cells[cells$gate == gate_var, ]
    ggplot(data_subset, aes(x = !!sym(F_var))) + 
      geom_density(alpha = .2, fill = "#FF6FF6") +
      ggtitle(paste(gate_var, "-", F_var))
  })
}) %>% purrr::flatten()

# Combine all plots in a grid
density_plots <- wrap_plots(plots_list, ncol = length(gate_vars))  # adjust ncol and nrow as needed
print(density_plots)

# FOR THE WHOLE COUNT SET OF DATA .................................
filtered_cells <- cells %>%
  filter(!(gate %in% c("all_phyto", "bacteria")))

# Generate a list of plots
plots_list <- purrr::map(Fsize_vars, ~ {
  F_var <- .x
  ggplot(filtered_cells, aes(x = !!sym(F_var))) + 
    geom_density(alpha = .2, fill = "#FF6FF6") +
    ggtitle(paste("Density of", F_var))
})

# Combine all plots in a grid
density_plot_all <- wrap_plots(plots_list, ncol = 1)  # adjust ncol as needed
print(density_plot_all)

#--------- Reshape data ----------------------------------------------------------------------------------------
Fpop <- cells %>% 
  select(light, treatment, bottle, gate, Fsize_bv) %>%
  pivot_wider(names_from = gate, values_from = Fsize_bv)
str(Fpop)
Fpop <- Fpop %>% 
  select(-bacteria, -all_phyto)

# Calculate the sum of known gates
Fpop$known_sum <- rowSums(Fpop[, c("Cocco", "micro", "nanos", "peuks", "syn")], na.rm = TRUE)
# Calculate the 'other' category as 1 minus this sum
Fpop$other <- 1 - Fpop$known_sum
# You can then remove the 'known_sum' column if you don't need it
Fpop$known_sum <- NULL

#------FITTING MODELS composition:  -----------------------------------------------------------------------------------
# fitting a multivariate model for Fpop
gates <- Fpop[, c("Cocco", "micro", "nanos", "peuks", "syn", "other")]
Fpop_model <- lm(formula = clr(gates) ~ light * treatment, data = Fpop)

#------individual models still with centred log ratio transform-----------------------------------------------------------
Fpop_clr <- clr(Fpop[,c("Cocco", "micro", "nanos", "peuks", "syn", "other")])
Fpop_clr <- cbind(Fpop[, !names(Fpop) %in% c("Cocco", "micro", "nanos", "peuks", "syn", "other")], Fpop_clr)

model_Cocco <- lm(Fpop_clr$Cocco ~ light * treatment, data = Fpop_clr)
model_micro <- lm(Fpop_clr$micro ~ light * treatment, data = Fpop_clr)
model_nanos <- lm(Fpop_clr$nanos ~ light * treatment, data = Fpop_clr)
model_peuks <- lm(Fpop_clr$peuks ~ light * treatment, data = Fpop_clr)
model_syn <- lm(Fpop_clr$syn ~ light * treatment, data = Fpop_clr)

visreg(model_peuks, "treatment", by="light", overlay = TRUE, data = Fpop_clr,
       gg = TRUE, points = list(shape = 19), 
       scale = "linear", partial = TRUE) 

par(mfrow = c(3, 2))
plot(model_peuks, which = 1:6)
plot(model_nanos, which = 1:6)
plot(model_micro, which = 1:6)
plot(model_Cocco, which = 1:6)

# Convert the clr-transformed data to a data frame
Fpop_clr_df <- as.data.frame(Fpop_clr)

# Melt the data frame to long format for easier plotting with ggplot2
Fpop_clr_long <- reshape2::melt(Fpop_clr_df)

# Create a density plot for each clr-transformed variable
ggplot(Fpop_clr_long, aes(value)) + 
  geom_density(alpha = .2, fill = "#FF6FF6") + 
  facet_wrap(~variable, scales = "free") + 
  theme_bw() +
  xlab("clr-transformed value") +
  ylab("Frequency")

# ---- Model diagnostics for the response variables in the multivariate model -------------------
# Create a list to store the diagnostic grids
diagnostic_grids <- list()

# Loop through each response variable
for (response_var in colnames(gates)) {
  # Extract residuals for the current response variable
  residuals_var <- residuals(Fpop_model)[, response_var]
  
  # Fitted values for the current response variable
  fitted_values_var <- fitted(Fpop_model)[, response_var]
  
  # Calculate leverage values for the current response variable
  X <- model.matrix(Fpop_model)
  leverage_var <- diag(X %*% solve(t(X) %*% X) %*% t(X))
  
  # Create the diagnostic plots for the current response variable
  plots <- list()
  
  # Residuals vs. Fitted
  plot1 <- ggplot(data.frame(Fitted = fitted_values_var, Residuals = residuals_var),
                  aes(x = Fitted, y = Residuals)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    labs(x = "Fitted Values", y = "Residuals") +
    ggtitle(paste("Residuals vs Fitted for", response_var))
  
  # Scale-Location plot
  plot2 <- ggplot(data.frame(Fitted = fitted_values_var,
                             Sqrt_Residuals = sqrt(abs(residuals_var))),
                  aes(x = Fitted, y = Sqrt_Residuals)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    labs(x = "Fitted Values", y = "Square Root of Residuals") +
    ggtitle(paste("Scale-Location plot for", response_var))
  
  # QQ plot
  plot3 <- ggplot(data.frame(Standardized_Residuals = residuals_var),
                  aes(sample = Standardized_Residuals)) +
    stat_qq() +
    stat_qq_line() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    ggtitle(paste("QQ Plot of Residuals for", response_var))
  
  # Create the diagnostic grid for the current response variable
  grid <- wrap_plots(plot1, plot2, plot3)
  
  # Add the diagnostic grid for the current response variable to the main list
  diagnostic_grids[[response_var]] <- grid
}

# Print the diagnostic grids for each response variable
for (response_var in colnames(gates)) {
  grid <- diagnostic_grids[[response_var]]
  print(grid)
}

# --------------------------emmeans for the multivariate model -----------------------------------------------------
Fpop_emmeans <- summary(emmeans(Fpop_model, specs = pairwise ~ treatment:light |rep.meas))
Fpop_emmeans_TRTinL <- emmeans(Fpop_model, ~ treatment|light |rep.meas)
TRTinL_contrasts <- (mvcontrast(Fpop_emmeans_TRTinL, method = "pairwise", show.ests = TRUE))
Fpop_emmeans_LinTRT <- emmeans(Fpop_model, ~ light|treatment |rep.meas)
LinTRT_contrasts <- (mvcontrast(Fpop_emmeans_LinTRT, method = "pairwise", show.ests = TRUE))

Fpop_contrasts <- bind_rows(emmeans = Fpop_emmeans$emmeans, 
                          light = as.data.frame(TRTinL_contrasts$estimates), 
                          trt = as.data.frame(LinTRT_contrasts$estimates))

# Export the data frame to an Excel file
write_xlsx(Fpop_contrasts, "Fpop_emmeans_multi.xlsx")


#---- Plots with just data on them ---------------------------------------------------------
# select columns that are required
vars <- c("peuks", "nanos", "micro", "Cocco", "syn")
#create summary
Fpop_summary <- Fpop %>%
  group_by(light, treatment) %>%
  dplyr::summarise(across(all_of(vars), 
                          list(mean = ~mean(.x * 100, na.rm = TRUE), 
                               se = ~sd(.x * 100, na.rm = TRUE) / sqrt(sum(!is.na(.x))))))


# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, 
                                          show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, legend_position = "top",
                                          hline_y) {
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

plot_peuks <- generate_plots_with_errorbars(Fpop_summary, "peuks_mean", "peuks_se", "peuks", show_xlab = FALSE, y_min = 0, y_max = 1.2,
                                            y_lab = expression(paste("Fpop"[pico], " (%)")), legend_position = "none", hline_y = 0.21) + 
                                            scale_y_continuous(breaks = seq(0, 1.2, by = 0.4))
plot_nanos <- generate_plots_with_errorbars(Fpop_summary, "nanos_mean", "nanos_se", "nanos", show_xlab = FALSE, y_min = 0, y_max = 40,
                                            y_lab = expression(paste("Fpop"[nano], " (%)")), legend_position = "none", hline_y = 11.45) + 
                                            scale_y_continuous(breaks = seq(0, 40, by = 10))
plot_micro <- generate_plots_with_errorbars(Fpop_summary, "micro_mean", "micro_se", "micro", show_xlab = FALSE, y_min = 0, y_max = 60,
                                            y_lab = expression(paste("Fpop"[micro], " (%)")), legend_position = "none", hline_y = 2.79) + 
                                            scale_y_continuous(breaks = seq(0, 60, by = 15))
plot_cocco <- generate_plots_with_errorbars(Fpop_summary, "Cocco_mean", "Cocco_se", "Cocco", show_xlab = FALSE, y_min = 0, y_max = 15,
                                            y_lab = expression(paste("Fpop"[cocco], " (%)")), legend_position = "none", hline_y = 10.99) + 
                                            scale_y_continuous(breaks = seq(0, 15, by = 5))
plot_cyano <- generate_plots_with_errorbars(Fpop_summary, "syn_mean", "syn_se", "syn", show_xlab = FALSE, y_min = 0, y_max = 0.1,
                                            y_lab = expression(paste("Fpop (%)")), legend_position = "none", hline_y = 0.03) + 
                                            scale_y_continuous(breaks = seq(0, 0.1, by = 0.02))

#---adding asterisks for significance ----------------------------------------------------------------
asterisks_peuks <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.38), y_pos = c(-1.34, -1.19, -2.12, -2.52, -2.74, -3.5, -1.19), 
  label = c("a", "a", "b", "a", "a", "b", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_nano <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.022, 0.738), y_pos = c(2.46, 2.34, 2.84, 1.37, 1.49, 1.56, 2.46, 2.84), 
  label = c("a", "a", "b", "a", "a", "a","\u263C", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black", "black"))

asterisks_micro <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.738), y_pos = c(0.22, 1.7, 3.46, 1.61, 0.31, 0.83, 3.46), 
  label = c("a", "a", "b", "a", "a", "a", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

asterisks_cocco <- data.frame(
  x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.022), y_pos = c(1.44, 0.1, -0.5, -0.14, 1.23, 0.8, 1.44), 
  label = c("a", "a", "b", "a", "a", "a", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))

# -----------------the code to generate each visreg plot ---------------------------------------------------------------------
visreg_Fpop <- visreg(Fpop_model, "treatment", by="light", overlay = TRUE, data = Fpop,
       gg = TRUE, points = list(shape = 19), scale = "linear", partial = TRUE) 

visreg_peuks <-  visreg_Fpop[[4]] +
                  theme_bw() + 
                  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
                  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
                  ylab("") +  
                  xlab("") +
                  scale_y_continuous(limits = c(-4, -1), breaks = seq(-4, -1, by = 1),
                                     sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop[~"pico"])))) +
                  theme(axis.text = element_text(size = 10, colour = "black"),
                        panel.grid.minor = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
                        axis.text.y.right = element_text(size = 10, colour = "black"),
                        legend.title = element_blank(), legend.text = element_text(size = 9), 
                        legend.direction = "horizontal", legend.position = "top") +
                  geom_text(data = asterisks_peuks, aes(x = x_pos, y = y_pos, label = label), 
                            colour = asterisks_peuks$colours, size = 4, inherit.aes = FALSE)  +
                  guides(fill = "none")

visreg_nanos <-  visreg_Fpop[[3]] +
                  theme_bw() + 
                  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
                  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
                  ylab("") +  
                  xlab("") +
                  scale_y_continuous(limits = c(1, 3), breaks = seq(1, 3, by = 0.5),
                                     sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop[~"nano"])))) +
                  theme(axis.text = element_text(size = 10, colour = "black"),
                        panel.grid.minor = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
                        axis.text.y.right = element_text(size = 10, colour = "black"),
                        legend.title = element_blank(), legend.text = element_text(size = 9), 
                        legend.direction = "horizontal", legend.position = "top") +
                  geom_text(data = asterisks_nano, aes(x = x_pos, y = y_pos, label = label), 
                            colour = asterisks_nano$colours, size = 4, inherit.aes = FALSE)  +
                  guides(fill = "none")

visreg_micro <-  visreg_Fpop[[2]] +
                  theme_bw() + 
                  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
                  scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
                  ylab("") +  
                  xlab("") +
                  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, by = 1),
                                     sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop[~"micro"])))) +
                  theme(axis.text = element_text(size = 10, colour = "black"),
                        panel.grid.minor = element_blank(),
                        axis.title.x = element_blank(),
                        axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
                        axis.text.y.right = element_text(size = 10, colour = "black"),
                        legend.title = element_blank(), legend.text = element_text(size = 9), 
                        legend.direction = "horizontal", legend.position = "top") + 
                  geom_text(data = asterisks_micro, aes(x = x_pos, y = y_pos, label = label), 
                            colour = asterisks_micro$colours, size = 4, inherit.aes = FALSE)  +
                  guides(fill = "none")

visreg_cocco <-  visreg_Fpop[[1]] +
                        theme_bw() + 
                        scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
                        scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
                        ylab("") +  
                        xlab("") +
                        scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, by = 1),
                                           sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop[~"cocco"])))) +
                        theme(axis.text = element_text(size = 10, colour = "black"),
                              panel.grid.minor = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
                              axis.text.y.right = element_text(size = 10, colour = "black"),
                              legend.title = element_blank(), legend.text = element_text(size = 9), 
                              legend.direction = "horizontal", legend.position = "top") +
                        geom_text(data = asterisks_cocco, aes(x = x_pos, y = y_pos, label = label), 
                                  colour = asterisks_cocco$colours, size = 4, inherit.aes = FALSE)  +
                        guides(fill = "none")

visreg_cyano <-  visreg_Fpop[[5]] +
                        theme_bw() + 
                        scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
                        scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
                        ylab("") +  
                        xlab("") +
                        scale_y_continuous(limits = c(-5, -3), breaks = seq(-5, -3, by = 0.5),
                                            sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop[~"cyano"])))) +
                        theme(axis.text = element_text(size = 10, colour = "black"),
                              panel.grid.minor = element_blank(),
                              axis.title.x = element_blank(),
                              axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
                              axis.text.y.right = element_text(size = 10, colour = "black"),
                              legend.title = element_blank(), legend.text = element_text(size = 9), 
                              legend.direction = "horizontal", legend.position = "top") +
                        guides(fill = "none")

# -----------putting plots together------------------------------------------------------------------------------------------------------
Fpop_model_plot <- guide_area() / (plot_peuks | visreg_peuks) / (plot_nanos | visreg_nanos)/ 
  (plot_micro | visreg_micro) / (plot_cocco | visreg_cocco) +
  plot_layout(guides = 'collect', nrow = (5), heights = c(1,10,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C','D','E','F','G','H'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = "none")

ggsave("Fpop_model_plot.svg", Fpop_model_plot, width = 7.2, height = 9.8) 

#---- Plot with all of the gates on it ---------------------------------------------------------
# summarise the Fsize from cells df
Fsize_summary <- cells %>%
  dplyr::filter(!(gate %in% c("all_phyto", "bacteria"))) %>%
  dplyr::group_by(light, treatment, gate) %>%
  dplyr::summarize(
    mean_Fpop = mean(Fsize_bv, na.rm = TRUE) * 100,
    se_Fpop = (sd(Fsize_bv, na.rm = TRUE) / sqrt(n())) * 100) %>%
  mutate(light = as.factor(light), 
         treatment = as.factor(treatment),
         gate = as.factor(gate) %>% fct_recode("Coccolithophores" = "Cocco", "Microeukaryotes" = "micro",
                                               "Nanoeukaryotes" = "nanos", "Picoeukaryotes" = "peuks", "Cyanobacteria" = "syn"))

Fsize_summary$gate <- factor(Fsize_summary$gate, levels = c("Cyanobacteria","Picoeukaryotes", "Nanoeukaryotes", 
                                                            "Microeukaryotes", "Coccolithophores"), ordered = TRUE)
Fsize_summary$light <- factor(Fsize_summary$light, levels = c("Low Light", "High Light"), ordered = TRUE)

# -----Create the ggplot object--- FOR A BAR CHART-------------------------------------------------------------------------------------------
Fpop_plot <- ggplot(Fsize_summary, aes(x = treatment, y = mean_Fpop, fill = gate)) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "#F39200FF") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 28.5, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 25.2, linetype = "dashed", color = "#007B3DFF") +
  geom_hline(yintercept = 3.6, linetype = "dashed", color = "#31B7BCFF") +
  geom_errorbar(aes(ymin = mean_Fpop - se_Fpop, ymax = mean_Fpop + se_Fpop), width = 0.2) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_x_discrete(breaks = unique(Fsize_summary$treatment)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 80, by = 10)) +
  labs(x = "Treatment", y = "Fpop (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),legend.direction = "horizontal",legend.position = "none") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  facet_wrap(~light)


Fpop_bar <- ggplot(Fsize_summary, aes(x = treatment, y = mean_Fpop, fill = gate)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_x_discrete(breaks = unique(Fsize_summary$treatment)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Treatment", y = expression(paste("Fpop (%)"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"), 
        panel.grid.minor = element_blank(), legend.title = element_blank(),
        legend.direction = "horizontal", legend.position = "top") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  facet_wrap(~light)

# ---- PCA analysis --------------------------------------------------------------------------------------------
gates <- Fpop[, c("syn", "peuks", "nanos", "micro", "Cocco")]
Fpop$clr_transformed <- clr(gates)

pca_result <- prcomp(Fpop$clr_transformed, center = TRUE, scale. = TRUE)
scores <- as.data.frame(pca_result$x)
loadings <- as.data.frame(pca_result$rotation)
biplot(pca_result, scale = 0)
rownames(loadings) <- c("syn" = "Cyanobacteria", 
                        "peuks" = "Picoeukaryotes", 
                        "nanos" = "Nanoeukaryotes", 
                        "micro" = "Microeukaryotes", 
                        "Cocco" = "Coccolithophores")[rownames(loadings)]

# Plot scores
biplot_gg <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Fpop$treatment, colour = Fpop$light), size = 5) + 
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2),
               arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "black") +
  scale_x_continuous(limits = c(-4, 2), breaks = seq(-4, 2, by = 1)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  geom_text(data = loadings, aes(x = PC1*2, y = PC2*2, label = rownames(loadings)), 
            vjust = -0.5, color = "black") +
  labs(x = paste("PC1:", round(100 * summary(pca_result)$importance[2, 1], 1), "%"), 
       y = paste("PC2:", round(100 * summary(pca_result)$importance[2, 2], 1), "%")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))
biplot_gg2 <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Fpop$treatment, colour = Fpop$light), size = 5) + 
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2),
               arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "black") +
  scale_x_continuous(limits = c(-4, 2), breaks = seq(-4, 2, by = 1)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  geom_text(data = loadings, aes(x = PC1*2, y = PC2*2, label = rownames(loadings)), 
            vjust = -0.5, color = "black") +
  labs(x = paste("PC1:", round(100 * summary(pca_result)$importance[2, 1], 1), "%"), 
       y = paste("PC2:", round(100 * summary(pca_result)$importance[2, 2], 1), "%")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"), axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))


# plotting the biplot and results together....
Fpop_model_plot <- guide_area() / Fpop_plot / (biplot_gg | biplot_gg2) +
  plot_layout(guides = 'collect', nrow = 3, heights = c(1,10.5,11)) +
  plot_annotation(tag_levels = list(c('a','b'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 9),
        legend.direction = "horizontal", 
        legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) &
  guides(fill = "none")

Fpop_model_plot

ggsave("Fpop_plot_DCM.svg", Fpop_model_plot, width = 7.2, height = 7) 
