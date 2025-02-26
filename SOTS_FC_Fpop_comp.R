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
library(compositions)

#------load the flow cytometry data  into R---------------------------------------------------------------------------------------
cells <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/sots_flowcyt.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
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
  group_by(treatment, bottle, time, gate = ifelse(gate %in% nano_gates, "nanos", gate)) %>%
  summarise(across(c(counts, Fsize, Fsize_bv, FV12, FB3), ~sum(.x, na.rm = TRUE)), .groups = "drop")

# grab data from initial conditions
time_0 <- cells %>% filter(time == "0") %>% mutate_at(vars(counts:FB3), as.numeric)
cells <- cells %>% filter(time == 5)

#Check the structure
str(cells)

#-------------- create a normalised version of the Fsize_bv data ---------------------------------------------------------------
# Get the normalization factor
norm_factor <- cells %>%
  filter(gate == "all_phyto") %>%
  select(treatment, bottle, time, norm_Fsize_bv = Fsize_bv)
# Join the normalization factor to the original data
cells_norm <- cells %>%
  left_join(norm_factor, by = c("treatment", "bottle", "time"))
# Create normalized Fsize_bv
cells_norm <- cells_norm %>%
  mutate(Fsize_bv = Fsize_bv / norm_Fsize_bv)
# Select desired columns
cells_norm <- cells_norm %>%
  filter(gate != "all_phyto") %>%
  select(treatment, bottle, time, gate, Fsize_bv)

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
density_plots <- wrap_plots(plots_list, ncol = length(Fsize_vars))  # adjust ncol and nrow as needed
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
  select(treatment, bottle, gate, Fsize_bv) %>%
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

# also a version where these are normalised to 1
Fpop <- Fpop %>%
  rowwise() %>%
  mutate(sum_val = sum(Cocco, micro, nanos, peuks, syn)) %>%
  mutate(
    Cocco = Cocco / sum_val,
    micro = micro / sum_val,
    nanos = nanos / sum_val,
    peuks = peuks / sum_val,
    syn = syn / sum_val
  ) %>%
  select(-other, -sum_val)  # Remove the 'other' and temporary 'sum_val' columns


#------FITTING MODELS composition:  -----------------------------------------------------------------------------------
# #running a MANOVA
gates <- Fpop[, c("syn", "peuks", "nanos", "micro", "Cocco")]
Fpop$clr_transformed <- clr(gates)

Fpop_manova_clr <- manova(Fpop$clr_transformed[,2:5] ~ treatment, data = Fpop) # drop a term for d.f. in Wilk's test
summary(Fpop_manova_clr, test = "Wilks")
Fpop_manova_clr <- manova(Fpop$clr_transformed ~ treatment, data = Fpop)

# --------------------------emmeans for the multivariate model -----------------------------------------------------
Fpop_emmeans <- summary(emmeans(Fpop_manova_clr, specs = pairwise ~ treatment |rep.meas))
Fpop_emmeans_TRT <- emmeans(Fpop_manova_clr, ~ treatment |rep.meas)
TRT_contrasts <- (mvcontrast(Fpop_emmeans_TRT, method = "pairwise", adjust = "sidak", show.ests = TRUE))

Fpop_contrasts <- bind_rows(emmeans = Fpop_emmeans$emmeans, 
                            trt = as.data.frame(TRT_contrasts$estimates))

# Export the data frame to an Excel file
write_xlsx(Fpop_contrasts, "SOTS_Fpop_emmeans_multi.xlsx")

# ---- PCA analysis --------------------------------------------------------------------------------------------
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
  geom_point(aes(shape = Fpop$treatment), size = 5) + 
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*2, yend = PC2*2),
               arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#D51317FF") +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  geom_text(data = loadings, aes(x = PC1*2, y = PC2*2, label = rownames(loadings)), 
            vjust = -0.5, color = "#D51317FF") +
  labs(x = paste("PC1:", round(100 * summary(pca_result)$importance[2, 1], 1), "%"), 
       y = paste("PC2:", round(100 * summary(pca_result)$importance[2, 2], 1), "%")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))

ggsave("Fpop_biplot.svg", biplot_gg, width = 5.6, height = 4.6)

# ---- Model diagnostics for the response variables in the multivariate model -------------------
# Create a list to store the diagnostic grids
diagnostic_grids <- list()

# Loop through each response variable
for (response_var in colnames(gates)) {
  # Extract residuals for the current response variable
  residuals_var <- residuals(Fpop_manova_clr)[, response_var]
  
  # Fitted values for the current response variable
  fitted_values_var <- fitted(Fpop_manova_clr)[, response_var]
  
  # Calculate leverage values for the current response variable
  X <- model.matrix(Fpop_manova_clr)
  leverage_var <- diag(X %*% solve(t(X) %*% X) %*% t(X))
  
  # Create the diagnostic plots for the current response variable
  plots <- list()
  
  # Residuals vs. Fitted
  plot1 <- ggplot(data.frame(Fitted = fitted_values_var, Residuals = residuals_var),
                  aes(x = Fitted, y = Residuals)) +
    geom_point() +
    labs(x = "Fitted Values", y = "Residuals") +
    ggtitle(paste("Residuals vs Fitted for", response_var))
  
  # Scale-Location plot
  plot2 <- ggplot(data.frame(Fitted = fitted_values_var,
                             Sqrt_Residuals = sqrt(abs(residuals_var))),
                  aes(x = Fitted, y = Sqrt_Residuals)) +
    geom_point() +
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

# ------ plotting up the Fpop data ----------------------------------------------------------------------------------------
# summarise the Fsize from cells df
Fsize_summary <- cells_norm %>%
  dplyr::filter(!(gate %in% c("all_phyto", "bacteria"))) %>%
  dplyr::group_by(treatment, gate) %>%
  dplyr::summarize(
    mean_Fsize = mean(Fsize_bv*100, na.rm = TRUE),
    se_Fsize = (sd(Fsize_bv*100, na.rm = TRUE) / sqrt(n()))) %>%
  mutate(treatment = as.factor(treatment),
    gate = as.factor(gate) %>% fct_recode("Coccolithophores" = "Cocco", "Microeukaryotes" = "micro",
                                          "Nanoeukaryotes" = "nanos", "Picoeukaryotes" = "peuks", 
                                          "Cyanobacteria" = "syn"))

Fsize_summary$gate <- factor(Fsize_summary$gate, levels = c("Cyanobacteria", "Picoeukaryotes", "Nanoeukaryotes", 
                                                              "Microeukaryotes", "Coccolithophores"), ordered = TRUE)
Fsize_summary <- na.omit(Fsize_summary)# Remove rows with NA values

Fpop_plot <- ggplot(Fsize_summary, aes(x = treatment, y = mean_Fsize, fill = gate)) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#F39200FF") +
  geom_hline(yintercept = 6.4, linetype = "dashed", color = "#EFD500FF") +
  geom_hline(yintercept = 59.9, linetype = "dashed", color = "#95C11FFF") +
  geom_hline(yintercept = 12.3, linetype = "dashed", color = "#007B3DFF") +
  geom_errorbar(aes(ymin = mean_Fsize - se_Fsize, ymax = mean_Fsize + se_Fsize), width = 0.2) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_x_discrete(breaks = unique(Fsize_summary$treatment)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)) +
  labs(x = "Treatment", y = "Fpop (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
    panel.grid.minor = element_blank(), 
    legend.title = element_blank(),legend.direction = "horizontal",legend.position = "none") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

#---adding asterisks for significance ----------------------------------------------------------------
asterisks <- data.frame(
  x_pos = c(0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928, 0.212, 0.57, 0.928), 
  y_pos = c(-2.77, -2.28, -2.51, -0.87, -0.509, -1.18, 2.64, 2.36, 2.77, 1.78, 1.31, 1.705, 0.72, 0.61, 0.72), 
  label = c("a", "b", "c", "a", "b", "c", "a", "b", "a", "a", "b", "a", "a", "a", "a"), 
  colours = c("#F39200", "#F39200", "#F39200","#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D", "#007B3D","#31B7BC","#31B7BC","#31B7BC"))

# --- plotting VISREG model objects ------------------------------------------------------------------------------------
visreg_Fpop <- visreg(Fpop_manova_clr, "treatment", data = Fpop, plot = FALSE) 
# combining the visreg objects for each gate
visreg_Fpop_all <- visregList(visreg_Fpop[[1]],visreg_Fpop[[2]],visreg_Fpop[[3]],visreg_Fpop[[4]],visreg_Fpop[[5]], 
           collapse = TRUE, labels=c("Cyanobacteria", "Picoeukaryotes", "Nanoeukaryotes", "Microeukaryotes", "Coccolithophores"))
#defining the gate factor.
visreg_Fpop_all$res$gate <- as.factor(visreg_Fpop_all$res$visregCollapse) 
visreg_Fpop_all$fit$gate <- as.factor(visreg_Fpop_all$fit$visregCollapse) 
visreg_Fpop_all$res$gate <- factor(visreg_Fpop_all$res$gate, levels = c("Cyanobacteria", "Picoeukaryotes", "Nanoeukaryotes", 
                                                                        "Microeukaryotes", "Coccolithophores"), order = TRUE)
visreg_Fpop_all$fit$gate <- factor(visreg_Fpop_all$fit$gate, levels = c("Cyanobacteria", "Picoeukaryotes", "Nanoeukaryotes", 
                                                                        "Microeukaryotes", "Coccolithophores"), order = TRUE)

visreg_Fpop <- plot(visreg_Fpop_all, by="gate", gg = TRUE, scale = "linear", partial = TRUE, top = 'points', overlay = TRUE) + 
  theme_bw() + 
  ylab(expression(paste("ln(cells mL"^-1,")"))) +  
  ylab("")+
  xlab("Treatment") +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1),
                     sec.axis = sec_axis(~., name = expression(italic(clr)(Fpop)), breaks = seq(-4, 4, by = 1))) +
  scale_colour_manual(values = c("#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_fill_manual(values = c("#F3920033", "#EFD50033","#95C11F33","#007B3D33","#31B7BC33")) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
        axis.text.y.right = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  geom_text(data = asterisks, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks$colours, size = 4, inherit.aes = FALSE)  +
  guides(element_blank()) 

# -----------putting plots together------------------------------------------------------------------------------------------------------
Fpop_model_plot <-  guide_area() / (Fpop_plot | visreg_Fpop) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(1,10)) +
  plot_layout(nrow = (2), heights = c(0.5,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top", 
        legend.spacing.x = unit(0.1, "cm"), legend.box = "vertical",
        legend.spacing.y = unit(-0.3, "cm"), 
        plot.tag.position = c(0.01, 0.97)) &
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("SOTS_Fpop_model_plot.svg", Fpop_model_plot, width = 7.2, height = 3.97) 

# -----------Fpop bar plot?? ------------------------------------------------------------------------------------------------------
Fpop_bar_plot <- ggplot(Fsize_summary, aes(x = treatment, y = mean_Fsize, fill = gate)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Cyanobacteria"    = "#F39200FF",
                               "Picoeukaryotes"   = "#EFD500FF",
                               "Nanoeukaryotes"   = "#95C11FFF",
                               "Microeukaryotes"  = "#007B3DFF",
                               "Coccolithophores" = "#31B7BCFF")) +
  scale_x_discrete(breaks = unique(Fsize_summary$treatment)) +
  scale_y_continuous(limits = c(-0.01, 100.01), breaks = seq(0, 100, by = 20)) + 
  labs(x = "Treatment", y = "Fpop (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "top") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

Fpop_model_plot_bar <-  guide_area() / (Fpop_bar_plot | visreg_Fpop) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(1,10)) +
  plot_layout(nrow = (2), heights = c(0.5,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top", 
        legend.spacing.x = unit(0.1, "cm"), legend.box = "vertical",
        legend.spacing.y = unit(-0.3, "cm"), 
        plot.tag.position = c(0.01, 0.97)) &
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("SOTS_Fpop_model_plot_bar.svg", Fpop_model_plot_bar, width = 7.2, height = 3.97)

#---- plot with Fpop data and PCR... -------------------------------------------------
asterisks <- data.frame(
  x_pos = c(1.3, 2.3, 3.3, 1.3, 2.3, 3.3, 1.3, 2.3, 3.3, 1.3, 2.3, 3.3, 1.3, 2.3, 3.3), 
  y_pos = c(0, 0, 0, 2, 3, 2, 62, 62, 65, 25, 22, 22, 10, 12, 10),  
  label = c("a", "b", "c", "a", "b", "c", "a", "b", "a", "a", "b", "a", "a", "a", "a"), 
  colours = c("#F39200", "#F39200", "#F39200","#EFD500", "#EFD500", "#EFD500", "#95C11F", "#95C11F", "#95C11F", "#007B3D", "#007B3D", "#007B3D","#31B7BC","#31B7BC","#31B7BC"))

Fpop_plot <- ggplot(Fsize_summary, aes(x = treatment, y = mean_Fsize, color = gate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_Fsize - se_Fsize, ymax = mean_Fsize + se_Fsize), width = 0.2) +
  scale_color_manual(values = c("#F39200FF", "#EFD500FF","#95C11FFF","#007B3DFF","#31B7BCFF")) +
  scale_x_discrete(breaks = unique(Fsize_summary$treatment)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20)) +
  labs(x = "Treatment", y = "Fpop (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),legend.direction = "horizontal",legend.position = "none") +
  geom_text(data = asterisks, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks$colours, size = 4, inherit.aes = FALSE)  +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

Fpop_model_plot_bar <-  guide_area() / (Fpop_plot | biplot_gg) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(1,10)) +
  plot_layout(nrow = (2), heights = c(0.5,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 8), 
        legend.direction = "horizontal", legend.position = "top", 
        legend.spacing.x = unit(0.1, "cm"), legend.box = "horizontal",
        legend.spacing.y = unit(-0.3, "cm"), 
        plot.tag.position = c(0.01, 0.97)) &
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


#--- playing around with correlations -------------------------
library(GGally)
ggpairs(Fpop, mapping = aes(fill = treatment, colour = treatment), columns = c("treatment", "Cocco", "micro", "nanos", "peuks", "syn")) +
  theme_bw()