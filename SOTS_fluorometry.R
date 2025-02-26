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

#------load the fluorometry data  into R---------------------------------------------------------------------------------------
LIFT <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/SOTS experiment/sots_lift.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(treatment = as.factor(treatment) %>% fct_recode("+DFB" = "DFB", "+Fe" = "Fe"),
         bottle = as.factor(bottle),
         time = as.factor(time),
         Sig = ifelse(!is.na(Sig), Sig/100, NA)) %>%
         arrange(desc(time))

# Reorder  levels
LIFT$treatment <- factor(LIFT$treatment, levels = c("initial","+DFB", "Control", "+Fe"), ordered = TRUE)
#Check the structure
str(LIFT)

# Divide the data into two data frames, one for time = "5" and one for time = "10"
time_0 <- LIFT %>% filter(time == "0") %>% select(-XG, -Fm, -Fv, -FvFo) %>% mutate_at(vars(Fo:Tau3), as.numeric)
time_2 <- LIFT %>% filter(time == "2") %>% select(-XG, -Fm, -Fv, -FvFo) %>% mutate_at(vars(Fo:Tau3), as.numeric)
time_4 <- LIFT %>% filter(time == "4") %>% select(-XG,-Fm, -Fv, -FvFo) %>% mutate_at(vars(Fo:Tau3), as.numeric) 
time_6 <- LIFT %>% filter(time == "6") %>% select(-XG,-Fm, -Fv, -FvFo) %>% mutate_at(vars(Fo:Tau3), as.numeric) 

# check the structure of the two data frames
str(time_2)
str(time_4)
str(time_6)

# reshape the data frame and create a density plot, then bar plot
time_6_long <- time_6 %>%
  gather(variable, value, Fo:Tau3) %>%
  mutate(value = as.numeric(value))

# ----------- brief look at data distributions ----------------------------------------------------------------------------------       
ggplot(time_6_long, aes(x = value)) +
  geom_density() +
  facet_wrap(~ variable, ncol = 3, scales = "free")

ggplot(time_6_long, aes(y=value)) +
  geom_boxplot() +
  facet_wrap(variable ~ ., ncol = 3, scales = 'free')

#----make a naughty little plot of Fo with time ------------------------------------------------------
Fo_time <- ggplot(LIFT, aes(y = Fo, x = time, shape = treatment)) +
  stat_summary(fun.y = mean, geom = "point", size = 5, aes(group = treatment), color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(group = treatment)) +
  ylab("Fo") +
  xlab("Time (days)") +
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

#----------scatterplot of Sig vs Fv/Fm ---------------------------------------------------------------------------------------------------
fvfmVSsigma <- ggplot(LIFT, aes(y = FvFm, x = Sig, colour = time, shape = treatment)) +
  geom_point(size = 5) +
  ylab("Fv/Fm") +
  xlab(expression(sigma[PSII])) +
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  scale_colour_manual(values = c("#95C11FFF","#DCDCDCFF" ,"#A9A9A9FF", "#000000FF")) +
  scale_x_continuous(limits = c(6, 10), breaks = seq(6, 10, by = 1)) +
  scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, by = 0.1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))

ggsave("fvfmVSsigma_plot.svg", fvfmVSsigma, width = 5.6, height = 4.6)

# Filtering out the data for linear fits where time is not equal to 0
LIFT_filtered_for_fit <- subset(LIFT, time != 0)

# Function to compute R squared for each group
compute_r2_base <- function(time_val, data) {
  subset_data <- data[data$time == time_val, ]
  model <- lm(FvFm ~ Sig, data = subset_data)
  r_squared <- summary(model)$r.squared
  return(data.frame(time = time_val, r2 = r_squared))
}

times <- unique(LIFT_filtered_for_fit$time)
r2_values <- do.call(rbind, lapply(times, compute_r2_base, data = LIFT_filtered_for_fit))


fvfmVSsigma <- ggplot(LIFT, aes(y = FvFm, x = Sig, colour = time, shape = treatment)) +
  geom_point(size = 5) +
  geom_smooth(data = LIFT_filtered_for_fit, method = "lm", se = FALSE, aes(group = time), size = 0.5) +  # Make lines thinner with size = 0.5
  ylab("Fv/Fm") +
  xlab(expression(sigma[PSII])) +
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  scale_colour_manual(values = c("#95C11FFF","#DCDCDCFF" ,"#A9A9A9FF", "#000000FF")) +
  scale_x_continuous(limits = c(6, 10), breaks = seq(6, 10, by = 1)) +
  scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, by = 0.1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

ggsave("fvfmVSsigma_plot.svg", fvfmVSsigma, width = 5.6, height = 4.6)


#---------------- Fit all models at once----------------------------------------------------------------------------------------------------
# List of variables
variables_list <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP")

#fitting the models    
models_list <- variables_list %>%
  map(~{
    if (.x == "Tau3") {formula_str <- paste0("I((", .x, "))", " ~ treatment")} 
    else {formula_str <- paste0(.x, " ~ treatment")}
    glm(as.formula(formula_str), data = time_6, family = gaussian, na.action = na.exclude)}) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Extract model summaries into a tidy format
model_summaries <- models_list %>% map(broom::tidy)
model_details <- models_list %>% map(broom::glance)
# View model summaries and details
model_summaries
model_details

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
std_residuals_list <- models_list %>% 
  map(~rstandard(.x)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Create a list of QQ plots for each model in models_list and models_list_reduced
qq_plot_list <- map(names(std_residuals_list), 
                    ~create_qq_plot(std_residuals_list[.x], paste0(.x, " Full")))

# Convert the list of plots to a patchwork object
qq_plot_patchwork <- wrap_plots(qq_plot_list, ncol = 2)

# Print the plot
print(qq_plot_patchwork)

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

# Create a list of Residuals vs Fitted plots for each model in models_list and models_list_reduced
res_fit_plot_list <- lapply(names(models_list), function(x) {
  res_fit_plot(models_list[[x]], paste0(x, " Full"))
})

# Convert the list of plots to a patchwork object
residual_plot_patchwork <- wrap_plots(res_fit_plot_list, ncol = 2)

# Print the plot
print(residual_plot_patchwork)

# ---------------Plot the raw data alongside the model comparisons using visreg.--------------------------------------------------------------
# Calculate the mean and standard error for each unique combination
LIFT_summary <- LIFT %>%
  filter(treatment != "initial") %>%
  group_by(treatment, time) %>%
  summarise(
    FvFm_mean = mean(FvFm),
    FvFm_se = sd(FvFm) / sqrt(n()),
    Sig_mean = mean(Sig),
    Sig_se = sd(Sig) / sqrt(n()),
    Tau1_mean = mean(Tau1),
    Tau1_se = sd(Tau1) / sqrt(n()),
    Tau2_mean = mean(Tau2),
    Tau2_se = sd(Tau2) / sqrt(n()),
    Tau3_mean = mean(Tau3)/1000,
    Tau3_se = (sd(Tau3) / sqrt(n()))/1000,
    PQP_mean = mean(PQP),
    PQP_se = sd(PQP) / sqrt(n())
  )

# Generate function for facet grid with point and error bars
generate_plots_with_errorbars <- function(data, variable_mean, variable_se, variable_name, show_legend = FALSE, 
                                          show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, legend_position = "top", hline_y) {
  p <- ggplot(data, aes(x = treatment, y = !!sym(variable_mean), color = time)) +
    geom_hline(yintercept = hline_y, linetype = "dashed", color = "#95C11FFF") +
    geom_point(size = 3, shape = 19) +
    geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), width = 0.15) +
    scale_colour_manual(values = c("#C6C5C5FF" ,"#706F6FFF", "#000000FF")) +
    scale_fill_manual(values = c("#C6C5C5FF" ,"#706F6FFF", "#000000FF"))  +
    scale_x_discrete(breaks = unique(data$treatment)) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = if (show_xlab) "Treatment" else NULL,
         y = y_lab,
         color = "Time (days)") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, colour = "black"), 
          panel.grid.minor = element_blank(), 
          legend.direction = "horizontal", legend.position = legend_position) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)) 
  return(p)
}

# plots for each different photophysiological outcome.
tau <- "\u03C4" #tau unicode
plot_FvFm <- generate_plots_with_errorbars(LIFT_summary, "FvFm_mean", "FvFm_se", "FvFm", show_xlab = FALSE, y_min = 0.3, y_max = 0.6,
                                           y_lab = "Fv/Fm", legend_position = "top", hline_y = 0.5335) + 
                                          scale_y_continuous(breaks = seq(0.3, 0.6, by = 0.05))
plot_Sig <- generate_plots_with_errorbars(LIFT_summary, "Sig_mean", "Sig_se", "Sig", show_xlab = TRUE, y_min = 6, y_max = 10, 
                                          y_lab = bquote(sigma[PSII]), legend_position = "none", hline_y = 7.628375)+ 
  scale_y_continuous(breaks = seq(6, 10, by = 0.5))
plot_Tau1 <- generate_plots_with_errorbars(LIFT_summary, "Tau1_mean", "Tau1_se", "Tau1", show_xlab = FALSE, y_min = 0, y_max = 3000,
                                           y_lab = expression(tau[1]~ (µs)), legend_position = "none", hline_y = 674.25) + 
  scale_y_continuous(breaks = seq(0, 3000, by = 500))
plot_Tau2 <- generate_plots_with_errorbars(LIFT_summary, "Tau2_mean", "Tau2_se", "Tau2", show_xlab = FALSE, y_min = 0, y_max = 5000,
                                           y_lab = expression(tau[2]~ (µs)), legend_position = "none", hline_y = 922.5) + 
  scale_y_continuous(breaks = seq(0, 5000, by = 500))
plot_Tau3 <- generate_plots_with_errorbars(LIFT_summary, "Tau3_mean", "Tau3_se", "Tau3", show_xlab = TRUE, y_min = 0, y_max = 30,
                                           y_lab = expression(tau[3]~ (ms)), legend_position = "none", hline_y = 8.05675) + 
  scale_y_continuous(breaks = seq(0, 30, by = 5))
plot_PQP <- generate_plots_with_errorbars(LIFT_summary, "PQP_mean", "PQP_se", "PQP", show_xlab = FALSE, y_min = 0, y_max = 9,
                                           y_lab = expression("PQP size (PQ:PSII)"), legend_position = "none", hline_y = 4.506) + 
  scale_y_continuous(breaks = seq(0, 9, by = 1.5))

#---------- Investigate group comparisons using emmeans-------------------------------------------------------------------------------------------
emm_lift <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
  map(~emmeans(models_list[[.x]], specs = pairwise ~ treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Combine emmeans and contrasts for each data frame
emmeans_FvFm <- bind_rows(emmeans = emm_lift$FvFm$emmeans, contrast = emm_lift$FvFm$contrasts)
emmeans_Sig <- bind_rows(emmeans = emm_lift$Sig$emmeans, contrast = emm_lift$Sig$contrasts)
emmeans_Tau1 <- bind_rows(emmeans = emm_lift$Tau1$emmeans, contrast = emm_lift$Tau1$contrasts)
emmeans_Tau2 <- bind_rows(emmeans = emm_lift$Tau2$emmeans, contrast = emm_lift$Tau2$contrasts)
emmeans_Tau3 <- bind_rows(emmeans = emm_lift$Tau3$emmeans, contrast = emm_lift$Tau3$contrasts)
emmeans_PQP <- bind_rows(emmeans = emm_lift$PQP$emmeans, contrast = emm_lift$PQP$contrasts)

# Create a list of combined data frames
combined_lift_emm <- list(
  emmeans_FvFm = emmeans_FvFm,
  emmeans_Sig = emmeans_Sig,
  emmeans_Tau1 = emmeans_Tau1,
  emmeans_Tau2 = emmeans_Tau2,
  emmeans_Tau3 = emmeans_Tau3,
  emmeans_PQP = emmeans_PQP
)

# Export the combined data frames to an Excel file
write_xlsx(combined_lift_emm, "LIFT_emmeanss.xlsx")

# -------------------------------Create visreg plots -----------------------------------------------------------------------
#---adding asterisks for significance ----------------------------------------------------------------
asterisks_fvfm <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(0.395, 0.385, 0.5), 
  label = c("a", "a", "b"), colours = c("black", "black", "black"))

asterisks_sig <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(8.80, 9.40, 8.55), 
  label = c("ab", "a", "b"), colours = c("black", "black", "black"))

asterisks_tau1 <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(2900, 2500, 2300), 
  label = c("a", "a", "a"), colours = c("black", "black", "black"))

asterisks_tau2 <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(4300, 4000, 4500), 
  label = c("a", "a", "a"), colours = c("black", "black", "black"))

asterisks_tau3 <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(20000, 23000, 28000), 
  label = c("a", "a", "b"), colours = c("black", "black", "black"))

asterisks_PQP <- data.frame(
  x_pos = c(0.212, 0.57, 0.928), y_pos = c(6.1, 5.7, 6.8), 
  label = c("a", "a", "a"), colours = c("black", "black", "black"))

# Remove y-axis label for visreg plots
visreg_FvFm <- visreg(models_list$FvFm, "treatment", top = 'points', partial = TRUE, gg = TRUE,
         line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black",shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(0.3, 0.6), breaks = seq(0.3, 0.6, by = 0.05)) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top") +
  geom_text(data = asterisks_fvfm, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_fvfm$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_Sig <- visreg(models_list$Sig, "treatment", top = 'points', partial = TRUE, gg = TRUE,
                     line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black",shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(6, 10), breaks = seq(4, 12, by = 1)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  geom_text(data = asterisks_sig, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_sig$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_Tau1 <- visreg(models_list$Tau1, "treatment", top = 'points', partial = TRUE, gg = TRUE,
                      line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black", shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  geom_text(data = asterisks_tau1, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_tau1$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_Tau2 <- visreg(models_list$Tau2, "treatment", top = 'points', partial = TRUE, gg = TRUE,
                      line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black",shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, by = 500)) +
  ylab("") +  
  xlab("") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(), axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  geom_text(data = asterisks_tau2, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_tau2$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_Tau3 <- visreg(models_list$Tau3, "treatment", top = 'points', partial = TRUE, gg = TRUE,
                      line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black",shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  geom_text(data = asterisks_tau3, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_tau3$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

visreg_PQP <- visreg(models_list$PQP, "treatment", top = 'points', partial = TRUE, gg = TRUE,
                      line.par = list(col="#706F6FFF"), fill.par = list(fill = "#706F6F", alpha = 0.3), points.par = list(col="black",shape = 1, size = 2)) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 1.5)) +
  ylab("") +  
  xlab("Treatment") +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  geom_text(data = asterisks_PQP, aes(x = x_pos, y = y_pos, label = label), 
            colour = asterisks_PQP$colours, size = 4, inherit.aes = FALSE)  +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


#---------------------- combined plots and printing ------------------------------------------
fvfmsig_plot <- guide_area() / (plot_FvFm | visreg_FvFm) / (plot_Sig | visreg_Sig) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

tau_plot <- guide_area() / (plot_Tau1 | visreg_Tau1) / (plot_Tau2 | visreg_Tau2) / (plot_PQP | visreg_PQP) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

PQP_plot <- guide_area() / (plot_PQP | visreg_PQP) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("fvfmsig_model_plot.svg", fvfmsig_plot, width = 7.2, height = 7.4)
ggsave("tau_model_plot2.svg", tau_plot, width = 7.2, height = 7.4)   
ggsave("PQP_model_plot.svg", PQP_plot, width = 7.2, height = 3.52)

#----- Making a PCA and biplot -------------------------------------------------------------------------
# Exclude non-numeric columns like factors
numeric_data <- time_6[, sapply(time_6, is.numeric)]
numeric_data$Fo <- NULL

# Scale the data
scaled_data <- scale(numeric_data)
#make the PCA
pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
biplot(pca_result, scale = 0)
# Extract scores and loadings from PCA result
scores <- as.data.frame(pca_result$x)
loadings <- as.data.frame(pca_result$rotation)

# Combine treatment information with scores
scores$treatment <- time_6$treatment

biplot_gg <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = treatment), size = 5) + 
  scale_shape_manual(values = c("initial" = 4, "+DFB" = 13, "Control" = 1, "+Fe" = 19)) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*3, yend = PC2*3),
               arrow = arrow(type = "open", length = unit(0.1, "inches")), color = "#D51317FF") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  geom_text(data = loadings, aes(x = PC1*2.5, y = PC2*1, label = rownames(loadings)), 
            vjust = -0.5, color = "#D51317FF") +
  labs(x = paste("PC1:", round(100 * summary(pca_result)$importance[2, 1], 1), "%"), 
       y = paste("PC2:", round(100 * summary(pca_result)$importance[2, 5], 1), "%")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))

ggsave("Fluorescence_biplot.svg", biplot_gg, width = 5.6, height = 4.6)

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
FvFm_stan_full <- stan_glm(FvFm ~  treatment, data = time_6, family = gaussian, iter = 5000)
Sig_stan_full <- stan_glm(Sig ~ treatment, data = time_6, family = gaussian, iter = 5000)
Tau1_stan_full <- stan_glm(Tau1 ~ treatment, data = time_6, family = gaussian, iter = 5000)
Tau2_stan_full <- stan_glm(Tau2 ~ treatment, data = time_6, family = gaussian, iter = 5000)
Tau2_stan_light <- stan_glm(Tau2 ~ treatment, data = time_6, family = gaussian, iter = 5000)
Tau3_stan_full <- stan_glm(Tau3 ~ treatment, data = time_6, family = gaussian, iter = 5000)
PQP_stan_full <- stan_glm(PQP ~ treatment, data = time_6, family = gaussian, iter = 5000)