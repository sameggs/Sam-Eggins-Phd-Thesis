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

#------load the fluorometry data  into R---------------------------------------------------------------------------------------
    LIFT <- read.csv("/Users/eggboy/Dropbox/Science/Data/Voyages/DCM experiment/dcm_lift.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
      mutate(light = as.factor(light), 
             inhibitor = as.factor(inhibitor),
             treatment = as.factor(treatment) %>% fct_recode("DCM" = "initial","+DFB" = "DFB", "+Fe" = "Fe"),
             bottle = as.factor(bottle),
             time = as.factor(time))
# Reorder  levels
    LIFT$treatment <- factor(LIFT$treatment, levels = c("DCM","+DFB", "Control", "+Fe"), ordered = TRUE)
#Check the structure
    str(LIFT)

# Divide the data into two data frames, one for time = "5" and one for time = "10"
    time_0 <- LIFT %>% filter(time == "0") %>% select(-XG, -Fm, -Fv, -FvFo) %>% mutate_at(vars(FvFm:Tau3), as.numeric)
    time_5 <- LIFT %>% filter(time %in% c("5", "0")) %>% select(-XG, -Fm, -Fv, -FvFo) %>% mutate_at(vars(FvFm:Tau3), as.numeric)
    time_10 <- LIFT %>% filter(time %in% c("10", "0")) %>% select(-XG,-Fm, -Fv, -FvFo) %>% mutate_at(vars(FvFm:Tau3), as.numeric) 

# check the structure of the two data frames
    str(time_5)
    str(time_10)

# reshape the data frame and create a density plot, then bar plot
    time_10_long <- time_10 %>%
      gather(variable, value, Fo:Tau3) %>%
      mutate(value = as.numeric(value))
 
# ----------- brief look at data distributions ----------------------------------------------------------------------------------       
    ggplot(time_10_long, aes(x = value)) +
      geom_density() +
      facet_wrap(~ variable, ncol = 3, scales = "free")
    
    ggplot(time_10_long, aes(y=value)) +
      geom_boxplot() +
      facet_wrap(variable ~ ., ncol = 3, scales = 'free')

#----------scatterplot of Sig vs Fv/Fm ---------------------------------------------------------------------------------------------------
fvfmVSsigma_5 <- ggplot(time_5, aes(y = FvFm, x = Sig, color = light, fill = light, shape = treatment)) +
      geom_point(size = 5) +
      ylab("") +
      xlab(bquote(sigma^PSII)) +
      scale_x_continuous(limits = c(400, 1200), breaks = seq(400, 1200, by = 200)) +
      scale_y_continuous(limits = c(0.3, 0.9), breaks = seq(0.3, 0.9, by = 0.1)) +
      scale_shape_manual(values = c("DCM" = 4 ,"+DFB" = 1, "Control" = 21, "+Fe" = 19)) +
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      theme_bw() +
      theme(legend.title = element_blank(), panel.grid.minor = element_blank())
    
fvfmVSsigma_10 <- ggplot(time_10, aes(y = FvFm, x = Sig, color = light, fill = light, shape = treatment)) +
      geom_point(size = 5) +
      ylab("Fv/Fm") +
      xlab("") +
      scale_x_continuous(limits = c(400, 1200), breaks = seq(400, 1200, by = 200)) +
      scale_y_continuous(limits = c(0.3, 0.9), breaks = seq(0.3, 0.9, by = 0.1)) +
      scale_shape_manual(values = c("DCM" = 4 ,"+DFB" = 1, "Control" = 21, "+Fe" = 19)) +
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      theme_bw() +
      theme(legend.title = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_blank())

fvfmsig_X <- guide_area() / (fvfmVSsigma_5 | fvfmVSsigma_10) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("fvfmXsigma_plot.svg", fvfmsig_X, width = 7.2, height = 4.2)

#---------------- Fit all models at once----------------------------------------------------------------------------------------------------
time_10 <- time_10 %>% filter(time %in% c("10"))
# List of variables
variables_list <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP")

#fitting the models    
models_list <- variables_list %>% 
      map(~glm(as.formula(paste0(.x, " ~ light + treatment + light*treatment")), data = time_10, family = gaussian, na.action = na.exclude)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
    
models_list_reduced <- variables_list %>% 
      map(~glm(as.formula(paste0(.x, " ~ light + treatment")), data = time_10, family = gaussian, na.action = na.exclude)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
    
# Extract model summaries into a tidy format
    model_summaries <- models_list %>% map(broom::tidy)
    model_details <- models_list %>% map(broom::glance)
# View model summaries and details
    model_summaries
    model_details
    
#---- FItting MODELS using Bayesian methods --------------------------------------------------------
FvFm_stan_full <- stan_glm(FvFm ~ light * treatment, data = time_10, family = gaussian, iter = 5000)
FvFm_stan_reduced <- stan_glm(FvFm ~ light + treatment, data = time_10, family = gaussian, iter = 5000)
Sig_stan_full <- stan_glm(Sig ~ light * treatment, data = time_10, family = gaussian, iter = 15000)
Sig_stan_reduced <- stan_glm(Sig ~ light + treatment, data = time_10, family = gaussian, iter = 15000)
Tau1_stan_full <- stan_glm(Tau1 ~ light * treatment, data = time_10, family = gaussian, iter = 10000)
Tau1_stan_reduced <- stan_glm(Tau1 ~ light + treatment, data = time_10, family = gaussian, iter = 10000)
Tau2_stan_full <- stan_glm(Tau2 ~ light * treatment, data = time_10, family = gaussian, iter = 5000)
Tau2_stan_reduced <- stan_glm(Tau2 ~ light + treatment, data = time_10, family = gaussian, iter = 5000)
Tau2_stan_light <- stan_glm(Tau2 ~ light + treatment, data = time_10, family = gaussian, iter = 5000)
Tau3_stan_full <- stan_glm(Tau3 ~ light * treatment, data = time_10, family = gaussian, iter = 5000)
Tau3_stan_reduced <- stan_glm(Tau3 ~ light + treatment, data = time_10, family = gaussian, iter = 5000)
PQP_stan_full <- stan_glm(PQP ~ light * treatment, data = time_10, family = gaussian, iter = 5000)
PQP_stan_reduced <- stan_glm(PQP ~ light + treatment, data = time_10, family = gaussian, iter = 5000)

# Compute LOO estimates for each model    
FvFm_loo_full <- loo(FvFm_stan_full, k_threshold = 0.7)
FvFm_loo_reduced <- loo(FvFm_stan_reduced, k_threshold = 0.7)
Sig_loo_full <- loo(Sig_stan_full, k_threshold = 0.7)
Sig_loo_reduced <- loo(Sig_stan_reduced, k_threshold = 0.7)
Tau1_loo_full <- loo(Tau1_stan_full, k_threshold = 0.7)
Tau1_loo_reduced <- loo(Tau1_stan_reduced, k_threshold = 0.7)
Tau2_loo_full <- loo(Tau2_stan_full, k_threshold = 0.7)
Tau2_loo_reduced <- loo(Tau2_stan_reduced, k_threshold = 0.7)
Tau3_loo_full <- loo(Tau3_stan_full, k_threshold = 0.7)
Tau3_loo_reduced <- loo(Tau3_stan_reduced, k_threshold = 0.7)
PQP_loo_full <- loo(PQP_stan_full, k_threshold = 0.7)
PQP_loo_reduced <- loo(PQP_stan_reduced, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(FvFm_loo_full, FvFm_loo_reduced)
loo_compare(Sig_loo_full, Sig_loo_reduced)
loo_compare(Tau1_loo_full, Tau1_loo_reduced)
loo_compare(Tau2_loo_full, Tau2_loo_reduced)
loo_compare(Tau3_loo_full, Tau3_loo_reduced)
loo_compare(PQP_loo_full, PQP_loo_reduced)

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
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3","PQP"))
    
    std_residuals_list_reduced <- models_list_reduced %>% 
      map(~rstandard(.x)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
    
    # Create a list of QQ plots for each model in models_list and models_list_reduced
    qq_plot_list <- map(names(std_residuals_list), 
                        ~create_qq_plot(std_residuals_list[.x], paste0(.x, " Full")))
    
    qq_plot_list_reduced <- map(names(std_residuals_list_reduced), 
                                ~create_qq_plot(std_residuals_list_reduced[.x], paste0(.x, " Reduced")))
    
    # Combine the two lists of plots
    combined_qq_plots <- map2(qq_plot_list, qq_plot_list_reduced, ~(.x | .y))
    
    # Convert the list of plots to a patchwork object
    qq_plot_patchwork <- wrap_plots(combined_qq_plots, ncol = 2)
    
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
    
    res_fit_plot_list_reduced <- lapply(names(models_list_reduced), function(x) {
      res_fit_plot(models_list_reduced[[x]], paste0(x, " Reduced"))
    })
    
    # Combine the two lists of plots
    combined_res_fit_plots <- mapply(function(x, y) {x | y}, 
                                     res_fit_plot_list, 
                                     res_fit_plot_list_reduced, 
                                     SIMPLIFY = FALSE)
    
    # Convert the list of plots to a patchwork object
    residual_plot_patchwork <- wrap_plots(combined_res_fit_plots, ncol = 2)
    
    # Print the plot
    print(residual_plot_patchwork)
      
#---------- Investigate group comparisons using emmeans-------------------------------------------------------------------------------------------
emm_lift <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
      map(~emmeans(models_list[[.x]], specs = pairwise ~ light:treatment)) %>%
      map(~summary(., infer = TRUE)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
emm_lift_light <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
      map(~emmeans(models_list[[.x]], specs = pairwise ~ light|treatment)) %>%
      map(~summary(., infer = TRUE)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
emm_lift_trt <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
      map(~emmeans(models_list[[.x]], specs = pairwise ~ treatment|light)) %>%
      map(~summary(., infer = TRUE)) %>%
      set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
    
# Combine emmeans and contrasts for each data frame
emmeans_FvFm <- bind_rows(emmeans = emm_lift$FvFm$emmeans, light = emm_lift_light$FvFm$contrasts, trt = emm_lift_trt$FvFm$contrasts)
emmeans_Sig <- bind_rows(emmeans = emm_lift$Sig$emmeans, light = emm_lift_light$Sig$contrasts, trt = emm_lift_trt$Sig$contrasts)
emmeans_Tau1 <- bind_rows(emmeans = emm_lift$Tau1$emmeans, light = emm_lift_light$Tau1$contrasts, trt = emm_lift_trt$Tau1$contrasts)
emmeans_Tau2 <- bind_rows(emmeans = emm_lift$Tau2$emmeans, light = emm_lift_light$Tau2$contrasts, trt = emm_lift_trt$Tau2$contrasts)
emmeans_Tau3 <- bind_rows(emmeans = emm_lift$Tau3$emmeans, light = emm_lift_light$Tau3$contrasts, trt = emm_lift_trt$Tau3$contrasts)
emmeans_PQP <- bind_rows(emmeans = emm_lift$PQP$emmeans, light = emm_lift_light$PQP$contrasts, trt = emm_lift_trt$PQP$contrasts)

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

# ---------------Plot the raw data alongside the model comparisons using visreg.--------------------------------------------------------------
 # Calculate the mean and standard error for each unique combination
    time_10_summary <- time_10 %>%
      group_by(treatment, light) %>%
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
        guides(
          color = guide_legend(nrow = 1, byrow = TRUE),
          shape = guide_legend(nrow = 1, byrow = TRUE)
        ) 
      return(p)
    }

 # plots for each different photophysiological outcome.
    tau <- "\u03C4" #tau unicode
    plot_FvFm <- generate_plots_with_errorbars(time_10_summary, "FvFm_mean", "FvFm_se", "FvFm", show_xlab = FALSE, y_min = 0.25, y_max = 0.85,
                                               y_lab = "Fv/Fm", legend_position = "none", hline_y = 0.501) + 
                                               scale_y_continuous(breaks = seq(0.3, 0.8, by = 0.1))
    plot_Sig <- generate_plots_with_errorbars(time_10_summary, "Sig_mean", "Sig_se", "Sig", show_xlab = TRUE, y_min = 500, y_max = 1300, 
                                              y_lab = bquote(sigma[PSII]), legend_position = "none", hline_y = 790.63)
    plot_Tau1 <- generate_plots_with_errorbars(time_10_summary, "Tau1_mean", "Tau1_se", "Tau1", show_xlab = FALSE, y_min = -50, y_max = 1550,
                                               y_lab = expression(tau[1]~ (µs)), legend_position = "none", hline_y = 257.25)
    plot_Tau2 <- generate_plots_with_errorbars(time_10_summary, "Tau2_mean", "Tau2_se", "Tau2", show_xlab = FALSE, y_min = 0, y_max = 3000,
                                               y_lab = expression(tau[2]~ (µs)), legend_position = "none", hline_y = 581)
    plot_Tau3 <- generate_plots_with_errorbars(time_10_summary, "Tau3_mean", "Tau3_se", "Tau3", show_xlab = TRUE, y_min = 0, y_max = 25,
                                               y_lab = expression(tau[3]~ (ms)), legend_position = "none", hline_y = 10.088)   
    plot_PQP <- generate_plots_with_errorbars(time_10_summary, "PQP_mean", "PQP_se", "Tau1", show_xlab = FALSE, y_min = 0, y_max = 9,
                                               y_lab = expression("PQP size (PQ:PSII)"), legend_position = "none", hline_y = 5.21225) + 
                                                scale_y_continuous(breaks = seq(0, 9, by = 3))
    
#------ Now an attempt to put both the times on each plot -------------------------------------------------------
    # Calculate the mean and standard error for each unique combination
    LIFT_summary <- LIFT %>%
      filter(time %in% c("5", "10")) %>%
      group_by(treatment, light, time) %>%
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
    
    #make function to generate plots...
    generate_plots_with_errorbars <- function(data, variable_mean, variable_se, show_legend = FALSE, 
                                              show_xlab = TRUE, y_min = NULL, y_max = NULL, y_lab = NULL, legend_position = "top", hline_y) {
      
      # Create the x_pos column for positioning
      data$x_pos <- as.numeric(as.factor(data$treatment)) + ifelse(data$time == "5", -0.1, 0.1)
      
      p <- ggplot(data, aes(x = x_pos, y = !!sym(variable_mean), color = light, fill = time, shape = time)) +
        geom_hline(yintercept = hline_y, linetype = "dashed", color = "#95C11FFF") +
        geom_line(aes(group = interaction(light, treatment)), linewidth = 0.5, show.legend = FALSE) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = !!sym(variable_mean) - !!sym(variable_se), ymax = !!sym(variable_mean) + !!sym(variable_se)), 
                      width = 0.15) +
        scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
        scale_fill_manual(values = c("white", "black")) +
        scale_shape_manual(values = c("5" = 21, "10" = 21)) +
        scale_x_continuous(breaks = 2:4, limits = c(1.55,4.45), labels = unique(data$treatment)) +
        coord_cartesian(ylim = c(y_min, y_max)) +
        labs(x = if (show_xlab) "Treatment" else NULL,
             y = y_lab) +
        theme_bw() +
        theme(axis.text = element_text(size = 10, colour = "black"), 
              panel.grid.minor = element_blank(), legend.title = element_blank(),
              legend.direction = "horizontal", legend.position = legend_position) +
        guides(
          color = guide_legend(nrow = 1, byrow = TRUE),
          shape = guide_legend(nrow = 1, byrow = TRUE)
        ) 
      return(p)
    }
    
    plot_FvFm_time <- generate_plots_with_errorbars(LIFT_summary, "FvFm_mean", "FvFm_se", "FvFm", show_xlab = FALSE, y_min = 0.25, y_max = 0.85,
                                               y_lab = "Fv/Fm", legend_position = "none", hline_y = 0.501) + 
      scale_y_continuous(breaks = seq(0.3, 0.8, by = 0.1))
    plot_Sig_time <- generate_plots_with_errorbars(LIFT_summary, "Sig_mean", "Sig_se", "Sig", show_xlab = TRUE, y_min = 500, y_max = 1300, 
                                              y_lab = bquote(sigma[PSII]), legend_position = "none", hline_y = 790.63)
    plot_Tau1_time <- generate_plots_with_errorbars(LIFT_summary, "Tau1_mean", "Tau1_se", "Tau1", show_xlab = FALSE, y_min = -50, y_max = 1550,
                                               y_lab = expression(tau[1]~ (µs)), legend_position = "none", hline_y = 257.25)
    plot_Tau2_time <- generate_plots_with_errorbars(LIFT_summary, "Tau2_mean", "Tau2_se", "Tau2", show_xlab = FALSE, y_min = 0, y_max = 3000,
                                               y_lab = expression(tau[2]~ (µs)), legend_position = "none", hline_y = 581)
    plot_Tau3_time <- generate_plots_with_errorbars(LIFT_summary, "Tau3_mean", "Tau3_se", "Tau3", show_xlab = TRUE, y_min = 0, y_max = 25,
                                               y_lab = expression(tau[3]~ (ms)), legend_position = "none", hline_y = 10.088)   
    plot_PQP_time <- generate_plots_with_errorbars(LIFT_summary, "PQP_mean", "PQP_se", "Tau1", show_xlab = FALSE, y_min = 0, y_max = 9,
                                              y_lab = expression("PQP size (PQ:PSII)"), legend_position = "none", hline_y = 5.21225) + 
      scale_y_continuous(breaks = seq(0, 9, by = 3))    
    
# -------------------------------Create visreg plots -----------------------------------------------------------------------
#---adding asterisks for significance ----------------------------------------------------------------
    asterisks_fvfm <- data.frame(
      x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.072), y_pos = c(0.361, 0.302, 0.482, 0.624, 0.658, 0.733, 0.83), 
      label = c("a", "b", "c", "a", "a", "b", "all:\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black"))
    
    asterisks_sig <- data.frame(
      x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858), y_pos = c(1213, 1185, 504, 844, 799, 813), 
      label = c("a", "a", "b", "a", "a", "b"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD"))

    asterisks_tau1 <- data.frame(
      x_pos = c(0.38), y_pos = c(1260), 
      label = c("\u263C"), colours = c("black"))
    
    asterisks_PQP <- data.frame(
      x_pos = c(0.142, 0.5, 0.858, 0.142, 0.5, 0.858, 0.022, 0.738), y_pos = c(7.91, 7.57, 6.16, 2.66, 3.47, 1.45, 7.91, 6.16), 
      label = c("a", "a", "a", "ab", "aa", "ab", "\u263C", "\u263C"), colours = c("#D51317", "#D51317", "#D51317", "#0094CD", "#0094CD", "#0094CD", "black", "black"))
    
# Remove y-axis label for visreg plots
    visreg_FvFm <- visreg(models_list$FvFm, "treatment", by="light", overlay = TRUE, 
                          gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(0.25, 0.85), breaks = seq(0.3, 0.8, by = 0.1)) +
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
    
    visreg_Sig <- visreg(models_list$Sig, "treatment", by="light", overlay = TRUE, 
                         gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(500, 1300), breaks = seq(600, 1200, by = 200)) +
      ylab("") +  
      xlab("Treatment") +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.text.y = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "none") +
      geom_text(data = asterisks_sig, aes(x = x_pos, y = y_pos, label = label), 
                colour = asterisks_sig$colours, size = 4, inherit.aes = FALSE)  +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    
    visreg_Tau1 <- visreg(models_list$Tau1, "treatment", by="light", overlay = TRUE, 
                         gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(-50, 1550), breaks = seq(0, 1500, by = 500)) +
      ylab("") +  
      xlab("") +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.text.y = element_blank(), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "none") +
      geom_text(data = asterisks_tau1, aes(x = x_pos, y = y_pos, label = label), 
                colour = asterisks_tau1$colours, size = 4, inherit.aes = FALSE)  +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    
    visreg_Tau2 <- visreg(models_list$Tau2, "treatment", by="light", overlay = TRUE, 
                          gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 1000)) +
      ylab("") +  
      xlab("") +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.text.y = element_blank(), axis.title.x = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "none") +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    
    visreg_Tau3 <- visreg(models_list$Tau3, "treatment", by="light", overlay = TRUE, 
                          gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 10000)) +
      ylab("") +  
      xlab("Treatment") +
      theme(axis.text = element_text(size = 10, colour = "black"),
            axis.text.y = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "none") +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    
    visreg_PQP <- visreg(models_list$PQP, "treatment", by="light", overlay = TRUE, 
                          gg = TRUE, points.par = list(shape = 1, size = 2)) + 
      theme_bw() + 
      scale_colour_manual(values = c("#D51317FF", "#0094CDFF")) +
      scale_fill_manual(values = c("#D5131733", "#0094CD33"))  +
      scale_y_continuous(limits = c(0, 9), breaks = seq(0, 9, by = 3)) +
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
            legend.title = element_blank(), legend.text = element_text(size = 9), 
            legend.direction = "horizontal", legend.position = "top",
            plot.tag.position = c(0.01, 0.97)) & 
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))
    
tau_plot <- guide_area() / (plot_Tau1 | visreg_Tau1) / (plot_Tau2 | visreg_Tau2) / (plot_Tau3 | visreg_Tau3) +
      plot_layout(guides = 'collect', nrow = (4), heights = c(1.5,10,10,10)) +
      plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
      theme(axis.text = element_text(size = 10, colour = "black"),
            legend.title = element_blank(), legend.text = element_text(size = 9), 
            legend.direction = "horizontal", legend.position = "top",
            plot.tag.position = c(0.01, 0.97)) & 
      guides(fill = guide_legend(nrow = 1, byrow = TRUE))

PQP_plot <- guide_area() / (plot_PQP | visreg_PQP) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("fvfmsig_model_plot_DCM.svg", fvfmsig_plot, width = 7.2, height = 7.4)
ggsave("tau_model_plot_DCM.svg", tau_plot, width = 7.2, height = 7.4)   
ggsave("PQP_model_plot_DCM.svg", PQP_plot, width = 7.2, height = 3.52)

#---- combined plots but with the T5 data on the left -----------------------------
fvfmsig_plot_time <- guide_area() / (plot_FvFm_time | visreg_FvFm) / (plot_Sig_time | visreg_Sig) / (plot_PQP_time | visreg_PQP)+
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

tau_plot_time <- guide_area() / (plot_Tau1_time | visreg_Tau1) / (plot_Tau2_time | visreg_Tau2) / (plot_Tau3_time | visreg_Tau3) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1.5,10,10,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

PQP_plot_time <- guide_area() / (plot_PQP_time | visreg_PQP) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10)) +
  plot_annotation(tag_levels = list(c('A','B','C', 'D', 'E', 'F'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("fvfmsig_model_plot_DCM_time.svg", fvfmsig_plot_time, width = 7.2, height = 7.4)
ggsave("tau_model_plot_DCM_time.svg", tau_plot_time, width = 7.2, height = 7.4)   
ggsave("PQP_model_plot_DCM_time.svg", PQP_plot_time, width = 7.2, height = 3.52)
    
#------ OKAy ----------- time for mid point examination ----------------------------
#---------------- Fit all models at once----------------------------------------------------------------------------------------------------
time_5 <- time_5 %>% filter(time %in% c("5"))
# List of variables
variables_list <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP")

#fitting the models    
models_list <- variables_list %>% 
  map(~glm(as.formula(paste0(.x, " ~ light + treatment + light*treatment")), data = time_5, family = gaussian, na.action = na.exclude)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

models_list_reduced <- variables_list %>% 
  map(~glm(as.formula(paste0(.x, " ~ light + treatment")), data = time_5, family = gaussian, na.action = na.exclude)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Extract model summaries into a tidy format
model_summaries <- models_list %>% map(broom::tidy)
model_details <- models_list %>% map(broom::glance)
# View model summaries and details
model_summaries
model_details

#---- FItting MODELS using Bayesian methods --------------------------------------------------------
FvFm_stan_full <- stan_glm(FvFm ~ light * treatment, data = time_5, family = gaussian, iter = 5000)
FvFm_stan_reduced <- stan_glm(FvFm ~ light + treatment, data = time_5, family = gaussian, iter = 5000)
Sig_stan_full <- stan_glm(Sig ~ light * treatment, data = time_5, family = gaussian, iter = 15000)
Sig_stan_reduced <- stan_glm(Sig ~ light + treatment, data = time_5, family = gaussian, iter = 15000)
Tau1_stan_full <- stan_glm(Tau1 ~ light * treatment, data = time_5, family = gaussian, iter = 10000)
Tau1_stan_reduced <- stan_glm(Tau1 ~ light + treatment, data = time_5, family = gaussian, iter = 10000)
Tau2_stan_full <- stan_glm(Tau2 ~ light * treatment, data = time_5, family = gaussian, iter = 5000)
Tau2_stan_reduced <- stan_glm(Tau2 ~ light + treatment, data = time_5, family = gaussian, iter = 5000)
Tau2_stan_light <- stan_glm(Tau2 ~ light + treatment, data = time_5, family = gaussian, iter = 5000)
Tau3_stan_full <- stan_glm(Tau3 ~ light * treatment, data = time_5, family = gaussian, iter = 5000)
Tau3_stan_reduced <- stan_glm(Tau3 ~ light + treatment, data = time_5, family = gaussian, iter = 5000)
PQP_stan_full <- stan_glm(PQP ~ light * treatment, data = time_5, family = gaussian, iter = 5000)
PQP_stan_reduced <- stan_glm(PQP ~ light + treatment, data = time_5, family = gaussian, iter = 5000)

# Compute LOO estimates for each model    
FvFm_loo_full <- loo(FvFm_stan_full, k_threshold = 0.7)
FvFm_loo_reduced <- loo(FvFm_stan_reduced, k_threshold = 0.7)
Sig_loo_full <- loo(Sig_stan_full, k_threshold = 0.7)
Sig_loo_reduced <- loo(Sig_stan_reduced, k_threshold = 0.7)
Tau1_loo_full <- loo(Tau1_stan_full, k_threshold = 0.7)
Tau1_loo_reduced <- loo(Tau1_stan_reduced, k_threshold = 0.7)
Tau2_loo_full <- loo(Tau2_stan_full, k_threshold = 0.7)
Tau2_loo_reduced <- loo(Tau2_stan_reduced, k_threshold = 0.7)
Tau3_loo_full <- loo(Tau3_stan_full, k_threshold = 0.7)
Tau3_loo_reduced <- loo(Tau3_stan_reduced, k_threshold = 0.7)
PQP_loo_full <- loo(PQP_stan_full, k_threshold = 0.7)
PQP_loo_reduced <- loo(PQP_stan_reduced, k_threshold = 0.7)

# Compare the models using loo_compare()
loo_compare(FvFm_loo_full, FvFm_loo_reduced)
loo_compare(Sig_loo_full, Sig_loo_reduced)
loo_compare(Tau1_loo_full, Tau1_loo_reduced)
loo_compare(Tau2_loo_full, Tau2_loo_reduced)
loo_compare(Tau3_loo_full, Tau3_loo_reduced)
loo_compare(PQP_loo_full, PQP_loo_reduced)

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
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3","PQP"))

std_residuals_list_reduced <- models_list_reduced %>% 
  map(~rstandard(.x)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Create a list of QQ plots for each model in models_list and models_list_reduced
qq_plot_list <- map(names(std_residuals_list), 
                    ~create_qq_plot(std_residuals_list[.x], paste0(.x, " Full")))

qq_plot_list_reduced <- map(names(std_residuals_list_reduced), 
                            ~create_qq_plot(std_residuals_list_reduced[.x], paste0(.x, " Reduced")))

# Combine the two lists of plots
combined_qq_plots <- map2(qq_plot_list, qq_plot_list_reduced, ~(.x | .y))

# Convert the list of plots to a patchwork object
qq_plot_patchwork <- wrap_plots(combined_qq_plots, ncol = 2)

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

res_fit_plot_list_reduced <- lapply(names(models_list_reduced), function(x) {
  res_fit_plot(models_list_reduced[[x]], paste0(x, " Reduced"))
})

# Combine the two lists of plots
combined_res_fit_plots <- mapply(function(x, y) {x | y}, 
                                 res_fit_plot_list, 
                                 res_fit_plot_list_reduced, 
                                 SIMPLIFY = FALSE)

# Convert the list of plots to a patchwork object
residual_plot_patchwork <- wrap_plots(combined_res_fit_plots, ncol = 2)

# Print the plot
print(residual_plot_patchwork)

#--------------emmeans table ----------------------------------------
emm_lift <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
  map(~emmeans(models_list[[.x]], specs = pairwise ~ light:treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
emm_lift_light <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
  map(~emmeans(models_list[[.x]], specs = pairwise ~ light|treatment)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))
emm_lift_trt <- c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP") %>% 
  map(~emmeans(models_list[[.x]], specs = pairwise ~ treatment|light)) %>%
  map(~summary(., infer = TRUE)) %>%
  set_names(c("FvFm", "Sig", "Tau1", "Tau2", "Tau3", "PQP"))

# Combine emmeans and contrasts for each data frame
emmeans_FvFm <- bind_rows(emmeans = emm_lift$FvFm$emmeans, light = emm_lift_light$FvFm$contrasts, trt = emm_lift_trt$FvFm$contrasts)
emmeans_Sig <- bind_rows(emmeans = emm_lift$Sig$emmeans, light = emm_lift_light$Sig$contrasts, trt = emm_lift_trt$Sig$contrasts)
emmeans_Tau1 <- bind_rows(emmeans = emm_lift$Tau1$emmeans, light = emm_lift_light$Tau1$contrasts, trt = emm_lift_trt$Tau1$contrasts)
emmeans_Tau2 <- bind_rows(emmeans = emm_lift$Tau2$emmeans, light = emm_lift_light$Tau2$contrasts, trt = emm_lift_trt$Tau2$contrasts)
emmeans_Tau3 <- bind_rows(emmeans = emm_lift$Tau3$emmeans, light = emm_lift_light$Tau3$contrasts, trt = emm_lift_trt$Tau3$contrasts)
emmeans_PQP <- bind_rows(emmeans = emm_lift$PQP$emmeans, light = emm_lift_light$PQP$contrasts, trt = emm_lift_trt$PQP$contrasts)

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
write_xlsx(combined_lift_emm, "LIFT_emmeans_T5.xlsx")