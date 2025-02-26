library(tidyverse)
library(ggplot2)
library(ggsci)
library(scales)
library(multcomp)
library(broom)
library(ggstar)

setwd("/Users/eggboy/Dropbox/Science/Data/Growth Rates") #setwd

# import the growth rates
growth_rates <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/Zn_Growth.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
mutate(species = as.factor(species),
       EDTA = as.factor(EDTA), Zn = as.factor(Zn), Co = as.factor(Co), Cd = as.factor(Cd),
       tube = as.numeric(tube),
       transfer = as.numeric(transfer),
       growth = as.numeric(growth))

# make it so that there are numeric values for the concentrations of metals
growth_rates$Zn_conc <- ifelse(growth_rates$Zn == "-", 0.4,
                          ifelse(growth_rates$EDTA == 10 & growth_rates$Zn == "+", 7.9,
                          ifelse(growth_rates$EDTA == 10 & growth_rates$Zn == "h", 3.95,
                          ifelse(growth_rates$EDTA == 10 & growth_rates$Zn == "m", 1.2,       
                          ifelse(growth_rates$EDTA == 10 & growth_rates$Zn == "x5", 39.5,
                          ifelse(growth_rates$EDTA == 100 & growth_rates$Zn == "m", 1.2,
                          ifelse(growth_rates$EDTA == 100 & growth_rates$Zn == "+", 79, NA)))))))

growth_rates$Co_conc <- ifelse(growth_rates$EDTA == 10 & growth_rates$Co == "-", 0.01,
                          ifelse(growth_rates$EDTA == 100 & growth_rates$Co == "-", 0.1,
                          ifelse(growth_rates$EDTA == 100 & growth_rates$Co == "--", 0.01,
                          ifelse(growth_rates$EDTA == 10 & growth_rates$Co == "+", 5,
                          ifelse(growth_rates$EDTA == 100 & growth_rates$Co == "+", 50, NA)))))

growth_rates$Cd_conc <- ifelse(growth_rates$Cd == "-", 0.002,
                          ifelse(growth_rates$Cd == "0.2", 0.2,
                          ifelse(growth_rates$Cd == "1", 1,
                          ifelse(growth_rates$Cd == "2", 2,
                          ifelse(growth_rates$Cd == "10", 10,
                          ifelse(growth_rates$Cd == "100", 100, NA))))))

# now for the free ion concentrations:
growth_rates <- growth_rates %>% mutate(freeCd = case_when(
      EDTA == 10 & Cd_conc == 0.002 ~ 0.055,
      EDTA == 10 & Cd_conc == 0.2 ~ 5.55,
      EDTA == 10 & Cd_conc == 1 ~ 27.8,
      EDTA == 10 & Cd_conc == 10 ~ 257,
      EDTA == 100 & Cd_conc == 0.002 ~ 0.0051,
      EDTA == 100 & Cd_conc == 2 ~ 5.10,
      EDTA == 100 & Cd_conc == 10 ~ 25.5,
      EDTA == 100 & Cd_conc == 100 ~ 255,
      TRUE ~ NA_real_),
      Cd2 = case_when(
        EDTA == 10 & Cd_conc == 0.002 ~ 0.0021,
        EDTA == 10 & Cd_conc == 0.2 ~ 0.21,
        EDTA == 10 & Cd_conc == 1 ~ 1.06,
        EDTA == 10 & Cd_conc == 10 ~ 10.6,
        EDTA == 100 & Cd_conc == 0.002 ~ 0.00019,
        EDTA == 100 & Cd_conc == 2 ~ 0.20,
        EDTA == 100 & Cd_conc == 10 ~ 0.98,
        EDTA == 100 & Cd_conc == 100 ~ 9.8,
        TRUE ~ NA_real_),
      trtCd = case_when(
        EDTA == 10 & Cd_conc == 0.002 ~ "-Cd",
        EDTA == 10 & Cd_conc %in% c(0.2, 1, 10) ~ paste0(Cd, "Cd"),
        EDTA == 100 & Cd_conc %in% c(0.002, 2, 10, 100) ~ paste0(Cd, "Cd"),
        TRUE ~ NA_character_),
      freeZn = case_when(
        EDTA == 10 & Zn_conc == 0.4 ~ 1.8,
        EDTA == 10 & Zn_conc == 1.2 ~ 5.37,
        EDTA == 10 & Zn_conc == 3.95 ~ 17.7,
        EDTA == 10 & Zn_conc == 7.9 ~ 35.4,
        EDTA == 10 & Zn_conc == 39.5 ~ 178,
        EDTA == 100 & Zn_conc == 0.4 ~ 0.16,
        EDTA == 100 & Zn_conc == 1.2 ~ 0.49,
        EDTA == 100 & Zn_conc == 79 ~ 35.4, #actually 32.4
        TRUE ~ NA_real_),
      Zn2 = case_when(
        EDTA == 10 & Zn_conc == 0.4 ~ 0.94,
        EDTA == 10 & Zn_conc == 1.2 ~ 2.83,
        EDTA == 10 & Zn_conc == 3.95 ~ 9.3,
        EDTA == 10 & Zn_conc == 7.9 ~ 18.6,
        EDTA == 10 & Zn_conc == 39.5 ~ 93.4,
        EDTA == 100 & Zn_conc == 0.4 ~ 0.085,
        EDTA == 100 & Zn_conc == 1.2 ~ 0.26,
        EDTA == 100 & Zn_conc == 79 ~ 18.6, #actually 16.8
        TRUE ~ NA_real_),
      trtZn = case_when(
        EDTA == 10 & Zn_conc == 0.4 ~ "-Zn",
        EDTA == 10 & Zn_conc == 1.2 ~ "MZn",
        EDTA == 10 & Zn_conc == 3.95 ~ "Half Zn",
        EDTA == 10 & Zn_conc == 7.9 ~ "+Zn",
        EDTA == 10 & Zn_conc == 39.5 ~ "x5 Zn",
        EDTA == 100 & Zn_conc == 0.4 ~ "--Zn",
        EDTA == 100 & Zn_conc == 1.2 ~ "--mZn",
        EDTA == 100 & Zn_conc == 79 ~ "+Zn",
        TRUE ~ NA_character_),
      freeCo = case_when(
        EDTA == 10 & Co_conc == 0.01 ~ 0.024,
        EDTA == 10 & Co_conc == 5 ~ 12.2,
        EDTA == 100 & Co_conc == 0.1 ~ 0.022,
        EDTA == 100 & Co_conc == 0.01 ~ 0.0022,
        EDTA == 100 & Co_conc == 50 ~ 12.2, #actually 11.0
        TRUE ~ NA_real_),
      Co2 = case_when(
        EDTA == 10 & Co_conc == 0.01 ~ 0.019,
        EDTA == 10 & Co_conc == 5 ~ 9.3,
        EDTA == 100 & Co_conc == 0.1 ~ 0.017,
        EDTA == 100 & Co_conc == 0.01 ~ 0.0017,
        EDTA == 100 & Co_conc == 50 ~ 9.3, #actually 8.34
        TRUE ~ NA_real_),
      trtCo = case_when(
        EDTA == 10 & Co_conc == 0.01 ~ "-Co",
        EDTA == 10 & Co_conc == 5 ~ "+Co",
        EDTA == 100 & Co_conc == 0.01 ~ "--Co",
        EDTA == 100 & Co_conc == 0.1 ~ "-Co",
        EDTA == 100 & Co_conc == 50 ~ "+Co",
        TRUE ~ NA_character_))

growth_rates <- growth_rates %>%
  mutate(species = factor(species, levels = c("Lennoxia faveolata", "Fragilariopsis yeti", 
                                              "Chaetoceros flexuosus", "Chaetoceros neogracile", 
                                              "Proboscia inermis", "Thalassiosira antarctica", 
                                              "Eucampia antarctica", "Synedra sp.", 
                                              "Phaeocystis antarctica (SX9)", 
                                              "Phaeocystis antarctica (PFZ3)")), 
         Zn_trt = factor(trtZn, levels = c("--Zn", "--mZn", "-Zn", "MZn", "+Zn", "Half Zn", "x5 Zn")))
         
Co_plot <- ggplot(growth_rates %>% filter(trtZn == "+Zn", trtCd == "-Cd"), 
       aes(x = as.factor(freeCo), y = growth)) +
  stat_summary(fun = mean, geom = "bar", fill = "grey", color = "black", size = 0.6) +  
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
               geom = "errorbar", width = 0.2) +  
  facet_wrap(~species) +
  ylab(expression(mu~"(d"^-1*")")) +
  xlab(expression("[Co'] (pM)")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), 
        legend.text = element_text(size = 9))

ggsave("SO_Co_plot.svg", Co_plot, width = 6, height = 4.5)

# Create the summary dataframe
growth_summary <- growth_rates %>%
  group_by(species, EDTA, trtCd, trtCo, trtZn, Zn2, Co2, Cd2, freeCd, freeZn, freeCo) %>%
  summarise(n = n(),  # Count the number of observations in each group
            growth_mean = mean(growth, na.rm = TRUE),  # Mean of growth rates
            growth_sd = sd(growth, na.rm = TRUE)) %>%  # Standard deviation of growth rates
  ungroup()

# Add relative growth rate and relative SD to the summary
growth_summary <- growth_summary %>%
  group_by(species) %>%
  mutate(max_growth = max(growth_mean, na.rm = TRUE),  # Calculate max growth for each species
         relative_growth = growth_mean / max_growth,  # Calculate relative growth rate
         relative_sd = growth_sd / max_growth) %>%    # Calculate relative standard deviation
  ungroup()

growth_summary_Zn <- growth_rates %>%
  filter(trtCd == "-Cd") %>%  # Filter for rows where trtCd is "-Cd"
  group_by(species, trtCo, trtZn, freeZn, freeCo, Zn2, Co2) %>% 
  summarise(
    n = n(),  # Count the number of observations in each group
    growth_mean = mean(growth, na.rm = TRUE),  # Mean of growth rates
    growth_sd = sd(growth, na.rm = TRUE)       # Standard deviation of growth rates
  ) %>%
  ungroup()

growth_summary_Zn <- growth_summary_Zn %>%
  group_by(species) %>%
  mutate(max_growth = max(growth_mean, na.rm = TRUE),  # Calculate max growth for each species
         relative_growth = growth_mean / max_growth,  # Calculate relative growth rate
         relative_sd = growth_sd / max_growth) %>%    # Calculate relative standard deviation
  ungroup()


Zn_plot <- ggplot(growth_summary_Zn, aes(x = freeZn, y = relative_growth, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = relative_growth - relative_sd, ymax = relative_growth + relative_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  ylab(expression(mu / mu[max])) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))


ggsave("SO_ZnCo_plot2.svg", Zn_plot, width = 7.2, height = 6)

# Add max observed growth directly to the existing data frame
growth_summary_Zn <- growth_summary_Zn %>%
  group_by(species, trtCo) %>%
  mutate(max_growth_Co = max(relative_growth, na.rm = TRUE)) %>%
  ungroup() 

# Define the Monod model fitting function with fixed max_growth_Co
fit_monod <- function(data) {
  tryCatch({
    nls(
      relative_growth ~ max_growth_Co * freeZn / (Ks + freeZn),  # max_growth_Co is fixed
      data = data,
      start = list(Ks = 1)
    )
  }, error = function(e) NULL)  # Handle fitting errors gracefully
}

# Group data, fit models, and extract results
model_results <- growth_summary_Zn %>%
  group_by(species, trtCo) %>%
  nest() %>%  # Nest the data for each group
  mutate(
    model = map(data, fit_monod),  # Fit the Monod model
    params = map(model, ~ if (!is.null(.)) broom::tidy(.) else tibble()),  # Extract parameters
    predictions = map2(model, data, ~ if (!is.null(.x)) {
      .y %>% mutate(predicted = predict(.x, newdata = .y))
    } else tibble())  # Generate predictions
  )

# Unnest parameters and predictions separately
params_results <- model_results %>%
  dplyr::select(species, trtCo, params) %>%
  unnest(cols = params)

prediction_results <- model_results %>%
  dplyr::select(species, trtCo, predictions) %>%
  unnest(cols = predictions)


# Extract parameters for mu_max
mu_max_results <- params_results %>%
  filter(term == "mu_max") %>%
  dplyr::select(species, trtCo, estimate, std.error, statistic, p.value) %>%
  rename(mu_max = estimate)

# Extract parameters for Ks
Ks_results <- params_results %>%
  filter(term == "Ks") %>%
  dplyr::select(species, trtCo, estimate, std.error, statistic, p.value) %>%
  rename(Ks = estimate)

# Calculate the difference and propagated standard error
Ks_comparison <- Ks_results %>%
  dplyr::select(species, trtCo, Ks, std.error) %>%  # Drop statistic and p.value
  pivot_wider(
    names_from = trtCo, 
    values_from = c(Ks, std.error), 
    names_sep = "_"
  )

Ks_comparison <- Ks_comparison %>%
  mutate(
    delta_Ks = `Ks_+Co` - `Ks_-Co`,  # Difference in Ks
    se_delta_Ks = sqrt(`std.error_+Co`^2 + `std.error_-Co`^2),  # Propagate SE
    z_value = delta_Ks / se_delta_Ks,  # z-value for the difference
    p_value = 2 * pnorm(-abs(z_value))  # Two-tailed p-value
  )

#---- trying to do some ANOVAs -------------
library(rstatix)

games_howell_Zn <- growth_rates %>%
  filter(trtCd == "-Cd", trtCo != "--Co", growth != 0) %>%  # Filter for relevant data and exclude growth = 0
  group_by(species, trtCo) %>%
  filter(n_distinct(trtZn) > 1) %>%  # Ensure at least 2 levels of trtZn
  games_howell_test(growth ~ trtZn)

games_howell_Zn <- growth_rates %>%
  filter(trtCd == "-Cd", trtCo != "--Co", growth != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(growth ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co <- growth_rates %>%
  filter(trtCd == "-Cd", growth != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(growth ~ trtCo)

# View results
games_howell_Co

# ----Cd results ....----------
Cd_plot <- ggplot(growth_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                  aes(x = freeCd, y = `relative_growth`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `relative_growth` - `relative_sd`, ymax = `relative_growth` + `relative_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression(mu / mu[max])) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))


ggsave("SO_Cd_plotXXXX.svg", Cd_plot, width = 7.2, height = 6)

#---- trying to do some ANOVAs -------------
library(rstatix)

anova_results <- growth_rates %>%
  filter(trtCo != "--Co", growth != 0) %>%
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Exclude groups with only 1 level of trtCd
  welch_anova_test(growth ~ trtCd)

anova_results

games_howell_Cd <- growth_rates %>%
  filter(trtCo != "--Co", growth != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>%  # Set hierarchical order
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Ensure at least 2 levels of trtCd
  games_howell_test(growth ~ trtCd)  # Perform Games-Howell test


games_howell_Cd <- growth_rates %>%
  filter(trtCo != "--Co", growth != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(freeCd) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(growth ~ freeCd)  # Perform Games-Howell test

# View results
games_howell_Cd

one_sample_tests <- growth_rates %>%
  filter(trtCo != "--Co") %>%  # Exclude unwanted treatments
  mutate(
    trtCd = factor(trtCd,
                   levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"),
                   ordered = TRUE),
    trtZn = factor(trtZn,
                   levels = c("+Zn", "MZn", "-Zn", "--Zn"),
                   ordered = TRUE)
  ) %>%
  group_by(species, trtCo, trtZn) %>%  
  filter(n_distinct(trtCd) > 1) %>%  # Only consider groups with at least 2 treatments
  group_modify(~ {
    # Within the subgroup (.x), first identify which trtCd levels have all growth = 0
    zero_levels <- .x %>%
      group_by(trtCd) %>%
      summarize(all_zero = all(growth == 0), .groups = "drop") %>%
      filter(all_zero) %>%
      pull(trtCd)
    
    # If there is no treatment that is exclusively zero, do nothing (return an empty tibble)
    if (length(zero_levels) == 0) {
      return(tibble())
    }
    
    # Otherwise, for every treatment that is not exclusively zero, run a one-sample t-test
    # comparing its growth values to 0.
    nonzero_results <- .x %>%
      filter(!(trtCd %in% zero_levels)) %>%
      group_by(trtCd) %>%
      do({
        cur_data <- .
        # Require at least two observations to run a t-test
        if(nrow(cur_data) < 2) {
          tibble(statistic = NA_real_, p_value = NA_real_)
        } else {
          test_res <- t.test(cur_data$growth, mu = 0)
          tibble(statistic = test_res$statistic, p_value = test_res$p.value)
        }
      }) %>%
      ungroup() %>%
      mutate(zero_treatment = paste(zero_levels, collapse = ", "))
    
    return(nonzero_results)
  }) %>%
  ungroup()

# View the one-sample t-test results:
one_sample_tests

growth_summary_Cd <- growth_rates %>%
  filter(!trtZn %in% c("Half Zn", "x5 Zn", "--mZn")) %>% 
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>% # Filter for rows where trtCd is "-Cd"
  group_by(species, trtCo, trtZn, trtCd, freeZn, freeCo, freeCd, Zn2, Co2, Cd2) %>% 
  summarise(
    n = n(),  # Count the number of observations in each group
    growth_mean = mean(growth, na.rm = TRUE),  # Mean of growth rates
    growth_sd = sd(growth, na.rm = TRUE)       # Standard deviation of growth rates
  ) %>%
  ungroup()

growth_summary_Cd <- growth_summary_Cd %>%
  group_by(species) %>%
  mutate(max_growth = max(growth_mean, na.rm = TRUE),  # Calculate max growth for each species
         relative_growth = growth_mean / max_growth,  # Calculate relative growth rate
         relative_sd = growth_sd / max_growth) %>%    # Calculate relative standard deviation
  ungroup()

#------Temperate bugs ---------------
# import the growth rates
growth_rates_temp <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/Zn_Growth_temperate.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(species = as.factor(species),
         EDTA = as.factor(EDTA), Zn = as.factor(Zn), Co = as.factor(Co), Cd = as.factor(Cd),
         well = as.numeric(well),
         transfer = as.numeric(transfer),
         growth = as.numeric(growth))

# make it so that there are numeric values for the concentrations of metals
growth_rates_temp$Zn_conc <- ifelse(growth_rates_temp$Zn == "-", 0.5,
                               ifelse(growth_rates_temp$EDTA == 100 & growth_rates_temp$Zn == "+", 80.2,
                                      ifelse(growth_rates_temp$EDTA == 100 & growth_rates_temp$Zn == "m", 4.5, NA)))

growth_rates_temp$Co_conc <- ifelse(growth_rates_temp$EDTA == 100 & growth_rates_temp$Co == "-", 0.1,
                               ifelse(growth_rates_temp$EDTA == 100 & growth_rates_temp$Co == "m", 1.1,
                                      ifelse(growth_rates_temp$EDTA == 100 & growth_rates_temp$Co == "+", 50.4, NA)))

growth_rates_temp$Cd_conc <- ifelse(growth_rates_temp$Cd == "-", 0.02,
                               ifelse(growth_rates_temp$Cd == "2", 2,
                                      ifelse(growth_rates_temp$Cd == "10", 10,
                                             ifelse(growth_rates_temp$Cd == "50", 50,
                                                    ifelse(growth_rates_temp$Cd == "250", 250,
                                                           ifelse(growth_rates_temp$Cd == "1250", 1250, 
                                                                  ifelse(growth_rates_temp$Cd == "2500", 2500, 
                                                                         ifelse(growth_rates_temp$Cd == "6250", 6250,
                                                                                ifelse(growth_rates_temp$Cd == "12500", 12500,NA)))))))))

# now for the free ion concentrations:
growth_rates_temp <- growth_rates_temp %>% mutate(freeCd = case_when(
  EDTA == 100 & Cd_conc == 0.02 ~ 0.0776,
  EDTA == 100 & Cd_conc == 2 ~ 7.76,
  EDTA == 100 & Cd_conc == 10 ~ 38.8,
  EDTA == 100 & Cd_conc == 50 ~ 194,
  EDTA == 100 & Cd_conc == 250 ~ 973,
  EDTA == 100 & Cd_conc == 1250 ~ 4915,
  EDTA == 100 & Cd_conc == 2500 ~ 9956,
  EDTA == 100 & Cd_conc == 6250 ~ 25900,
  EDTA == 100 & Cd_conc == 12500 ~ 55550,
  TRUE ~ NA_real_),
  Cd2 = case_when(
    EDTA == 100 & Cd_conc == 0.02 ~ 0.0028,
    EDTA == 100 & Cd_conc == 2 ~ 0.28,
    EDTA == 100 & Cd_conc == 10 ~ 1.38,
    EDTA == 100 & Cd_conc == 50 ~ 6.9,
    EDTA == 100 & Cd_conc == 250 ~ 34.6,
    EDTA == 100 & Cd_conc == 1250 ~ 174.5,
    EDTA == 100 & Cd_conc == 2500 ~ 353.6,
    EDTA == 100 & Cd_conc == 6250 ~ 920,
    EDTA == 100 & Cd_conc == 12500 ~ 1973,
    TRUE ~ NA_real_),
  trtCd = case_when(
    EDTA == 100 & Cd_conc == 0.02 ~ "-Cd",
    EDTA == 100 & Cd_conc %in% c(2, 10, 50, 250, 1250, 2500, 6250, 12500) ~ paste0(Cd, "Cd"),
    TRUE ~ NA_character_),
  freeZn = case_when(
    EDTA == 100 & Zn_conc == 0.5 ~ 0.205,
    EDTA == 100 & Zn_conc == 4.5 ~ 1.84,
    EDTA == 100 & Zn_conc == 80.2 ~ 32.9,
    TRUE ~ NA_real_),
  Zn2 = case_when(
    EDTA == 100 & Zn_conc == 0.5 ~ 0.098,
    EDTA == 100 & Zn_conc == 4.5 ~ 0.88,
    EDTA == 100 & Zn_conc == 80.2 ~ 15.7,
    TRUE ~ NA_real_),
  trtZn = case_when(
    EDTA == 100 & Zn_conc == 0.5 ~ "-Zn",
    EDTA == 100 & Zn_conc == 4.5 ~ "MZn",
    EDTA == 100 & Zn_conc == 80.2 ~ "+Zn",
    TRUE ~ NA_character_),
  freeCo = case_when(
    EDTA == 100 & Co_conc == 0.1 ~ 0.019,
    EDTA == 100 & Co_conc == 1.1 ~ 0.207,
    EDTA == 100 & Co_conc == 50.4 ~ 9.5,
    TRUE ~ NA_real_),
  Co2 = case_when(
    EDTA == 100 & Co_conc == 0.1 ~ 0.014,
    EDTA == 100 & Co_conc == 1.1 ~ 0.15,
    EDTA == 100 & Co_conc == 50.4 ~ 7.0,
    TRUE ~ NA_real_),
  trtCo = case_when(
    EDTA == 100 & Co_conc == 0.1 ~ "-Co",
    EDTA == 100 & Co_conc == 1.1 ~ "MCo",
    EDTA == 100 & Co_conc == 50.4 ~ "+Co",
    TRUE ~ NA_character_))

growth_rates_temp <- growth_rates_temp %>%
  mutate(species = factor(species, levels = c("Thalassiosira pseudonana", "Thalassiosira oceanica", 
                                              "Phaeodactylum tricornutum", "Coccolithus braarudii", 
                                              "Gephyrocapsa huxleyi", "Diacronema lutheri", 
                                              "Isochrysis galbana", "Chlamydomonas concordia", 
                                              "Ostreococcus tauri", "Bathycoccus prasinos",
                                              "Thoracosphaera heimii")), 
         Zn_trt = factor(trtZn, levels = c("-Zn", "MZn", "+Zn")),
         Co_trt = factor(trtZn, levels = c("-Co", "MCo", "+Co")))


# Create the summary dataframe
growth_summary_temp <- growth_rates_temp %>%
  group_by(species, EDTA, trtCd, trtCo, trtZn, freeCd, freeZn, freeCo) %>%
  summarise(n = n(),  # Count the number of observations in each group
            growth_mean = mean(growth, na.rm = TRUE),  # Mean of growth rates
            growth_sd = sd(growth, na.rm = TRUE)) %>%  # Standard deviation of growth rates
  ungroup()

# Add relative growth rate and relative SD to the summary
growth_summary_temp <- growth_summary_temp %>%
  group_by(species) %>%
  mutate(max_growth = max(growth_mean, na.rm = TRUE),  # Calculate max growth for each species
         relative_growth = growth_mean / max_growth,  # Calculate relative growth rate
         relative_sd = growth_sd / max_growth) %>%    # Calculate relative standard deviation
  ungroup()

Cd_plot_temp <- ggplot(growth_summary_temp, 
                  aes(x = freeCd, y = `relative_growth`, color = trtZn, fill = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `relative_growth` - `relative_sd`, ymax = `relative_growth` + `relative_sd`), width = 0.3) +
  geom_star(aes(starshape = trtCo, fill = trtCo), size = 3) +
  scale_shape_manual(values = c(22,23,1)) +
  scale_starshape_manual(values = c("-Co" = 14, "MCo" = 15, "+Co" = 28)) +
  scale_colour_manual(values = c("-Zn" = "#D51317FF", "MZn" = "#6A5372FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("-Co" = "white", "MCo" = "grey", "+Co" = "#0094CDFF")) +
  scale_starshape_manual(values = c("lowestCo" = 8, "lowCo" = 5, "midCo" = 7, "hiCo" = 10)) +  # Define star shapes
  scale_x_log10(labels = label_number(), breaks = scales::trans_breaks("log10", function(x) 10^x),) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression(mu / mu[max])) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("Temp_Cd_plot.svg", Cd_plot_temp, width = 7.2, height = 6)


#------Literature bugs ---------------
# import the growth rates
growth_rates_lit <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/Lit_Zn_growth.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Species = as.factor(Species),
         Clone = as.factor(Clone), 
         Location = as.factor(Location), 
         Group = as.factor(Group), 
         Habitat = as.factor(Habitat),
         pZn = as.numeric(pZn), pCo = as.numeric(pCo),
         pmolZn2 = as.numeric(pmolZn2), pmolCo2 = as.numeric(pmolCo2),
         growth = as.numeric(growth),
         sd = as.numeric(sd),
         freeZn = 10^(-pZn)*(10^12), freeCo = 10^(-pZn)*(10^12))

growth_summary_lit <- growth_rates_lit %>%
  group_by(Species, Ref, Clone) %>%
  mutate(max_growth = max(growth, na.rm = TRUE),  # Calculate max growth for each species
         relative_growth = growth / max_growth) %>%    # Calculate relative growth rate
  ungroup()

growth_summary_lit_Zn <- growth_summary_lit %>% filter(Zn_Co_exp == "Zn") %>%
  filter(pmolZn2 <= 20000) %>%
  mutate(trtCo = case_when(
      pmolCo2 >= 0.0003 & pmolCo2 <= 0.003 ~ "lowestCo",
      pmolCo2 >= 0.017 & pmolCo2 <= 0.06 ~ "lowCo",
      pmolCo2 >= 0.98 & pmolCo2 <= 2.52 ~ "midCo",
      pmolCo2 >= 8.9 & pmolCo2 <= 70.8 ~ "hiCo",
      TRUE ~ NA_character_))

str(growth_summary_lit_Zn)

ggplot(growth_summary_lit_Zn %>% filter(Ref == "Saito and Goepfert 2008"), aes(x = pmolZn2, y = relative_growth, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23,5)) +
  #facet_wrap(~Species) +  
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), limits = c(0,1)) +
  ylab(expression(mu / mu[max])) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

# Add max_growth calculation based on Species, Ref, Clone, and trtCo
growth_summary_lit_Zn <- growth_summary_lit_Zn %>%
  group_by(Species, Ref, Clone, trtCo) %>%
  mutate(max_growth = max(relative_growth, na.rm = TRUE)) %>%
  ungroup()  # Ungroup to prevent unintended behavior in later operations

# Define the Monod model fitting function with fixed mu_max
fit_monod <- function(data) {
  tryCatch({
    nls(
      relative_growth ~ max_growth * pmolZn2 / (Ks + pmolZn2),  # max_growth is fixed
      data = data,
      start = list(Ks = 1)
    )
  }, error = function(e) NULL)  # Handle fitting errors gracefully
}

# Group data, fit models, and extract results
model_results <- growth_summary_lit_Zn %>%
  group_by(Species, Ref, Clone, trtCo, Group, Ocean, Habitat) %>%
  nest() %>%  # Nest the data for each group
  mutate(
    model = map(data, fit_monod),  # Fit the Monod model
    params = map(model, ~ if (!is.null(.)) broom::tidy(.) else tibble()),  # Extract parameters
    predictions = map2(model, data, ~ if (!is.null(.x)) {
      .y %>% mutate(predicted = predict(.x, newdata = .y))
    } else tibble())  # Generate predictions
  )

# Unnest parameters and predictions separately
params_results <- model_results %>%
  dplyr::select(Species,Ref,Clone,trtCo,Group,Ocean,Habitat, params) %>%
  unnest(cols = params)

prediction_results <- model_results %>%
  dplyr::select(Species,Ref,Clone,trtCo,Group,Ocean,Habitat, predictions) %>%
  unnest(cols = predictions)

# Extract parameters for mu_max
mu_max_results <- params_results %>%
  filter(term == "mu_max") %>%
  dplyr::select(Species,Ref,Clone,trtCo,Group,Ocean,Habitat, estimate, std.error, statistic, p.value) %>%
  rename(mu_max = estimate)

# Extract parameters for Ks
Ks_results <- params_results %>%
  filter(term == "Ks") %>%
  dplyr::select(Species,Ref,Clone,trtCo,Group,Ocean,Habitat, estimate, std.error, statistic, p.value) %>%
  rename(Ks = estimate)


# edited some data and made some estimates where the model didn't work
Ks_data_all <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/Ks_Zn_data_setmumax.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Species = as.factor(Species),
         Clone = as.factor(Clone), 
         Ref = as.factor(Ref), 
         Group = as.factor(Group), 
         Ks = as.numeric(Ks), trtCo = as.factor(trtCo),
         error = as.numeric(error), Less = as.factor(Less))

# Step 1: Define the desired order for Environment and Groups
environment_order <- c("Temperate neritic", "Gulf of Mexico", "NW Atlantic", 
                       "NE Atlantic", "NE Pacific", "Tasman Sea", 
                       "Subantarctic", "Southern Ocean")

group_order <- c("diatom", "haptophyte", "green algae", "dinoflagellate", "cyanobacteria")

# Step 2: Explicitly order Environment
Ks_data_all <- Ks_data_all %>%
  mutate(
    Environment = factor(Environment, levels = environment_order)  # Explicitly set order
  )

# Step 3: Create Species_Clone and sort data
Ks_data_all <- Ks_data_all %>%
  mutate(
    Species_Clone = paste0(Species, " (", Clone, ")")  # Combine Species and Clone
  ) %>%
  arrange(Environment, factor(Group, levels = group_order), Species_Clone) %>%
  mutate(
    Species_Clone = factor(Species_Clone, levels = unique(Species_Clone))  # Preserve order
  )

# Step 4: Combine Environment with Species_Clone for labels
Ks_data_all <- Ks_data_all %>%
  mutate(
    Full_Label = factor(paste(Environment, Species_Clone, sep = ": "), levels = unique(paste(Environment, Species_Clone, sep = ": ")))
  )

# Reverse the levels for Full_Label
Ks_data_all <- Ks_data_all %>%
  mutate(Full_Label = factor(Full_Label, levels = rev(levels(Full_Label))))

#ensure errors don't drop below zero
Ks_data_all <- Ks_data_all %>%
  mutate(xmin = pmax(Ks - error, 0.0001),  
    xmax = Ks + error)

# Step 5: Create the plot
Ks_plot <- ggplot(Ks_data_all, aes(x = Ks, y = Full_Label, color = Group, fill = trtCo)) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 1) +  # Use corrected xmin and xmax
  geom_star(aes(starshape = trtCo, fill = trtCo), size = 3) +  # Use stars instead of points
  scale_x_log10(labels = scales::label_number()) + 
  coord_cartesian(xlim = c(0.01, 100)) + 
  scale_y_discrete(position = "right") +  # Respect the reversed y-axis order
  scale_colour_manual(values = c("haptophyte" = "#D51317FF", "dinoflagellate" = "#6A5372FF", "diatom" = "#0094CDFF", "green algae" = "#007B3DFF", "cyanobacteria" = "#31B7BCFF")) +
  scale_fill_manual(values = c("lowestCo" = "white", "lowCo" = "lightgrey", "midCo" = "darkgrey", "hiCo" = "black")) +
  scale_starshape_manual(values = c("lowestCo" = 13, "lowCo" = 15, "midCo" = 28, "hiCo" = 1)) +  # Define star shapes
  ylab("Species (Clone) by Environment") +
  xlab(expression("Ks (pM Zn')")) +  # Updated label for x-axis
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9, colour = "black"), axis.title.y = element_blank(),
    legend.title = element_text(size = 9),
    legend.position = "top",
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 9))

Ks_plot

ggsave("KsZn_plot_labs3.svg", Ks_plot, width = 7.2, height = 8.38)


Ks_plot <- ggplot(Ks_data_all, aes(x = Ks, y = Full_Label, color = Group, fill = trtCo)) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 1) +  # Use corrected xmin and xmax
  geom_star(aes(starshape = trtCo, fill = trtCo), size = 3) +  # Use stars instead of points
  scale_x_log10(labels = scales::label_number()) + 
  coord_cartesian(xlim = c(0.01, 100)) + 
  scale_y_discrete(position = "right") +  # Respect the reversed y-axis order
  scale_colour_manual(values = c("haptophyte" = "#D51317FF", "dinoflagellate" = "#6A5372FF", "diatom" = "#0094CDFF", "green algae" = "#007B3DFF", "cyanobacteria" = "#31B7BCFF")) +
  scale_fill_manual(values = c("lowestCo" = "white", "lowCo" = "lightgrey", "midCo" = "darkgrey", "hiCo" = "black")) +
  scale_starshape_manual(values = c("lowestCo" = 13, "lowCo" = 15, "midCo" = 28, "hiCo" = 1)) +  # Define star shapes
  ylab("Species (Clone) by Environment") +
  xlab(expression("Ks (pM Zn')")) +  # Updated label for x-axis
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "blue") +  # Add vertical line at x = 0.85
  geom_vline(xintercept = 6.7, linetype = "dashed", color = "red") +    # Add vertical line at x = 6.7
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9, colour = "black"), axis.title.y = element_blank(),
    legend.title = element_text(size = 9),
    legend.position = "top",
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 9)
  )

Ks_plot


