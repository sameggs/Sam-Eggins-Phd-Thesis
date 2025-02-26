library(tidyverse)
library(ggplot2)
library(ggsci)
library(scales)
library(multcomp)
library(broom)
library(ggstar)
library(emmeans)

setwd("/Users/eggboy/Dropbox/Science/Data/LIFT/Southern Ocean CdZnCo/R") #setwd

# import the growth rates
LIFT <- read.csv("/Users/eggboy/Dropbox/Science/Data/LIFT/Southern Ocean CdZnCo/R/LIFT_Data.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
mutate(species = as.factor(species),
       EDTA = as.factor(EDTA), Zn = as.factor(Zn), Co = as.factor(Co), Cd = as.factor(Cd),
       rep = as.numeric(rep), wavelength = as.factor(wavelength), generation = as.numeric(generation)) %>%
       mutate(across(Fo:Tau3QA, as.numeric))

# make it so that there are numeric values for the concentrations of metals
LIFT$Zn_conc <- ifelse(LIFT$Zn == "-", 0.4,
                          ifelse(LIFT$EDTA == 10 & LIFT$Zn == "+", 7.9,
                          ifelse(LIFT$EDTA == 10 & LIFT$Zn == "h", 3.95,
                          ifelse(LIFT$EDTA == 10 & LIFT$Zn == "m", 1.2,       
                          ifelse(LIFT$EDTA == 10 & LIFT$Zn == "x5", 39.5,
                          ifelse(LIFT$EDTA == 100 & LIFT$Zn == "m", 1.2,
                          ifelse(LIFT$EDTA == 100 & LIFT$Zn == "+", 79, NA)))))))

LIFT$Co_conc <- ifelse(LIFT$EDTA == 10 & LIFT$Co == "-", 0.01,
                          ifelse(LIFT$EDTA == 100 & LIFT$Co == "-", 0.1,
                          ifelse(LIFT$EDTA == 100 & LIFT$Co == "--", 0.01,
                          ifelse(LIFT$EDTA == 10 & LIFT$Co == "+", 5,
                          ifelse(LIFT$EDTA == 100 & LIFT$Co == "+", 50, NA)))))

LIFT$Cd_conc <- ifelse(LIFT$Cd == "-", 0.002,
                          ifelse(LIFT$Cd == "0.2", 0.2,
                          ifelse(LIFT$Cd == "1", 1,
                          ifelse(LIFT$Cd == "2", 2,
                          ifelse(LIFT$Cd == "10", 10,
                          ifelse(LIFT$Cd == "100", 100, NA))))))

# now for the free ion concentrations:
LIFT <- LIFT %>% mutate(freeCd = case_when(
      EDTA == 10 & Cd_conc == 0.002 ~ 0.05,
      EDTA == 10 & Cd_conc == 0.2 ~ 5,
      EDTA == 10 & Cd_conc == 1 ~ 24.9,
      EDTA == 10 & Cd_conc == 10 ~ 251,
      EDTA == 100 & Cd_conc == 0.002 ~ 0.0051,
      EDTA == 100 & Cd_conc == 2 ~ 5.10,
      EDTA == 100 & Cd_conc == 10 ~ 24.9,
      EDTA == 100 & Cd_conc == 100 ~ 251,
      TRUE ~ NA_real_),
      Cd2 = case_when(
        EDTA == 10 & Cd_conc == 0.002 ~ 0.002,
        EDTA == 10 & Cd_conc == 0.2 ~ 0.19,
        EDTA == 10 & Cd_conc == 1 ~ 0.95,
        EDTA == 10 & Cd_conc == 10 ~ 9.6,
        EDTA == 100 & Cd_conc == 0.002 ~ 0.0002,
        EDTA == 100 & Cd_conc == 2 ~ 0.20,
        EDTA == 100 & Cd_conc == 10 ~ 0.95,
        EDTA == 100 & Cd_conc == 100 ~ 9.6,
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
        EDTA == 10 & Co_conc == 0.01 ~ 0.02,
        EDTA == 10 & Co_conc == 5 ~ 11,
        EDTA == 100 & Co_conc == 0.1 ~ 0.02,
        EDTA == 100 & Co_conc == 0.01 ~ 0.002,
        EDTA == 100 & Co_conc == 50 ~ 11, #actually 11.0
        TRUE ~ NA_real_),
      Co2 = case_when(
        EDTA == 10 & Co_conc == 0.01 ~ 0.017,
        EDTA == 10 & Co_conc == 5 ~ 8.33,
        EDTA == 100 & Co_conc == 0.1 ~ 0.017,
        EDTA == 100 & Co_conc == 0.01 ~ 0.0017,
        EDTA == 100 & Co_conc == 50 ~ 8.33, #actually 8.34
        TRUE ~ NA_real_),
      trtCo = case_when(
        EDTA == 10 & Co_conc == 0.01 ~ "-Co",
        EDTA == 10 & Co_conc == 5 ~ "+Co",
        EDTA == 100 & Co_conc == 0.01 ~ "--Co",
        EDTA == 100 & Co_conc == 0.1 ~ "-Co",
        EDTA == 100 & Co_conc == 50 ~ "+Co",
        TRUE ~ NA_character_))

LIFT <- LIFT %>%
  mutate(species = factor(species, levels = c("Lennoxia faveolata", "Fragilariopsis yeti", 
                                              "Chaetoceros flexuosus", "Chaetoceros neogracile", 
                                              "Proboscia inermis", "Thalassiosira antarctica", 
                                              "Eucampia antarctica", "Synedra sp.", 
                                              "Phaeocystis antarctica (SX9)", 
                                              "Phaeocystis antarctica (PFZ3)")), 
         Zn_trt = factor(trtZn, levels = c("--Zn", "-Zn", "MZn", "+Zn")))

# Create the summary dataframe
LIFT_summary <- LIFT %>%
  group_by(species, EDTA, trtCd, trtCo, trtZn, Zn2, Co2, Cd2, freeCd, freeZn, freeCo) %>%
  summarise(n = n(),  # Count the number of observations in each group
            FvFm_mean = mean(FvFm, na.rm = TRUE),  # Mean of growth rates
            FvFm_sd = sd(FvFm, na.rm = TRUE),
            Sig_mean = mean(Sig, na.rm = TRUE),  # Mean of growth rates
            Sig_sd = sd(Sig, na.rm = TRUE),
            Tau1_mean = mean(Tau1QA, na.rm = TRUE),
            Tau1_sd = sd(Tau1QA, na.rm = TRUE),
            Tau2_mean = mean(Tau2QA, na.rm = TRUE),
            Tau2_sd = sd(Tau2QA, na.rm = TRUE),
            Tau3_mean = mean(Tau3QA, na.rm = TRUE),
            Tau3_sd = sd(Tau3QA, na.rm = TRUE)) %>%  # Standard deviation of growth rates
  ungroup()

LIFT_summary_Zn <- LIFT %>%
  filter(trtCd == "-Cd") %>%  # Filter for rows where trtCd is "-Cd"
  group_by(species, trtCo, trtZn, freeZn, freeCo, Zn2, Co2) %>% 
  summarise(
    n = n(), 
    FvFm_mean = mean(FvFm, na.rm = TRUE),  
    FvFm_sd = sd(FvFm, na.rm = TRUE),
    Sig_mean = mean(Sig, na.rm = TRUE),  
    Sig_sd = sd(Sig, na.rm = TRUE),
    Tau1_mean = mean(Tau1QA, na.rm = TRUE),
    Tau1_sd = sd(Tau1QA, na.rm = TRUE),
    Tau2_mean = mean(Tau2QA, na.rm = TRUE),
    Tau2_sd = sd(Tau2QA, na.rm = TRUE),
    Tau3_mean = mean(Tau3QA, na.rm = TRUE),
    Tau3_sd = sd(Tau3QA, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(
    FvFm_max = max(FvFm_mean, na.rm = TRUE), # Find the max FvFm for each species
    relative_FvFm = FvFm_mean / FvFm_max     # Calculate relative FvFm
  ) %>%
  ungroup()


Zn_plot_FvFm <- ggplot(LIFT_summary_Zn, aes(x = freeZn, y = FvFm_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = FvFm_mean - FvFm_sd, ymax = FvFm_mean + FvFm_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  #scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  ylab(expression("Fv/Fm")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Zn_plot_Sigma <- ggplot(LIFT_summary_Zn, aes(x = freeZn, y = Sig_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = Sig_mean - Sig_sd, ymax = Sig_mean + Sig_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  #scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  ylab(expression("Sigma PSII")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("SO_ZnCo_FvFm.svg", Zn_plot_FvFm, width = 7.2, height = 6)
ggsave("SO_ZnCo_Sigma.svg", Zn_plot_Sigma, width = 7.2, height = 6)

#ggsave("SO_ZnCo_FvFm.svg", Zn_plot, width = 7.2, height = 6)

# Taus ----------
Zn_plot_Tau1 <- ggplot(LIFT_summary_Zn, aes(x = freeZn, y = Tau1_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = Tau1_mean - Tau1_sd, ymax = Tau1_mean + Tau1_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  #scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  ylab(expression("Fv/Fm")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Zn_plot_Tau2 <- ggplot(LIFT_summary_Zn, aes(x = freeZn, y = Tau2_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = Tau2_mean - Tau2_sd, ymax = Tau2_mean + Tau2_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, by = 5000)) +
  ylab(expression("Fv/Fm")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Zn_plot_Tau3 <- ggplot(LIFT_summary_Zn, aes(x = freeZn, y = Tau3_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22,1,23)) +
  geom_errorbar(aes(ymin = Tau3_mean - Tau3_sd, ymax = Tau3_mean + Tau3_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 10000)) +
  ylab(expression("Fv/Fm")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("Zn_plot_Tau1.svg", Zn_plot_Tau1, width = 7.2, height = 6)
ggsave("Zn_plot_Tau2.svg", Zn_plot_Tau2, width = 7.2, height = 5.6)
ggsave("Zn_plot_Tau3.svg", Zn_plot_Tau3, width = 7.2, height = 5.6)


#---- trying to do some ANOVAs -------------
library(rstatix)

games_howell_Zn_FvFm <- LIFT %>%
  filter(trtCd == "-Cd", trtCo != "--Co", FvFm != 0) %>%  # Filter for relevant data and exclude 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(FvFm ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_FvFm

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_FvFm <- LIFT %>%
  filter(trtCd == "-Cd", FvFm != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(FvFm ~ trtCo)

# View results
games_howell_Co_FvFm

# now for Sigma
games_howell_Zn_Sigma <- LIFT %>%
  filter(trtCd == "-Cd", trtCo != "--Co", Sig != 0) %>%  # Filter for relevant data and exclude 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(Sig ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_Sigma

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_Sigma <- LIFT %>%
  filter(trtCd == "-Cd", Sig != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(Sig ~ trtCo)

# View results
games_howell_Co_Sigma

# now for Sigma
games_howell_Zn_Tau3 <- LIFT %>%
  filter(trtCd == "-Cd", trtCo != "--Co", Tau3QA != 0) %>%  # Filter for relevant data and exclude 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(Tau3QA ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_Tau3

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_Tau3 <- LIFT %>%
  filter(trtCd == "-Cd", trtCo != "--Co", !is.na(Tau1QA), Tau3QA != 0) %>%  # Remove NA values and exclude 0
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(Tau3QA ~ trtCo)

# View results
games_howell_Co_Tau3

# ----- Cd ---------------------------------
# Cd results ....
Cd_plot_FvFm <- ggplot(LIFT_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                       aes(x = freeCd, y = `FvFm_mean`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `FvFm_mean` - `FvFm_sd`, ymax = `FvFm_mean` + `FvFm_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression("Fv/Fm")) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Cd_plot_Sigma <- ggplot(LIFT_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                        aes(x = freeCd, y = `Sig_mean`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `Sig_mean` - `Sig_sd`, ymax = `Sig_mean` + `Sig_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression("Sigma PSII")) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("SO_Cd_plot_FvFm.svg", Cd_plot_FvFm, width = 7.2, height = 6)
ggsave("SO_Cd_plot_Sigma.svg", Cd_plot_Sigma, width = 7.2, height = 6)

Cd_plot_Tau1 <- ggplot(LIFT_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                       aes(x = freeCd, y = `Tau1_mean`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `Tau1_mean` - `Tau1_sd`, ymax = `Tau1_mean` + `Tau1_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression("Tau1")) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Cd_plot_Tau2 <- ggplot(LIFT_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                       aes(x = freeCd, y = `Tau2_mean`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `Tau2_mean` - `Tau2_sd`, ymax = `Tau2_mean` + `Tau2_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 14000, by = 1000)) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression("tau2")) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Cd_plot_Tau3 <- ggplot(LIFT_summary %>% filter(trtCo != "--Co", trtZn != "Half Zn", trtZn != "x5 Zn", !(trtZn == "+Zn" & trtCo == "-Co")), 
                       aes(x = freeCd, y = `Tau3_mean`, color = trtZn, fill = trtZn, shape = trtCo)) +
  geom_line(aes(group = interaction(trtZn, trtCo)), size = 0.6) +
  geom_errorbar(aes(ymin = `Tau3_mean` - `Tau3_sd`, ymax = `Tau3_mean` + `Tau3_sd`), width = 0.3) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,23)) +
  scale_colour_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_fill_manual(values = c("--Zn" = "#D51317FF", "-Zn" = "#8E3E54FF", "MZn" = "#476990FF", "+Zn" = "#0094CDFF")) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(limits = c(0, 60000), breaks = seq(0, 60000, by = 10000)) +
  facet_wrap(~species, nrow = 3, ncol = 4) +  
  ylab(expression("tau3")) +
  xlab(expression("[Cd'] (pmol L"^-1*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("Cd_plot_Tau1.svg", Cd_plot_Tau1, width = 7.2, height = 6)
ggsave("Cd_plot_Tau2.svg", Cd_plot_Tau2, width = 7.2, height = 5.6)
ggsave("Cd_plot_Tau3.svg", Cd_plot_Tau3, width = 7.2, height = 5.6)

#---- trying to do some ANOVAs -------------
library(rstatix)

anova_results_FvFm <- LIFT %>%
  filter(trtCo != "--Co", FvFm != 0) %>%
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Exclude groups with only 1 level of trtCd
  welch_anova_test(FvFm ~ trtCd)

anova_results_FvFm

games_howell_FvFm_Cd <- LIFT %>%
  filter(trtCo != "--Co", FvFm != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>%  # Set hierarchical order
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Ensure at least 2 levels of trtCd
  games_howell_test(FvFm ~ trtCd)  # Perform Games-Howell test

# View results
games_howell_FvFm_Cd

games_howell_Sig_Cd <- LIFT %>%
  filter(trtCo != "--Co", Sig != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>%  # Set hierarchical order
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Ensure at least 2 levels of trtCd
  games_howell_test(Sig ~ trtCd)  # Perform Games-Howell test

# View results
games_howell_Sig_Cd

games_howell_Tau3QA_Cd <- LIFT %>%
  filter(trtCo != "--Co", Tau3QA != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>%  # Set hierarchical order
  group_by(species, trtCo, trtZn) %>%
  filter(n_distinct(trtCd) > 1) %>%  # Ensure at least 2 levels of trtCd
  games_howell_test(Tau3QA ~ trtCd)  # Perform Games-Howell test

# View results
games_howell_Tau3QA_Cd

LIFT_summary_Cd <- LIFT %>%
  filter(!trtZn %in% c("Half Zn", "x5 Zn", "--mZn")) %>% 
  mutate(
    trtCd = factor(trtCd, levels = c("-Cd", "0.2Cd", "1Cd", "2Cd", "10Cd", "100Cd"), ordered = TRUE),
    trtZn = factor(trtZn, levels = c("+Zn", "MZn", "-Zn", "--Zn"), ordered = TRUE)
  ) %>% # Filter for rows where trtCd is "-Cd"
  group_by(species, trtCo, trtZn, trtCd, freeZn, freeCo, freeCd, Zn2, Co2, Cd2) %>% 
  summarise(
    n = n(), 
    FvFm_mean = mean(FvFm, na.rm = TRUE),  
    FvFm_sd = sd(FvFm, na.rm = TRUE),
    Sig_mean = mean(Sig, na.rm = TRUE),  
    Sig_sd = sd(Sig, na.rm = TRUE),
    Tau1_mean = mean(Tau1QA, na.rm = TRUE),
    Tau1_sd = sd(Tau1QA, na.rm = TRUE),
    Tau2_mean = mean(Tau2QA, na.rm = TRUE),
    Tau2_sd = sd(Tau2QA, na.rm = TRUE),
    Tau3_mean = mean(Tau3QA, na.rm = TRUE),
    Tau3_sd = sd(Tau3QA, na.rm = TRUE)) %>%
  ungroup()



# -------comparing Fv/Fm and growth rate...-------------------
combined_summary_Zn <- LIFT_summary_Zn %>%
  inner_join(growth_summary_Zn, 
             by = c("species", "trtCo", "trtZn"))

# Calculate line equations and R-squared values
line_eqn_data <- combined_summary_Zn %>%
  group_by(species) %>%
  summarise(
    slope = coef(lm(FvFm_mean ~ relative_growth))[2],
    intercept = coef(lm(FvFm_mean ~ relative_growth))[1],
    r_squared = summary(lm(FvFm_mean ~ relative_growth))$r.squared
  ) %>%
  mutate(
    label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2), 
                   "\nR² = ", round(r_squared, 2))
  )

# Scatter plot of relative Fv/Fm vs. relative growth
FvFm_plot_growthrate <- ggplot(combined_summary_Zn, aes(x = relative_growth, y = FvFm_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(22,1,23)) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
  geom_smooth(size = 0.6, method = "lm", se = FALSE, aes(group = species), linetype = "dashed", colour = "black") +  # Add a line for each species
  facet_wrap(~species, nrow = 3, ncol = 4) +
  ylab(expression("Fv/Fm")) +
  xlab(expression(mu / mu[max])) +
  geom_text(data = line_eqn_data, 
            aes(x = 0.3, y = 0.4, 
                label = label), 
            inherit.aes = FALSE, size = 3) + # Add line equations and R²
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ggsave("FvFm.svg", FvFm_plot_growthrate, width = 7.2, height = 6)

# now with literally all of the data...
combined_summary <- LIFT_summary %>%
  inner_join(growth_summary, 
             by = c("species", "trtCo", "trtZn", "trtCd", "EDTA")) %>%
  filter(relative_growth != 0) %>%
  mutate(freeCd.x = factor(freeCd.x))

# Calculate line equations and R-squared values
line_eqn_data <- combined_summary %>%
  group_by(species) %>%
  summarise(
    slope = coef(lm(FvFm_mean ~ relative_growth))[2],
    intercept = coef(lm(FvFm_mean ~ relative_growth))[1],
    r_squared = summary(lm(FvFm_mean ~ relative_growth))$r.squared
  ) %>%
  mutate(
    label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2), 
                   "\nR² = ", round(r_squared, 2))
  )

# Scatter plot of relative Fv/Fm vs. relative growth
FvFm_plot_growthrate_Cd <- ggplot(combined_summary, aes(x = relative_growth, y = FvFm_mean, color = trtCo, fill = trtZn, starshape = freeCd.x)) +
  geom_star(aes(starshape = freeCd.x), size = 3) + 
  scale_starshape_manual(values = c("0.0051" = 13, "0.05" = 15, "5.1" = 28, "24.9" = 14, "251" = 1)) + 
  scale_colour_manual(values = c("--Co" = "white","-Co" = "grey","+Co" = "black")) +
  scale_fill_manual(values = c("--Zn" = "#E69F00FF","--mZn" = "#D51317FF","-Zn" = "#8E3E54FF","MZn" = "#476990FF","Half Zn" = "#0094CDFF","+Zn" = "#1B9E77FF","x5 Zn" = "#00441BFF")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1)) +
  geom_smooth(size = 0.6, method = "lm", se = FALSE, aes(group = species), linetype = "dashed", colour = "black") +  # Add a line for each species
  facet_wrap(~species, nrow = 3, ncol = 4) +
  ylab(expression("Fv/Fm")) +
  xlab(expression(mu / mu[max])) +
  geom_text(data = line_eqn_data, 
            aes(x = 0.3, y = 0.4, 
                label = label), 
            inherit.aes = FALSE, size = 3) + # Add line equations and R²
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1), starshape = guide_legend(nrow = 1))

ggsave("FvFm_growth_plot.svg", FvFm_plot_growthrate_Cd, width = 7.2, height = 6)

# checking that its colour blind friendly
install.packages("colorspace")
library(colorspace)

# Define your color palette
colors <- c("#E69F00","#D51317","#8E3E54","#476990","#0094CD","#1B9E77","#00441B")

# Visualize the original palette
swatchplot(colors)
swatchplot(deutan(colors))
swatchplot(protan(colors))
swatchplot(tritan(colors))

FvFm_sigma <- ggplot(combined_summary, aes(x = Sig_mean, y = FvFm_mean, color = trtCo, fill = trtZn, starshape = freeCd.x)) +
  geom_star(aes(starshape = freeCd.x), size = 3) + 
  scale_starshape_manual(values = c("0.0051" = 13, "0.05" = 15, "5.1" = 28, "24.9" = 14, "251" = 1)) + 
  scale_colour_manual(values = c("--Co" = "white","-Co" = "grey","+Co" = "black")) +
  scale_fill_manual(values = c("--Zn" = "#E69F00FF","--mZn" = "#D51317FF","-Zn" = "#8E3E54FF","MZn" = "#476990FF","Half Zn" = "#0094CDFF","+Zn" = "#1B9E77FF","x5 Zn" = "#00441BFF")) +
  facet_wrap(~species, nrow = 3, ncol = 4) +
  geom_smooth(size = 0.6, method = "lm", se = FALSE, aes(group = species), linetype = "dashed", colour = "black") +  # Add a line for each species
  ylab(expression("Fv/Fm")) +
  xlab(expression(mu / mu[max])) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top", legend.box = "horizontal",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1), starshape = guide_legend(nrow = 1))

