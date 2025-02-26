library(tidyverse)
library(ggplot2)
library(ggsci)
library(scales)
library(multcomp)
library(broom)
library(ggstar)
library(patchwork)
library(multcomp)

setwd("/Users/eggboy/Dropbox/Science/Data/Cell Sizes") #setwd

# import the growth rates
cells <- read.csv("/Users/eggboy/Dropbox/Science/Data/Cell Sizes/Zn_Chla_CV_SA.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
mutate(species = as.factor(species),
       EDTA = as.factor(EDTA), Zn = as.factor(Zn), Co = as.factor(Co), Cd = as.factor(Cd),
       rep = as.numeric(rep),
       CV = as.numeric(CV),
       SAV = as.numeric(SAV),
       Chla_cell = as.numeric(Chla_cell),
       Chla_CV = as.numeric(Chla_CV),
       Colonies = as.numeric(Colonies),
       Col_size = as.numeric(Col_size))

# make it so that there are numeric values for the concentrations of metals
cells$Zn_conc <- ifelse(cells$Zn == "-", 0.4,
                          ifelse(cells$EDTA == 10 & cells$Zn == "+", 7.9,
                          ifelse(cells$EDTA == 10 & cells$Zn == "h", 3.95,
                          ifelse(cells$EDTA == 10 & cells$Zn == "m", 1.2,       
                          ifelse(cells$EDTA == 10 & cells$Zn == "x5", 39.5,
                          ifelse(cells$EDTA == 100 & cells$Zn == "m", 1.2,
                          ifelse(cells$EDTA == 100 & cells$Zn == "+", 79, NA)))))))

cells$Co_conc <- ifelse(cells$EDTA == 10 & cells$Co == "-", 0.01,
                          ifelse(cells$EDTA == 100 & cells$Co == "-", 0.1,
                          ifelse(cells$EDTA == 100 & cells$Co == "--", 0.01,
                          ifelse(cells$EDTA == 10 & cells$Co == "+", 5,
                          ifelse(cells$EDTA == 100 & cells$Co == "+", 50, NA)))))

cells$Cd_conc <- ifelse(cells$Cd == "-", 0.002)

# now for the free ion concentrations:
cells <- cells %>% mutate(freeCd = case_when(
      EDTA == 10 & Cd_conc == 0.002 ~ 0.0055,
      EDTA == 100 & Cd_conc == 0.002 ~ 0.0051,
      TRUE ~ NA_real_),
      Cd2 = case_when(
        EDTA == 10 & Cd_conc == 0.002 ~ 0.0021,
        EDTA == 100 & Cd_conc == 0.002 ~ 0.00019,
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

cells <- cells %>%
  mutate(species = factor(species, levels = c("Lennoxia faveolata", "Fragilariopsis yeti", 
                                              "Chaetoceros flexuosus", "Chaetoceros neogracilis", 
                                              "Proboscia inermis", "Thalassiosira antarctica", 
                                              "Eucampia antarctica", "Synedra sp.", 
                                              "Phaeocystis antarctica (SX9)", 
                                              "Phaeocystis antarctica (PFZ3)")), 
         Zn_trt = factor(trtZn, levels = c("--Zn", "--mZn","-Zn", "MZn", "+Zn")))
         

# Create the summary dataframe
cells_summary <- cells %>%
  group_by(species, EDTA, trtCd, trtCo, trtZn, Zn2, Co2, Cd2, freeCd, freeZn, freeCo) %>%
  summarise(n = n(),  # Count the number of observations in each group
            CV_mean = mean(CV, na.rm = TRUE),  
            CV_sd = sd(CV, na.rm = TRUE),
            SAV_mean = mean(SAV, na.rm = TRUE),  
            SAV_sd = sd(SAV, na.rm = TRUE),
            Chla_cell_mean = mean(Chla_cell, na.rm = TRUE),  
            Chla_cell_sd = sd(Chla_cell, na.rm = TRUE),
            Chla_CV_mean = mean(Chla_CV, na.rm = TRUE), 
            Chla_CV_sd = sd(Chla_CV, na.rm = TRUE),
            Colonies_mean = mean(Colonies, na.rm = TRUE),  
            Colonies_sd = sd(Colonies, na.rm = TRUE),
            Col_size_mean = mean(Col_size, na.rm = TRUE),  
            Col_size_sd = sd(Col_size, na.rm = TRUE),
            ) %>%  # Standard deviation of growth rates
  ungroup()

Zn_plot_SAV <- ggplot(cells_summary, aes(x = freeZn, y = SAV_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22, 1, 23)) +
  geom_errorbar(aes(ymin = SAV_mean - SAV_sd, ymax = SAV_mean + SAV_sd), width = 0.2) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4, scales = "free_y") +  
  scale_x_log10(labels = label_number()) +
  ylab(expression("SA:Volume (um^3)")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_y_continuous(expand = c(0, 0.05), labels = scales::label_number(accuracy = 0.01)) +  # 5% padding
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))


Zn_plot_CV <- ggplot(cells_summary, aes(x = freeZn, y = CV_mean, color = trtCo, fill = trtCo, starshape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  # Connect points within each group
  geom_star(aes(size = trtCo)) +  # Plot the mean as points
  scale_starshape_manual(values = c(21,15,14)) +
  scale_size_manual(values = c("--Co" = 1.5, "-Co" = 1.5, "+Co" = 2)) +
  geom_errorbar(aes(ymin = CV_mean - CV_sd, ymax = CV_mean + CV_sd), width = 0) +  # Add error bars
  facet_wrap(~species, nrow = 3, ncol = 4, scales = "free_y") +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), limits = c(0, NA)) +
  scale_x_log10(labels = label_number()) +
  annotation_logticks(sides = "t", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) + # Add downward log ticks
  ylab(expression("Cell Volume (um^3)")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "#D51317FF", "-Co" = "#D51317FF", "+Co" = "#D51317FF")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "#D51317FF")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))


Zn_plot_ChlaCV <- ggplot(cells_summary, aes(x = freeZn, y = Chla_CV_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  
  geom_errorbar(aes(ymin = Chla_CV_mean - Chla_CV_sd, ymax = Chla_CV_mean + Chla_CV_sd), width = 0.2) +  # Add error bars
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22, 1, 23)) +
  facet_wrap(~species, nrow = 3, ncol = 4, scales = "free_y") +  
  scale_x_log10(labels = label_number()) +
  ylab(expression("Chla")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(0, NA), labels = scales::label_number(accuracy = 0.1)) + # Set y-axis limits to start at 0
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

Zn_plot_Chla_cell <- ggplot(cells_summary, aes(x = freeZn, y = Chla_cell_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  
  geom_errorbar(aes(ymin = Chla_cell_mean - Chla_cell_sd, ymax = Chla_cell_mean + Chla_cell_sd), width = 0.2) +  # Add error bars
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(22, 1, 23)) +
  facet_wrap(~species, nrow = 3, ncol = 4, scales = "free_y") +  
  scale_x_log10(labels = label_number()) +
  ylab(expression("Chla per cell")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + # Set y-axis limits to start at 0
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

ZnCVChla <- Zn_plot_CV | Zn_plot_SAV | Zn_plot_ChlaCV 


ggsave("SO_ZnCo_CV.svg", Zn_plot_CV, width = 7.3, height = 5.6)
ggsave("SO_ZnCo_SAV.svg", Zn_plot_SAV, width = 7.4, height = 5.6)
ggsave("SO_ZnCo_ChlaCV.svg", Zn_plot_ChlaCV, width = 7.2, height = 5.6)


SX9_colonies <- ggplot(cells_summary %>% filter(species == "Phaeocystis antarctica (SX9)"), aes(x = Zn2, y = Colonies_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  
  geom_errorbar(aes(ymin = Colonies_mean - Colonies_sd, ymax = Colonies_mean + Colonies_sd), width = 0.2) +  # Add error bars
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(1, 23)) +
  scale_x_log10(labels = label_number()) +
  ylab(expression("Colonial fraction (%)")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_y_continuous(limits = c(0, 50)) + # Set y-axis limits to start at 0
  annotation_logticks(sides = "t", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) + # Add downward log ticks
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

SX9_colony_size <- ggplot(cells_summary %>% filter(species == "Phaeocystis antarctica (SX9)"), aes(x = Zn2, y = Col_size_mean, color = trtCo, fill = trtCo, shape = trtCo)) +
  geom_line(aes(group = trtCo), size = 0.6) +  
  geom_errorbar(aes(ymin = Col_size_mean - Col_size_sd, ymax = Col_size_mean + Col_size_sd), width = 0.2) +  # Add error bars
  geom_point(size = 3) +  # Plot the mean as points
  scale_shape_manual(values = c(1, 23)) +
  scale_x_log10(labels = label_number()) +
  ylab(expression("Colony Volume")) +
  xlab(expression("[Zn'] (pmol L"^-1*")")) +
  scale_colour_manual(values = c("--Co" = "black", "-Co" = "black", "+Co" = "black")) +
  scale_fill_manual(values = c("--Co" = "white", "-Co" = "white", "+Co" = "black")) +
  scale_y_continuous(limits = c(0, NA)) + # Set y-axis limits to start at 0
  annotation_logticks(sides = "t", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) + # Add downward log ticks
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(size = 9, colour = "black"),
        legend.title = element_text(size = 9), legend.position = "top",
        legend.text = element_text(size = 9), strip.text = element_text(size = 9))

SX9_plot <- SX9_colonies | SX9_colony_size

ggsave("SO_SX9_colonies.svg", SX9_plot, width = 4.18, height = 2.49)

#---- statistics ----------------
#---- trying to do some ANOVAs -------------
library(rstatix)

games_howell_Zn_CV <- cells %>%
  filter(trtCd == "-Cd", trtCo != "--Co", CV != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(CV ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_CV

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_CV <- cells %>%
  filter(trtCd == "-Cd", CV != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(CV ~ trtCo)

# View results
games_howell_Co_CV

# now for SA/V
games_howell_Zn_SAV <- cells %>%
  filter(trtCd == "-Cd", trtCo != "--Co", SAV != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(SAV ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_SAV

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_SAV <- cells %>%
  filter(trtCd == "-Cd", SAV != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(SAV ~ trtCo)

# View results
games_howell_Co_SAV

# Chlorophyll
games_howell_Zn_Chla_CV <- cells %>%
  filter(trtCd == "-Cd", trtCo != "--Co", CV != 0) %>%  # Filter for relevant data and exclude growth = 0
  mutate(freeZn = as.factor(freeZn)) %>%  # Convert freeZn to a factor
  group_by(species, trtCo) %>%
  filter(n_distinct(freeZn) > 1) %>%  # Ensure at least 2 levels of freeZn
  games_howell_test(Chla_CV ~ freeZn)  # Perform Games-Howell test

# View results
games_howell_Zn_Chla_CV

# Perform Games-Howell test for trtCo at fixed Zn levels
games_howell_Co_Chla_CV <- cells %>%
  filter(trtCd == "-Cd", CV != 0) %>%
  group_by(species, freeZn) %>%
  filter(n_distinct(trtCo) > 1) %>%  # Ensure there are at least two levels of trtCo
  games_howell_test(Chla_CV ~ trtCo)

# View results
games_howell_Co_Chla_CV

# ----- plotting Chla/cell as function of CV
CV_plot <- ggplot(cells_summary, aes(x = CV_mean, y = Chla_cell_mean, colour = species, fill = species, starshape = species)) +
  geom_errorbar(aes(ymin = Chla_cell_mean - Chla_cell_sd, ymax = Chla_cell_mean + Chla_cell_sd), width = 0.02) +  # Add error bars for y
  geom_errorbarh(aes(xmin = CV_mean - CV_sd, xmax = CV_mean + CV_sd), height = 0.02) +  # Add error bars for x
  geom_star(size = 3) +  # Add points for means
  scale_x_log10() +  # Log scale for x-axis (CV)
  scale_y_log10() +  # Log scale for y-axis (Chla_cell_mean)
  scale_starshape_manual(values = c(16,13,1,14,12,15,6,25,4,3)) +  # Different shapes for Co treatments
  scale_fill_frontiers() +   
  scale_colour_frontiers() +  # Color and fill for species
  annotation_logticks(sides = "lb", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +  # Log ticks on left and bottom
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 9, colour = "black"),
    legend.title = element_text(size = 9), 
    legend.text = element_text(size = 9), 
    legend.position = "top",  # Place legend at the top
    strip.text = element_text(size = 9)
  )

SAV_plot <- ggplot(cells_summary, aes(x = SAV_mean, y = Chla_cell_mean, colour = species, fill = species, starshape = species)) +
  geom_errorbar(aes(ymin = Chla_cell_mean - Chla_cell_sd, ymax = Chla_cell_mean + Chla_cell_sd), width = 0.02) +  # Add error bars for y
  geom_errorbarh(aes(xmin = SAV_mean - SAV_sd, xmax = SAV_mean + SAV_sd), height = 0.02) +  # Add error bars for x
  geom_star(size = 3) +  # Add points for means
  labs(x = "Surface Area-to-Volume Ratio (SA:V, log scale)", 
       y = "Chlorophyll a per Cell (Chla_cell mean, log scale)") +
  scale_x_log10() +  # Log scale for x-axis
  scale_y_log10() +  # Log scale for y-axis
  scale_starshape_manual(values = c(16,13,1,14,12,15,6,25,4,3)) +  # Different shapes for Co treatments
  scale_fill_frontiers() +  
  scale_colour_frontiers() +# Fill for Co treatments
  annotation_logticks(sides = "lb", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) +  # Log ticks on left and bottom
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(), 
    axis.text = element_text(size = 9, colour = "black"),
    legend.title = element_text(size = 9), 
    legend.text = element_text(size = 9), 
    legend.position = "none",  # Place legend at the top
    strip.text = element_text(size = 9)
  )

Chla_size_plot <- (CV_plot | SAV_plot) + 
  plot_layout(guides = "auto") + 
  theme(legend.position = "top")



