setwd("/Users/eggboy/Dropbox/Science/Data/Rubisco") #setwd

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(lmtest)

Rubisco_temp <- read.csv("/Users/eggboy/Dropbox/Science/Data/Rubisco/TempKinetics.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate_at(vars(kcat_mean:Sco_SD), as.numeric) %>%
  mutate(Species = factor(Species, levels = c("N. tabacum", "P. inermis", "P. tricornutum", "C. flexuosus")))
str(Rubisco_temp)

arrhenius <- read.csv("/Users/eggboy/Dropbox/Science/Data/Rubisco/tempparameters.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate_at(vars(deltaH:c), as.numeric) %>%
  mutate(Species = factor(Species, levels = c("N. tabacum", "P. inermis", "P. tricornutum", "C. flexuosus")))
str(arrhenius)


kcat_nt <- function(Temp) {return(exp(26.74 - (63510/ (8.314 * (Temp + 273.15)))))}
kcat_pi <- function(Temp) {return(exp(28.75 - (69060/ (8.314 * (Temp + 273.15)))))}
kcat_pt <- function(Temp) {return(exp(27.49 - (64990/ (8.314 * (Temp + 273.15)))))}
kcat_cf <- function(Temp) {return(exp(27.15 - (66150/ (8.314 * (Temp + 273.15)))))}

kcat_nt(4)  # For N. tabacum
kcat_pi(4)  # For P. inermis
kcat_pt(4)  # For P. tricornutum
kcat_cf(4)  # For C. flexuosus


O2rates <- data.frame(
  Temp = c(20, 4, 4, 4),
  kcat_mean = c(2.23, 0.461, 0.608, 0.341),
  kcat_SD = c(0.073, 0.030, 0.024, 0.023),
  Species = c("P. tricornutum", "P. tricornutum", "P. inermis", "C. flexuosus"))

kcat_plot <- ggplot(Rubisco_temp, aes(x = Temp, colour = Species)) +
  geom_errorbar(aes(ymin = kcat_mean - kcat_SD, ymax = kcat_mean + kcat_SD), width = 0.8) +
  geom_point(aes(y = kcat_mean), size = 2, shape = 5) +
  stat_function(fun = kcat_nt, color = "#A9A9A9FF") + 
  stat_function(fun = kcat_pi, color = "#000000FF") +
  stat_function(fun = kcat_pt, color = "#164194FF") +
  stat_function(fun = kcat_cf, color = "#D51317FF") +
  ylab(expression(k[cat]^C~(s^-1))) +
  xlab(expression("Temperature (°C)")) + 
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  scale_colour_manual(values = c("#A9A9A9FF","#000000FF","#164194FF", "#D51317FF")) +
  scale_fill_manual(values = c("#A9A9A9FF","#000000FF","#164194FF", "#D51317FF")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) +
  geom_errorbar(data = O2rates, aes(x = Temp, ymin = kcat_mean - kcat_SD, ymax = kcat_mean + kcat_SD, colour = Species), width = 0.9) +
  geom_point(data = O2rates, aes(x = Temp, y = kcat_mean, fill = Species, colour = Species), size = 2, shape = 23) 
  

Kc_nt <- function(Temp) {return(exp(16.66 - (33320/ (8.314 * (Temp + 273.15)))))}
Kc_pi <- function(Temp) {return(exp(21.66 - (46030/ (8.314 * (Temp + 273.15)))))}
Kc_pt <- function(Temp) {return(exp(22.4 - (46020/ (8.314 * (Temp + 273.15)))))}
Kc_cf <- function(Temp) {return(exp(21.67 - (44040/ (8.314 * (Temp + 273.15)))))}

Kc_nt(4)  # For N. tabacum
Kc_pi(4)  # For P. inermis
Kc_pt(4)  # For P. tricornutum
Kc_cf(4)  # For C. flexuosus

Kc_plot <- ggplot(Rubisco_temp, aes(x = Temp)) +
  geom_errorbar(aes(ymin = Kc_mean - Kc_SD, ymax = Kc_mean + Kc_SD, colour = Species), width = 0.9) +
  geom_point(aes(y = Kc_mean, colour = Species), size = 2, shape = 5) +
  stat_function(fun = Kc_nt, color = "#A9A9A9FF") + 
  stat_function(fun = Kc_pi, color = "#000000FF") +
  stat_function(fun = Kc_pt, color = "#164194FF") +
  stat_function(fun = Kc_cf, color = "#D51317FF") +
  ylab(expression(K[c]^"21%O2")) +
  xlab(expression("Temperature (°C)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  scale_colour_manual(values = c("#A9A9A9FF","#000000FF","#164194FF", "#D51317FF")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 
  

Sco_nt <- function(Temp) {return(exp(-6.94 - (-28130/ (8.314 * (Temp + 273.15)))))}
Sco_pi <- function(Temp) {return(exp(-6.89 - (-28940/ (8.314 * (Temp + 273.15)))))}
Sco_pt <- function(Temp) {return(exp(-6.89 - (-28970/ (8.314 * (Temp + 273.15)))))}
Sco_cf <- function(Temp) {return(exp(-6.77 - (-27770/ (8.314 * (Temp + 273.15)))))}

Sco_nt(25)  # For N. tabacum
Sco_pi(25)  # For P. inermis
Sco_pt(20)  # For P. tricornutum
Sco_cf(25)  # For C. flexuosus

Sco_plot <- ggplot(Rubisco_temp, aes(x = Temp, colour = Species)) +
  geom_point(aes(y = Sco_mean), size = 2, shape = 5) +
  geom_errorbar(aes(ymin = Sco_mean - Sco_SD, ymax = Sco_mean + Sco_SD), width = 0.9) +
  ylab(expression(S[C/O]~(mol.mol^-1))) +
  xlab(expression("Temperature (°C)")) + 
  scale_y_continuous(limits = c(50, 250), breaks = seq(0, 250, by = 25)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  scale_colour_manual(values = c("#A9A9A9FF","#000000FF","#164194FF", "#D51317FF")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) +
  stat_function(fun = Sco_nt, color = "#A9A9A9FF") + 
  stat_function(fun = Sco_pi, color = "#000000FF") +
  stat_function(fun = Sco_pt, color = "#164194FF") +
  stat_function(fun = Sco_cf, color = "#D51317FF")

KcatKc_nt <- function(Temp) {return(exp(17.08 - (30416/ (8.314 * (Temp + 273.15)))))}
KcatKc_pi <- function(Temp) {return(exp(14.18 - (23509/ (8.314 * (Temp + 273.15)))))}
KcatKc_pt <- function(Temp) {return(exp(12.01 - (19020/ (8.314 * (Temp + 273.15)))))}
KcatKc_cf <- function(Temp) {return(exp(12.25 - (21789/ (8.314 * (Temp + 273.15)))))}

KcatKc_pi(20)
KcatKc_cf(20)
KcatKc_pt(4)
KcatKc_pt(20)

KcatKc_plot <- ggplot(Rubisco_temp, aes(x = Temp, colour = Species)) +
  ylab(expression(k[cat]^C / K[c]^"21%O2"~(mM^-1 * "." * s^-1))) +
  xlab(expression("Temperature (°C)")) + 
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 25)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  scale_fill_manual(values = c("#A9A9A9FF","#000000FF","#164194FF", "#D51317FF")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) +
  stat_function(fun = KcatKc_nt, color = "#A9A9A9FF") + 
  stat_function(fun = KcatKc_pi, color = "#000000FF") +
  stat_function(fun = KcatKc_pt, color = "#164194FF") +
  stat_function(fun = KcatKc_cf, color = "#D51317FF")


# calculating the constants n that
# Convert Temp to Kelvin
Rubisco_temp$Temp_K <- Rubisco_temp$Temp + 273.15
Rubisco_temp$Inv_Temp <- 1 / Rubisco_temp$Temp_K

# Create a new column for kcat/Kc
Rubisco_temp$kcat_over_Kc <- 1000*Rubisco_temp$kcat_mean / Rubisco_temp$Kc_mean

# List to store results
results <- list()

# Loop over each species and fit the model
for(species in unique(Rubisco_temp$Species)) {
  
  # Filter data for the species
  species_data <- Rubisco_temp %>% filter(Species == species)
  
  # Fit linear model
  model <- lm(log(kcat_over_Kc) ~ Inv_Temp, data = species_data)
  
  # Get coefficients
  coefficients <- coef(model)
  
  # Calculate deltaH and c
  deltaH <- -coefficients["Inv_Temp"] * 8.314
  c <- coefficients["(Intercept)"]
  
  results[[species]] <- list(deltaH = deltaH, c = c, model = model)
}

# Display results
results

#---- putting  all of the plots together ----------------------------------

kinetics_all <- guide_area() / (kcat_plot | Kc_plot) / (Sco_plot | KcatKc_plot) +
  plot_layout(guides = 'collect', nrow = (3), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("Rubisco_kinetics2.svg", kinetics_all, width = 7.2, height = 7)

