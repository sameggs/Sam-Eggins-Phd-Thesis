setwd("/Users/eggboy/Dropbox/Science/Data/Rubisco") #setwd

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(lmtest)

Rubisco_percent <- read.csv("/Users/eggboy/Dropbox/Science/Data/Rubisco/Rubisco_percent.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Rubisco = as.numeric(Rubisco)) %>%
  mutate(Species = factor(Species, levels = c("P. tricornutum", "P. inermis", "C. flexuosus", "T. antarctica", "E. antarctica", "tobacco")),
         Temp = as.factor(Temp))
str(Rubisco_percent)

Percent <- ggplot(Rubisco_percent, aes(y = Rubisco, x = Species, colour = Temp)) +
  geom_boxplot() +
  ylab(expression("Rubisco (%TSP)")) +
  xlab(expression("Species")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5), minor_breaks = seq(0, 20, by = 2.5)) +
  scale_colour_manual(values = c("#000000FF","#164194FF")) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9),
        legend.position = "none")


Rubisco_mean <- Rubisco_percent %>%
  group_by(Species) %>%
  summarise(Mean_Rubisco = mean(Rubisco, na.rm = TRUE))

Rubisco_kinetics <- read.csv("/Users/eggboy/Dropbox/Science/Data/Rubisco/Rubisco_ref.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Kcat = as.numeric(Kcat),
         Kc = as.numeric(Kc),
         Sco = as.numeric(Sco)) 

Rubisco_kinetics <- tail(Rubisco_kinetics, 5)

Rubisco_kinetics <- Rubisco_kinetics %>%
  mutate(Species = case_when(
    Species == "Proboscia inermis" ~ "P. inermis",
    Species == "Phaeodactylum tricornutum" ~ "P. tricornutum",
    Species == "Chaeotoceros flexuosus" ~ "C. flexuosus",
    Species == "Eucampia antarctica" ~ "E. antarctica",
    Species == "Thalassiosira antarctica" ~ "T. antarctica",
    TRUE ~ Species  # Keeps other species names unchanged
  ))

print(Rubisco_kinetics)

Rubiscodata <- merge(Rubisco_kinetics, Rubisco_mean, by = "Species", all = TRUE)

KCvsPercent <- ggplot(Rubiscodata, aes(y = Mean_Rubisco, x = Kc)) +
  geom_point(size = 5) +
  xlab(expression(K[c]^"21%O2")) +
  ylab(expression("Rubisco (%TSP)")) +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))

ggsave("PercentRubiscoPlot2.svg", Percent, width = 4, height = 4.3)
ggsave("KCvsPercentRubiscoPlot2.svg", KCvsPercent, width = 4.6, height = 4.6)
