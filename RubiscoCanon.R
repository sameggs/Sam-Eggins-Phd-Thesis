setwd("/Users/eggboy/Dropbox/Science/Data/Rubisco") #setwd

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(patchwork)

Rubisco <- read.csv("/Users/eggboy/Dropbox/Science/Data/Rubisco/Rubisco_ref.csv", fileEncoding="UTF-8-BOM", header = TRUE, na.strings = "") %>%
  mutate(Kcat = as.numeric(Kcat),
         Kc_air = as.numeric(Kc_air),
         Sco = as.numeric(Sco),
         Content = as.numeric(Content)) %>%
  mutate(Kcat_over_Kc_air = 1000* Kcat / Kc_air) %>%
  mutate(Kcat_over_Kc = 1000* Kcat / Kc)

str(Rubisco)


# Linear model for Lineage B
model_B <- lm(Kc_air ~ Kcat, data = subset(Rubisco, Lineage == "B" & Aquatic == "terrestrial"))
summary_B <- summary(model_B)
r_squared_B <- summary_B$r.squared

model_Baq <- lm(Kc_air ~ Kcat, data = subset(Rubisco, Lineage == "B"))
summary_Baq <- summary(model_Baq)
r_squared_Baq <- summary_Baq$r.squared

# Linear model for Lineage D
model_D <- lm(Kc_air ~ Kcat, data = subset(Rubisco, Lineage == "D"))
summary_D <- summary(model_D)
r_squared_D <- summary_D$r.squared

print(paste("R squared for Lineage B:", r_squared_B))
print(paste("R squared for Lineage Baq:", r_squared_Baq))
print(paste("R squared for Lineage D:", r_squared_D))

# Spearman's correlation for Lineage B (terrestrial only)
cor_B <- cor.test(~ Kcat + Kc_air, data = filter(Rubisco, Lineage == "B" & Aquatic == "terrestrial"), method = "spearman")
print(paste("Spearman's rho for Lineage B (terrestrial):", cor_B$estimate))
print(paste("P-value for Lineage B (terrestrial):", cor_B$p.value))

# Spearman's correlation for Lineage B (all)
cor_Baq <- cor.test(~ Kcat + Kc_air, data = filter(Rubisco, Lineage == "B"), method = "spearman")
print(paste("Spearman's rho for Lineage B (all):", cor_Baq$estimate))
print(paste("P-value for Lineage B (all):", cor_Baq$p.value))

# Spearman's correlation for Lineage D
cor_D <- cor.test(~ Kcat + Kc_air, data = filter(Rubisco, Lineage == "D"), method = "spearman")
print(paste("Spearman's rho for Lineage D:", cor_D$estimate))
print(paste("P-value for Lineage D:", cor_D$p.value))

KcatKc_air <- ggplot(Rubisco, aes(y = Kc_air, x = Kcat, fill = Organism, shape = Organism, colour = Organism)) +
  geom_point(size = 3) +
  ylab(expression(K[c]^"21%O2")) +
  xlab(expression(k[cat]^C~(s^-1))) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3_aquatic_new" = 21,"C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3_aquatic_new" = "pink","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3_aquatic_new" = "pink","C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 25), minor_breaks = seq(0, 150, by = 25)) +
  scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, by = 1), minor_breaks = seq(1, 7, by = 1)) +
  geom_smooth(data = subset(Rubisco, Lineage == "B" & Aquatic == "terrestrial"), aes(group = Lineage), method = "lm", se = FALSE, color = "#007B3DFF", size = 0.5) +  
  geom_smooth(data = subset(Rubisco, Lineage == "B"), aes(group = Lineage), method = "lm", se = FALSE, color = "#31B7BCFF", size = 0.5) +# Line for Lineage B
  geom_smooth(data = subset(Rubisco, Lineage == "D"), aes(group = Lineage), method = "lm", se = FALSE, color = "#D51317FF", size = 0.5) +   # Line for Lineage D
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 


KcatSco <- ggplot(Rubisco, aes(y = Sco, x = Kcat, fill = Organism, shape = Organism, colour = Organism)) +
  geom_point(size = 3) +
  ylab(expression(S[C/O]~(mol.mol^-1)))+
  xlab(expression(k[cat]^C~(s^-1))) + 
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3_aquatic_new" = 21,"C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3_aquatic_new" = "pink","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3_aquatic_new" = "pink","C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_x_continuous(limits = c(1, 7), breaks = seq(1, 7, by = 1)) +
  scale_y_continuous(limits = c(50, 250), breaks = seq(50, 250, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9))

# combine plots
Rubisco_canon <-  (KcatKc_air / KcatSco) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(10,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "vertical", legend.position = c(0.9,1),
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("Rubisco_canon_neww.svg", Rubisco_canon, width = 5.4, height = 7.4)



Rubisco_content_Kc <- ggplot(Rubisco, aes(x = Content, y = Kc_air)) +
  geom_point(aes(fill = Organism, shape = Organism, colour = Organism), size = 3) +
  ylab(expression(K[c])) +
  xlab(expression("Rubisco content (% TSP)")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  geom_smooth(data = Rubisco, method = "lm", se = FALSE, color = "black", size = 0.5) +  # Line for Lineage B
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

Rubisco_content_eff <- ggplot(Rubisco, aes(x = Content, y = Kcat_over_Kc_air)) +
  geom_point(aes(fill = Organism, shape = Organism, colour = Organism), size = 3) +
  ylab(expression(k[cat]^C / K[c]^"21%O2"~(mM^-1 * "." * s^-1))) +
  xlab(expression("Rubisco content (% TSP)")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23,"haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  scale_x_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5)) +
  geom_smooth(data = Rubisco, method = "lm", se = FALSE, color = "black", size = 0.5) +  # Line for Lineage B
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

# combine plots
Rubisco_content <-  (Rubisco_content_Kc / Rubisco_content_eff) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(10,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "vertical", legend.position = c(0.9,1),
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("Rubisco_content_new2.svg", Rubisco_content, width = 5.4, height = 7.4)

# a bit more muckin around
Rubisco_content_Kc_temp <- ggplot(Rubisco, aes(x = Content, y = Kc)) +
  geom_point(aes(fill = Growtemp, shape = Aquatic), colour = "black", size = 3) +
  ylab(expression(K[c])) +
  xlab(expression("Rubisco content (% TSP)")) +
  scale_fill_gradient2(
    low = "#31B7BC", 
    mid = "white", midpoint = 12.5,
    high = "#D51317", limits = c(0,25),
    name = "Temperature") +
  scale_shape_manual(values = c("terrestrial" = 21, "semi-aquatic" = 21, "aquatic" = 23)) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, by = 12.5)) +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5)) +
  geom_smooth(data = Rubisco, method = "lm", se = FALSE, color = "black", size = 0.5) +  # Line for Lineage B
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) +
  guides(fill = guide_colorbar(frame.colour = "black",  frame.linewidth = 0.2, ticks.colour = "black"))

ggsave("Rubisco_content_temp.svg", Rubisco_content_Kc_temp, width = 4.5, height = 2.7)



kcatKcair <- ggplot(Rubisco) +
  geom_jitter(aes(x = Lineage, y = Kcat_over_Kc_air, fill = Organism, colour = Organism, shape = Organism), size = 3) +
  ylab(expression(k[cat]^C / K[c]^"21%O2"~(mM^-1 * "." * s^-1))) +
  xlab(expression("Rubisco content (% TSP)")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 275), breaks = seq(0, 275, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

Kc <- ggplot(Rubisco) +
  geom_jitter(aes(x = Lineage, y = Kc_air, fill = Organism, colour = Organism, shape = Organism), size = 3) +
  ylab(expression(K[C])) +
  xlab(expression("Lineage")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

Ko <- ggplot(Rubisco) +
  geom_jitter(aes(x = Lineage, y = Ko, fill = Organism, colour = Organism, shape = Organism), size = 3) +
  ylab(expression(K[O])) +
  xlab(expression("Lineage")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(0, 2500), breaks = seq(0, 2500, by = 250)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 

Sco <- ggplot(Rubisco) +
  geom_jitter(aes(x = Lineage, y = Sco, fill = Organism, colour = Organism, shape = Organism), size = 3) +
  ylab(expression(S[C/O]~(mol.mol^-1))) +
  xlab(expression("Lineage")) +
  scale_shape_manual(values = c("C3" = 21, "C3_aquatic" = 21, "C3-C4" = 21, "CAM" = 21,"C4" = 21, "chlorophyte" = 23, "ochrophyte" = 23, "ochrophyte_d" = 23, "ochrophyte_dts" = 23, "haptophyte" = 23, "rhodophyte" = 22)) +
  scale_colour_manual(values = c("C3" = "#007B3DFF", "C3_aquatic" = "#007B3DFF","C3-C4" = "#007B3DFF", "CAM" = "#95C11FFF","C4" = "black", "chlorophyte" = "#007B3DFF", "ochrophyte" = "black", "ochrophyte_d" = "black", "ochrophyte_dts" = "black","haptophyte" = "#D51317FF", "rhodophyte" = "#D51317FF")) +
  scale_fill_manual(values = c("C3" = "white", "C3_aquatic" = "#31B7BCFF", "C3-C4" = "#95C11FFF", "CAM" = "white", "C4" = "#007B3DFF", "chlorophyte" = "#31B7BCFF", "ochrophyte" = "#8B4513", "ochrophyte_d" = "#D51317FF", "ochrophyte_dts" = "black", "haptophyte" = "white", "rhodophyte" = "white")) +
  scale_y_continuous(limits = c(50, 250), breaks = seq(50, 250, by = 25)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9)) 


# combine plots
Rubisco_canon2 <-  (kcatKcair / Sco) | (Kc / Ko) +
  plot_layout(guides = 'collect', nrow = (2), heights = c(10,10)) +
  plot_annotation(tag_levels = list(c('A','B'))) &
  theme(axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "vertical", legend.position = "none",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plots 
ggsave("Rubisco_canon_jitter2.svg", Rubisco_canon2, width = 7.2, height = 6)


