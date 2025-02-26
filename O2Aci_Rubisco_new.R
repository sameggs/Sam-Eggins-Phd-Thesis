#' **Fitting Michaelis-Menton curve to ACi data**

# set working directory
setwd("/Users/eggboy/Dropbox/Science/Data/Oxygen Evolution/R Carbon Uptake")

# loading the tidyverse and ggplot
library(ggplot2)
library(ggpmisc)
library(ggthemes)
library(tidyverse)
library(broom)
library(ggsci)
library(patchwork)

#import data - remember to check the directory username
ACi <- read.csv("/Users/eggboy/Dropbox/Science/Data/Oxygen Evolution/csvdatafiles/Ci_curves_all_species_respincluded.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
mutate(temp = as.factor(temp),
       species = as.factor(species),
       assay = as.factor(assay),
       inhibitor = as.factor(inhibitor),
       DICadj = as.numeric(DICadj),
       O2rate_rubisco = as.numeric(O2rate_rubisco),
       O2rate_cell = as.numeric(O2rate_cell))%>%
  filter(DIC != 0, !is.na(O2rate_rubisco))

#cleaning up the data frame a bit    
ACi <- subset(ACi, DIC != 0)

# Subset the data
Pt20 <- subset(ACi, species == "P. tricornutum" & temp == 20)
Pt4 <- subset(ACi, species == "P. tricornutum" & temp == 4)
Pi <- subset(ACi, species == "P. inermis" & temp == 4)
Cf <- subset(ACi, species == "C. flexuosus" & temp == 4)


# ---- P tricornutum 20 C -----------------------------------------------------------------------------
toprub = 200 #setting initial conditions for O2max
kco2 = 0.2 #for half sat CO2
kdic = 100 #for half sat DIC

#fitting MM kinetics to whole data set
Pt20_nls_DIC <-  nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pt20, start=list(O2max=toprub,Km=kdic))
Pt20_nls_CO2 <-  nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pt20, start=list(O2max=toprub,Km=kco2))
# confidence intervals for the model using the profiling method
summary(Pt20_nls_DIC) # the standard errors produced by the summary which assumes normal distributions for model parameters
# profiling confint function..
confint(Pt20_nls_DIC)
confint(Pt20_nls_CO2)
plot(profile(Pt20_nls_DIC, "Km"))
plot(profile(Pt20_nls_DIC, "O2max")) #the distribution of the model parameters are plotted - and they're the same distributions for either the DIC or the CO2 model of course

# Extract and combine coefficients and standard errors
coefficients_Pt20 <- bind_rows(
  tidy(Pt20_nls_DIC) %>% mutate(species = "DIC"),
  tidy(Pt20_nls_CO2) %>% mutate(species = "CO2"))

# Reshape the table to include estimates and standard errors
table_Pt20 <- coefficients_Pt20 %>%
  pivot_wider(
    id_cols = species,
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pt20

# Now fitting MM curves to each individual curve instead ---
# Split the data by 'assay'
split_Pt20 <- Pt20 %>%
  filter(!is.na(O2rate_rubisco)) %>%
  split(.$assay)
split_Pt20 <- Filter(function(df) nrow(df) > 0, split_Pt20) # Filtering (base R - based on a list) data frames that are empty

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=180000,Km=100))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=180000,Km=1))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Pt20 <- map_df(split_Pt20, fit_nls_dic)
runco2efficients_Pt20 <- map_df(split_Pt20, fit_nls_co2)

# Assuming runcoefficients_Pt20 and runco2efficients_Pt20 have the same structure
combined_coefficients_Pt20 <- bind_rows(
  runcoefficients_Pt20 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt20 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt20 <- combined_coefficients_Pt20 %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# plotting up in ggplot
#plotting results for DIC
CO2_Rubisco_plot_Pt20 <- ggplot(na.omit(Pt20), aes(x = DICadj, y = O2rate_rubisco, shape = assay)) +
  geom_point(colour = "black") +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  theme_bw() +
  scale_shape_manual(values = c(19,17,15,1,2,0)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 50)) +
  scale_x_continuous(limits = c(0, 14.07), breaks = seq(0, 15, by = 2.5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  
  # plot the line of best fit    
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "2"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "3"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "4"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "5"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "6"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "7"), colour = "#333333") 

DIC_Rubisco_plot_Pt20 <- ggplot(na.omit(Pt20), aes(x = DICadj, y = O2rate_rubisco, shape = assay)) +
  geom_point(colour = "black") +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  theme_bw() +
  scale_shape_manual(values = c(19,17,15,1,2,0)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 50)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  
  # plot the line of best fit    
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "2"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "3"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "4"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "5"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "6"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt20, assay == "7"), colour = "#333333") 


# ---- P tricornutum 4 C -----------------------------------------------------------------------------
toprub = 20 #setting initial conditions for O2max
kco2 = 0.5 #for half sat CO2
kdic = 50 #for half sat DIC

#fitting MM kinetics to whole data set
Pt4_nls_DIC <-  nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pt4, start=list(O2max=toprub,Km=kdic))
Pt4_nls_CO2 <-  nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pt4, start=list(O2max=toprub,Km=kco2))
# confidence intervals for the model using the profiling method
summary(Pt4_nls_DIC) # the standard errors produced by the summary which assumes normal distributions for model parameters
# profiling confint function..
confint(Pt4_nls_DIC)
confint(Pt4_nls_CO2)
plot(profile(Pt4_nls_DIC, "Km"))
plot(profile(Pt4_nls_DIC, "O2max")) #the distribution of the model parameters are plotted - and they're the same distributions for either the DIC or the CO2 model of course

# Extract and combine coefficients and standard errors
coefficients_Pt4 <- bind_rows(
  tidy(Pt4_nls_DIC) %>% mutate(species = "DIC"),
  tidy(Pt4_nls_CO2) %>% mutate(species = "CO2"))

# Reshape the table to include estimates and standard errors
table_Pt4 <- coefficients_Pt4 %>%
  pivot_wider(
    id_cols = species,
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pt4

# Now fitting MM curves to each individual curve instead ---
# Split the data by 'assay'
split_Pt4 <- Pt4 %>%
  filter(!is.na(O2rate_rubisco)) %>%
  split(.$assay)
split_Pt4 <- Filter(function(df) nrow(df) > 0, split_Pt4) # Filtering (base R - based on a list) data frames that are empty

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=180000,Km=100))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=180000,Km=1))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Pt4 <- map_df(split_Pt4, fit_nls_dic)
runco2efficients_Pt4 <- map_df(split_Pt4, fit_nls_co2)

# Assuming runcoefficients_Pt4 and runco2efficients_Pt4 have the same structure
combined_coefficients_Pt4 <- bind_rows(
  runcoefficients_Pt4 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt4 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt4 <- combined_coefficients_Pt4 %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# plotting up in ggplot
#plotting results for DIC
CO2_Rubisco_plot_Pt4 <- ggplot(na.omit(Pt4), aes(x = CO2adj, y = O2rate_rubisco, shape = assay)) +
  geom_point(colour = "black") +
  theme_bw() +
  scale_shape_manual(values = c(19,17,15,1,2,0)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 25)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  
  # plot the line of best fit    
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "2"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "3"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "4"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "5"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "6"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kco2)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "7"), colour = "#333333") 

DIC_Rubisco_plot_Pt4 <- ggplot(na.omit(Pt4), aes(x = DICadj, y = O2rate_rubisco, shape = assay)) +
  geom_point(colour = "black") +
  theme_bw() +
  scale_shape_manual(values = c(19,17,15,1,2,0)) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 25)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  
  # plot the line of best fit    
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "2"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "3"), colour = "#333333") +
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "4"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "5"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "6"), colour = "#333333") + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ Vmax * x / (Km + x), start = list(Vmax = toprub, Km = kdic)),
              se = FALSE, size = 0.5, data = subset(Pt4, assay == "7"), colour = "#333333") 

# ---- P inermis 4 C -----------------------------------------------------------------------------
toprub = 1 #setting initial conditions for O2max
kco2 = 0.5 #for half sat CO2
kdic = 50 #for half sat DIC

#fitting MM kinetics to whole data set
Pi_nls_DIC <-  nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pi, start=list(O2max=toprub,Km=kdic))
Pi_nls_CO2 <-  nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pi, start=list(O2max=toprub,Km=kco2))
# confidence intervals for the model using the profiling method
summary(Pi_nls_DIC) # the standard errors produced by the summary which assumes normal distributions for model parameters
# profiling confint function..
confint(Pi_nls_DIC)
confint(Pi_nls_CO2)
plot(profile(Pi_nls_DIC, "Km"))
plot(profile(Pi_nls_DIC, "O2max")) #the distribution of the model parameters are plotted - and they're the same distributions for either the DIC or the CO2 model of course

# Extract and combine coefficients and standard errors
coefficients_Pi <- bind_rows(
  tidy(Pi_nls_DIC) %>% mutate(species = "DIC"),
  tidy(Pi_nls_CO2) %>% mutate(species = "CO2"))

# Reshape the table to include estimates and standard errors
table_Pi <- coefficients_Pi %>%
  pivot_wider(
    id_cols = species,
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pi

# Now fitting MM curves to each individual curve instead ---
# Split the data by 'assay'
split_Pi <- Pi %>%
  filter(!is.na(O2rate_rubisco)) %>%
  split(.$assay)
split_Pi <- Filter(function(df) nrow(df) > 0, split_Pi) # Filtering (base R - based on a list) data frames that are empty

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=180000,Km=100))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=180000,Km=1))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Pi <- map_df(split_Pi, fit_nls_dic)
runco2efficients_Pi <- map_df(split_Pi, fit_nls_co2)

# Assuming runcoefficients_Pi and runco2efficients_Pi have the same structure
combined_coefficients_Pi <- bind_rows(
  runcoefficients_Pi %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pi %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pi <- combined_coefficients_Pi %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Bootstrapping all data
nls_DIC_Pi <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pi, start=list(O2max=toprub,Km=1))
nls_CO2_Pi <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pi, start=list(O2max=toprub,Km=1))
boot_DIC_Pi <- nlsBoot(nls_DIC_Pi)
boot_CO2_Pi <- nlsBoot(nls_CO2_Pi)

# Extract bootCI data
bootCI_DIC_Pi <- as_tibble(boot_DIC_Pi$bootCI, rownames = "coef_name") %>%
  mutate(species = "DIC")  
bootCI_CO2_Pi <- as_tibble(boot_CO2_Pi$bootCI, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_bootCI_Pi <- rbind(bootCI_DIC_Pi, bootCI_CO2_Pi)

# Extract estiboot data
estiboot_DIC_Pi <- as_tibble(boot_DIC_Pi$estiboot, rownames = "coef_name") %>%
  mutate(species = "DIC")  # Add species column
estiboot_CO2_Pi <- as_tibble(boot_CO2_Pi$estiboot, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_estiboot_Pi <- rbind(estiboot_DIC_Pi, estiboot_CO2_Pi)  

# Combine bootCI and estiboot data
boot_data_Pi <- left_join(combined_bootCI_Pi, combined_estiboot_Pi, by = c("coef_name", "species"))

# Reshape the table to include estimates and standard errors
boot_table_Pi <- boot_data_Pi %>%
  pivot_wider(
    id_cols = c(species),
    names_from = coef_name,
    values_from = c(Estimate, "Std. error", "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
boot_table_Pi

# Generic function for Michaelis-Menten fit
MM_fit_Pi <- function(x, species, inhibitor) {
  coefs <- filter(boot_table_Pi, species == !!species)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MMPi_DIC <- function(x) MM_fit_Pi(x, "DIC")
MMPi_CO2 <- function(x) MM_fit_Pi(x, "CO2")


# MM fits over confidence intervals 
upPi_DIC <- sapply(Pi$DICadj, FUN = function(x) boot_DIC_Pi$bootCI[1,3] * x / (boot_DIC_Pi$bootCI[2,2] + x))
lwPi_DIC <- sapply(Pi$DICadj, FUN = function(x) boot_DIC_Pi$bootCI[1,2] * x / (boot_DIC_Pi$bootCI[2,3] + x))
upPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) boot_CO2_Pi$bootCI[1,3] * x / (boot_CO2_Pi$bootCI[2,2] + x))
lwPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) boot_CO2_Pi$bootCI[1,2] * x / (boot_CO2_Pi$bootCI[2,3] + x))

# plotting up in ggplot
CO2_Rubisco_plot_Pi <- ggplot(Pi, aes(x = CO2adj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwPi_CO2, ymax = upPi_CO2), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MMPi_CO2, colour = "black", size = 0.7)

DIC_Rubisco_plot_Pi <- ggplot(Pi, aes(x = DICadj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwPi_DIC, ymax = upPi_DIC), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MMPi_DIC, colour = "black", size = 0.7)


# ---- C flexuosus 4 C -----------------------------------------------------------------------------
toprub = 0.5 #setting initial conditions for O2max
kco2 = 0.5 #for half sat CO2
kdic = 50 #for half sat DIC

#fitting MM kinetics to whole data set
Cf_nls_DIC <-  nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Cf, start=list(O2max=toprub,Km=kdic))
Cf_nls_CO2 <-  nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Cf, start=list(O2max=toprub,Km=kco2))
# confidence intervals for the model using the profiling method
summary(Cf_nls_DIC) # the standard errors produced by the summary which assumes normal distributions for model parameters
# profiling confint function..
confint(Cf_nls_DIC)
confint(Cf_nls_CO2)
plot(profile(Cf_nls_DIC, "Km"))
plot(profile(Cf_nls_DIC, "O2max")) #the distribution of the model parameters are plotted - and they're the same distributions for either the DIC or the CO2 model of course

# Extract and combine coefficients and standard errors
coefficients_Cf <- bind_rows(
  tidy(Cf_nls_DIC) %>% mutate(species = "DIC"),
  tidy(Cf_nls_CO2) %>% mutate(species = "CO2"))

# Reshape the table to include estimates and standard errors
table_Cf <- coefficients_Cf %>%
  pivot_wider(
    id_cols = species,
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Cf

# Now fitting MM curves to each individual curve instead ---
# Split the data by 'assay'
split_Cf <- Cf %>%
  filter(!is.na(O2rate_rubisco)) %>%
  split(.$assay)
split_Cf <- Filter(function(df) nrow(df) > 0, split_Cf) # Filtering (base R - based on a list) data frames that are empty

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=0.5,Km=50))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=0.5,Km=0.5))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Cf <- map_df(split_Cf, fit_nls_dic)
runco2efficients_Cf <- map_df(split_Cf, fit_nls_co2)

# Assuming runcoefficients_Cf and runco2efficients_Cf have the same structure
combined_coefficients_Cf <- bind_rows(
  runcoefficients_Cf %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Cf %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Cf <- combined_coefficients_Cf %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Bootstrapping all data
nls_DIC_Cf <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Cf, start=list(O2max=toprub,Km=1))
nls_CO2_Cf <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Cf, start=list(O2max=toprub,Km=1))
boot_DIC_Cf <- nlsBoot(nls_DIC_Cf)
boot_CO2_Cf <- nlsBoot(nls_CO2_Cf)

# Extract bootCI data
bootCI_DIC_Cf <- as_tibble(boot_DIC_Cf$bootCI, rownames = "coef_name") %>%
  mutate(species = "DIC")  
bootCI_CO2_Cf <- as_tibble(boot_CO2_Cf$bootCI, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_bootCI_Cf <- rbind(bootCI_DIC_Cf, bootCI_CO2_Cf)

# Extract estiboot data
estiboot_DIC_Cf <- as_tibble(boot_DIC_Cf$estiboot, rownames = "coef_name") %>%
  mutate(species = "DIC")  # Add species column
estiboot_CO2_Cf <- as_tibble(boot_CO2_Cf$estiboot, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_estiboot_Cf <- rbind(estiboot_DIC_Cf, estiboot_CO2_Cf)  

# Combine bootCI and estiboot data
boot_data_Cf <- left_join(combined_bootCI_Cf, combined_estiboot_Cf, by = c("coef_name", "species"))

# Reshape the table to include estimates and standard errors
boot_table_Cf <- boot_data_Cf %>%
  pivot_wider(
    id_cols = c(species),
    names_from = coef_name,
    values_from = c(Estimate, "Std. error", "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
boot_table_Cf

# Generic function for Michaelis-Menten fit
MM_fit_Cf <- function(x, species, inhibitor) {
  coefs <- filter(boot_table_Cf, species == !!species)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MMCf_DIC <- function(x) MM_fit_Cf(x, "DIC")
MMCf_CO2 <- function(x) MM_fit_Cf(x, "CO2")


# MM fits over confidence intervals 
upCf_DIC <- sapply(Cf$DICadj, FUN = function(x) boot_DIC_Cf$bootCI[1,3] * x / (boot_DIC_Cf$bootCI[2,2] + x))
lwCf_DIC <- sapply(Cf$DICadj, FUN = function(x) boot_DIC_Cf$bootCI[1,2] * x / (boot_DIC_Cf$bootCI[2,3] + x))
upCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) boot_CO2_Cf$bootCI[1,3] * x / (boot_CO2_Cf$bootCI[2,2] + x))
lwCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) boot_CO2_Cf$bootCI[1,2] * x / (boot_CO2_Cf$bootCI[2,3] + x))

# plotting the ACi curves
CO2_Rubisco_plot_Cf <- ggplot(Cf, aes(x = CO2adj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwCf_CO2, ymax = upCf_CO2), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MMCf_CO2, colour = "black", size = 0.7)

DIC_Rubisco_plot_Cf <- ggplot(Cf, aes(x = DICadj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwCf_DIC, ymax = upCf_DIC), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.1)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MMCf_DIC, colour = "black", size = 0.7)


# ---- P tricornutum 20 C removing dud values -----------------------------------------------------------------------------
toprub = 2 #setting initial conditions for O2max

# Selecting the 2 highest Rubisco assays
selected_Pt20 <- split_Pt20[c("4", "5")]
Pt20_filtered <- Pt20 %>%
  filter(assay %in% c("4", "5"))

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=toprub,Km=100))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=toprub,Km=1))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Pt20 <- map_df(selected_Pt20, fit_nls_dic)
runco2efficients_Pt20 <- map_df(selected_Pt20, fit_nls_co2)

# Assuming runcoefficients_Pt20 and runco2efficients_Pt20 have the same structure
combined_coefficients_Pt20 <- bind_rows(
  runcoefficients_Pt20 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt20 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt20_selected <- combined_coefficients_Pt20 %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Bootstrapping all data
nls_select_DIC_Pt20 <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pt20_filtered, start=list(O2max=toprub,Km=1))
nls_select_CO2_Pt20 <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pt20_filtered, start=list(O2max=toprub,Km=1))
boot_DIC_Pt20 <- nlsBoot(nls_select_DIC_Pt20)
boot_CO2_Pt20 <- nlsBoot(nls_select_CO2_Pt20)

# Extract bootCI data
bootCI_DIC_Pt20 <- as_tibble(boot_DIC_Pt20$bootCI, rownames = "coef_name") %>%
  mutate(species = "DIC")  
bootCI_CO2_Pt20 <- as_tibble(boot_CO2_Pt20$bootCI, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_bootCI_Pt20 <- rbind(bootCI_DIC_Pt20, bootCI_CO2_Pt20)

# Extract estiboot data
estiboot_DIC_Pt20 <- as_tibble(boot_DIC_Pt20$estiboot, rownames = "coef_name") %>%
  mutate(species = "DIC")  # Add species column
estiboot_CO2_Pt20 <- as_tibble(boot_CO2_Pt20$estiboot, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_estiboot_Pt20 <- rbind(estiboot_DIC_Pt20, estiboot_CO2_Pt20)  

# Combine bootCI and estiboot data
boot_data_Pt20 <- left_join(combined_bootCI_Pt20, combined_estiboot_Pt20, by = c("coef_name", "species"))

# Reshape the table to include estimates and standard errors
boot_table_Pt20 <- boot_data_Pt20 %>%
  pivot_wider(
    id_cols = c(species),
    names_from = coef_name,
    values_from = c(Estimate, "Std. error", "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
boot_table_Pt20

# Generic function for Michaelis-Menten fit
MM_fit_Pt20 <- function(x, species, inhibitor) {
  coefs <- filter(boot_table_Pt20, species == !!species)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MM20_DIC <- function(x) MM_fit_Pt20(x, "DIC")
MM20_CO2 <- function(x) MM_fit_Pt20(x, "CO2")


# MM fits over confidence intervals 
up20_DIC <- sapply(Pt20_filtered$DICadj, FUN = function(x) boot_DIC_Pt20$bootCI[1,3] * x / (boot_DIC_Pt20$bootCI[2,2] + x))
lw20_DIC <- sapply(Pt20_filtered$DICadj, FUN = function(x) boot_DIC_Pt20$bootCI[1,2] * x / (boot_DIC_Pt20$bootCI[2,3] + x))
up20_CO2 <- sapply(Pt20_filtered$CO2adj, FUN = function(x) boot_CO2_Pt20$bootCI[1,3] * x / (boot_CO2_Pt20$bootCI[2,2] + x))
lw20_CO2 <- sapply(Pt20_filtered$CO2adj, FUN = function(x) boot_CO2_Pt20$bootCI[1,2] * x / (boot_CO2_Pt20$bootCI[2,3] + x))

# plotting up in ggplot
#plotting results for DIC
CO2_Rubisco_plot_Pt20_selected <- ggplot(Pt20_filtered, aes(x = CO2adj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw20_CO2, ymax = up20_CO2), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5)) +
  scale_x_continuous(limits = c(0, 14.07), breaks = seq(0, 15, by = 2.5)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MM20_CO2, colour = "black", size = 0.7)
  
DIC_Rubisco_plot_Pt20_selected <- ggplot(Pt20_filtered, aes(x = DICadj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw20_DIC, ymax = up20_DIC), fill = "#706F6F", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MM20_DIC, colour = "black", size = 0.7)
  

# ---- P tricornutum 4 C removing dud values -----------------------------------------------------------------------------
toprub = 0.5 #setting initial conditions for O2max

# Selecting the 2 highest Rubisco assays
selected_Pt4 <- split_Pt4[c("5","6")]
Pt4_filtered <- Pt4 %>%
  filter(assay %in% c("5", "6")) 

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=toprub,Km=100))
  return(coef(bestfit))
}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=toprub,Km=1))
  return(coef(bestfit))
}

# Apply the functions and combine results
runcoefficients_Pt4 <- map_df(selected_Pt4, fit_nls_dic)
runco2efficients_Pt4 <- map_df(selected_Pt4, fit_nls_co2)

# Assuming runcoefficients_Pt4 and runco2efficients_Pt4 have the same structure
combined_coefficients_Pt4 <- bind_rows(
  runcoefficients_Pt4 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt4 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt4_selected <- combined_coefficients_Pt4 %>%
  group_by(CoefficientType) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

# Bootstrapping all data
nls_select_DIC_Pt4 <- nls(O2rate_rubisco ~ O2max*DICadj/(Km+DICadj), Pt4_filtered, start=list(O2max=toprub,Km=1))
nls_select_CO2_Pt4 <- nls(O2rate_rubisco ~ O2max*CO2adj/(Km+CO2adj), Pt4_filtered, start=list(O2max=toprub,Km=1))
boot_DIC_Pt4 <- nlsBoot(nls_select_DIC_Pt4)
boot_CO2_Pt4 <- nlsBoot(nls_select_CO2_Pt4)

# Extract bootCI data
bootCI_DIC_Pt4 <- as_tibble(boot_DIC_Pt4$bootCI, rownames = "coef_name") %>%
  mutate(species = "DIC")  
bootCI_CO2_Pt4 <- as_tibble(boot_CO2_Pt4$bootCI, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_bootCI_Pt4 <- rbind(bootCI_DIC_Pt4, bootCI_CO2_Pt4)

# Extract estiboot data
estiboot_DIC_Pt4 <- as_tibble(boot_DIC_Pt4$estiboot, rownames = "coef_name") %>%
  mutate(species = "DIC")  # Add species column
estiboot_CO2_Pt4 <- as_tibble(boot_CO2_Pt4$estiboot, rownames = "coef_name") %>%
  mutate(species = "CO2")  
combined_estiboot_Pt4 <- rbind(estiboot_DIC_Pt4, estiboot_CO2_Pt4) 

# Combine bootCI and estiboot data
boot_data_Pt4 <- left_join(combined_bootCI_Pt4, combined_estiboot_Pt4, by = c("coef_name", "species"))

# Reshape the table to include estimates and standard errors
boot_table_Pt4 <- boot_data_Pt4 %>%
  pivot_wider(
    id_cols = c(species),
    names_from = coef_name,
    values_from = c(Estimate, "Std. error", "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
boot_table_Pt4

# Generic function for Michaelis-Menten fit
MM_fit_Pt4 <- function(x, species, inhibitor) {
  coefs <- filter(boot_table_Pt4, species == !!species)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MM4_DIC <- function(x) MM_fit_Pt4(x, "DIC")
MM4_CO2 <- function(x) MM_fit_Pt4(x, "CO2")

# MM fits over confidence intervals 
up4_DIC <- sapply(Pt4_filtered$DICadj, FUN = function(x) boot_DIC_Pt4$bootCI[1,3] * x / (boot_DIC_Pt4$bootCI[2,2] + x))
lw4_DIC <- sapply(Pt4_filtered$DICadj, FUN = function(x) boot_DIC_Pt4$bootCI[1,2] * x / (boot_DIC_Pt4$bootCI[2,3] + x))
up4_CO2 <- sapply(Pt4_filtered$CO2adj, FUN = function(x) boot_CO2_Pt4$bootCI[1,3] * x / (boot_CO2_Pt4$bootCI[2,2] + x))
lw4_CO2 <- sapply(Pt4_filtered$CO2adj, FUN = function(x) boot_CO2_Pt4$bootCI[1,2] * x / (boot_CO2_Pt4$bootCI[2,3] + x))

CO2_Rubisco_plot_Pt4_selected <- ggplot(Pt4_filtered, aes(x = CO2adj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw4_CO2, ymax = up4_CO2), fill = "#706F6F", colour= NA, alpha = 0.2) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.2)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MM4_CO2, colour = "black", size = 0.7)

DIC_Rubisco_plot_Pt4_selected <- ggplot(Pt4_filtered, aes(x = DICadj, y = O2rate_rubisco)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw4_DIC, ymax = up4_DIC), fill = "#706F6F", colour= NA, alpha = 0.2) +
  geom_point(aes(fill = assay), colour = "black", shape = 23) +
  scale_fill_frontiers() +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.1)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  ylab(bquote('mol'~O[2]~.~mol~Rubisco^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcM)')) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  stat_function(fun = MM4_DIC, colour = "black", size = 0.7)

# --- plotting all together but without the duds ---------------

ACi_CO2_rubisco_selected <- guide_area() / (CO2_Rubisco_plot_Pt20_selected | CO2_Rubisco_plot_Pt4_selected) / (CO2_Rubisco_plot_Pi | CO2_Rubisco_plot_Cf) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ACi_DIC_rubisco_selected <- guide_area() / (DIC_Rubisco_plot_Pt20_selected | DIC_Rubisco_plot_Pt4_selected) / (DIC_Rubisco_plot_Pi | DIC_Rubisco_plot_Cf) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("ACi_CO2_Rubisco_plot_2.svg", ACi_CO2_rubisco_selected, width = 7.2, height = 7.4)
ggsave("ACi_DIC_Rubisco_plot_2.svg", ACi_DIC_rubisco_selected, width = 7.2, height = 7.4)

#----- residuals vs fitted -------------------------
plot(residuals(Pt20_nls_DIC) ~ fitted(Pt20_nls_DIC), 
     xlab = "Fitted values", ylab = "Residuals", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)

plot(residuals(Pt4_nls_DIC) ~ fitted(Pt4_nls_DIC), 
     xlab = "Fitted values", ylab = "Residuals", 
     main = "Residuals vs Fitted")
abline(h = 0, col = "red", lty = 2)
