#' **Fitting Michaelis-Menton curve to ACi data**

# set working directory
setwd("/Users/eggboy/Dropbox/Science/Data/Oxygen Evolution/R Carbon Uptake")

# loading the tidyverse and ggplot
library(ggplot2)
library(ggpmisc)
library(ggthemes)
library(tidyverse)
library(broom)
library(ggsci) # for colour schemes on plots
library(patchwork) # for arranging plots
library(nlstools) # for non-linear least squares fitting
library(car) # for levene's test

#import data - remember to check the directory username
ACi <- read.csv("/Users/eggboy/Dropbox/Science/Data/Oxygen Evolution/csvdatafiles/Ci_curves_all_species_respincluded.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(temp = as.factor(temp),
         species = as.factor(species),
         assay = as.factor(assay),
         inhibitor = as.factor(inhibitor),
         DICadj = as.numeric(DICadj),
         O2rate_cell = as.numeric(O2rate_cell),
         O2rate_cell = as.numeric(O2rate_cell)) %>%
  filter(!is.na(O2rate_cell))

#ordering factor levels
ACi$inhibitor <- factor(ACi$inhibitor, levels=c("none", "AZ", "EZ"))
# applying a scaling factor so the nls function can handle the data
ACi$O2rate_cell <- ACi$O2rate_cell * 3600 * 1e15

# Subset the data
Pt20 <- subset(ACi, species == "P. tricornutum" & temp == 20)
Pt4 <- subset(ACi, species == "P. tricornutum" & temp == 4 & assay != 20) # removing assay 20 for the per cell data
Pi <- subset(ACi, species == "P. inermis" & temp == 4)
Cf <- subset(ACi, species == "C. flexuosus" & temp == 4)

# ------------------------------- P tricornutum 20 C -----------------------------------------------------------------------------------
top = 600 #setting initial conditions for O2max
kco2 = 0.2 #for half sat CO2
kdic = 100 #for half sat DIC

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=top,Km=100)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=top,Km=1)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Fit models and retain 'inhibitor' information
runcoefficients_Pt20 <- Pt20 %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_dic = fit_nls_dic(.)) %>%
  unnest_wider(coefs_dic)

runco2efficients_Pt20 <- Pt20 %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_co2 = fit_nls_co2(.)) %>%
  unnest_wider(coefs_co2)

# Assuming runcoefficients_Pt4 and runco2efficients_Pt4 have the same structure
combined_coefficients_Pt20 <- bind_rows(
  runcoefficients_Pt20 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt20 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt20 <- combined_coefficients_Pt20 %>%
  group_by(CoefficientType, inhibitor) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

#fitting MM kinetics to whole data set
Pt20_nls_DIC <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt20[which(Pt20$inhibitor == "none"),], start=list(O2max=top,Km=kdic))
Pt20_nls_DIC_AZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt20[which(Pt20$inhibitor == "AZ"),], start=list(O2max=top,Km=kdic))
Pt20_nls_DIC_EZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt20[which(Pt20$inhibitor == "EZ"),], start=list(O2max=top,Km=kdic))
Pt20_nls_CO2 <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt20[which(Pt20$inhibitor == "none"),], start=list(O2max=top,Km=kco2))
Pt20_nls_CO2_AZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt20[which(Pt20$inhibitor == "AZ"),], start=list(O2max=top,Km=kco2))
Pt20_nls_CO2_EZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt20[which(Pt20$inhibitor == "EZ"),], start=list(O2max=top,Km=kco2))
#list the models
models_Pt20 <- list(Pt20_nls_DIC = Pt20_nls_DIC,
                    Pt20_nls_DIC_AZ = Pt20_nls_DIC_AZ,
                    Pt20_nls_DIC_EZ = Pt20_nls_DIC_EZ,
                    Pt20_nls_CO2 = Pt20_nls_CO2,
                    Pt20_nls_CO2_AZ = Pt20_nls_CO2_AZ,
                    Pt20_nls_CO2_EZ = Pt20_nls_CO2_EZ)

# Apply summary and confidence intervals using profiling to each model in the list
model_summaries_Pt20 <- map(models_Pt20, summary)
model_confints_Pt20 <- map(models_Pt20, confint)

# Setting up the plotting area for a grid (3 models x 2 parameters)
par(mfrow = c(3, 2))
# Loop over each model and plot profiles for Km and O2max
for(model in models_Pt20) {
  plot(profile(model, "Km"))
  plot(profile(model, "O2max"))
}
par(mfrow = c(1, 1)) # Resetting to default plotting layout

# Extract and combine coefficients and standard errors
coefficients_Pt20 <- map_df(names(models_Pt20), ~{
  model_data <- tidy(models_Pt20[[.x]])
  species <- ifelse(grepl("DIC", .x), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  mutate(model_data, species = species, inhibitor = inhibitor)
})

# Transforming confidence interval data
confint_data_Pt20 <- map_df(names(model_confints_Pt20), ~{
  confint_data <- as_tibble(model_confints_Pt20[[.x]], rownames = "term") %>%
    mutate(species = ifelse(grepl("DIC", .x), "DIC", "CO2"),
           inhibitor = case_when(
             grepl("AZ", .x) ~ "AZ",
             grepl("EZ", .x) ~ "EZ",
             TRUE ~ "none")) 
  confint_data
}, .id = "model")
# joining the data together
coefficients_Pt20 <- left_join(coefficients_Pt20, select(confint_data_Pt20, -model), 
                               by = c("species", "inhibitor", "term"))

# Reshape the table to include estimates and standard errors
table_Pt20 <- coefficients_Pt20 %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
    names_from = term,
    values_from = c(estimate, std.error, "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pt20

#bootstrapping
boot_models_Pt20 <- map(models_Pt20, nlsBoot)

# Function to extract bootCI and estiboot, and add species and inhibitor
extract_boot_data <- function(model_name, model_boot) {
  # Extract and convert bootCI and estiboot data to tibble, capturing row names
  bootCI_data <- as_tibble(model_boot$bootCI, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  estiboot_data <- as_tibble(model_boot$estiboot, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  
  # Combine bootCI and estiboot data
  combined_data <- cbind(bootCI_data, estiboot_data[,-1]) # Exclude repeated coef_name column from estiboot_data
  
  # Add species and inhibitor columns
  species <- ifelse(grepl("DIC", model_name), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", model_name) ~ "AZ",
    grepl("EZ", model_name) ~ "EZ",
    TRUE ~ "none"
  )
  combined_data <- cbind(species = species, inhibitor = inhibitor, combined_data)
  
  return(combined_data)
}

# Apply the function to each bootstrapped model
boot_data_Pt20 <- map2_df(names(boot_models_Pt20), boot_models_Pt20, extract_boot_data)

# View the table
boot_data_Pt20

# Reshape the table to include estimates and standard errors
boot_table_Pt20 <- boot_data_Pt20 %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
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
  coefs <- filter(boot_table_Pt20, species == !!species, inhibitor == !!inhibitor)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MM20_DIC <- function(x) MM_fit_Pt20(x, "DIC", "none")
MMAZ20_DIC <- function(x) MM_fit_Pt20(x, "DIC", "AZ")
MMEZ20_DIC <- function(x) MM_fit_Pt20(x, "DIC", "EZ")
MM20_CO2 <- function(x) MM_fit_Pt20(x, "CO2", "none")
MMAZ20_CO2 <- function(x) MM_fit_Pt20(x, "CO2", "AZ")
MMEZ20_CO2 <- function(x) MM_fit_Pt20(x, "CO2", "EZ")

#functions for confidence intervals
# Helper function to get coefficients and confidence intervals
get_coef_conf_Pt20 <- function(species, inhibitor, coef_name, conf_level) {
  filtered_data <- boot_data_Pt20 %>%
    filter(species == !!species, inhibitor == !!inhibitor, coef_name == !!coef_name)
  result <- filtered_data %>%
    select(!!sym(conf_level)) %>%
    unlist()
  return(result)
}

# Applying the helper function to created MM fits over confidence intervals 
up20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "none", "O2max", "97.5%") * x / (get_coef_conf_Pt20("DIC", "none", "Km", "2.5%") + x))
lw20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "none", "O2max", "2.5%") * x / (get_coef_conf_Pt20("DIC", "none", "Km", "97.5%") + x))
upAZ20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pt20("DIC", "AZ", "Km", "2.5%") + x))
lwAZ20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pt20("DIC", "AZ", "Km", "97.5%") + x))
upEZ20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pt20("DIC", "EZ", "Km", "2.5%") + x))
lwEZ20_DIC <- sapply(Pt20$DICadj, FUN = function(x) get_coef_conf_Pt20("DIC", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pt20("DIC", "EZ", "Km", "97.5%") + x))
up20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "none", "O2max", "97.5%") * x / (get_coef_conf_Pt20("CO2", "none", "Km", "2.5%") + x))
lw20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "none", "O2max", "2.5%") * x / (get_coef_conf_Pt20("CO2", "none", "Km", "97.5%") + x))
upAZ20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pt20("CO2", "AZ", "Km", "2.5%") + x))
lwAZ20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pt20("CO2", "AZ", "Km", "97.5%") + x))
upEZ20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pt20("CO2", "EZ", "Km", "2.5%") + x))
lwEZ20_CO2 <- sapply(Pt20$CO2adj, FUN = function(x) get_coef_conf_Pt20("CO2", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pt20("CO2", "EZ", "Km", "97.5%") + x))


# plotting up in ggplot
#plotting results for DIC
CO2_cell_plot_Pt20 <- ggplot(Pt20, aes(x = CO2adj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw20_CO2, ymax = up20_CO2), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZ20_CO2, ymax = upAZ20_CO2), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZ20_CO2, ymax = upEZ20_CO2), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('fmol'~O[2]~.~hr^-1~cell^-1)) +  
  xlab(bquote('CO2 (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-400, 1200), breaks = seq(-400, 1200, by = 200)) +
  scale_x_continuous(limits = c(0, 14.07), breaks = seq(0, 15, by = 2.5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MM20_CO2, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZ20_CO2, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZ20_CO2, colour= "#0094CDFF", size = 0.7) 

DIC_cell_plot_Pt20 <- ggplot(Pt20, aes(x = DICadj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw20_DIC, ymax = up20_DIC), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZ20_DIC, ymax = upAZ20_DIC), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZ20_DIC, ymax = upEZ20_DIC), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('fmol'~O[2]~.~hr^-1~cell^-1)) +  
  xlab(bquote('DIC (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-400, 1200), breaks = seq(-400, 1200, by = 200)) +
  scale_x_continuous(limits = c(0, 3000), breaks = seq(0, 3000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MM20_DIC, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZ20_DIC, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZ20_DIC, colour= "#0094CDFF", size = 0.7)

# ------ conducting model comparisons ----------------------------------------------------
# Create a function to apply nlsConfRegions to a single model
apply_nlsConfRegions <- function(model) {
  nlsConfRegions(model, length = 500, exp = 2)
}

# Apply the function to all models in models_Pt20
Bealeconf_regions_Pt20 <- lapply(models_Pt20, apply_nlsConfRegions)

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
Bealeconf_Pt20_DIC <- list(
  Beale_DIC = data.frame(x = Bealeconf_regions_Pt20$Pt20_nls_DIC$cr[, 1], y = Bealeconf_regions_Pt20$Pt20_nls_DIC$cr[, 2]),
  Beale_DIC_AZ = data.frame(x = Bealeconf_regions_Pt20$Pt20_nls_DIC_AZ$cr[, 1], y = Bealeconf_regions_Pt20$Pt20_nls_DIC_AZ$cr[, 2]),
  Beale_DIC_EZ = data.frame(x = Bealeconf_regions_Pt20$Pt20_nls_DIC_EZ$cr[, 1], y = Bealeconf_regions_Pt20$Pt20_nls_DIC_EZ$cr[, 2])
)

Beale_Pt20 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 350, by = 50)) +
  scale_x_continuous(limits = c(400, 800), breaks = seq(400, 800, by = 50)) +
  geom_point(data = Bealeconf_Pt20_DIC$Beale_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = Bealeconf_Pt20_DIC$Beale_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = Bealeconf_Pt20_DIC$Beale_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
boot_Pt20_DIC <- list(
  boots_DIC = data.frame(x = boot_models_Pt20$Pt20_nls_DIC$coefboot[, 1], y = boot_models_Pt20$Pt20_nls_DIC$coefboot[, 2]),
  boots_DIC_AZ = data.frame(x = boot_models_Pt20$Pt20_nls_DIC_AZ$coefboot[, 1], y = boot_models_Pt20$Pt20_nls_DIC_AZ$coefboot[, 2]),
  boots_DIC_EZ = data.frame(x = boot_models_Pt20$Pt20_nls_DIC_EZ$coefboot[, 1], y = boot_models_Pt20$Pt20_nls_DIC_EZ$coefboot[, 2])
)

Bootplot_Pt20 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 350), breaks = seq(0, 350, by = 50)) +
  scale_x_continuous(limits = c(400, 800), breaks = seq(400, 800, by = 50)) +
  geom_point(data = boot_Pt20_DIC$boots_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = boot_Pt20_DIC$boots_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = boot_Pt20_DIC$boots_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# -- fitting models to individual assays for hypothesis testing -------------------
fit_nls_DIC <- function(data) {
  nls_model <- nls(O2rate_cell ~ O2max * DICadj / (Km + DICadj), 
                   data = data,
                   start = list(O2max = top, Km = kdic))
  return(nls_model)
}

# Fitting the model for each assay and recording the associated inhibitor
assay_models_Pt20 <- list()
for (assay in unique(Pt20$assay)) {
  assay_data <- Pt20[Pt20$assay == assay, ]
  inhibitor_used <- unique(assay_data$inhibitor) # Record the inhibitor used
  # Fit the model
  fitted_model <- fit_nls_DIC(assay_data) 
  # Store the results 
  model_key <- paste("Assay", assay, "Inhibitor", inhibitor_used, sep = "_")
  assay_models_Pt20[[model_key]] <- fitted_model
}

#extract summaries for the models 
assay_summaries_Pt20 <- map(assay_models_Pt20, summary)

# Extract and combine coefficients and standard errors
assay_coefficients_Pt20 <- map_df(names(assay_models_Pt20), ~{
  model_data <- tidy(assay_models_Pt20[[.x]])
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  assay_number <- str_extract(.x, "(?<=Assay_)[0-9]+") # Extract the assay number
  mutate(model_data, inhibitor = inhibitor, assay = as.numeric(assay_number))
})

# Reshape the table to include estimates and standard errors
assay_table_Pt20 <- assay_coefficients_Pt20 %>%
  pivot_wider(
    id_cols = c(inhibitor, assay),
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
assay_table_Pt20

Boot_plusdata_Pt20 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 50)) +  
  scale_y_continuous(limits = c(300, 1100), breaks = seq(300, 1100, by = 100)) +  
  geom_point(data = boot_Pt20_DIC$boots_DIC, aes(y = x, x = y), color = "#D51317FF", shape = 3) +  
  geom_point(data = boot_Pt20_DIC$boots_DIC_AZ, aes(y = x, x = y), color = "#95C11FFF", shape = 3) +  
  geom_point(data = boot_Pt20_DIC$boots_DIC_EZ, aes(y = x, x = y), color = "#0094CDFF", shape = 3) +  
  geom_point(data = assay_table_Pt20, aes(y = estimate_O2max, x = estimate_Km, fill = inhibitor), colour = "black", shape = 23, size = 2) +  
  geom_errorbarh(data = avgcoef_Pt20 %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, xmin = K - Kse, xmax = K + Kse, colour = inhibitor), colour = "black", height = 20, size = 0.4) +
  geom_errorbar(data = avgcoef_Pt20 %>% filter(CoefficientType == "DIC"), aes(x = K, ymin = VmaxO2 - O2se, ymax = VmaxO2 + O2se, colour = inhibitor), colour = "black", width = 15, size = 0.4) +
  geom_point(data = avgcoef_Pt20 %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, x = K, fill = inhibitor), colour = "black", shape = 22, size = 4) +  
  scale_fill_manual(values = c("none" = "#D51317FF", "AZ" = "#95C11FFF", "EZ" = "#0094CDFF")) + 
  scale_colour_manual(values = c("none" = "#D5131755", "AZ" = "#95C11F55", "EZ" = "#0094CD55")) +
  xlab(bquote(K[0.5~~DIC]~'(\u03bcM)')) + 
  ylab(bquote(VO[2]^'max'~'(fmol'~O[2]~.~hr^-1~cell^-1~')')) 

# ---- ANOVA ---------
Km_anova_Pt20 <- aov(estimate_Km ~ inhibitor, data = assay_table_Pt20)
summary(Km_anova_Pt20)
TukeyHSD(Km_anova_Pt20)
qqnorm(residuals(Km_anova_Pt20))
qqline(residuals(Km_anova_Pt20))
shapiro.test(residuals(Km_anova_Pt20))
leveneTest(estimate_Km ~ inhibitor, data = assay_table_Pt20)
kruskal.test(estimate_Km ~ inhibitor, data = assay_table_Pt20)
dunn.test(assay_table_Pt20$estimate_Km, assay_table_Pt20$inhibitor, method="bonferroni")

O2_anova_Pt20 <- aov(estimate_O2max ~ inhibitor, data = assay_table_Pt20)
summary(O2_anova_Pt20)
TukeyHSD(O2_anova_Pt20)
qqnorm(residuals(O2_anova_Pt20))
qqline(residuals(O2_anova_Pt20))
shapiro.test(residuals(O2_anova_Pt20))
leveneTest(estimate_O2max ~ inhibitor, data = assay_table_Pt20)
kruskal.test(estimate_O2max ~ inhibitor, data = assay_table_Pt20)
dunn.test(assay_table_Pt20$estimate_O2max, assay_table_Pt20$inhibitor, method="bonferroni")

# ---- P tricornutum 4 C -----------------------------------------------------------------------------
top = 150 #setting initial conditions for O2max
kco2 = 2 #for half sat CO2
kdic = 50 #for half sat DIC

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=top,Km=kdic)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=top,Km=kco2)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Fit models and retain 'inhibitor' information
runcoefficients_Pt4 <- Pt4 %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_dic = fit_nls_dic(.)) %>%
  unnest_wider(coefs_dic)

runco2efficients_Pt4 <- Pt4 %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_co2 = fit_nls_co2(.)) %>%
  unnest_wider(coefs_co2)

# Assuming runcoefficients_Pt4 and runco2efficients_Pt4 have the same structure
combined_coefficients_Pt4 <- bind_rows(
  runcoefficients_Pt4 %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pt4 %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pt4 <- combined_coefficients_Pt4 %>%
  group_by(CoefficientType, inhibitor) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n()),) %>%
  mutate(across(where(is.numeric), ~round(., 4)))


#fitting MM kinetics to whole data set
Pt4_nls_DIC <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt4[which(Pt4$inhibitor == "none"),], start=list(O2max=top,Km=kdic))
Pt4_nls_DIC_AZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt4[which(Pt4$inhibitor == "AZ"),], start=list(O2max=top,Km=kdic))
Pt4_nls_DIC_EZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pt4[which(Pt4$inhibitor == "EZ"),], start=list(O2max=top,Km=kdic))
Pt4_nls_CO2 <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt4[which(Pt4$inhibitor == "none"),], start=list(O2max=top,Km=kco2))
Pt4_nls_CO2_AZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt4[which(Pt4$inhibitor == "AZ"),], start=list(O2max=top,Km=kco2))
Pt4_nls_CO2_EZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pt4[which(Pt4$inhibitor == "EZ"),], start=list(O2max=top,Km=kco2))
#list the models
models_Pt4 <- list(Pt4_nls_DIC = Pt4_nls_DIC,
                    Pt4_nls_DIC_AZ = Pt4_nls_DIC_AZ,
                    Pt4_nls_DIC_EZ = Pt4_nls_DIC_EZ,
                    Pt4_nls_CO2 = Pt4_nls_CO2,
                    Pt4_nls_CO2_AZ = Pt4_nls_CO2_AZ,
                    Pt4_nls_CO2_EZ = Pt4_nls_CO2_EZ)

# Apply summary and confidence intervals using profiling to each model in the list
model_summaries_Pt4 <- map(models_Pt4, summary)
model_confints_Pt4 <- map(models_Pt4, confint)

# Setting up the plotting area for a grid (3 models x 2 parameters)
par(mfrow = c(3, 2))
# Loop over each model and plot profiles for Km and O2max
for(model in models_Pt4) {
  plot(profile(model, "Km"))
  plot(profile(model, "O2max"))
}
par(mfrow = c(1, 1)) # Resetting to default plotting layout

# Extract and combine coefficients and standard errors
coefficients_Pt4 <- map_df(names(models_Pt4), ~{
  model_data <- tidy(models_Pt4[[.x]])
  species <- ifelse(grepl("DIC", .x), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  mutate(model_data, species = species, inhibitor = inhibitor)
})

# Transforming confidence interval data
confint_data_Pt4 <- map_df(names(model_confints_Pt4), ~{
  confint_data <- as_tibble(model_confints_Pt4[[.x]], rownames = "term") %>%
    mutate(species = ifelse(grepl("DIC", .x), "DIC", "CO2"),
           inhibitor = case_when(
             grepl("AZ", .x) ~ "AZ",
             grepl("EZ", .x) ~ "EZ",
             TRUE ~ "none"))
  confint_data
}, .id = "model")
# joining the data together
coefficients_Pt4 <- left_join(coefficients_Pt4, select(confint_data_Pt4, -model), 
                               by = c("species", "inhibitor", "term"))

# Reshape the table to include estimates and standard errors
table_Pt4 <- coefficients_Pt4 %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
    names_from = term,
    values_from = c(estimate, std.error, "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pt4

#bootstrapping
boot_models_Pt4 <- map(models_Pt4, nlsBoot)

# Function to extract bootCI and estiboot, and add species and inhibitor
extract_boot_data <- function(model_name, model_boot) {
  # Extract and convert bootCI and estiboot data to tibble, capturing row names
  bootCI_data <- as_tibble(model_boot$bootCI, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  estiboot_data <- as_tibble(model_boot$estiboot, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  
  # Combine bootCI and estiboot data
  combined_data <- cbind(bootCI_data, estiboot_data[,-1]) # Exclude repeated coef_name column from estiboot_data
  
  # Add species and inhibitor columns
  species <- ifelse(grepl("DIC", model_name), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", model_name) ~ "AZ",
    grepl("EZ", model_name) ~ "EZ",
    TRUE ~ "none"
  )
  combined_data <- cbind(species = species, inhibitor = inhibitor, combined_data)
  
  return(combined_data)
}

# Apply the function to each bootstrapped model
boot_data_Pt4 <- map2_df(names(boot_models_Pt4), boot_models_Pt4, extract_boot_data)

# View the table
boot_data_Pt4

# Reshape the table to include estimates and standard errors
boot_table_Pt4 <- boot_data_Pt4 %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
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
  coefs <- filter(boot_table_Pt4, species == !!species, inhibitor == !!inhibitor)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MM4_DIC <- function(x) MM_fit_Pt4(x, "DIC", "none")
MMAZ4_DIC <- function(x) MM_fit_Pt4(x, "DIC", "AZ")
MMEZ4_DIC <- function(x) MM_fit_Pt4(x, "DIC", "EZ")
MM4_CO2 <- function(x) MM_fit_Pt4(x, "CO2", "none")
MMAZ4_CO2 <- function(x) MM_fit_Pt4(x, "CO2", "AZ")
MMEZ4_CO2 <- function(x) MM_fit_Pt4(x, "CO2", "EZ")

#functions for confidence intervals
# Helper function to get coefficients and confidence intervals
get_coef_conf_Pt4 <- function(species, inhibitor, coef_name, conf_level) {
  filtered_data <- boot_data_Pt4 %>%
    filter(species == !!species, inhibitor == !!inhibitor, coef_name == !!coef_name)
  result <- filtered_data %>%
    select(!!sym(conf_level)) %>%
    unlist()
  return(result)
}

# Applying the helper function to created MM fits over confidence intervals 
up4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "none", "O2max", "97.5%") * x / (get_coef_conf_Pt4("DIC", "none", "Km", "2.5%") + x))
lw4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "none", "O2max", "2.5%") * x / (get_coef_conf_Pt4("DIC", "none", "Km", "97.5%") + x))
upAZ4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pt4("DIC", "AZ", "Km", "2.5%") + x))
lwAZ4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pt4("DIC", "AZ", "Km", "97.5%") + x))
upEZ4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pt4("DIC", "EZ", "Km", "2.5%") + x))
lwEZ4_DIC <- sapply(Pt4$DICadj, FUN = function(x) get_coef_conf_Pt4("DIC", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pt4("DIC", "EZ", "Km", "97.5%") + x))
up4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "none", "O2max", "97.5%") * x / (get_coef_conf_Pt4("CO2", "none", "Km", "2.5%") + x))
lw4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "none", "O2max", "2.5%") * x / (get_coef_conf_Pt4("CO2", "none", "Km", "97.5%") + x))
upAZ4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pt4("CO2", "AZ", "Km", "2.5%") + x))
lwAZ4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pt4("CO2", "AZ", "Km", "97.5%") + x))
upEZ4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pt4("CO2", "EZ", "Km", "2.5%") + x))
lwEZ4_CO2 <- sapply(Pt4$CO2adj, FUN = function(x) get_coef_conf_Pt4("CO2", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pt4("CO2", "EZ", "Km", "97.5%") + x))

# plotting up in ggplot
CO2_cell_plot_Pt4 <- ggplot(Pt4, aes(x = CO2adj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw4_CO2, ymax = up4_CO2), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZ4_CO2, ymax = upAZ4_CO2), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZ4_CO2, ymax = upEZ4_CO2), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-6, 18), breaks = seq(-6, 18, by = 3)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MM4_CO2, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZ4_CO2, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZ4_CO2, colour= "#0094CDFF", size = 0.7) 


DIC_cell_plot_Pt4 <- ggplot(Pt4, aes(x = DICadj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lw4_DIC, ymax = up4_DIC), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZ4_DIC, ymax = upAZ4_DIC), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZ4_DIC, ymax = upEZ4_DIC), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-6, 18), breaks = seq(-6, 18, by = 3)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MM4_DIC, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZ4_DIC, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZ4_DIC, colour= "#0094CDFF", size = 0.7) 

# ------ conducting model comparisons ----------------------------------------------------
# Create a function to apply nlsConfRegions to a single model
apply_nlsConfRegions <- function(model) {
  nlsConfRegions(model, length = 500, exp = 2)
}

# Apply the function to all models in models_Pt4
Bealeconf_regions_Pt4 <- lapply(models_Pt4, apply_nlsConfRegions)

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
Bealeconf_Pt4_DIC <- list(
  Beale_DIC = data.frame(x = Bealeconf_regions_Pt4$Pt4_nls_DIC$cr[, 1], y = Bealeconf_regions_Pt4$Pt4_nls_DIC$cr[, 2]),
  Beale_DIC_AZ = data.frame(x = Bealeconf_regions_Pt4$Pt4_nls_DIC_AZ$cr[, 1], y = Bealeconf_regions_Pt4$Pt4_nls_DIC_AZ$cr[, 2]),
  Beale_DIC_EZ = data.frame(x = Bealeconf_regions_Pt4$Pt4_nls_DIC_EZ$cr[, 1], y = Bealeconf_regions_Pt4$Pt4_nls_DIC_EZ$cr[, 2])
)

Beale_Pt4 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 2500), breaks = seq(0, 2500, by = 500)) +
  scale_x_continuous(limits = c(60, 180), breaks = seq(60, 180, by = 30)) +
  geom_point(data = Bealeconf_Pt4_DIC$Beale_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = Bealeconf_Pt4_DIC$Beale_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = Bealeconf_Pt4_DIC$Beale_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
boot_Pt4_DIC <- list(
  boots_DIC = data.frame(x = boot_models_Pt4$Pt4_nls_DIC$coefboot[, 1], y = boot_models_Pt4$Pt4_nls_DIC$coefboot[, 2]),
  boots_DIC_AZ = data.frame(x = boot_models_Pt4$Pt4_nls_DIC_AZ$coefboot[, 1], y = boot_models_Pt4$Pt4_nls_DIC_AZ$coefboot[, 2]),
  boots_DIC_EZ = data.frame(x = boot_models_Pt4$Pt4_nls_DIC_EZ$coefboot[, 1], y = boot_models_Pt4$Pt4_nls_DIC_EZ$coefboot[, 2])
)

Bootplot_Pt4 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 2500), breaks = seq(0, 2500, by = 500)) +
  scale_x_continuous(limits = c(60, 180), breaks = seq(60, 180, by = 30)) +
  geom_point(data = boot_Pt4_DIC$boots_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = boot_Pt4_DIC$boots_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = boot_Pt4_DIC$boots_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(mol'~O[2]~.~cell^-1~s^-1~')')) 

# -- fitting models to individual assays for hypothesis testing -------------------
fit_nls_DIC <- function(data) {
  nls_model <- nls(O2rate_cell ~ O2max * DICadj / (Km + DICadj), 
                   data = data,
                   start = list(O2max = top, Km = kdic))
  return(nls_model)
}

# Fitting the model for each assay and recording the associated inhibitor
assay_models_Pt4 <- list()
for (assay in unique(Pt4$assay)) {
  assay_data <- Pt4[Pt4$assay == assay, ]
  inhibitor_used <- unique(assay_data$inhibitor) # Record the inhibitor used
  # Fit the model
  fitted_model <- fit_nls_DIC(assay_data) 
  # Store the results 
  model_key <- paste("Assay", assay, "Inhibitor", inhibitor_used, sep = "_")
  assay_models_Pt4[[model_key]] <- fitted_model
}

#extract summaries for the models 
assay_summaries_Pt4 <- map(assay_models_Pt4, summary)

# Extract and combine coefficients and standard errors
assay_coefficients_Pt4 <- map_df(names(assay_models_Pt4), ~{
  model_data <- tidy(assay_models_Pt4[[.x]])
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  assay_number <- str_extract(.x, "(?<=Assay_)[0-9]+") # Extract the assay number
  mutate(model_data, inhibitor = inhibitor, assay = as.numeric(assay_number))
})

# Reshape the table to include estimates and standard errors
assay_table_Pt4 <- assay_coefficients_Pt4 %>%
  pivot_wider(
    id_cols = c(inhibitor, assay),
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
assay_table_Pt4


Boot_plusdata_Pt4 <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_x_continuous(limits = c(0, 2500), breaks = seq(0, 2500, by = 250)) +
  scale_y_continuous(limits = c(8, 20), breaks = seq(8, 20, by = 2)) +
  geom_point(data = boot_Pt4_DIC$boots_DIC, aes(y = x, x = y), color = "#D51317FF", shape = 3) +
  geom_point(data = boot_Pt4_DIC$boots_DIC_AZ, aes(y = x, x = y), color = "#95C11FFF", shape = 3) +
  geom_point(data = boot_Pt4_DIC$boots_DIC_EZ, aes(y = x, x = y), color = "#0094CDFF", shape = 3) +
  geom_point(data = assay_table_Pt4, aes(y = estimate_O2max, x = estimate_Km, fill = inhibitor), colour = "black", shape = 23, size = 2) +
  geom_errorbarh(data = avgcoef_Pt4 %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, xmin = K - Kse, xmax = K + Kse, colour = inhibitor), colour = "black", height = 0.3, size = 0.4) +
  geom_errorbar(data = avgcoef_Pt4 %>% filter(CoefficientType == "DIC"), aes(x = K, ymin = VmaxO2 - O2se, ymax = VmaxO2 + O2se, colour = inhibitor), colour = "black", width = 62.5, size = 0.4) +
  geom_point(data = avgcoef_Pt4 %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, x = K, fill = inhibitor), colour = "black", shape = 22, size = 4) +
  scale_fill_manual(values = c("none" = "#D51317FF", "AZ" = "#95C11FFF", "EZ" = "#0094CDFF")) + 
  scale_colour_manual(values = c("none" = "#D5131755", "AZ" = "#95C11F55", "EZ" = "#0094CD55")) +
  xlab(bquote(K[0.5~~DIC]~'(\u03bcM)')) + 
  ylab(bquote(VO[2]^'max'~'(fmol'~O[2]~.~hr^-1~cell^-1~')')) 

# ---- ANOVA ---------
Km_anova_Pt4 <- aov(estimate_Km ~ inhibitor, data = assay_table_Pt4)
summary(Km_anova_Pt4)
TukeyHSD(Km_anova_Pt4)
qqnorm(residuals(Km_anova_Pt4))
qqline(residuals(Km_anova_Pt4))
shapiro.test(residuals(Km_anova_Pt4))
leveneTest(estimate_Km ~ inhibitor, data = assay_table_Pt4)
kruskal.test(estimate_Km ~ inhibitor, data = assay_table_Pt4)
dunn.test(assay_table_Pt4$estimate_Km, assay_table_Pt4$inhibitor, method="bonferroni")

O2_anova_Pt4 <- aov(estimate_O2max ~ inhibitor, data = assay_table_Pt4)
summary(O2_anova_Pt4)
TukeyHSD(O2_anova_Pt4)
qqnorm(residuals(O2_anova_Pt4))
qqline(residuals(O2_anova_Pt4))
shapiro.test(residuals(O2_anova_Pt4))
leveneTest(estimate_O2max ~ inhibitor, data = assay_table_Pt4)
kruskal.test(estimate_O2max ~ inhibitor, data = assay_table_Pt4)
dunn.test(assay_table_Pt4$estimate_O2max, assay_table_Pt4$inhibitor, method="bonferroni")


# ---- P inermis 4 C -----------------------------------------------------------------------------
top = 150 #setting initial conditions for O2max
kco2 = 0.5 #for half sat CO2
kdic = 50 #for half sat DIC

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=top,Km=100)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=top,Km=1)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Fit models and retain 'inhibitor' information
runcoefficients_Pi <- Pi %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_dic = fit_nls_dic(.)) %>%
  unnest_wider(coefs_dic)

runco2efficients_Pi <- Pi %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_co2 = fit_nls_co2(.)) %>%
  unnest_wider(coefs_co2)

# Assuming runcoefficients_Pi and runco2efficients_Pi have the same structure
combined_coefficients_Pi <- bind_rows(
  runcoefficients_Pi %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Pi %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Pi <- combined_coefficients_Pi %>%
  group_by(CoefficientType, inhibitor) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n())) %>%
  mutate(across(where(is.numeric), ~round(., 4)))


#fitting MM kinetics to whole data set
Pi_nls_DIC <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pi[which(Pi$inhibitor == "none"),], start=list(O2max=top,Km=kdic))
Pi_nls_DIC_AZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pi[which(Pi$inhibitor == "AZ"),], start=list(O2max=top,Km=kdic))
Pi_nls_DIC_EZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Pi[which(Pi$inhibitor == "EZ"),], start=list(O2max=top,Km=kdic))
Pi_nls_CO2 <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pi[which(Pi$inhibitor == "none"),], start=list(O2max=top,Km=kco2))
Pi_nls_CO2_AZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pi[which(Pi$inhibitor == "AZ"),], start=list(O2max=top,Km=kco2))
Pi_nls_CO2_EZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Pi[which(Pi$inhibitor == "EZ"),], start=list(O2max=top,Km=kco2))
#list the models
models_Pi <- list(Pi_nls_DIC = Pi_nls_DIC,
                   Pi_nls_DIC_AZ = Pi_nls_DIC_AZ,
                   Pi_nls_DIC_EZ = Pi_nls_DIC_EZ,
                   Pi_nls_CO2 = Pi_nls_CO2,
                   Pi_nls_CO2_AZ = Pi_nls_CO2_AZ,
                   Pi_nls_CO2_EZ = Pi_nls_CO2_EZ)

# Apply summary and confidence intervals using profiling to each model in the list
model_summaries_Pi <- map(models_Pi, summary)
model_confints_Pi <- map(models_Pi, confint)

# Setting up the plotting area for a grid (3 models x 2 parameters)
par(mfrow = c(3, 2))
# Loop over each model and plot profiles for Km and O2max
for(model in models_Pi) {
  plot(profile(model, "Km"))
  plot(profile(model, "O2max"))
}
par(mfrow = c(1, 1)) # Resetting to default plotting layout

# Extract and combine coefficients and standard errors
coefficients_Pi <- map_df(names(models_Pi), ~{
  model_data <- tidy(models_Pi[[.x]])
  species <- ifelse(grepl("DIC", .x), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  mutate(model_data, species = species, inhibitor = inhibitor)
})

# Transforming confidence interval data
confint_data_Pi <- map_df(names(model_confints_Pi), ~{
  confint_data <- as_tibble(model_confints_Pi[[.x]], rownames = "term") %>%
    mutate(species = ifelse(grepl("DIC", .x), "DIC", "CO2"),
           inhibitor = case_when(
             grepl("AZ", .x) ~ "AZ",
             grepl("EZ", .x) ~ "EZ",
             TRUE ~ "none"))
  confint_data
}, .id = "model")
# joining the data together
coefficients_Pi <- left_join(coefficients_Pi, select(confint_data_Pi, -model), 
                              by = c("species", "inhibitor", "term"))

# Reshape the table to include estimates and standard errors
table_Pi <- coefficients_Pi %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
    names_from = term,
    values_from = c(estimate, std.error, "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Pi

#bootstrapping
boot_models_Pi <- map(models_Pi, nlsBoot)

# Function to extract bootCI and estiboot, and add species and inhibitor
extract_boot_data <- function(model_name, model_boot) {
  # Extract and convert bootCI and estiboot data to tibble, capturing row names
  bootCI_data <- as_tibble(model_boot$bootCI, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  estiboot_data <- as_tibble(model_boot$estiboot, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  
  # Combine bootCI and estiboot data
  combined_data <- cbind(bootCI_data, estiboot_data[,-1]) # Exclude repeated coef_name column from estiboot_data
  
  # Add species and inhibitor columns
  species <- ifelse(grepl("DIC", model_name), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", model_name) ~ "AZ",
    grepl("EZ", model_name) ~ "EZ",
    TRUE ~ "none"
  )
  combined_data <- cbind(species = species, inhibitor = inhibitor, combined_data)
  
  return(combined_data)
}

# Apply the function to each bootstrapped model
boot_data_Pi <- map2_df(names(boot_models_Pi), boot_models_Pi, extract_boot_data)

# View the table
boot_data_Pi

# Reshape the table to include estimates and standard errors
boot_table_Pi <- boot_data_Pi %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
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
  coefs <- filter(boot_table_Pi, species == !!species, inhibitor == !!inhibitor)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MMPi_DIC <- function(x) MM_fit_Pi(x, "DIC", "none")
MMAZPi_DIC <- function(x) MM_fit_Pi(x, "DIC", "AZ")
MMEZPi_DIC <- function(x) MM_fit_Pi(x, "DIC", "EZ")
MMPi_CO2 <- function(x) MM_fit_Pi(x, "CO2", "none")
MMAZPi_CO2 <- function(x) MM_fit_Pi(x, "CO2", "AZ")
MMEZPi_CO2 <- function(x) MM_fit_Pi(x, "CO2", "EZ")

#functions for confidence intervals
# Helper function to get coefficients and confidence intervals
get_coef_conf_Pi <- function(species, inhibitor, coef_name, conf_level) {
  filtered_data <- boot_data_Pi %>%
    filter(species == !!species, inhibitor == !!inhibitor, coef_name == !!coef_name)
  result <- filtered_data %>%
    select(!!sym(conf_level)) %>%
    unlist()
  return(result)
}

# Applying the helper function to created MM fits over confidence intervals 
upPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "none", "O2max", "97.5%") * x / (get_coef_conf_Pi("DIC", "none", "Km", "2.5%") + x))
lwPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "none", "O2max", "2.5%") * x / (get_coef_conf_Pi("DIC", "none", "Km", "97.5%") + x))
upAZPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pi("DIC", "AZ", "Km", "2.5%") + x))
lwAZPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pi("DIC", "AZ", "Km", "97.5%") + x))
upEZPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pi("DIC", "EZ", "Km", "2.5%") + x))
lwEZPi_DIC <- sapply(Pi$DICadj, FUN = function(x) get_coef_conf_Pi("DIC", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pi("DIC", "EZ", "Km", "97.5%") + x))
upPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "none", "O2max", "97.5%") * x / (get_coef_conf_Pi("CO2", "none", "Km", "2.5%") + x))
lwPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "none", "O2max", "2.5%") * x / (get_coef_conf_Pi("CO2", "none", "Km", "97.5%") + x))
upAZPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Pi("CO2", "AZ", "Km", "2.5%") + x))
lwAZPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Pi("CO2", "AZ", "Km", "97.5%") + x))
upEZPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Pi("CO2", "EZ", "Km", "2.5%") + x))
lwEZPi_CO2 <- sapply(Pi$CO2adj, FUN = function(x) get_coef_conf_Pi("CO2", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Pi("CO2", "EZ", "Km", "97.5%") + x))


# plotting up in ggplot
CO2_cell_plot_Pi <- ggplot(Pi, aes(x = CO2adj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwPi_CO2, ymax = upPi_CO2), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZPi_CO2, ymax = upAZPi_CO2), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZPi_CO2, ymax = upEZPi_CO2), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-300, 900), breaks = seq(-300, 900, by = 150)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MMPi_CO2, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZPi_CO2, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZPi_CO2, colour= "#0094CDFF", size = 0.7) 


DIC_cell_plot_Pi <- ggplot(Pi, aes(x = DICadj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwPi_DIC, ymax = upPi_DIC), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZPi_DIC, ymax = upAZPi_DIC), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZPi_DIC, ymax = upEZPi_DIC), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-300, 900), breaks = seq(-300, 900, by = 150)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MMPi_DIC, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZPi_DIC, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZPi_DIC, colour= "#0094CDFF", size = 0.7) 
  

# ------ conducting model comparisons ----------------------------------------------------
# Create a function to apply nlsConfRegions to a single model
apply_nlsConfRegions <- function(model) {
  nlsConfRegions(model, length = 500, exp = 2)
}

# Apply the function to all models in models_Pi
Bealeconf_regions_Pi <- lapply(models_Pi, apply_nlsConfRegions)

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
Bealeconf_Pi_DIC <- list(
  Beale_DIC = data.frame(x = Bealeconf_regions_Pi$Pi_nls_DIC$cr[, 1], y = Bealeconf_regions_Pi$Pi_nls_DIC$cr[, 2]),
  Beale_DIC_AZ = data.frame(x = Bealeconf_regions_Pi$Pi_nls_DIC_AZ$cr[, 1], y = Bealeconf_regions_Pi$Pi_nls_DIC_AZ$cr[, 2]),
  Beale_DIC_EZ = data.frame(x = Bealeconf_regions_Pi$Pi_nls_DIC_EZ$cr[, 1], y = Bealeconf_regions_Pi$Pi_nls_DIC_EZ$cr[, 2])
)

Beale_Pi <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, by = 500)) +
  scale_x_continuous(limits = c(100, 300), breaks = seq(90, 180, by = 30)) +
  geom_point(data = Bealeconf_Pi_DIC$Beale_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = Bealeconf_Pi_DIC$Beale_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = Bealeconf_Pi_DIC$Beale_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
boot_Pi_DIC <- list(
  boots_DIC = data.frame(x = boot_models_Pi$Pi_nls_DIC$coefboot[, 1], y = boot_models_Pi$Pi_nls_DIC$coefboot[, 2]),
  boots_DIC_AZ = data.frame(x = boot_models_Pi$Pi_nls_DIC_AZ$coefboot[, 1], y = boot_models_Pi$Pi_nls_DIC_AZ$coefboot[, 2]),
  boots_DIC_EZ = data.frame(x = boot_models_Pi$Pi_nls_DIC_EZ$coefboot[, 1], y = boot_models_Pi$Pi_nls_DIC_EZ$coefboot[, 2])
)

Bootplot_Pi <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, by = 500)) +
  scale_x_continuous(limits = c(90, 180), breaks = seq(90, 180, by = 30)) +
  geom_point(data = boot_Pi_DIC$boots_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = boot_Pi_DIC$boots_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = boot_Pi_DIC$boots_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# -- fitting models to individual assays for hypothesis testing -------------------
fit_nls_DIC <- function(data) {
  nls_model <- nls(O2rate_cell ~ O2max * DICadj / (Km + DICadj), 
                   data = data,
                   start = list(O2max = top, Km = kdic))
  return(nls_model)
}

# Fitting the model for each assay and recording the associated inhibitor
assay_models_Pi <- list()
for (assay in unique(Pi$assay)) {
  assay_data <- Pi[Pi$assay == assay, ]
  inhibitor_used <- unique(assay_data$inhibitor) # Record the inhibitor used
  # Fit the model
  fitted_model <- fit_nls_DIC(assay_data) 
  # Store the results 
  model_key <- paste("Assay", assay, "Inhibitor", inhibitor_used, sep = "_")
  assay_models_Pi[[model_key]] <- fitted_model
}

#extract summaries for the models 
assay_summaries_Pi <- map(assay_models_Pi, summary)

# Extract and combine coefficients and standard errors
assay_coefficients_Pi <- map_df(names(assay_models_Pi), ~{
  model_data <- tidy(assay_models_Pi[[.x]])
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  assay_number <- str_extract(.x, "(?<=Assay_)[0-9]+") # Extract the assay number
  mutate(model_data, inhibitor = inhibitor, assay = as.numeric(assay_number))
})

# Reshape the table to include estimates and standard errors
assay_table_Pi <- assay_coefficients_Pi %>%
  pivot_wider(
    id_cols = c(inhibitor, assay),
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
assay_table_Pi


Boot_plusdata_Pi <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, by = 200)) +  
  scale_y_continuous(limits = c(300, 1100), breaks = seq(300, 1100, by = 100)) +  
  geom_point(data = boot_Pi_DIC$boots_DIC, aes(y = x, x = y), color = "#D51317FF", shape = 3) +  
  geom_point(data = boot_Pi_DIC$boots_DIC_AZ, aes(y = x, x = y), color = "#95C11FFF", shape = 3) +  
  geom_point(data = boot_Pi_DIC$boots_DIC_EZ, aes(y = x, x = y), color = "#0094CDFF", shape = 3) +  
  geom_point(data = assay_table_Pi, aes(y = estimate_O2max, x = estimate_Km, fill = inhibitor), colour = "black", shape = 23, size = 2) +  
  geom_errorbarh(data = avgcoef_Pi %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, xmin = K - Kse, xmax = K + Kse, colour = inhibitor), colour = "black", height = 20, size = 0.4) +
  geom_errorbar(data = avgcoef_Pi %>% filter(CoefficientType == "DIC"), aes(x = K, ymin = VmaxO2 - O2se, ymax = VmaxO2 + O2se, colour = inhibitor), colour = "black", width = 40, size = 0.4) +
  geom_point(data = avgcoef_Pi %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, x = K, fill = inhibitor), colour = "black", shape = 22, size = 4) +  # Swapped x and y
  scale_fill_manual(values = c("none" = "#D51317FF", "AZ" = "#95C11FFF", "EZ" = "#0094CDFF")) + 
  scale_colour_manual(values = c("none" = "#D5131755", "AZ" = "#95C11F55", "EZ" = "#0094CD55")) +
  xlab(bquote(K[0.5~~DIC]~'(\u03bcM)')) + 
  ylab(bquote(VO[2]^'max'~'(fmol'~O[2]~.~hr^-1~cell^-1~')')) 


# ---- ANOVA ---------
Km_anova_Pi <- aov(estimate_Km ~ inhibitor, data = assay_table_Pi)
summary(Km_anova_Pi)
TukeyHSD(Km_anova_Pi)
qqnorm(residuals(Km_anova_Pi))
qqline(residuals(Km_anova_Pi))
shapiro.test(residuals(Km_anova_Pi))
leveneTest(estimate_Km ~ inhibitor, data = assay_table_Pi)
kruskal.test(estimate_Km ~ inhibitor, data = assay_table_Pi)
dunn.test(assay_table_Pi$estimate_Km, assay_table_Pi$inhibitor, method="bonferroni")

O2_anova_Pi <- aov(estimate_O2max ~ inhibitor, data = assay_table_Pi)
summary(O2_anova_Pi)
TukeyHSD(O2_anova_Pi)
qqnorm(residuals(O2_anova_Pi))
qqline(residuals(O2_anova_Pi))
shapiro.test(residuals(O2_anova_Pi))
leveneTest(estimate_O2max ~ inhibitor, data = assay_table_Pi)
kruskal.test(estimate_O2max ~ inhibitor, data = assay_table_Pi)
dunn.test(assay_table_Pi$estimate_O2max, assay_table_Pi$inhibitor, method="bonferroni")

# ---- C flexuosus 4 C -----------------------------------------------------------------------------
top = 90 #setting initial conditions for O2max
kco2 = 0.5 #for half sat CO2
kdic = 50 #for half sat DIC

# Function to fit nls model for DIC
fit_nls_dic <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), df, start=list(O2max=top,Km=100)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Function to fit nls model for CO2
fit_nls_co2 <- function(df) {
  bestfit <- tryCatch(
    nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), df, start=list(O2max=top,Km=1)),
    error = function(e) return(NULL))
  
  if (!is.null(bestfit) && inherits(bestfit, "nls")) {
    return(coef(bestfit))
  } else {return(rep(NA, length(coef(bestfit)))) }}

# Fit models and retain 'inhibitor' information
runcoefficients_Cf <- Cf %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_dic = fit_nls_dic(.)) %>%
  unnest_wider(coefs_dic)

runco2efficients_Cf <- Cf %>%
  filter(!is.na(O2rate_cell)) %>%
  mutate(inhibitor = as.character(inhibitor)) %>%  # Temporarily convert to character
  group_by(assay, inhibitor) %>%
  do(coefs_co2 = fit_nls_co2(.)) %>%
  unnest_wider(coefs_co2)

# Assuming runcoefficients_Cf and runco2efficients_Cf have the same structure
combined_coefficients_Cf <- bind_rows(
  runcoefficients_Cf %>% mutate(CoefficientType = "DIC"),
  runco2efficients_Cf %>% mutate(CoefficientType = "CO2"))

# Calculate the averages and standard deviations
avgcoef_Cf <- combined_coefficients_Cf %>%
  group_by(CoefficientType, inhibitor) %>%
  summarize(
    VmaxO2 = mean(O2max, na.rm = TRUE),
    O2sd = sd(O2max, na.rm = TRUE),
    O2se = sd(O2max, na.rm = TRUE) / sqrt(n()),
    K = mean(Km, na.rm = TRUE),
    Ksd = sd(Km, na.rm = TRUE),
    Kse = sd(Km, na.rm = TRUE) / sqrt(n()),) %>%
  mutate(across(where(is.numeric), ~round(., 4)))


#fitting MM kinetics to whole data set
Cf_nls_DIC <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Cf[which(Cf$inhibitor == "none"),], start=list(O2max=top,Km=kdic))
Cf_nls_DIC_AZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Cf[which(Cf$inhibitor == "AZ"),], start=list(O2max=top,Km=kdic))
Cf_nls_DIC_EZ <-  nls(O2rate_cell ~ O2max*DICadj/(Km+DICadj), Cf[which(Cf$inhibitor == "EZ"),], start=list(O2max=top,Km=kdic))
Cf_nls_CO2 <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Cf[which(Cf$inhibitor == "none"),], start=list(O2max=top,Km=kco2))
Cf_nls_CO2_AZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Cf[which(Cf$inhibitor == "AZ"),], start=list(O2max=top,Km=kco2))
Cf_nls_CO2_EZ <-  nls(O2rate_cell ~ O2max*CO2adj/(Km+CO2adj), Cf[which(Cf$inhibitor == "EZ"),], start=list(O2max=top,Km=kco2))
#list the models
models_Cf <- list(Cf_nls_DIC = Cf_nls_DIC,
                   Cf_nls_DIC_AZ = Cf_nls_DIC_AZ,
                   Cf_nls_DIC_EZ = Cf_nls_DIC_EZ,
                   Cf_nls_CO2 = Cf_nls_CO2,
                   Cf_nls_CO2_AZ = Cf_nls_CO2_AZ,
                   Cf_nls_CO2_EZ = Cf_nls_CO2_EZ)

# Apply summary and confidence intervals using profiling to each model in the list
model_summaries_Cf <- map(models_Cf, summary)
model_confints_Cf <- map(models_Cf, confint)

# Setting up the plotting area for a grid (3 models x 2 parameters)
par(mfrow = c(3, 2))
# Loop over each model and plot profiles for Km and O2max
for(model in models_Cf) {
  plot(profile(model, "Km"))
  plot(profile(model, "O2max"))
}
par(mfrow = c(1, 1)) # Resetting to default plotting layout

# Extract and combine coefficients and standard errors
coefficients_Cf <- map_df(names(models_Cf), ~{
  model_data <- tidy(models_Cf[[.x]])
  species <- ifelse(grepl("DIC", .x), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  mutate(model_data, species = species, inhibitor = inhibitor)
})

# Transforming confidence interval data
confint_data_Cf <- map_df(names(model_confints_Cf), ~{
  confint_data <- as_tibble(model_confints_Cf[[.x]], rownames = "term") %>%
    mutate(species = ifelse(grepl("DIC", .x), "DIC", "CO2"),
           inhibitor = case_when(
             grepl("AZ", .x) ~ "AZ",
             grepl("EZ", .x) ~ "EZ",
             TRUE ~ "none"))
  confint_data
}, .id = "model")
# joining the data together
coefficients_Cf <- left_join(coefficients_Cf, select(confint_data_Cf, -model), 
                              by = c("species", "inhibitor", "term"))

# Reshape the table to include estimates and standard errors
table_Cf <- coefficients_Cf %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
    names_from = term,
    values_from = c(estimate, std.error, "2.5%", "97.5%"),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
table_Cf

#bootstrapping
boot_models_Cf <- map(models_Cf, nlsBoot)

# Function to extract bootCI and estiboot, and add species and inhibitor
extract_boot_data <- function(model_name, model_boot) {
  # Extract and convert bootCI and estiboot data to tibble, capturing row names
  bootCI_data <- as_tibble(model_boot$bootCI, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  estiboot_data <- as_tibble(model_boot$estiboot, rownames = "coef_name") %>%
    rename_with(~ ifelse(. == "coef_name", "coef_name", .))
  
  # Combine bootCI and estiboot data
  combined_data <- cbind(bootCI_data, estiboot_data[,-1]) # Exclude repeated coef_name column from estiboot_data
  
  # Add species and inhibitor columns
  species <- ifelse(grepl("DIC", model_name), "DIC", "CO2")
  inhibitor <- case_when(
    grepl("AZ", model_name) ~ "AZ",
    grepl("EZ", model_name) ~ "EZ",
    TRUE ~ "none"
  )
  combined_data <- cbind(species = species, inhibitor = inhibitor, combined_data)
  
  return(combined_data)
}

# Apply the function to each bootstrapped model
boot_data_Cf <- map2_df(names(boot_models_Cf), boot_models_Cf, extract_boot_data)

# View the table
boot_data_Cf

# Reshape the table to include estimates and standard errors
boot_table_Cf <- boot_data_Cf %>%
  pivot_wider(
    id_cols = c(species, inhibitor),
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
  coefs <- filter(boot_table_Cf, species == !!species, inhibitor == !!inhibitor)
  coefs$Estimate_O2max * x / (coefs$Estimate_Km + x)
}

# Now, applying this to each C-species and inhibitor combo
MMCf_DIC <- function(x) MM_fit_Cf(x, "DIC", "none")
MMAZCf_DIC <- function(x) MM_fit_Cf(x, "DIC", "AZ")
MMEZCf_DIC <- function(x) MM_fit_Cf(x, "DIC", "EZ")
MMCf_CO2 <- function(x) MM_fit_Cf(x, "CO2", "none")
MMAZCf_CO2 <- function(x) MM_fit_Cf(x, "CO2", "AZ")
MMEZCf_CO2 <- function(x) MM_fit_Cf(x, "CO2", "EZ")

#functions for confidence intervals
# Helper function to get coefficients and confidence intervals
get_coef_conf_Cf <- function(species, inhibitor, coef_name, conf_level) {
  filtered_data <- boot_data_Cf %>%
    filter(species == !!species, inhibitor == !!inhibitor, coef_name == !!coef_name)
  result <- filtered_data %>%
    select(!!sym(conf_level)) %>%
    unlist()
  return(result)
}

# Applying the helper function to created MM fits over confidence intervals 
upCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "none", "O2max", "97.5%") * x / (get_coef_conf_Cf("DIC", "none", "Km", "2.5%") + x))
lwCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "none", "O2max", "2.5%") * x / (get_coef_conf_Cf("DIC", "none", "Km", "97.5%") + x))
upAZCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Cf("DIC", "AZ", "Km", "2.5%") + x))
lwAZCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Cf("DIC", "AZ", "Km", "97.5%") + x))
upEZCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Cf("DIC", "EZ", "Km", "2.5%") + x))
lwEZCf_DIC <- sapply(Cf$DICadj, FUN = function(x) get_coef_conf_Cf("DIC", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Cf("DIC", "EZ", "Km", "97.5%") + x))
upCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "none", "O2max", "97.5%") * x / (get_coef_conf_Cf("CO2", "none", "Km", "2.5%") + x))
lwCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "none", "O2max", "2.5%") * x / (get_coef_conf_Cf("CO2", "none", "Km", "97.5%") + x))
upAZCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "AZ", "O2max", "97.5%") * x / (get_coef_conf_Cf("CO2", "AZ", "Km", "2.5%") + x))
lwAZCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "AZ", "O2max", "2.5%") * x / (get_coef_conf_Cf("CO2", "AZ", "Km", "97.5%") + x))
upEZCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "EZ", "O2max", "97.5%") * x / (get_coef_conf_Cf("CO2", "EZ", "Km", "2.5%") + x))
lwEZCf_CO2 <- sapply(Cf$CO2adj, FUN = function(x) get_coef_conf_Cf("CO2", "EZ", "O2max", "2.5%") * x / (get_coef_conf_Cf("CO2", "EZ", "Km", "97.5%") + x))


# plotting up in ggplot
CO2_cell_plot_Cf <- ggplot(Cf, aes(x = CO2adj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwCf_DIC, ymax = upCf_DIC), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZCf_DIC, ymax = upAZCf_DIC), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZCf_DIC, ymax = upEZCf_DIC), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('CO2 (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-20, 60), breaks = seq(-20, 60, by = 10)) +
  scale_x_continuous(limits = c(0, 29.513), breaks = seq(0, 30, by = 5)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MMCf_CO2, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZCf_CO2, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZCf_CO2, colour= "#0094CDFF", size = 0.7) 
 

DIC_cell_plot_Cf <- ggplot(Cf, aes(x = DICadj, y = O2rate_cell)) +
  theme_bw() +
  geom_ribbon(aes(ymin = lwCf_DIC, ymax = upCf_DIC), fill = "#D51317", colour= NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwAZCf_DIC, ymax = upAZCf_DIC), fill = "#95C11F", colour = NA, alpha = 0.1) +
  geom_ribbon(aes(ymin = lwEZCf_DIC, ymax = upEZCf_DIC), fill = "#0094CD", colour= NA, alpha = 0.1) +
  geom_point(aes(fill = inhibitor, shape = inhibitor), size = 2,colour = "black", ) +
  scale_fill_manual(values = c("#D51317FF" ,"#95C11FFF","#0094CDFF")) +
  scale_shape_manual(values = c(23,23,23)) +
  ylab(bquote('mol'~O[2]~.~cell^-1~s^-1)) + 
  xlab(bquote('DIC (\u03bcmol'~L^-1~')')) +
  scale_y_continuous(limits = c(-20, 60), breaks = seq(-20, 60, by = 10)) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by = 500)) +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  # plot the lines    
  stat_function(fun = MMCf_DIC, colour = "#D51317FF", size = 0.7) +
  stat_function(fun = MMAZCf_DIC, colour = "#95C11FFF", size = 0.7) +
  stat_function(fun = MMEZCf_DIC, colour= "#0094CDFF", size = 0.7) 

# ------ conducting model comparisons ----------------------------------------------------
# Create a function to apply nlsConfRegions to a single model
apply_nlsConfRegions <- function(model) {
  nlsConfRegions(model, length = 500, exp = 2)
}

# Apply the function to all models in models_Cf
Bealeconf_regions_Cf <- lapply(models_Cf, apply_nlsConfRegions)

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
Bealeconf_Cf_DIC <- list(
  Beale_DIC = data.frame(x = Bealeconf_regions_Cf$Cf_nls_DIC$cr[, 1], y = Bealeconf_regions_Cf$Cf_nls_DIC$cr[, 2]),
  Beale_DIC_AZ = data.frame(x = Bealeconf_regions_Cf$Cf_nls_DIC_AZ$cr[, 1], y = Bealeconf_regions_Cf$Cf_nls_DIC_AZ$cr[, 2]),
  Beale_DIC_EZ = data.frame(x = Bealeconf_regions_Cf$Cf_nls_DIC_EZ$cr[, 1], y = Bealeconf_regions_Cf$Cf_nls_DIC_EZ$cr[, 2])
)

Beale_Cf <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 1000)) +
  scale_x_continuous(limits = c(4, 24), breaks = seq(4, 24, by = 2)) +
  geom_point(data = Bealeconf_Cf_DIC$Beale_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = Bealeconf_Cf_DIC$Beale_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = Bealeconf_Cf_DIC$Beale_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# Extract the data from the confidence regions for DIC, DIC_AZ, and DIC_EZ
boot_Cf_DIC <- list(
  boots_DIC = data.frame(x = boot_models_Cf$Cf_nls_DIC$coefboot[, 1], y = boot_models_Cf$Cf_nls_DIC$coefboot[, 2]),
  boots_DIC_AZ = data.frame(x = boot_models_Cf$Cf_nls_DIC_AZ$coefboot[, 1], y = boot_models_Cf$Cf_nls_DIC_AZ$coefboot[, 2]),
  boots_DIC_EZ = data.frame(x = boot_models_Cf$Cf_nls_DIC_EZ$coefboot[, 1], y = boot_models_Cf$Cf_nls_DIC_EZ$coefboot[, 2])
)

Bootplot_Cf <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, by = 1000)) +
  scale_x_continuous(limits = c(6, 22), breaks = seq(6, 22, by = 2)) +
  geom_point(data = boot_Cf_DIC$boots_DIC, aes(x, y), color = "#D51317FF", shape = 1) +
  geom_point(data = boot_Cf_DIC$boots_DIC_AZ, aes(x, y), color = "#95C11FFF", shape = 2) +
  geom_point(data = boot_Cf_DIC$boots_DIC_EZ, aes(x, y), color = "#0094CDFF", shape = 3) +
  ylab(bquote(K[0.5~~DIC])) + 
  xlab(bquote(VO[2]^'max'~'(\u03bcmol'~O[2]~hr^-1~mg~Chl~italic(a)^-1~')')) 

# -- fitting models to individual assays for hypothesis testing -------------------
fit_nls_DIC <- function(data) {
  nls_model <- nls(O2rate_cell ~ O2max * DICadj / (Km + DICadj), 
                   data = data,
                   start = list(O2max = top, Km = kdic))
  return(nls_model)
}

# Fitting the model for each assay and recording the associated inhibitor
assay_models_Cf <- list()
for (assay in unique(Cf$assay)) {
  assay_data <- Cf[Cf$assay == assay, ]
  inhibitor_used <- unique(assay_data$inhibitor) # Record the inhibitor used
  # Fit the model
  fitted_model <- fit_nls_DIC(assay_data) 
  # Store the results 
  model_key <- paste("Assay", assay, "Inhibitor", inhibitor_used, sep = "_")
  assay_models_Cf[[model_key]] <- fitted_model
}

#extract summaries for the models 
assay_summaries_Cf <- map(assay_models_Cf, summary)

# Extract and combine coefficients and standard errors
assay_coefficients_Cf <- map_df(names(assay_models_Cf), ~{
  model_data <- tidy(assay_models_Cf[[.x]])
  inhibitor <- case_when(
    grepl("AZ", .x) ~ "AZ",
    grepl("EZ", .x) ~ "EZ",
    TRUE ~ "none")
  assay_number <- str_extract(.x, "(?<=Assay_)[0-9]+") # Extract the assay number
  mutate(model_data, inhibitor = inhibitor, assay = as.numeric(assay_number))
})

# Reshape the table to include estimates and standard errors
assay_table_Cf <- assay_coefficients_Cf %>%
  pivot_wider(
    id_cols = c(inhibitor, assay),
    names_from = term,
    values_from = c(estimate, std.error),
    names_sep = "_"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  as_tibble()
# View the table
assay_table_Cf

Boot_plusdata_Cf <- ggplot() +
  theme_bw() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),  
        panel.grid.minor = element_blank(), 
        legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "none") +
  scale_x_continuous(limits = c(0, 8000), breaks = seq(0, 8000, by = 1000)) +  
  scale_y_continuous(limits = c(20, 80), breaks = seq(20, 80, by = 10)) +  
  geom_point(data = boot_Cf_DIC$boots_DIC, aes(y = x, x = y), color = "#D51317FF", shape = 3) +  
  geom_point(data = boot_Cf_DIC$boots_DIC_AZ, aes(y = x, x = y), color = "#95C11FFF", shape = 3) +  
  geom_point(data = boot_Cf_DIC$boots_DIC_EZ, aes(y = x, x = y), color = "#0094CDFF", shape = 3) +  
  geom_point(data = assay_table_Cf, aes(y = estimate_O2max, x = estimate_Km, fill = inhibitor), colour = "black", shape = 23, size = 2) +  
  geom_errorbarh(data = avgcoef_Cf %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, xmin = K - Kse, xmax = K + Kse, colour = inhibitor), colour = "black", height = 1.5, size = 0.4) +
  geom_errorbar(data = avgcoef_Cf %>% filter(CoefficientType == "DIC"), aes(x = K, ymin = VmaxO2 - O2se, ymax = VmaxO2 + O2se, colour = inhibitor), colour = "black", width = 200, size = 0.4) + geom_point(data = avgcoef_Cf %>% filter(CoefficientType == "DIC"), aes(y = VmaxO2, x = K, fill = inhibitor), colour = "black", shape = 22, size = 4) +  
  scale_fill_manual(values = c("none" = "#D51317FF", "AZ" = "#95C11FFF", "EZ" = "#0094CDFF")) + 
  scale_colour_manual(values = c("none" = "#D5131755", "AZ" = "#95C11F55", "EZ" = "#0094CD55")) +
  xlab(bquote(K[0.5~~DIC]~'(\u03bcM)')) + 
  ylab(bquote(VO[2]^'max'~'(fmol'~O[2]~.~hr^-1~cell^-1~')')) 

# ---- ANOVA ---------
Km_anova_Cf <- aov(estimate_Km ~ inhibitor, data = assay_table_Cf)
summary(Km_anova_Cf)
TukeyHSD(Km_anova_Cf)
qqnorm(residuals(Km_anova_Cf))
qqline(residuals(Km_anova_Cf))
shapiro.test(residuals(Km_anova_Cf))
leveneTest(estimate_Km ~ inhibitor, data = assay_table_Cf)
kruskal.test(estimate_Km ~ inhibitor, data = assay_table_Cf)
dunn.test(assay_table_Cf$estimate_Km, assay_table_Cf$inhibitor, method="bonferroni")

O2_anova_Cf <- aov(estimate_O2max ~ inhibitor, data = assay_table_Cf)
summary(O2_anova_Cf)
TukeyHSD(O2_anova_Cf)
qqnorm(residuals(O2_anova_Cf))
qqline(residuals(O2_anova_Cf))
shapiro.test(residuals(O2_anova_Cf))
leveneTest(estimate_O2max ~ inhibitor, data = assay_table_Cf)
kruskal.test(estimate_O2max ~ inhibitor, data = assay_table_Cf)
dunn.test(assay_table_Cf$estimate_O2max, assay_table_Cf$inhibitor, method="bonferroni")

# --- Plotting plots all together -----------------

ACi_DIC_cell <- guide_area() / (DIC_cell_plot_Pt20 | DIC_cell_plot_Pt4) / (DIC_cell_plot_Pi | DIC_cell_plot_Cf) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ACi_CO2_cell <- guide_area() / (CO2_cell_plot_Pt20 | CO2_cell_plot_Pt4) / (CO2_cell_plot_Pi | CO2_cell_plot_Cf) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("ACi_cell_DIC_plot_2.svg", ACi_DIC_cell, width = 7.2, height = 7.4)
ggsave("ACi_cell_CO2_plot_2.svg", ACi_CO2_cell, width = 7.2, height = 7.4)

MM_fit_plots_cell <- guide_area() / (Boot_plusdata_Pt20 | Boot_plusdata_Pt4) / (Boot_plusdata_Pi | Boot_plusdata_Cf) +
  plot_layout(guides = 'collect', nrow = (4), heights = c(1,10,10)) +
  plot_annotation(tag_levels = list(c('a','b','c', 'd'))) &
  theme(axis.text = element_text(size = 10, colour = "black"), 
        legend.title = element_text(size = 9), legend.text = element_text(size = 9), 
        legend.direction = "horizontal", legend.position = "top",
        plot.tag.position = c(0.01, 0.97)) & 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

ggsave("MM_fit_boot_plot_cell.svg", MM_fit_plots_cell, width = 7.2, height = 7.4)
