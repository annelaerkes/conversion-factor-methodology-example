# Calculate tissue contaminant conversion factors (k) ----
#
# Example on the calculation of tissue conversion factors used for the Integrative 
# Chemical and Biological Effects Monitoring project based on data from
# Faxneld et al. 2015, Danielsson et al. 2018 and Larsen 2025.
#
# Pipeline: read eelpout, perch and herring liver/muscle data -> harmonize columns
# and filter to a shared set of parameters -> bind into one table -> calculate the
# liver:muscle conversion factor (k) per parameter and species, following the
# geometric-mean method of Soerensen et al. 2023.


library(tidyverse)
library(readxl)
library(stringr)
library(ggplot2)

## read data ----
###  perch data, Faxneld et al. 2015 ----
perch <- read_excel("Data/perch_liver_muscle_long.xlsx", skip = 2) |>
  rename(Sample = 'specimen_no',
         Parameter = 'Metal',
         Liver = 'L',
         Muscle = 'M') |>
  mutate(Parameter = str_to_sentence(Parameter),  
         Species = 'Perch',
         Sample = as.character(Sample)) |>
  filter(Parameter == 'As' | Parameter == 'Se' | Parameter == 'Cu' | Parameter == 'Zn') |>
  # only include samples with both liver and muscle data above LOQ
  filter(Muscle > 0, Liver > 0)

### herring data, Danielsson et al. 2018 ----
herring <- read_excel("Data/herring_liver_muscle.xlsx", skip = 2) |>
  select(Parameter, Sample, Muscle, Liver) |>
  # only include samples with both liver and muscle data above LOQ
  drop_na() |>
  filter(Parameter == 'As' | Parameter == 'Se' | Parameter == 'Cu' | Parameter == 'Zn') |>
  mutate(Site = 'Lilla Vartan',
         Species = 'Herring')

### eelpout and perch data, Larsen 2025 ----
perch_eelpout <- read_excel("Data/perch_eelpout_liver_muscle.xlsx", skip = 2) |>
  rename(Parameter = 'contaminant',
         Site = 'samplingsite',
         Species = 'species') |>
  mutate(across(where(is.character), ~ str_replace_all(.x, fixed("<"), "-"))) |>
  filter(Parameter == 'As' | Parameter == 'Se' | Parameter == 'Cu' | Parameter == 'Zn',
         # only include samples with both liver and muscle data above LOQ
         muscle > 0, liver > 0) |>
  mutate(Muscle = (as.numeric(muscle))*DW_muscle/100, # from ww to dw 
         Liver = (as.numeric(liver)*DW_liver/100),
         Sample = as.character(ID)) |> # from ww to dw 
  select(Species, Parameter, Sample, Site, Muscle, Liver) 


## bind data and calculate k----
ind_k <- bind_rows(perch, herring, perch_eelpout) |>
  select(-c(Sample, Site)) |>
  mutate(log_k = log(Liver/Muscle), # equivalent to: (log_k = log(Liver) - log(Muscle)
         k_l_m_ind = Liver/Muscle)
  
## calculate average k per species and parameter ----
avg_k <- ind_k |> 
  group_by(Parameter,Species) |>
  summarise(
        across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
        n = n(),
        .groups = "drop"
  ) |>
  # producing the average k (liver:muscle ratio) 
  # multiply muscle by k to estimate liver, or divide liver by k to estimate muscle
  mutate(k_l_m =exp(log_k)) |>
  mutate(across(where(is.numeric), \(x) round(x, 1))) |>
  relocate(n, .after = last_col()) |>
  select(-log_k, -k_l_m_ind)

## write table with average k values to csv ----
write_csv(avg_k, "Results/fish_metal_conversion_factors.csv")

## plot individual and average k values ----
k <- ggplot(ind_k, aes(k_l_m_ind, Species))+
  geom_point(color = 'grey') +
  geom_point(avg_k, mapping = aes(k_l_m, Species), color = "black", size = 3) +
  geom_vline(xintercept = 1) +
  scale_x_log10() +
  facet_wrap(~Parameter, scales = "free") +
  theme_bw() +
  labs(x = "Liver:Muscle ratio (k)", y = "Species") +
  theme(legend.position = "bottom") +
  ggtitle("Liver:Muscle ratio (k) for As, Se, Cu and Zn in perch, herring and eelpout")

ggsave(plot = k,"Results/fish_metal_conversion_factors_plot.png", width = 8, height = 6, dpi = 300)
