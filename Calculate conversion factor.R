# Calculate tissue contaminant conversion factors (k) ----
#
# Example on the calculation of tissue conversion factors used for the Integrative 
# Chemical and Biological Effects Monitoring project based on data from
# Faxneld et al. 2015 and Danielsson et al. 2018
#
# Pipeline: read perch and herring liver/muscle data -> harmonize columns and
# filter to a shared set of parameters -> bind into one table -> calculate the
# liver:muscle conversion factor (k) per parameter and species, following the
# geometric-mean method of Soerensen et al. 2023.


library(tidyverse)
library(readxl)
library(stringr)

## read data ----
###  perch data ----
perch <- read_excel("data/perch_liver_muscle_long.xlsx", skip = 2) |>
  rename(Sample = 'specimen_no',
         Parameter = 'Metal',
         Liver = 'L',
         Muscle = 'M') |>
  mutate(Parameter = str_to_sentence(Parameter),  
         Species = 'Perch',
         Sample = as.character(Sample)) |>
  filter(Parameter == 'As' | Parameter == 'Se' | Parameter == 'Cu' | Parameter == 'Zn') |>
  # only include samples with both liver and muscle data above LOQ
  filter(Muscle > 0,
         Liver > 0)

### herring data ----
herring <- read_excel("data/herring_liver_muscle.xlsx", skip = 2) |>
  select(Parameter, Sample, Muscle, Liver) |>
  # only include samples with both liver and muscle data above LOQ
  drop_na() |>
  filter(Parameter == 'As' | Parameter == 'Se' | Parameter == 'Cu' | Parameter == 'Zn') |>
  mutate(Site = 'Lilla Vartan',
         Species = 'Herring')

## bind data and calculate k----
fish <- bind_rows(perch, herring) |>
  select(-c(Sample, Site)) |>
  mutate(log_k = log(Liver/Muscle)) |> 
  # alternative way of writing equation
  #mutate(log_k = log(Liver) - log(Muscle)) |>
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
  relocate(n, .after = last_col())
