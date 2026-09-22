# Calculate tissue conversion factors (k) ----

# Tissue conversion factors for Biological Effects project based on data from
# Faxneld et al. 2015 and Danielsson et al. 2018

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
  filter(Muscle > 0,
         Liver > 0)

### herring data ----
herring <- read_excel("data/herring_liver_muscle.xlsx", skip = 2) |>
  select(Parameter, Sample, Muscle, Liver) |>
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
  mutate(k_l_m =exp(log_k)) |>
  mutate(across(where(is.numeric), \(x) round(x, 1))) |>
  relocate(n, .after = last_col())
