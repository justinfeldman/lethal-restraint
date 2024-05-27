# Description: Create a new dataframe from the AP Lethal Restraint dataset
# standardizing variables, merging county/death investigation agency-level
# variables, and removing some based on exclusion criteria

library(tidyverse)
library(lubridate)

file_path <- file.path("..", "data", "01_case_level_data_lethal_restraint.csv")

# Load the original AP data file
lethal.restraint <- read.csv(file_path)

# Copy or recode variables for the analysis
lethal.restraint <- lethal.restraint %>%
  mutate(key = paste(i_state, d_year_death, d_last_name, d_first_name,
                     sep = "-"),
         year = d_year_death,
         month = month(mdy(d_dod)),
         age = d_age,
         gender = d_gender,
         race.ethnicity = case_when(
           d_ethnicity == "Hispanic" ~ "Hispanic",
           d_race %in% c("Other", "American Indian or Alaska Native", "Asian",
                         "Native Hawaiian or Pacific Islander") ~
             "Other persons of color",
           TRUE ~ d_race),
         manner = me_manner_death,
         outcome_homicide = if_else(me_manner_death == "homicide", 1, 0))

# Creating indicator variables from a collapsed set of force types

# Extract unique force types reported by AP
lethal.restraint <- lethal.restraint %>%
  mutate(list_uof  = gsub("[^A-Za-z0-9;]", "", i_uof_type))
         
unique_uof_types <- unique(unlist(strsplit(lethal.restraint$list_uof, ";")))

# Create indicator variables for force types as reported by AP
for (forcetype in unique_uof_types) {
  lethal.restraint[paste0("original_uof_", forcetype)] <- 
    as.integer(grepl(forcetype, lethal.restraint$list_uof))
}

# New indicator variables collapsing some of the original force categories
lethal.restraint <- lethal.restraint %>%
  mutate(uof_baton =
          if_else(original_uof_Struckwithbatonorotherobject == 1, 1, 0),
         uof_cew_darts =
          if_else(original_uof_TaserorConductedEnergyWeaponCEWdarts
            == 1, 1, 0),
         uof_cew_other = if_else(
           original_uof_TaserorConductedEnergyWeaponCEWdrivestun == 1 |
           original_uof_TaserorConductedEnergyWeaponCEWmodeunknown == 1, 1, 0),
         uof_chem_spray = if_else(original_uof_Pepperorchemicalspray == 1, 1, 0),
         uof_handcuffs = if_else(original_uof_Handcuffs == 1, 1, 0),
         uof_hobble_hogtie = if_else(original_uof_Hogtie == 1 |
                                      original_uof_Hobblerestraint == 1, 1, 0),
         uof_leg_restraint = if_else(original_uof_Legrestraint == 1, 1, 0),
         uof_neck = if_else(original_uof_Chokeholdorcarotidrestraint == 1, 1, 0),
         uof_other = if_else(original_uof_Other == 1 |
                                 original_uof_WRAPdevice == 1 |
                                 original_uof_Restraintchair == 1, 1, 0),
         uof_other_highly_hazardous = if_else(original_uof_Kneetoneck == 1 |
                                       original_uof_Beanbagrounds == 1 |
                                       original_uof_Caninebite == 1, 1, 0),
         uof_prone = if_else(
                      original_uof_Pronepositionnopressureappliedorunclear == 1 |
                      original_uof_Pronepositionwithofficerbodyweightorpressure == 1,
                      1, 0),
         uof_sedative = if_else(
                          original_uof_Chemicalrestraintsedative == 1, 1, 0),
         uof_spit_hood = if_else(original_uof_Spithood == 1, 1, 0),
         uof_struck_with_body = if_else(original_uof_Punched == 1 |
                                            original_uof_Kneestrike == 1 |
                                            original_uof_Elbowstrike == 1 |
                                            original_uof_Kicked == 1, 1, 0),
         uof_take_down_tackle = if_else(original_uof_Takedown == 1 |
                                            original_uof_Tackled == 1, 1, 0)
  )


# Delete original UoF indicator variables
lethal.restraint <- lethal.restraint %>%
  select(-starts_with("original_uof_"))

# Merge in covariates not in the main dataset
file_path <- file.path("..", "data", "additional covariates.csv")
additional_covariates <- read_csv(file_path)
lethal.restraint <- left_join(lethal.restraint,
                              additional_covariates, by = "key")

# Create quartiles for County % Republican
lethal.restraint <- lethal.restraint %>%
  mutate(county_pct_repub_q = ntile(county_pct_repub, 4))
  

# Generate autopsy agency previous determination variables
# Remove rows with no date, replace determinations with zero if missing
prev <- lethal.restraint %>% 
  filter(exclude == 0) %>%
  mutate(d_dod = mdy(d_dod)) %>%
  filter(!is.na(d_dod), autopsy_agency_stdized != '') %>%
  mutate_at(vars(outcome_homicide, outcome_cause_inj, outcome_cause_either), 
            ~coalesce(na_if(., NA), 0))

# Sort by autopsy agency, date
# Assign values for a set of previous determination variables based on whether
# 1 or more prior deaths were classified as homicides, mentioned injuries
# or mentioned force/injuries
prev <- prev %>%
  arrange(autopsy_agency_stdized, d_dod) %>%
  group_by(autopsy_agency_stdized) %>%
  mutate(prev_homicide = if_else(row_number() == 1, NA_integer_,
                                 lag(cumsum(outcome_homicide >= 1),
                                     default = 0)),
         prev_cause_inj = if_else(row_number() == 1, NA_integer_,
                                  lag(cumsum(outcome_cause_inj >= 1),
                                      default = 0)),
         prev_cause_either = if_else(row_number() == 1, NA_integer_,
                                     lag(cumsum(outcome_cause_either >= 1),
                                         default = 0)),
         first_death = if_else(row_number() == 1, T, F)
  ) %>%
  ungroup()

# Set prev_ variables to 1, rather than the cumulative sum
prev <- prev %>%
  mutate_at(vars(prev_homicide, prev_cause_inj, prev_cause_either), 
            ~if_else(. > 1, 1, .))

# Drop all variables already in lethal.restraint then merge to lethal.restraint
prev <- prev %>%
  select(key, prev_homicide, prev_cause_inj, prev_cause_either, first_death)

lethal.restraint <- left_join(lethal.restraint,
                              prev, by = "key")

# Save file
write_csv(lethal.restraint, "../data/lethal restraint analytic dataset.csv")
