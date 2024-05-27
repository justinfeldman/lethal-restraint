# Description: This code estimates IPTWs using WeightIt,
# then fits weighted logistic regression models using glm(),
# then uses MarginalEffects to estimate prevalence differences

# Note that the margin effect estimation is set up for the ATT estimand,
# for categorical variables, we specify a subset based on level of the variable
# then should only use the estimate for the subsetted variable value
# vs. reference group for a valid estimate,
# even though other contrasts appear too


library(WeightIt)
library(marginaleffects)
library(tidyverse)
library(rms)

set.seed(123)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data <- read.csv(file.path("..", "data",
                           "lethal restraint analytic dataset.csv"))

# Drop any observations with missing data for variables with any missingness
# (WeightIt doesn't do this automatically)
# N = 887 complete cases

data <- data %>%
  filter(!is.na(race.ethnicity), autopsy_agency_stdized != '', !is.na(month),
         jurisdiction_type != '', !is.na(outcome_cause_inj),
         !is.na(outcome_cause_either), !is.na(outcome_homicide), !is.na(age),
         exclude == 0)

data <- data %>%
  mutate(race.ethnicity = as.factor(race.ethnicity),
  county_pct_repub_q = relevel(as.factor(county_pct_repub_q), ref = 4),
  jurisdiction_type = relevel(as.factor(jurisdiction_type), ref = "me"),
  d_state_abbr = as.factor(d_state_abbr),
  county_pov = as.vector(scale(county_pov)),
  county_pct_repub = as.vector(scale(county_pct_repub)))

data$race.ethnicity <- relevel(data$race.ethnicity, ref = "White")

# Create a set of covariates used in all IPTW models
common_covariates <- "age + month + year + county_pov + county_rucc + d_state_abbr + "

# Add all UOF variables to the list of covariates based on prefix
common_covariates <- paste0(common_covariates,
                            paste(grep("^uof_", names(data),
                            value = TRUE), collapse = " + ")) 

# Estimate IPTWs for each exposure
# Seed needs to be set individually due to a dbatrts multithreading issue


# Exposure 1: CME Type
cme.weight <- weightit(as.formula(paste0("jurisdiction_type ~ race.ethnicity + 
                                          county_pct_repub + ",
                                         common_covariates)),
                       data = data, method = "bart", estimand = "att",
                       focal = "me", seed = 123)

# Exposure 2: race/ethnicity
race.weight <- weightit(as.formula(paste0("race.ethnicity ~ jurisdiction_type + 
                                          county_pct_repub + ",
                                         common_covariates)),
                       data = data, method = "bart", estimand = "att",
                       focal = "White")

# Exposure 3: % Republican (reference group = high Republican)
# Note: % Republican treated as quartiles here, continuous elsewhere

repub.weight <- weightit(as.formula(paste0("county_pct_repub_q ~
                                            jurisdiction_type +
                                           race.ethnicity + ", common_covariates)),
                        data = data, method = "bart", estimand = "att",
                        focal = 4)

# Exposure 4: Previous outcome values (3 sets of weights)

# Remove those without previous values (e.g. because they're first in agency)
data.prev <- data %>%
  filter(!is.na(prev_cause_either), !is.na(prev_cause_inj),
         !is.na(prev_homicide))

prev_homicide.weight <- weightit(as.formula(paste0("prev_homicide ~
                              county_pct_repub +
                              jurisdiction_type + race.ethnicity + ",
                           common_covariates)),
                           data = data.prev, method = "bart",
                           estimand = "att")

prev_inj.weight <- weightit(as.formula(paste0("prev_cause_inj ~
                              county_pct_repub +
                              jurisdiction_type + race.ethnicity + ",
                                              common_covariates)),
                            data = data.prev, method = "bart",
                            estimand = "att")

prev_either.weight <- weightit(as.formula(paste0("prev_cause_either ~
                              county_pct_repub +
                              jurisdiction_type + race.ethnicity + ",
                                                 common_covariates)),
                               data = data.prev, method = "bart")

# Add weights to dataframe (except previous outcome weights)
data$cme.weights <- cme.weight$weights
data$race.weights <- race.weight$weights
data$repub.weights <- repub.weight$weights

# Add previous outcome weights to subset
data.prev$homicide.weights <- prev_homicide.weight$weights
data.prev$inj.weights <- prev_inj.weight$weights
data.prev$either.weights <- prev_either.weight$weights

# Outcome models
# All continuous variables modeled as 4-knot restricted cubic splines
# except rural-urban continuum codes (treated as linear bc too few values)

# % Republican treated differently because it's continuous in models that
# use it as a covariate
# But quartiles in the models where it is the exposure of interest

# Covariates used in all outcome models
common_covariates_outcome <- paste0("race.ethnicity + jurisdiction_type + ",
                                    "rcs(age, 4) + rcs(month, 4) + ",
                                    "rcs(year, 4) + county_rucc +  d_state_abbr + ",
                                    paste(grep("^uof_", names(data),
                                                        value = TRUE),
                                          collapse = " + "))


# Models 1a, 1b, 1c: CME Type
model.cme.homicide <- glm(as.formula(paste0("outcome_homicide ~ ",
                                        common_covariates_outcome,
                                        " + rcs(county_pct_repub, 4)")),
                           data = data, family = quasibinomial,
                           weights = data$cme.weights)


model.cme.inj <- glm(as.formula(paste0("outcome_cause_inj ~ ",
                                            common_covariates_outcome,
                                            " + rcs(county_pct_repub, 4)")),
                          data = data, family = quasibinomial,
                          weights = data$cme.weights)     

model.cme.either <- glm(as.formula(paste0("outcome_cause_either ~ ",
                                       common_covariates_outcome,
                                       " + rcs(county_pct_repub, 4)")),
                          data = data, family = quasibinomial,
                          weights = data$cme.weights)  


# Model 1: ATT estimates for risk difference, SEs clustered for autopsy agency

avg_comparisons(model.cme.homicide,
                variables = "jurisdiction_type",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, jurisdiction_type == "sheriff-coroner"))

avg_comparisons(model.cme.homicide,
                variables = "jurisdiction_type",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, jurisdiction_type == "coroner"))

avg_comparisons(model.cme.inj,
               variables = "jurisdiction_type",
               vcov = ~autopsy_agency_stdized,
               newdata = subset(data, jurisdiction_type == "sheriff-coroner"))

avg_comparisons(model.cme.inj,
                variables = "jurisdiction_type",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, jurisdiction_type == "coroner"))

avg_comparisons(model.cme.either,
                variables = "jurisdiction_type",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, jurisdiction_type == "sheriff-coroner"))

avg_comparisons(model.cme.either,
                variables = "jurisdiction_type",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, jurisdiction_type == "coroner"))


# Model 2. Race
model.race.homicide <- glm(as.formula(paste0("outcome_homicide ~ ",
                                             common_covariates_outcome,
                                             " + rcs(county_pct_repub, 4)")),
                           data = data, family = quasibinomial,
                           weights = data$race.weights)

model.race.inj <- glm(as.formula(paste0("outcome_cause_inj ~ ",
                                        common_covariates_outcome,
                                        " + rcs(county_pct_repub, 4)")),
                      data = data, family = quasibinomial,
                      weights = data$race.weights)

model.race.either <- glm(as.formula(paste0("outcome_cause_either ~ ",
                                           common_covariates_outcome,
                                           " + rcs(county_pct_repub, 4)")),
                         data = data, family = quasibinomial,
                         weights = data$race.weights)

# ATT estimates for risk difference, SEs clustered for autopsy agency

avg_comparisons(model.race.homicide,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Black"))

avg_comparisons(model.race.homicide,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Hispanic"))

avg_comparisons(model.race.homicide,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Other persons of color"))

avg_comparisons(model.race.inj,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Black"))

avg_comparisons(model.race.inj,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Hispanic"))

avg_comparisons(model.race.inj,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Other persons of color"))

avg_comparisons(model.race.either,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Black"))

avg_comparisons(model.race.either,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Hispanic"))

avg_comparisons(model.race.either,
                variables = "race.ethnicity",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, race.ethnicity == "Other persons of color"))


# Model 3. Republican %
model.repub.homicide <- glm(as.formula(paste0("outcome_homicide ~ ",
                                             common_covariates_outcome,
                                             " + county_pct_repub_q")),
                           data = data, family = quasibinomial,
                           weights = data$repub.weights)

model.repub.inj <- glm(as.formula(paste0("outcome_cause_inj ~ ",
                                              common_covariates_outcome,
                                              " + county_pct_repub_q")),
                            data = data, family = quasibinomial,
                            weights = data$repub.weights)

model.repub.either <- glm(as.formula(paste0("outcome_cause_either ~ ",
                                         common_covariates_outcome,
                                         " + county_pct_repub_q")),
                       data = data, family = quasibinomial,
                       weights = data$repub.weights)


avg_comparisons(model.repub.homicide,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 1))

avg_comparisons(model.repub.homicide,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 2))

avg_comparisons(model.repub.homicide,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 3))


avg_comparisons(model.repub.inj,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 1))

avg_comparisons(model.repub.inj,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 2))

avg_comparisons(model.repub.inj,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 3))

avg_comparisons(model.repub.either,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 1))

avg_comparisons(model.repub.either,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 2))

avg_comparisons(model.repub.either,
                variables = "county_pct_repub_q",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, county_pct_repub_q == 3))

# Model 4. Previous manner/cause determinations
# Note: both the outcome and exposure differ in each model for these

model.prev.homicide <- glm(as.formula(paste0("outcome_homicide ~ ",
                                             common_covariates_outcome,
                                " + rcs(county_pct_repub, 4) + prev_homicide")),
                           data = data.prev, family = quasibinomial,
                           weights = data.prev$homicide.weights)

model.prev.inj <- glm(as.formula(paste0("outcome_cause_inj ~ ",
                                             common_covariates_outcome,
                                " + rcs(county_pct_repub, 4) + prev_cause_inj")),
                           data = data.prev, family = quasibinomial,
                           weights = data.prev$inj.weights)

model.prev.either <- glm(as.formula(paste0("outcome_cause_either ~ ",
                                             common_covariates_outcome,
                                " + rcs(county_pct_repub, 4) + prev_cause_either")),
                           data = data.prev, family = quasibinomial,
                           weights = data.prev$either.weights)

avg_comparisons(model.prev.homicide,
                variables = "prev_homicide",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, prev_homicide == 1))

avg_comparisons(model.prev.inj,
                variables = "prev_cause_inj",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, prev_cause_inj == 1))

avg_comparisons(model.prev.either,
                variables = "prev_cause_either",
                vcov = ~autopsy_agency_stdized,
                newdata = subset(data, prev_cause_either == 1))
