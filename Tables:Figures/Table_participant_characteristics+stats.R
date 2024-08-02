library(tidyverse)
library(table1)
library(htmltools)

# load the original covaraite data 

pheno_file <- "/Volumes/Sofer Lab/HCHS_SOL/Datasets/ms968_covariates_20200220.csv"
dat <- read.csv(pheno_file)

bmi_file <- "/Volumes/Sofer Lab/HCHS_SOL/Datasets/ms691_covariates_20191125.csv"
dat1 <- read.csv(bmi_file)
dat1 <- dat1 |> dplyr::select(ID, BMI)

dat <- dat |> inner_join(dat1,  by = ("ID" = "ID"))

# We also want the the mean and SD from global ancestry? TBD
# The covariates we want in addition to the ones fitted in the admixmap: BMI, hypertension, Type II disbetes, 
# maybe alcohol and cigarette use? 


# Load the covariates for batch 1(x1 from the association R script used to fit the model)
load("./admix_map_all/covariates_b1_model.Rdata")

# Load the covariates for batch 2(x2 from the association R script used to fit the model)
load("./admix_map_all/covariates_b2_model.Rdata")

# Participants from the two batches are not overlapping == 0
# sum(covariate_b1$SUBJECT_ID %in% covariate_b2$SUBJECT_ID)

# Match the indices and join other covariates that need to be reported
# Add a new column indicating which batch it is
# Batch1
col_report <- dat |> dplyr::select(ID, HYPERTENSION, DIABETES2_INDICATOR, DIABETES2, BMI) 
covariate_b1 <- covariate_b1 |> 
  inner_join(col_report, by = ("ID" = "ID")) |>
  mutate(batch = rep("Discovery Batch", nrow(covariate_b1)))

# Batch2
covariate_b2 <- covariate_b2 |> 
  inner_join(col_report, by = ("ID" = "ID")) |>
  mutate(batch = rep("Replication Batch", nrow(covariate_b2)))

# rbind the two dfs to obtain all covariates
covariate_total <- rbind(covariate_b1, covariate_b2) |>
  dplyr::select(AGE, BMI, GFRSCYS, GENDER, CENTER, gengrp6, HYPERTENSION, DIABETES2_INDICATOR, batch) |>
  mutate(HYPERTENSION = if_else(HYPERTENSION == 0, "No", "Yes"),
         DIABETES2_INDICATOR = if_else(DIABETES2_INDICATOR == 0, "No", "Yes"),
         GENDER = recode(GENDER, F = "Female", "M" = "Male"),
         GENDER = factor(GENDER),
         CENTER = recode(CENTER, "B" = "Bronx", "C" = "Chicago",
                         "M" = "Miami",  "S" = "San Diego"),
         CENTER = factor(CENTER),
         gengrp6 = factor(gengrp6),
         HYPERTENSION = factor(HYPERTENSION))
  # 5505 individuals in total

colnames(covariate_total) <-  c("Age (years)", "BMI", "eGFR", "Gender", "Recruitment Center",  "Genetic Analysis Group",
                                "Hypertension", "Diabetes", "batch")

# p-value comparisions in hypertension and typeII diabetes (both as significant)
hypertension <- table(covariate_total$Hypertension, covariate_total$batch)
chisq.test(hypertension)
diabetes <- table(covariate_total$`Type II Diabetes`, covariate_total$batch)
chisq.test(diabetes)
age <- t.test(`Age (years)`~covariate_total$batch, data = covariate_total)
print(age)

# compare eGFR levels, BMI, gender, center, and genetic analysis group 
# test mean eGFR levels difference
t.test(covariate_total$eGFR~covariate_total$batch, data = covariate_total)

# test mean BMI difference 
t.test(covariate_total$BMI~covariate_total$batch, data = covariate_total)

# test gender proporitons difference
chisq.test(table(covariate_total$Gender,  covariate_total$batch))

# test recruitment center distribution difference
center_p <- chisq.test(table(covariate_total$`Recruitment Center`,  covariate_total$batch))$p.value
bonferroni_p_value <- center_p * 6
min(bonferroni_p_value, 1)

# test genetic analysis group difference
ga_group_p <- chisq.test(table(covariate_total$`Genetic Analysis Group`,  covariate_total$batch))$p.value
bonferroni_p_value <- ga_group_p * 15
bonferroni_p_value <- min(bonferroni_p_value, 1)


### Start to reate the table for participant characteristics ####
# Create participant feature table 
units(covariate_total$BMI)   <- HTML("kg/m<sup>2</sup>")
units(covariate_total$eGFR)   <- HTML("mL/min/1.73m<sup>2</sup>")

# More display setting for what to render
render.mean_sd <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x)) 

  what <- switch(name,
                  `Age (years)` = "Mean (SD)",
                  BMI  = "Mean (SD)",
                  eGFR = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
  
}

# run table1() to generate tables based on mean and sd for continuous variables 
# and the total number of individuals in categorical variables.
table1(~ `Age (years)` + BMI + eGFR + Gender + `Recruitment Center`+ `Genetic Analysis Group`
      + Hypertension + `Diabetes`| batch, data = covariate_total,
      render = render.mean_sd) 

