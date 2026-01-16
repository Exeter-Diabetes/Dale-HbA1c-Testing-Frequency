.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

library(RMariaDB)
library(tidyverse)
library(glmmTMB)
library(parallel)
library(splines)

#Read the data
con2 <- dbConnect(MariaDB())

#Load in the data
dat <- dbReadTable(con2, "dh_hu_hba1c_visits_dataset_cleaned_v13")


#Remove all values where the start is >=10 - everyone is truncated at first 10 years.
dat_short <- dat %>% filter(start < 10)


#Calculate the offset
dat_short <- dat_short %>%
  mutate(
    frac_year = stop - floor(stop),
    frac_year = if_else(frac_year == 0, 1, frac_year),   # full year if integer
    offset_log = log(frac_year)                          # offset in years
  )

dat_short <- dat_short %>% 
  mutate(ethnicity_bin = as.factor(ethnicity_bin), gender = as.factor(gender), 
         patid = as.factor(patid), stop = as.numeric(stop), imd_decile = as.numeric(imd_decile),
         dm_diag_age = round(as.numeric(dm_diag_age), 1), 
         age_c = dm_diag_age - mean(dm_diag_age, na.rm = TRUE), count = as.numeric(count), 
         incident_comorbidity_score_timevar = as.numeric(incident_comorbidity_score_timevar), 
         preexisting_comorb_score = as.numeric(preexisting_comorb_score), 
         incident_depression_timevar  = factor(incident_depression_timevar, levels = c(0,1)), 
         pret2d_dep = factor(pret2d_dep, levels = c(0,1)))



n_cores = 15
null_model <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                        gender + age_c + ethnicity_bin + imd_decile + offset(offset_log),
                      data = postdep_test,   family = "compois", 
                      control = glmmTMBControl(parallel = n_cores), 
                      verbose = T)

saveRDS(null_model, "compois_postdep_null.RDS")


#Re-run model without comorbidity scores
predep_model <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                          gender + age_c + ethnicity_bin + imd_decile + offset(offset_log) + incident_depression_timevar, 
                        data = postdep_test,   family = "compois", 
                        control = glmmTMBControl(parallel = n_cores), 
                        verbose = T)

saveRDS(predep_model, "compois_postdep.RDS")

#And with
predep_model_comorb <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                   gender + age_c + ethnicity_bin + imd_decile + offset(offset_log) + preexisting_comorb_score + incident_comorbidity_score_timevar + incident_depression_timevar, 
                 data = postdep_test,   family = "compois", 
                 control = glmmTMBControl(parallel = n_cores), 
                 verbose = T)

saveRDS(predep_model_comorb, "compois_postdep_comorb.RDS")

#predep_models
n_cores = 15

#Without comorb
null_model <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                        gender + age_c + ethnicity_bin + imd_decile + offset(offset_log),
                      data = predep_df,   family = "compois", 
                      control = glmmTMBControl(parallel = n_cores), 
                      verbose = T)

saveRDS(null_model, "compois_predep_null.RDS")



predep_model <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                   gender + age_c + ethnicity_bin + imd_decile + offset(offset_log) + pret2d_dep, 
                 data = predep_df,   family = "compois", 
                 control = glmmTMBControl(parallel = n_cores), 
                 verbose = T)
saveRDS(predep_model, "compois_predep.RDS")

#With comorb
predep_comorb_model <- glmmTMB(count ~ ns(stop, 3) + (1 | patid) + 
                                 gender + age_c + ethnicity_bin + imd_decile + offset(offset_log) + preexisting_comorb_score + incident_comorbidity_score_timevar + pret2d_dep, 
                               data = predep_df,   family = "compois", 
                               control = glmmTMBControl(parallel = n_cores), 
                               verbose = T)

saveRDS(predep_comorb_model, "compois_predep_comorb.RDS")
