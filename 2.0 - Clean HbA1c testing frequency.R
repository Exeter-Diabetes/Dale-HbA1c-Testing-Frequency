#Create the final (for now) data frame
library(aurum)
library(dplyr)
library(data.table)
library(survival)
library(RMariaDB)
library(glmmTMB)
library(devtools)
library(DBI)
library(tidyverse)
rm(list = ls())
#install_github("Exeter-Diabetes/CPRD-analysis-package")

convert_int_to_numeric <- function(data) {
  data %>%
    mutate(across(where(is.numeric), as.numeric))
}

#####Script setup#####

#Steps:
# - Load in the Healthcare data, death data, pre- and post- depression data, age data,
#   and ethnicity. 
# - Filter the file for individuals with an unidentifiable diagnosis date. 
# - Add the pre- and post- depression fields with timing. 
# - Calculate the time to event for the all-cause mortality

#Script for all Qc steps that need to be included:
cprdenvname <- "CPRD_diabetes_data"
yaml <- ""

#Open the CPRD data
cprd <- CPRDData$new(cprdEnv = cprdenvname, cprdConf = yaml)

#glycaemic monitoring
analysis <- cprd$analysis("dh_hu")
health <- health %>% analysis$cached("hba1c_visits_clean") 

#Death
analysis <- cprd$analysis("all_patid")
death <- death %>% analysis$cached("death_end_dat")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")
deptype <- deptype %>% analysis$cached("first_depbroad_type")


#Pre-dep, post-dep, and type of depression
analysis <- cprd$analysis("at_diag")
predep <- predep %>% analysis$cached("dep_broad")
postdep <- postdep %>% analysis$cached("postt2d_dep_broad")

#Exclude people with a history code
deptype <- deptype %>% filter(preexisting == 1)

#Join the files together - this could get a bit complicated
analysis <- cprd$analysis("dh_hu")
dep_vars <- dep_vars %>% analysis$cached("pret2d_dep_severity_vars")
cohort <- cohort %>% analysis$cached("cohort_pracex_35_ins1yrex")




####Bind files together and clean####
full_table <- health %>% 
  left_join(cohort %>% select(patid, gender, dm_diag_age, dm_diag_date, cprd_ddate.x, gp_death_end_date)) %>%
  left_join(ethnicity %>% select(patid, ethnicity_5cat), by = "patid") %>% 
  left_join(predep %>% select(patid, pret2d_dep, pret2d_dep_diagdate, dep_t2d_diff), by = "patid") %>%
  left_join(postdep %>% select(patid, postt2d_dep, postt2d_dep_diagdate, dep_postt2d_diff), by = "patid") %>%
  anti_join(deptype, by = "patid") %>% left_join(dep_vars, by = "patid")

#full_table %>% count() #4030035
#full_table %>% summarise(n = n_distinct(patid)) #587990 individuals with full data - this is after removing individuals 
#With a pre-existing depression code that was historical - 56,095 individuals were lost.
#Save table
full_table %>% analysis$cached("hba1c_visits_clean_interim")

full_table <- full_table %>% analysis$cached("hba1c_visits_clean_interim")

#Define a time to death variable.
#if the value is NA, the new value is 0. 
#If the time is over 10 years, set to 0. 
#If the time is under 10 years, set to 1. 
full_table_mortality <- full_table %>%
  rename(cprd_ddate = cprd_ddate.x) %>% 
  mutate(mortality = ifelse(is.na(cprd_ddate), 0, 1)) %>% #If death date is NA, they didn't die.
  mutate(time_to_mortality = 
           ifelse(mortality == 0, as.numeric(datediff(gp_death_end_date, dm_diag_date))/365.25, 
                  as.numeric(datediff(cprd_ddate, dm_diag_date))/365.25)) %>% 
  mutate(time_to_mortality = round(time_to_mortality, 2)) %>%
  mutate(mortality_10yrs = ifelse(mortality == 1 & time_to_mortality > 10, 0, mortality), 
         time_to_mortality_10yrs = ifelse(time_to_mortality > 10, 10, time_to_mortality)) 
  

#Cache the table
full_table_mortality %>% analysis$cached("hba1c_visits_clean_interim_mort")


#Test whether there are any individuals with pre- and post-depression times
#Worked correctly - there is no-one with both pre-existing and incident depression. 
#dep_check <- full_table_mortality_local %>% filter(!is.na(dep_postt2d_diff) & !is.na(dep_t2d_diff)) %>% 
#  select(dep_postt2d_diff, dep_t2d_diff)
#Make sure no one dies after they're censored - this can't happen, as death date 
#Is baked into the censoring time..

full_table_mortality <- full_table_mortality %>% analysis$cached("hba1c_visits_clean_interim_mort")

#define a column - NA = pre-T2D dep, 0 = never depressed, 1 = post-t2d dep

#If pret2d_dep = 1, = NA.
#If pret2d_dep = 0 and postt2d_dep = 0, 0
#if pret2d_dep = 0 and postt2d_dep = 1, 1

full_table_mortality_local <- full_table_mortality %>% collect()


#Calculate incident depression.
full_table_mortality_local <- full_table_mortality_local %>%
  mutate(
    pret2d_dep = as.numeric(pret2d_dep),
    postt2d_dep = as.numeric(postt2d_dep)
  )

full_table_dep <- full_table_mortality_local %>% 
  mutate(incident_depression = case_when(
    pret2d_dep == 1 ~ NA,  # Use NA_real_ to ensure proper NA handling in numeric contexts
    pret2d_dep == 0 & postt2d_dep == 0 ~ 0, 
    pret2d_dep == 0 & postt2d_dep == 1 ~ 1))

#For individuals with pre-existing depression, calculate the number of codes before depression, 
#The number of total episodes prior to T2D, and the number of depression dates - can use Alex's code for this. 

#Calcu
full_table_dep_time <- full_table_dep %>%
   rename(time_to_incident_depression = dep_postt2d_diff)%>%
    mutate(incident_depression_year = ceiling(as.numeric(time_to_incident_depression/365.25)))


#Count depression occurring in the same year as the point at which people are affected.
#This actually accounts for people who develop depression later - good. 
full_table_dep_timevar <- full_table_dep_time %>% 
  mutate(incident_dep_timevar = case_when(is.na(incident_depression) ~ NA, 
                                          incident_depression == 0 ~ 0, 
                                          incident_depression == 1 ~ ifelse(as.numeric(incident_depression_year) <= as.numeric(hba1c_date_yrs), 1, 0)))


#######CHECK THIS CODE#######


####Remove Depression History####

#Remove individuals where the history of depression code.
#Pull the table locally
dep_local <- deptype %>% collect()
table(dep_local$dephist_cat_first)

#Anti-join the file.
full_table_final <- full_table_dep_timevar  %>% 
  anti_join(dep_local %>% filter(dephist_cat_first == 1) %>% 
              select(patid), by = "patid") #Just keep the column from here instead.



####Add ethnicity####
full_analysis <- full_table_final  %>%
  mutate(patid = as.factor(patid), hba1c_date_yrs = as.numeric(hba1c_date_yrs), 
         count = as.numeric(count), ethnicity_bin = case_when(
           ethnicity_5cat == 0 ~ 0, 
           ethnicity_5cat %in% c(1:4) ~ 1
         )) %>% filter(!is.na(ethnicity_bin)) %>% 
  mutate(ethnicity_bin = as.factor(ethnicity_bin)) 

full_analysis %>% distinct(patid) %>% count() #571852 

####Add TDI####

#Add Townsend Deprivation Index:
analysis = cprd$analysis("all_patid")
townsend <- townsend %>% analysis$cached("townsend_score")
townsend <- townsend %>% collect()

#Bind the table to the current df
full_analysis_te <- full_analysis %>% inner_join(townsend %>% mutate(patid = as.factor(patid)))

full_analysis_te <- convert_int_to_numeric(full_analysis_te)
full_analysis_te$incident_dep_timevar <- as.numeric(full_analysis_te$incident_dep_timevar)

full_analysis_te %>% distinct(patid) %>% count() #472205


#Smoking status - not included as not enough data.
#analysis = cprd$analysis("at_diag")
#smoking <- smoking %>% analysis$cached("smoking")
#smoking <- smoking %>% collect()
#sum(is.na(smoking$qrisk2_smoking_cat)) #2m are NA - don't include.
#rm(smoking)

#Load in the baseline BMI variable:
#analysis = cprd$analysis("all_patid")
#bmi <- bmi %>% analysis$cached("clean_bmi_medcodes")
#bmi <- bmi %>% collect()

#Filter the BMI table using the values in the hb table
#bmi_short <- bmi %>% mutate(patid = as.factor(patid)) %>% 
#  inner_join(full_analysis_te %>% select(patid, cprd_ddate, dm_diag_date) %>% 
#distinct(patid, .keep_all = T), by = "patid") #7.7m values

#Identify values within 3 months of diagnosis, or within 7 days after
#bmi_at_diag <- bmi_short %>% mutate(time_since_diagnosis = as.numeric(difftime(date, dm_diag_date, units = "days")))
#bmi_at_diag_bmi <- bmi_at_diag %>% 
#  filter(time_since_diagnosis < 30 & time_since_diagnosis > -30) 
#around 292k have baseline BMI - might be worth including or looking at it in 
#The future?


####Upload raw table####
con2 <- dbConnect(MariaDB())



#Create empty table
dbExecute(con2, "
CREATE TABLE dh_hu_hba1c_visits_interim_tds (
    patid BIGINT UNSIGNED, 
    hba1c_date_yrs BIGINT UNSIGNED, 
    count BIGINT UNSIGNED, 
    maximum_followup_year BIGINT UNSIGNED, 
    gender BIGINT UNSIGNED,
    dm_diag_age FLOAT, 
    dm_diag_date DATE, 
    cprd_ddate DATE, 
    gp_death_end_date DATE, 
    ethnicity_5cat BIGINT UNSIGNED, 
    pret2d_dep BIGINT UNSIGNED, 
    pret2d_dep_diagdate DATE, 
    dep_t2d_diff FLOAT, 
    postt2d_dep BIGINT UNSIGNED, 
    postt2d_dep_diagdate DATE, 
    dep_postt2d_diff FLOAT, 
    first_pret2d_dep_diff FLOAT, 
    n_pret2d_dep_codes BIGINT UNSIGNED, 
    pret2d_dep_recurrent BIGINT UNSIGNED, 
    last_pret2d_dep_diff FLOAT, 
    first_pret2d_dep_10yrs FLOAT, 
    last_pret2d_dep_10yrs FLOAT,
    mortality BIGINT UNSIGNED, 
    time_to_mortality FLOAT, 
    mortality_10yrs BIGINT UNSIGNED,
    time_to_mortality_10yrs FLOAT, 
    incident_depression BIGINT UNSIGNED, 
    incident_depression_year BIGINT UNSIGNED, 
    incident_dep_timevar BIGINT UNSIGNED, 
    ethnicity_bin BIGINT UNSIGNED, 
    imd_decile BIGINT UNSIGNED, 
    tds_2011 FLOAT
);
")

chunksize <- 10000
total_rows <- nrow(full_analysis_te)
nchunks <- ceiling(total_rows / chunksize)

# Start the transaction
dbExecute(con2, "START TRANSACTION;")

for (i in seq_len(nchunks)) {
  
  # Determine row indices for this chunk
  idx <- ((i - 1) * chunksize + 1):min(i * chunksize, total_rows)
  chunk <- full_analysis_te[idx, ]
  
  # Process each column vector-wise
  # The following 'lapply' returns a list where each element is a character vector
  # with each value already appropriately escaped/quoted
  processed <- lapply(chunk, function(col) {
    if (is.character(col) || is.factor(col)) {
      # Convert factors to character; escape single quotes and wrap in quotes
      paste0("'", gsub("'", "''", as.character(col)), "'")
    } else if (inherits(col, "Date")) {
      # Format dates or use NULL
      ifelse(is.na(col), "NULL", paste0("'", format(col, "%Y-%m-%d"), "'"))
    } else if (is.numeric(col)) {
      # For numerics, if NA/NaN/Inf then return NULL (without quotes)
      ifelse(is.na(col) | is.nan(col) | is.infinite(col),
             "NULL",
             format(col, scientific = FALSE))
    } else {
      # Fallback for other data types
      rep("NULL", length(col))
    }
  })
  
  # Now combine the processed columns into row strings.
  # do.call(..., sep=",") is equivalent to pasting together the elements by row.
  row_strings <- do.call(paste, c(processed, sep = ","))
  # Wrap each row in parentheses
  row_strings <- paste0("(", row_strings, ")")
  
  # Build the full multi-row INSERT query
  query <- paste0(
    "INSERT INTO dh_hu_hba1c_visits_interim_tds (",
    'patid, hba1c_date_yrs, count, maximum_followup_year, gender, dm_diag_age, dm_diag_date, ',
    'cprd_ddate, gp_death_end_date, ethnicity_5cat, pret2d_dep, pret2d_dep_diagdate, ',
    'dep_t2d_diff, postt2d_dep, postt2d_dep_diagdate, dep_postt2d_diff, first_pret2d_dep_diff, ',
    'n_pret2d_dep_codes, pret2d_dep_recurrent, last_pret2d_dep_diff, first_pret2d_dep_10yrs, ',
    'last_pret2d_dep_10yrs, mortality, time_to_mortality, mortality_10yrs, time_to_mortality_10yrs,  ',
    'incident_depression, incident_depression_year, incident_dep_timevar, ethnicity_bin, ',
    'imd_decile, tds_2011) VALUES ',
    paste(row_strings, collapse = ",")
  )
  
  # Execute the query for this chunk
  dbExecute(con2, query)
  
  print(paste("Inserted chunk", i, "of", nchunks))
}




######FROM THE CO-MORBIDITY SCORE SCRIPT########
#Calculate the co-morbidity score
analysis <- cprd$analysis("at_diag")
comorb <- comorb %>% analysis$cached("comorbidities_long")

analysis <- cprd$analysis("dh_hu")
hu <- hu %>% analysis$cached("hba1c_visits_interim_tds")
comorb_local <- comorb %>% collect()
hu <- hu %>% collect()


#The original list
#comorbidity_names <- c("frailty", "primary_hhf", "af", "angina", "bronchiectasis", 
#                       "cld", "copd", "dementia", "diabeticnephropathy", 
#                       "heartfailure", "ihd", "myocardialinfarction", 
#                       "neuropathy", "otherneuroconditions", "pad", 
#                       "pulmonaryfibrosis", "pulmonaryhypertension", 
#                       "retinopathy", "revasc", "rheumatoidarthritis", 
#                       "solid_cancer", "haem_cancer", "stroke", "tia")


#Select the more relevant comorbidites - CKD needs to be added to the score.
comorbidity_names <- c("frailty", "cld", "copd", "dementia", "diabeticnephropathy", 
                       "heartfailure", "myocardialinfarction", "ihd", "pad",  
                       "retinopathy", "stroke", "tia", "neuropathy")


#Collapse into grep statement and filter
grep_lines <- paste0(comorbidity_names, collapse = "|")

#Doesn't remove that many lines..
comorb_short <- comorb_local %>% select(patid, index_date, matches(grep_lines))

#get the variable names
index_variable_columns <- paste0("pre_index_date_", comorbidity_names)


#Get pre-indexes and select relevant ones
comorb_preindex <- comorb_short %>% select(patid, index_date, matches(index_variable_columns))
comorb_preindex[, -(1:2)] <- lapply(comorb_preindex[, -(1:2)], as.numeric)

#Load in the ckd data
analysis <- cprd$analysis("all_patid")
ckd_stage <- ckd_stage %>% analysis$cached("ckd_stages_from_algorithm")
ckd_stage_local <- ckd_stage %>% collect()

ckd_stage_local_filtered <- ckd_stage_local %>% filter(patid %in% hu$patid)

#Check out the ckd stage
ckd_chk <- ckd_stage_local %>%
  mutate(
    ckd_presence = ifelse(rowSums(!is.na(select(., stage_1:stage_5))) > 0, 1, 0),
    earliest_stage_date = apply(select(., stage_1:stage_5), 1, function(x) min(x, na.rm = TRUE))
  )

ckd_chk$earliest_stage_date[is.infinite(ckd_chk$earliest_stage_date)] <- NA
ckd_chk_bind <- ckd_chk %>% inner_join(hu %>% select(patid, dm_diag_date) %>% distinct(.keep_all = T), by = "patid")

ckd_stage_long <- ckd_stage %>%
  pivot_longer(
    cols = starts_with("stage"), 
    names_to = "stage", 
    values_to = "date"
  ) %>%
  filter(!is.na(date)) %>%  # Remove NA rows
  arrange(patid, date)  # Sort by patid and da

#Convert to the longer format that Katie used. This will allow this to be used with the current
#Code.
ckd_severity <- ckd_stage_long %>%
  mutate(stage_group = case_when(
    stage %in% c("stage_3a","stage_3b","stage_4", "stage_5") ~ "late_stage",
    TRUE ~ "early_stage"
  ))

ckd_severity_local <- ckd_severity %>% collect()


ckd_first <- ckd_severity_local %>%
  group_by(patid, stage_group) %>%
  summarise(first_entry = min(date), .groups = "drop")

head(ckd_first)

#check how many people are in late stage
ckd_chk <- ckd_first %>% filter(stage_group == "late_stage")  %>%
  count() #612082 - 25% - looks a lot more reasonable


#Next, add the diagnosis dates from the hu table, and calculate if these dates were before
#Diagnosis.
ckd_joined <- ckd_first %>% inner_join(hu %>% select(patid, dm_diag_date) %>% 
                                         distinct(.keep_all = T), by = "patid") %>%
  mutate(pre_t2d = ifelse(as.integer(difftime(first_entry, dm_diag_date, units = "days")) <= 7, 1, 0), 
         post_t2d = ifelse(as.integer(difftime(first_entry, dm_diag_date, units = "days")) > 7, 1, 0))

#Pivot wider to match format of previous table. 
#if CKD late is 1, CKD early must also be 1.
ckd_joined_wide <- ckd_joined %>%
  select(patid, stage_group, pre_t2d) %>%
  pivot_wider(names_from = stage_group, values_from = pre_t2d, values_fill = list(pre_t2d = 0)) %>%
  rename("pre_index_date_ckd_early" = "early_stage", "pre_index_date_ckd_late" = "late_stage") %>%
  mutate(pre_index_date_ckd_early = ifelse(pre_index_date_ckd_late == 1, 1,pre_index_date_ckd_early))

#Rbind with the previous complications table. If there are missing values, fill with 0.
#Remove the early stage score as it's not informative about disease.
comorb_preindex <- comorb_preindex %>% left_join(ckd_joined_wide) %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% select(-pre_index_date_ckd_early)


#Before the co-morbidity score can be calculated, some of these variables are
#Going to occur as the same variable, and thus should be collapsed into single
#Variables:

#This changes only multi-morbidity scores, not overall scores if n = 1.
comorb_prescore <- comorb_preindex %>% 
  mutate(pre_index_date_frailty = if_else(pre_index_date_frailty_mild == 1 | 
                                            pre_index_date_frailty_moderate == 1 | 
                                            pre_index_date_frailty_severe == 1, 1, 0)) %>%
  select(-c(pre_index_date_frailty_mild,pre_index_date_frailty_moderate, 
            pre_index_date_frailty_severe))


#Next, calculate the rowsums:
comorb_prescore <- comorb_prescore %>%
  mutate(preexisting_comorb_score = rowSums(select(., -(1:2))))

comorb_prescore <- comorb_prescore %>%
  mutate(row_number = row_number())

length(comorb_prescore$preexisting_comorb_score[comorb_prescore$preexisting_comorb_score > 0])/
  (nrow(comorb_prescore) + length(comorb_prescore$preexisting_comorb_score[comorb_prescore$preexisting_comorb_score > 0]))
#34.7% have a comorbidity at baseline... seems very reasonable.
#Only 24.7% now have a comorbidity at baseline...
#Seems reasonable i guess?

#this is done. Now the original table needs to be formatted into a long format, 
#And converted to match the count data from the HbA1c table.
#For each person... count the number of values after

#create a new data.frame.
#Limit to my list of comorbid conditions with the pre-index date and the post-index date.
#If the pre-index date is one, the post-index date should be set as NA, otherwise
#Keep the post index date.

#First, add the pre-index date and the post index first date columns to the file.

#Pre-index
comorb_local_ckd <- comorb_short %>% left_join(comorb_prescore %>% select(patid, pre_index_date_ckd_late))

#Create post-index first date.. which is the same as the post-index date.
#If post-T2D is 1, use the date available. 
#If it is 0, mark as NA.

ckd_post_wide <- ckd_joined %>%
  filter(stage_group == "late_stage" & post_t2d == 1) %>%  # Keep only late_stage and post_t2d == 1
  select(patid, first_entry, stage_group) %>%
  pivot_wider(names_from = stage_group, values_from = first_entry, values_fill = list(first_entry = NA)) %>%
  rename("post_index_date_first_ckd_late" = "late_stage")

comorb_local_ckd_full <- comorb_local_ckd %>% left_join(ckd_post_wide)


#Create a post-index column for the CKD_late variable.
comorbidity_names <- c("frailty", "cld", "copd", "dementia", "diabeticnephropathy", 
                       "heartfailure", "myocardialinfarction", "ihd", "pad",  
                       "retinopathy", "stroke", "tia", "neuropathy", "ckd_late")



#Keep the pre-index date columns, and the "post_index_date_first" columns. 
index_variable_columns <- paste0("pre_index_date_", comorbidity_names)
index_variable_columns <- c(paste0("post_index_date_first_", comorbidity_names), index_variable_columns)

#Get the columns
comorb_postindex <- comorb_local_ckd_full %>% select(patid, index_date, matches(index_variable_columns))


#Add the ckd after index to the file.
#for frailty cols - regenerate the new frailty column, and choose the earliest
#Date for the post-index frailty column 
comorb_postindex <- comorb_postindex %>%
  mutate(
    post_index_date_first_frailty = if_else(
      rowSums(!is.na(cbind(post_index_date_first_frailty_mild, 
                           post_index_date_first_frailty_moderate, 
                           post_index_date_first_frailty_severe))) == 0, # Check if all columns are NA
      NA, # If all are NA, return NA
      pmin(post_index_date_first_frailty_mild,
           post_index_date_first_frailty_moderate, 
           post_index_date_first_frailty_severe, na.rm = TRUE) # Otherwise, calculate the minimum
    )
  )

comorb_postindex <- comorb_postindex %>% 
  mutate(pre_index_date_frailty = if_else(pre_index_date_frailty_mild == 1 | 
                                            pre_index_date_frailty_moderate == 1 | 
                                            pre_index_date_frailty_severe == 1, 1, 0)) %>%
  select(-c(pre_index_date_frailty_mild,pre_index_date_frailty_moderate, 
            pre_index_date_frailty_severe, post_index_date_first_frailty_mild, 
            post_index_date_first_frailty_moderate, post_index_date_first_frailty_severe))




#Change co-morbidity names to match new cols
new_comorbidity_names <- c("frailty", "cld", "copd", "dementia", "diabeticnephropathy", 
                       "heartfailure", "myocardialinfarction", "ihd", "pad",  
                       "retinopathy", "stroke", "tia", "neuropathy", "ckd_late")

comorb_incident <- comorb_postindex
comorb_incident <- comorb_incident %>% filter(patid %in% hu$patid) #

#Set all individuals with a pre-existing value to NA.
for (comorbidity in new_comorbidity_names) {
  # Dynamically create column names
  pre_index_col <- paste0("pre_index_date_", comorbidity)
  post_index_col <- paste0("post_index_date_first_", comorbidity)
  new_col <- paste0("incident_", comorbidity)
  
  # Use mutate to create or modify the new column
  comorb_incident <- comorb_incident %>% 
    mutate(!!new_col := if_else(.data[[pre_index_col]] == 1, NA, .data[[post_index_col]]))
}

#All individuals accounted for..
#comorb_incident %>% distinct(patid, .keep_all=T) %>% count()
#hu %>% distinct(patid, .keep_all=T) %>% count()


#Now I've got the incident cols.. this should be fairly easy to calculate.
#Then, calculate the time to each event in year.
#As each individual can only have one date, each event is a new comorbidity.
#Therefore, we can calculate the culumative number of co-morbidities each year.
#From here, we can add the total number of pre-existing co-morbidities each year.

#Remove pre-existing cols, and convert to long format.
comorb_long <- comorb_incident %>%
  select(patid, index_date, matches("incident")) %>%
  pivot_longer(
    cols = starts_with("incident_"),       # Select columns that start with 'incident_'
    names_to = "incident_type",           # New column to store the type of incident
    values_to = "comorbidity_event_dt"           # New column to store the corresponding date
  )


#Find individual with most comorbidities:
most_non_na_group <- comorb_long %>%
  group_by(patid) %>%
  summarize(non_na_count = sum(!is.na(comorbidity_event_dt))) %>%
  arrange(desc(non_na_count)) %>%
  slice(1) # Someone has 12 comorbdities in total.


#Next, calculate the ceiling year: 
comorb_yrs <- comorb_long %>% mutate(comorbidity_event_yr = 
                                       ceiling(as.numeric(difftime(comorbidity_event_dt,index_date, units = "days")/365.25)))

#Check again

#Next, filter out any co-morbidities that didn't occur within 10yrs post-T2D.
comorb_filt %>% summarise(count = n_distinct(patid)) #Only keeps all individuals, even if NA.

most_non_na_group <- comorb_filt %>%
  group_by(patid) %>%
  summarize(non_na_count = sum(!is.na(comorbidity_event_dt))) %>%
  arrange(desc(non_na_count)) %>%
  slice(1) # Now 11.


#Re-calc maximum followup years in the table, and inner join with the Health use table.
hu_calced <- hu %>% mutate(maximum_followup_yrs = ceiling(as.numeric(difftime(gp_death_end_date, dm_diag_date, units = "days")/365.25)))


#Check number of rows in each file
#hu_calced %>% summarise(count = n_distinct(patid))

#Check how many rows would be available if using comorbidity data:
#hu_calced %>% inner_join(comorb_yrs %>% select(patid) %>% distinct(.keep_all = T)) %>%
#  summarise(count = n_distinct(patid)) #No one is lost - EXCELLENT!


comorb_hu <- comorb_yrs %>% inner_join(hu_calced %>% 
                                          select(patid, maximum_followup_yrs) %>% 
                                          distinct(patid, .keep_all=T),
                                        by = "patid")

#Filter all values which are higher than the maximum followup year.
comorb_max <- comorb_hu %>% filter(comorbidity_event_yr <= maximum_followup_yrs)

#Now?
most_non_na_group <- comorb_max %>%
  group_by(patid) %>%
  summarize(non_na_count = sum(!is.na(comorbidity_event_dt))) %>%
  arrange(desc(non_na_count)) %>%
  slice(1) # Now 14 in total. 


#Count all non-NA values in each year.
comorb_counted <- comorb_max %>%
  group_by(patid, comorbidity_event_yr) %>%
  summarize(incident_comorbidity_score = sum(!is.na(comorbidity_event_dt)))


#Next, create a mock data frame from the healthcare use file, 
#And add the comorbidity count data to it.
hu_template <- hu %>% select(patid, hba1c_date_yrs) %>% mutate(template_score = 0) %>%
  mutate(comorbidity_event_yr = as.numeric(hba1c_date_yrs)) %>% select(-hba1c_date_yrs)


#left join the co-morbidity file to the template file, then use the max
#Value of these columns as the post-index co-morbidity score, then
comorb_template <- hu_template %>% left_join(comorb_counted, by = c("patid","comorbidity_event_yr"))

#Find the max value of the columns - use NA.rm to get rid of NA values.
comorb_scores <- comorb_template %>% 
  mutate(incident_comorbidity_score = pmax(template_score, incident_comorbidity_score, na.rm = T)) %>%
  select(-template_score) %>% 
  left_join(comorb_prescore %>% select(patid, preexisting_comorb_score)) %>%
  group_by(patid) %>% 
  mutate(incident_comorbidity_score_timevar = cumsum(incident_comorbidity_score)) %>% 
  ungroup()

#Bind this score to the healthcare use file:
hu_comorb <- hu %>% mutate(hba1c_date_yrs = as.numeric(hba1c_date_yrs)) %>% 
             left_join(comorb_scores %>% 
             select(patid, comorbidity_event_yr, preexisting_comorb_score, incident_comorbidity_score_timevar), 
             by = c("patid" = "patid", "hba1c_date_yrs" = "comorbidity_event_yr")) %>% mutate(year_of_diagnosis = year(dm_diag_date))

names(hu_comorb)


#####Resave the table#####
con2 <- dbConnect(MariaDB())



#Create empty table
dbExecute(con2, "
CREATE TABLE dh_hu_hba1c_visits_dataset_cleaned_v10 (
    patid BIGINT UNSIGNED, 
    hba1c_date_yrs BIGINT UNSIGNED, 
    count BIGINT UNSIGNED, 
    maximum_followup_year BIGINT UNSIGNED, 
    gender BIGINT UNSIGNED,
    dm_diag_age FLOAT, 
    dm_diag_date DATE, 
    cprd_ddate DATE, 
    gp_death_end_date DATE, 
    ethnicity_5cat BIGINT UNSIGNED, 
    pret2d_dep BIGINT UNSIGNED, 
    pret2d_dep_diagdate DATE, 
    dep_t2d_diff FLOAT, 
    postt2d_dep BIGINT UNSIGNED, 
    postt2d_dep_diagdate DATE, 
    dep_postt2d_diff FLOAT, 
    first_pret2d_dep_diff FLOAT, 
    n_pret2d_dep_codes BIGINT UNSIGNED, 
    pret2d_dep_recurrent BIGINT UNSIGNED, 
    last_pret2d_dep_diff FLOAT, 
    first_pret2d_dep_10yrs FLOAT, 
    last_pret2d_dep_10yrs FLOAT,
    mortality BIGINT UNSIGNED, 
    time_to_mortality FLOAT, 
    mortality_10yrs BIGINT UNSIGNED,
    time_to_mortality_10yrs FLOAT, 
    incident_depression BIGINT UNSIGNED, 
    incident_depression_year FLOAT, 
    incident_dep_timevar BIGINT UNSIGNED, 
    ethnicity_bin BIGINT UNSIGNED, 
    imd_decile BIGINT UNSIGNED, 
    tds_2011 FLOAT,
    preexisting_comorb_score BIGINT UNSIGNED,
    incident_comorbidity_score_timevar BIGINT UNSIGNED,
    year_of_diagnosis BIGINT UNSIGNED
);
")

chunksize <- 10000
total_rows <- nrow(hu_comorb)
nchunks <- ceiling(total_rows / chunksize)

# Start the transaction
dbExecute(con2, "START TRANSACTION;")

for (i in seq_len(nchunks)) {
  
  # Determine row indices for this chunk
  idx <- ((i - 1) * chunksize + 1):min(i * chunksize, total_rows)
  chunk <- hu_comorb[idx, ]
  
  # Process each column vector-wise
  # The following 'lapply' returns a list where each element is a character vector
  # with each value already appropriately escaped/quoted
  processed <- lapply(chunk, function(col) {
    if (is.character(col) || is.factor(col)) {
      # Convert factors to character; escape single quotes and wrap in quotes
      paste0("'", gsub("'", "''", as.character(col)), "'")
    } else if (inherits(col, "Date")) {
      # Format dates or use NULL
      ifelse(is.na(col), "NULL", paste0("'", format(col, "%Y-%m-%d"), "'"))
    } else if (is.numeric(col)) {
      # For numerics, if NA/NaN/Inf then return NULL (without quotes)
      ifelse(is.na(col) | is.nan(col) | is.infinite(col),
             "NULL",
             format(col, scientific = FALSE))
    } else {
      # Fallback for other data types
      rep("NULL", length(col))
    }
  })
  
  # Now combine the processed columns into row strings.
  # do.call(..., sep=",") is equivalent to pasting together the elements by row.
  row_strings <- do.call(paste, c(processed, sep = ","))
  # Wrap each row in parentheses
  row_strings <- paste0("(", row_strings, ")")
  
  # Build the full multi-row INSERT query
  query <- paste0(
    "INSERT INTO dh_hu_hba1c_visits_dataset_cleaned_v10 (",
    'patid, hba1c_date_yrs, count, maximum_followup_year, gender, dm_diag_age, dm_diag_date, ',
    'cprd_ddate, gp_death_end_date, ethnicity_5cat, pret2d_dep, pret2d_dep_diagdate, ',
    'dep_t2d_diff, postt2d_dep, postt2d_dep_diagdate, dep_postt2d_diff, first_pret2d_dep_diff, ',
    'n_pret2d_dep_codes, pret2d_dep_recurrent, last_pret2d_dep_diff, first_pret2d_dep_10yrs, ',
    'last_pret2d_dep_10yrs, mortality, time_to_mortality, mortality_10yrs, time_to_mortality_10yrs,  ',
    'incident_depression, incident_depression_year, incident_dep_timevar, ethnicity_bin, ',
    'imd_decile, tds_2011, preexisting_comorb_score, incident_comorbidity_score_timevar,
    year_of_diagnosis) VALUES ',
    paste(row_strings, collapse = ",")
  )
  
  # Execute the query for this chunk
  dbExecute(con2, query)
  
  print(paste("Inserted chunk", i, "of", nchunks))
}

#Rename some of the cols:
analysis <- cprd$analysis("dh_hu")
hu <- hu %>% analysis$cached("hba1c_visits_dataset_cleaned_v10")
hb_abs <- hb_abs %>% analysis$cached("hba1c_index_measurement")
bmi_abs <- bmi_abs %>% analysis$cached("bmi_index_measurement")

#These files need to be regenerated.
hu_biomarkers <- hu %>% 
  left_join(hb_abs %>% select(patid, testvalue) %>% rename(hba1c_index = testvalue), by = "patid") %>%
  left_join(bmi_abs %>% select(patid, testvalue) %>% rename(bmi_index = testvalue), by = "patid")

hu_biomarkers %>% analysis$cached("hba1c_visits_dataset_cleaned_v11")

#Once all of this is done, we need to add values for the last "chunk" 
#That is missing between the year and the true censoring date, 
#For counts, depression and  comorbidity score. 
rm(list = ls())
cprdenvname <- "CPRD_diabetes_data"
yaml <- ""

#Open the CPRD data
cprd <- CPRDData$new(cprdEnv = cprdenvname, cprdConf = yaml)

#glycaemic monitoring
analysis <- cprd$analysis("dh_hu")
hu_biomarkers <- hu_biomarkers %>% analysis$cached("hba1c_visits_dataset_cleaned_v11")
hu_local <- hu_biomarkers %>% collect()


#Get a list of names for the df.
analysis <- cprd$analysis("at_diag")
hb <- hb %>% analysis$cached("full_hba1c_index_date_merge")

#Calculate the stop times
hu_local_start <- hu_local %>% convert_int_to_numeric() %>% 
  mutate(start = as.numeric(hba1c_date_yrs - 1), stop = hba1c_date_yrs, 
         time_to_censor = round(as.numeric(difftime(gp_death_end_date, dm_diag_date, units="days"))/365.25, 2))

#Get start times and the gap times
hu_names_times <- hu_local_start %>% group_by(patid) %>% 
  arrange(patid, hba1c_date_yrs) %>% 
  slice_tail() %>% select(patid, stop, time_to_censor) %>% 
  rename(start = stop, stop = time_to_censor)


#Keep the last row. This will be how the data new rows are generated.
last_row <- hu_local_start %>% group_by(patid) %>% 
  arrange(patid, hba1c_date_yrs) %>% 
  slice_tail()

#Count variable
#Remove people where the start and stop times are the same.
hu_names_times <- hu_names_times %>% as.data.frame %>%
  mutate(across(c(start, stop), as.numeric)) %>%
  filter(start < stop)

#Get hb values between these.
hb_max <- hb %>% filter(datediff >= 730) %>% collect()
hb_last <- hb_max %>% inner_join(hu_names_times) %>% 
  mutate(datediff = round(as.numeric(datediff)/365.25, 2)) %>%
  filter(datediff > start & datediff <= stop) %>% group_by(patid) %>% summarise(count = n())

#Last join the full list of names, and
hb_last <- hu_names_times %>% left_join(hb_last) %>%
  mutate(count = replace_na(count, 0))

#Change the last row
last_row_hb <- last_row %>% select(-count, -start, -stop) %>% inner_join(hb_last, by = "patid")

#Update the comorbidity score.
#remove unnecessary variables
rm("hb_max", "hu_local")

#load in the comorbidity file and convert to long.
analysis <- cprd$analysis("at_diag")
comorb <- comorb %>% analysis$cached("comorbidities_long")
comorb_local <- comorb %>% collect()

#Get the relevant names
comorbidity_names <- c("frailty", "cld", "copd", "dementia", "diabeticnephropathy", 
                       "heartfailure", "myocardialinfarction", "ihd", "pad",  
                       "retinopathy", "stroke", "tia", "neuropathy")



#Keep the pre-index date columns, and the "post_index_date_first" columns. 
index_variable_columns <- paste0("pre_index_date_", comorbidity_names)
index_variable_columns <- c(paste0("post_index_date_first_", comorbidity_names), index_variable_columns)


#Get the columns
comorb_postindex <- comorb_local %>% select(patid, index_date, matches(index_variable_columns))


#Add the ckd after index to the file.
#for frailty cols - regenerate the new frailty column, and choose the earliest
#Date for the post-index frailty column 
comorb_postindex <- comorb_postindex %>%
  mutate(
    post_index_date_first_frailty = if_else(
      rowSums(!is.na(cbind(post_index_date_first_frailty_mild, 
                           post_index_date_first_frailty_moderate, 
                           post_index_date_first_frailty_severe))) == 0, # Check if all columns are NA
      NA, # If all are NA, return NA
      pmin(post_index_date_first_frailty_mild,
           post_index_date_first_frailty_moderate, 
           post_index_date_first_frailty_severe, na.rm = TRUE) # Otherwise, calculate the minimum
    )
  )

comorb_postindex <- comorb_postindex %>% 
  mutate(pre_index_date_frailty = if_else(pre_index_date_frailty_mild == 1 | 
                                            pre_index_date_frailty_moderate == 1 | 
                                            pre_index_date_frailty_severe == 1, 1, 0)) %>%
  select(-c(pre_index_date_frailty_mild,pre_index_date_frailty_moderate, 
            pre_index_date_frailty_severe, post_index_date_first_frailty_mild, 
            post_index_date_first_frailty_moderate, post_index_date_first_frailty_severe))




#Change co-morbidity names to match new cols
new_comorbidity_names <- c("frailty", "cld", "copd", "dementia", "diabeticnephropathy", 
                           "heartfailure", "myocardialinfarction", "ihd", "pad",  
                           "retinopathy", "stroke", "tia", "neuropathy")

comorb_incident <- comorb_postindex
comorb_incident <- comorb_incident %>% filter(patid %in% last_row_hb$patid)

#Set all individuals with a pre-existing value to NA.
for (comorbidity in new_comorbidity_names) {
  # Dynamically create column names
  pre_index_col <- paste0("pre_index_date_", comorbidity)
  post_index_col <- paste0("post_index_date_first_", comorbidity)
  new_col <- paste0("incident_", comorbidity)
  
  # Use mutate to create or modify the new column
  comorb_incident <- comorb_incident %>% 
    mutate(!!new_col := if_else(.data[[pre_index_col]] == 1, NA, .data[[post_index_col]]))
}



#Now I've got the incident cols.. this should be fairly easy to calculate.
#Then, calculate the time to each event in year.
#As each individual can only have one date, each event is a new comorbidity.
#Therefore, we can calculate the culumative number of co-morbidities each year.
#From here, we can add the total number of pre-existing co-morbidities each year.

#Remove pre-existing cols, and convert to long format.
comorb_long <- comorb_incident %>%
  select(patid, index_date, matches("incident")) %>%
  pivot_longer(
    cols = starts_with("incident_"),       # Select columns that start with 'incident_'
    names_to = "incident_type",           # New column to store the type of incident
    values_to = "comorbidity_event_dt"           # New column to store the corresponding date
  ) %>% filter(!is.na(comorbidity_event_dt))

#remove comorbidities that aren't between the dates
comorb_filtered <- comorb_long %>% mutate(datediff = round(as.numeric(difftime(comorbidity_event_dt, index_date, units = "days"))/365.25,2)) %>% 
  inner_join(last_row_hb %>% select(patid, start, stop)) %>% 
  filter(datediff > start & datediff <= stop) # picked up 60k comorbs.. impressive.

#Count comorbidities
comorb_counts <- comorb_filtered %>% group_by(patid) %>% summarise(comorb_count = n())

#And CKD stage
analysis <- cprd$analysis("all_patid")
ckd_stage <- ckd_stage %>% analysis$cached("ckd_stages_from_algorithm")
ckd_stage_local <- ckd_stage %>% collect()

#Pivot to long.
ckd_long <- ckd_stage_local %>%
  pivot_longer(cols = starts_with("stage"), names_to = "stage", values_to = "date") %>%
  filter(!is.na(date)) %>%
  arrange(patid, stage)

#filter for severity
ckd_long_filtered <- ckd_long %>% filter(stage %in% c("stage_3a", "stage_3b", 
                                                      "stage_4", "stage_5"))

#Get Ids of people who develop CKD just prior to death.
ckd_dates <- ckd_long_filtered %>% 
  inner_join(last_row_hb %>% select(patid,dm_diag_date, start, stop)) %>%
  mutate(datediff = round(as.numeric(difftime(date, dm_diag_date, units = "days"))/365.25, 2)) %>%
  filter(datediff > start & datediff <= stop) %>% distinct(patid) %>% 
  mutate(ckd_score = 1)

#Add the individuals with no new comorbidities back in
last_row_final <- last_row_hb %>% 
  left_join(comorb_counts) %>% 
  left_join(ckd_dates %>% select(patid, ckd_score)) %>% as.data.frame() %>%
  mutate(comorb_count = replace_na(comorb_count, 0), 
         ckd_score = replace_na(ckd_score, 0)) %>% 
  mutate(incident_comorbidity_score_timevar = as.numeric(incident_comorbidity_score_timevar), 
         incident_comorbidity_score_timevar = incident_comorbidity_score_timevar + comorb_count + ckd_score)

#Arrange these rows like the previous rows.
last_row_formatted <- last_row_final %>% select(all_of(colnames(hu_local_start)))

#Format before joining
hu_local_formatted <- hu_local_start %>% convert_int_to_numeric()
last_row_formatted <- last_row_formatted %>% convert_int_to_numeric()

#Bind the dfs together
final_timepoints <- hu_local_formatted %>% bind_rows(last_row_formatted) %>%
  arrange(patid, hba1c_date_yrs)

#Calculate the "event" column.
final_mort <- final_timepoints %>% 
  mutate(death_time = round(as.numeric(difftime(cprd_ddate, dm_diag_date))/365.25, 2))
         
final_mort <- final_mort %>%         
    mutate(event = case_when(
    is.na(death_time) ~ 0,  # If date is missing, assign 0
    !is.na(death_time) & stop < death_time ~ 0,  # If within stop period, assign 0
    stop == death_time ~ 1  # If beyond stop period, assign 1
  ))

#Now finalise time-varying depression
final_dep <- final_mort %>% 
  mutate(time_to_incident_depression = 
         case_when(pret2d_dep == 1 ~ NA, 
                   pret2d_dep == 0 & postt2d_dep == 0 ~ NA,
                   pret2d_dep == 0 & postt2d_dep == 1 ~ round(as.numeric(difftime(postt2d_dep_diagdate, dm_diag_date, unit="days"))/365.25,2)))


#final_dep %>% filter(pret2d_dep == 0) %>% distinct(patid, .keep_all = T) %>% 
#filter(postt2d_dep == 1) %>% count() ##22875

final_dep_time <- final_dep %>% mutate(incident_depression_timevar = case_when(
  is.na(time_to_incident_depression) ~ 0,  # If date is missing, assign 0
  !is.na(time_to_incident_depression) & stop < time_to_incident_depression ~ 0,  # If within stop period, assign 0
  !is.na(time_to_incident_depression) & stop >= time_to_incident_depression ~ 1  # If beyond stop period, assign 1
))



#Check the df
for_survival <- final_dep_time %>% as.data.frame() %>%
  mutate(gender = as.factor(gender), pret2d_dep = as.factor(pret2d_dep), 
         ethnicity_bin = as.factor(ethnicity_bin), pret2d_dep_recurrent = as.factor(pret2d_dep_recurrent), 
         n_pret2d_dep_codes = as.numeric(n_pret2d_dep_codes), 
         last_pret2d_dep_10yrs = as.factor(last_pret2d_dep_10yrs), 
         first_pret2d_dep_10yrs = as.factor(first_pret2d_dep_10yrs), 
         first_pret2d_dep_diff = first_pret2d_dep_diff/365.25, 
         last_pret2d_dep_diff = last_pret2d_dep_diff/365.25)



#check the model
coxmod <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender, data = for_survival)

for_survival_incident <- for_survival %>% filter(pret2d_dep == 0)
for_survival_pre <- for_survival %>% filter(pret2d_dep == 1)

coxmod_recurrent <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + pret2d_dep_recurrent, data = for_survival_pre) #Recurrent is associated.
coxmod_depcodes <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + n_pret2d_dep_codes, data = for_survival_pre)
coxmod_first10 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + first_pret2d_dep_10yrs, data = for_survival_pre)
coxmod_last10 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + last_pret2d_dep_10yrs, data = for_survival_pre)
coxmod_timesincefirst <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + first_pret2d_dep_diff, data = for_survival_pre)
coxmod_timesincelast <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + last_pret2d_dep_diff, data = for_survival_pre)


for_upload <- for_survival %>% 
  mutate(maximum_followup_10yrs = pmin(maximum_followup_year, 10)) %>% 
  relocate(last_col(), .before=5)

for_upload <- for_upload %>% select(-incident_dep_timevar) %>% relocate(incident_depression_timevar, .before = "ethnicity_bin") %>%
  relocate(time_to_incident_depression, .before = "ethnicity_bin")
names(for_upload)




for_upload <- for_upload %>% mutate(pret2d_dep_recurrent = as.numeric(pret2d_dep_recurrent), 
                                    first_pret2d_dep_10yrs = as.numeric(first_pret2d_dep_10yrs), 
                                    last_pret2d_dep_10yrs = as.numeric(last_pret2d_dep_10yrs))





#Finally... upload the file
con2 <- dbConnect(MariaDB())



#Create empty table
dbExecute(con2, "
CREATE TABLE dh_hu_hba1c_visits_dataset_cleaned_v12 (
    patid BIGINT UNSIGNED, 
    hba1c_date_yrs BIGINT UNSIGNED, 
    count BIGINT UNSIGNED, 
    maximum_followup_year BIGINT UNSIGNED, 
    maximum_followup_10yrs BIGINT UNSIGNED,
    gender BIGINT UNSIGNED,
    dm_diag_age FLOAT, 
    dm_diag_date DATE, 
    cprd_ddate DATE, 
    gp_death_end_date DATE, 
    ethnicity_5cat BIGINT UNSIGNED, 
    pret2d_dep BIGINT UNSIGNED, 
    pret2d_dep_diagdate DATE, 
    dep_t2d_diff FLOAT, 
    postt2d_dep BIGINT UNSIGNED, 
    postt2d_dep_diagdate DATE, 
    dep_postt2d_diff FLOAT, 
    first_pret2d_dep_diff FLOAT, 
    n_pret2d_dep_codes BIGINT UNSIGNED, 
    pret2d_dep_recurrent BIGINT UNSIGNED, 
    last_pret2d_dep_diff FLOAT, 
    first_pret2d_dep_10yrs FLOAT, 
    last_pret2d_dep_10yrs FLOAT,
    mortality BIGINT UNSIGNED, 
    time_to_mortality FLOAT, 
    mortality_10yrs BIGINT UNSIGNED,
    time_to_mortality_10yrs FLOAT, 
    incident_depression BIGINT UNSIGNED, 
    incident_depression_year FLOAT, 
    incident_depression_timevar BIGINT UNSIGNED,
    time_to_incident_depression FLOAT,
    ethnicity_bin BIGINT UNSIGNED, 
    imd_decile BIGINT UNSIGNED, 
    tds_2011 FLOAT,
    preexisting_comorb_score BIGINT UNSIGNED,
    incident_comorbidity_score_timevar BIGINT UNSIGNED,
    year_of_diagnosis BIGINT UNSIGNED, 
    hba1c_index FLOAT, 
    bmi_index FLOAT,
    start FLOAT,
    stop FLOAT, 
    time_to_censor FLOAT, 
    death_time FLOAT,
    event BIGINT UNSIGNED 
);
")

chunksize <- 10000
total_rows <- nrow(for_upload)
nchunks <- ceiling(total_rows / chunksize)

# Start the transaction
dbExecute(con2, "START TRANSACTION;")

for (i in seq_len(nchunks)) {
  
  # Determine row indices for this chunk
  idx <- ((i - 1) * chunksize + 1):min(i * chunksize, total_rows)
  chunk <- for_upload[idx, ]
  
  # Process each column vector-wise
  # The following 'lapply' returns a list where each element is a character vector
  # with each value already appropriately escaped/quoted
  processed <- lapply(chunk, function(col) {
    if (is.character(col) || is.factor(col)) {
      # Convert factors to character; escape single quotes and wrap in quotes
      paste0("'", gsub("'", "''", as.character(col)), "'")
    } else if (inherits(col, "Date")) {
      # Format dates or use NULL
      ifelse(is.na(col), "NULL", paste0("'", format(col, "%Y-%m-%d"), "'"))
    } else if (is.numeric(col)) {
      # For numerics, if NA/NaN/Inf then return NULL (without quotes)
      ifelse(is.na(col) | is.nan(col) | is.infinite(col),
             "NULL",
             format(col, scientific = FALSE))
    } else {
      # Fallback for other data types
      rep("NULL", length(col))
    }
  })
  
  # Now combine the processed columns into row strings.
  # do.call(..., sep=",") is equivalent to pasting together the elements by row.
  row_strings <- do.call(paste, c(processed, sep = ","))
  # Wrap each row in parentheses
  row_strings <- paste0("(", row_strings, ")")
  
  # Build the full multi-row INSERT query
  query <- paste0(
    "INSERT INTO dh_hu_hba1c_visits_dataset_cleaned_v12 (",
    'patid, hba1c_date_yrs, count, maximum_followup_year, maximum_followup_10yrs, gender, dm_diag_age, dm_diag_date, ',
    'cprd_ddate, gp_death_end_date, ethnicity_5cat, pret2d_dep, pret2d_dep_diagdate, ',
    'dep_t2d_diff, postt2d_dep, postt2d_dep_diagdate, dep_postt2d_diff, first_pret2d_dep_diff, ',
    'n_pret2d_dep_codes, pret2d_dep_recurrent, last_pret2d_dep_diff, first_pret2d_dep_10yrs, ',
    'last_pret2d_dep_10yrs, mortality, time_to_mortality, mortality_10yrs, time_to_mortality_10yrs,  ',
    'incident_depression, incident_depression_year, incident_depression_timevar, time_to_incident_depression, ethnicity_bin, ',
    'imd_decile, tds_2011, preexisting_comorb_score, incident_comorbidity_score_timevar,
    year_of_diagnosis, hba1c_index, bmi_index, start, stop, time_to_censor, death_time, event
    ) VALUES ',
    paste(row_strings, collapse = ",")
  )
  
  # Execute the query for this chunk
  dbExecute(con2, query)
  
  print(paste("Inserted chunk", i, "of", nchunks))
}

