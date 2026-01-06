####Script preparation####
rm(list = ls())
library(aurum)
library(tidyverse)
library(ggplot2)
library(data.table)
library(DBI)
library(RMariaDB)
library(dbplyr)

#Script for all Qc steps that need to be included:
cprdenvname <- ""
yaml <- ""

#Open the CPRD data
cprd <- CPRDData$new(cprdEnv = cprdenvname, cprdConf = yaml)


#Load in tables
#HbA1c table
analysis <- cprd$analysis("at_diag")
hb <- hb %>% analysis$cached("full_hba1c_index_date_merge")


#Cohort table
analysis <- cprd$analysis("dh_t2")
cohort <- cohort %>% analysis$cached("cohort_pracex_35_ins1yrex")


####HbA1c table QC####
#Join the cohort table and the ids from the cohort table
hb_filtered <- hb %>% inner_join(cohort, by = "patid")

#Check if any individuals have a date on the same day.
#hb_filtered %>% group_by(patid, datediff) %>% filter(n() > 1) %>% count() 
#No individuals with multiple HbA1c measurements on the same day - yay!
#In here, I should probably include measurements that occur on the same day.
#Limit data to values within 10 years. 
hb_10 <- hb_filtered %>% filter(datediff >= 0)  #8003002 values left

#then filter individuals with only 3+ Hba1c measurements
hb_values <- hb_10 %>% group_by(patid) %>% 
  filter(n() >= 3) %>%
  ungroup()

#Save the table.
analysis <- cprd$analysis("dh_hu")
hb_values %>% analysis$cached("hba1c_interim_3hba1c") 

#Load the table back in
hb_values <- hb_values %>% analysis$cached("hba1c_interim_3hba1c")

#Count number of individuals remaining.
#hb_values %>% distinct(patid) %>% count() #714185 - good filter.

#Calculate the difference in years - as.integer rounds the values up -> the measure
#Calculates all values that happened in the previous year.
hb_years <- hb_values %>% mutate(hba1c_date_yrs = ceiling(datediff/365.25)) 

#Save the table.
hb_years %>% analysis$cached("hba1c_interim_10yrs_3hba1c_yrrounded")

#Reload the table
hb_years <- hb_years %>% analysis$cached("hba1c_interim_3hba1c_yrrounded")

#Calculate the maximum number of years of follow up.
hb_time <- hb_years %>%
  mutate(maximum_followup_year = sql("CEIL(TIMESTAMPDIFF(YEAR, `index_date`, `gp_death_end_date`))"))


#Filter out values which occur in the same year or after the maximum followup_year
#Values which are in the maximum follow up year are allowed.
hb_time_filt <- hb_time %>% filter(hba1c_date_yrs <= maximum_followup_year) #this cuts only 5% of values?


#From where, we can save the data again
hb_time_filt %>% analysis$cached("hba1c_interim_3hb_maxfollowup")

#Read the data back in
hb_filt <- hb_filter %>% analysis$cached("hba1c_interim_3hb_maxfollowup")


####Generate the counts table####

#Generate the counts table
hba1c_visits <-  hb_filt %>%
  mutate(hba1c_date_yrs = ifelse(hba1c_date_yrs == 0, 1, hba1c_date_yrs)) %>% #Values on index date are set to 1 year.
  count(patid, hba1c_date_yrs, name = "count") %>% 
  left_join(hb_filt %>% select(patid, maximum_followup_year) %>% #Add the maximum follow up years back in. 
  distinct(patid, .keep_all = T)) %>%
  analysis$cached("hba1c_visits_raw") #Save as a new table.

#hba1c_visits_raw %>% count() #4678935 rows
#hba1c_visits_raw %>% distinct(patid) %>% count() #For 712825

hba1c_visits_raw <- hba1c_visits_raw %>% analysis$cached("hba1c_visits_raw")

#Limit to individuals with at least two values in the first 10 years.
hba1c_visits_filtered <- hba1c_visits_raw %>% 
   mutate(max_values = pmin(maximum_followup_year, 10)) %>% #takes into account people who have more than 10yrs follow up.
   group_by(patid) %>%
   filter((n() >= 2) | (n() == 1 &  max_values > 1)) %>% ungroup() %>% #2 years follow up with at least 1 Hb measurement.
  analysis$cached("hba1c_visits_min2")


#Load in the file locally:
analysis <- cprd$analysis("dh_hu")
hba1c_visits_local <- hba1c_visits_filtered %>% 
  analysis$cached("hba1c_visits_min2") %>%
  collect()

#check how many people are left
hba1c_visits_local %>% head() 
hba1c_visits_local %>% summarize(unique_count = n_distinct(patid)) #644084


#order
hba1c_visits_local <- hba1c_visits_local %>% arrange(patid, hba1c_date_yrs)



#Maybe I could generate the DF, then filter by the date, then 
names_hb <- hba1c_visits_local %>% group_by(patid) %>% slice(1) %>% 
  ungroup() %>% 
  select(patid, maximum_followup_year) 

death_local <- cohort  %>% select(patid, dm_diag_date, gp_death_end_date) %>% collect()

hb_death <- hba1c_visits_local %>% left_join(death_local) 
#Checking the code above has worked correctly.
#hb_death <- hb_death %>% mutate(death_year = round(as.numeric(difftime(gp_death_end_date, dm_diag_date, units = "days")/365.25),1))

#Create a file.
new_df <-  names_hb %>%
  tidyr::uncount(maximum_followup_year) %>%
  group_by(patid) %>%
  mutate(hba1c_date_yrs = row_number(), count = 0) %>%
  ungroup() %>% inner_join(names_hb %>% select(patid, maximum_followup_year)) 


#Merge the files together
new_merge <- new_df %>%  
  bind_rows(hba1c_visits_local %>% mutate(count = as.numeric(count))) %>% arrange(patid, hba1c_date_yrs) %>%
  mutate(count = as.numeric(count)) %>%
  select(-max_values)

  

#Change to a dt, as this is faster when using a lot of data.
setDT(new_merge)  # Convert to data.table in place - takes about 5 mins to run.

#Keep only the maximum values - unnecessary 0s generated above will be removed if
#There is another value present. 
#This is a lot faster than using the complete() dplyr function. 
new_formatted <- new_merge[, .SD[which.max(count)[1]], by = .(patid, hba1c_date_yrs)]

#Check if this is correct - the df should be maximum followup * number of individuals
#For the number of rows. 

#All the numbers are adding up.
#new_formatted %>% distinct(patid) %>% count() #should be 644084 - 
#new_formatted_chk <- new_formatted %>% mutate(maximum_total = ifelse(as.numeric(maximum_followup_year) > 10, 10, as.numeric(maximum_followup_year)))
#mean_followup <- new_formatted_chk %>% distinct(patid, .keep_all=T) %>% summarise(mean_fu = mean(as.numeric(maximum_total)))


####Descriptive plots####

#The new data frame should now be the same length as new_df_filtered, 
#Except it should have populated values with 0s where no values were available. 
#filtered_patids <- new_formatted %>%
#  group_by(patid) %>%
#  filter(n() >= 10) %>%          # Only keep patid groups with 10 or more data points
#  ungroup()

#first_10_patids <- filtered_patids %>% 
#  filter(patid %in% unique(patid)[1:20])

#Plot out the IDs.
#ggplot(first_10_patids, aes(x = hba1c_date_yrs, y = count)) +
#  geom_line() +                           # Line for each trajectory
#  geom_point() +                          # Points to highlight data points
#  labs(
#       x = "Years since T2D diagnosis",
#       y = "Number of HbA1c measurements per year") +
#  facet_wrap(~ patid, nrow = 4, ncol = 5) + # Facet by patid, allowing y-axis to vary per panel
#  theme_minimal() + scale_x_continuous(breaks = 1:10) + scale_y_continuous(breaks = 0:8) + 
#  theme(strip.text = element_blank(), 
#  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))


#summary_data <- filtered_patids %>%
#  mutate(count = as.numeric(count)) %>% group_by(hba1c_date_yrs) %>%
#  summarize(
#    mean_count = mean(count, na.rm = TRUE),
#    n = n(),
#    # Calculate the confidence intervals using Poisson approximation
#    ci_lower = qpois(0.025, lambda = mean_count * n) / n,
#    ci_upper = qpois(0.975, lambda = mean_count * n) / n
#  )

#ggplot(summary_data, aes(x = hba1c_date_yrs, y = mean_count)) +
#  geom_col(fill = "black", alpha = 0.5) +                   # Bar plot for each mean count
#  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +  # Error bars for Poisson CI
#  labs(
#    x = "Years Since T2D Diagnosis",
#    y = "Number of HbA1c measurements"
#  ) +
#  theme_minimal() + scale_x_continuous(breaks = 1:10) + 
#  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5, 2.0)) + expand_limits(y = 2.0)

####Upload output of script####

con2 <- dbConnect(MariaDB())



#Create empty table
dbExecute(con2, "
CREATE TABLE dh_hu_hba1c_visits_clean (
    patid BIGINT UNSIGNED, 
    hba1c_date_yrs BIGINT,
    count BIGINT,
    maximum_followup_year BIGINT
);
")

# Define the chunk size and get row count
chunksize <- 10000
total_rows <- nrow(new_formatted)
nchunks <- ceiling(total_rows / chunksize)

# Start the transaction
dbExecute(con2, "START TRANSACTION;")

for (i in seq_len(nchunks)) {
  
  # Determine row indices for this chunk
  idx <- ((i - 1) * chunksize + 1):min(i * chunksize, total_rows)
  chunk <- new_formatted[idx, ]
  
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
    "INSERT INTO dh_hu_hba1c_visits_clean (",
    "patid, hba1c_date_yrs, count, maximum_followup_year) VALUES ",
    paste(row_strings, collapse = ",")
  )
  
  # Execute the query for this chunk
  dbExecute(con2, query)
  
  print(paste("Inserted chunk", i, "of", nchunks))
}


