# Introduction

This is the GitHub to provide the code for the manuscript entitled "Depression increases both glycated haemoglobin testing frequency and risk of all-cause mortality in type 2 diabetes"
(https://www.medrxiv.org/content/10.1101/2025.10.06.25337242v1).

Briefly, this manuscript aims to: 
- Determine the effect of depression on HbA1c testing frequency in individuals with type 2 diabetes
- Determine the effects of depression and HbA1c testing frequency on risk of all-cause mortality
- Determine whether the increased testing in individuals with depression offsets the increased risk of mortality in individuals with co-morbid t2d and depression.

These scripts are designed to reproduce the dataset as defined in this manuscript, using the CPRD Aurum 2023 December release. Please see the supplementary methods/results of the manuscript
for all code lists used in these scripts. Alternatively, most of the code lists are also available at: 
https://github.com/Exeter-Diabetes/CPRD-Codelists

Prior to running these scripts, the data set should be curated according to previous work:
https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts

# Script overview

There are four main scripts which comprise this project:

- **1.0 generate_testing_frequency.R**: This script generates the HbA1c testing frequency variable which is used in all further analyses. This variable summarises the number of HbA1c measurements an individual had on each unique date since index per year, and summarises these values into a more friendly table. Where individuals did not have any tests in a year, this value is imputed to be zero.

- **2.0 clean_testing_frequency.R**: This script adds the main testing frequency variable to the main table, and harmonises all data sources used in the manuscript into a single long-format table which is model ready.

- **3.0 generate_CPMMEM_models.R**: This scripts generates the Conway-Maxwell-Poisson Mixed effects models on HbA1c testing over the first ten years of diagnosis to determine the effect of pre-existing and incident depression on HbA1c testing frequency. This script also adds the interval variable to account for individuals who are censored mid-way through a year.

- **4.0 - time_to_event_models.R**: This script models the various time-to-event analyses we perform in the manuscript. First, we determine the effect of pre-existing and incident depression on all-cause mortality. Next, we determine the effect of HbA1c testing frequency on mortality. Finally, we adjust the model looking at the effect of depression on all-cause mortality for HbA1c testing frequency and look at the change in Hazard ratio. This script also includes all sensitivity analyses performed in the manuscript, including the Therneau-Grambsch tests to determine violation of the proportional hazards assumptions and the subsequent time-stratified cox models and relevant plots. 
