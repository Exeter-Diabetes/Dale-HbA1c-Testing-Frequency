# Introduction

This is the GitHub to provide the code for the manuscript entitled "Depression increases both glycated haemoglobin testing frequency and risk of all-cause mortality in type 2 diabetes"
(https://www.medrxiv.org/content/10.1101/2025.10.06.25337242v1).

These scripts are designed to reproduce the dataset as defined in this manuscript, using the CPRD Aurum 2023 December release. Please see the supplementary methods/results of the manuscript
for all code lists used in these scripts. Alternatively, most of the code lists are also available at: 
https://github.com/Exeter-Diabetes/CPRD-Codelists

Prior to running these scripts, the data set should be curated according to previous work:
https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts

# Script overview

There are four main scripts which comprise this project:

- **1.0 generate_testing_frequency.R**: This script generates the HbA1c testing frequency variable which is used in all further analyses. This variable summarises the number of HbA1c measurements an individual had on each unique date since index per year, and summarises these values into a more friendly table. Where individuals did not have any tests in a year, this value is imputed to be zero.

- **2.0 clean_testing_frequency.R**: This script adds the main testing frequency variable to the main table, and harmonises all data sources used in the manuscript into a single long-format table which is model ready.

- 3.0
