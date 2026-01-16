####Required####
rm(list = ls())
library(tidyverse)
library(survival)
library(splines)
library(aurum)
library(survminer)
library(rstpm2)
library(timereg)
library(gridExtra)
library(rstpm2)
library(openxlsx)
library(broom)
library(ggfortify)
install.packages("ggh4x")
library(ggh4x)
install.packages("rms")

#Helper function for extracting statistics from survival models into a pretty table
extract_hr_coxph <- function(models, model_names, variable) {
  if (length(models) != length(model_names)) {
    stop("Length of models and model_names must be the same.")
  }
  
  results <- data.frame(name = character(), HR = numeric(), lowerCI = numeric(),
                        upperCI = numeric(), CI_95 = character(), P_Value = numeric(), 
                        variable = character(),  
                        stringsAsFactors = FALSE)
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    name <- model_names[i]
    current_var <- variable[i]
    
    # Extract coefficients
    coefs <- summary(model)$coefficients
    conf_int <- summary(model)$conf.int
    
    if (!(current_var %in% rownames(coefs))) {
      print("current variable not in model")
      print(current_var)
      print(rownames(coefs))
      stop()  # Skip if the variable is not in the model
    }
    
    hr <- conf_int[current_var, "exp(coef)"]
    ci_lower <- conf_int[current_var, "lower .95"]
    ci_upper <- conf_int[current_var, "upper .95"]
    p_value <- coefs[current_var, "Pr(>|z|)"]
    
    # Create formatted 95% CI column
    ci_95 <- paste0("(", signif(ci_lower, 3), " - ", signif(ci_upper, 3), ")")
    
    # Append results to dataframe
    results <- rbind(results, data.frame(Name = name, HR = hr, lowerCI = ci_lower,
                                         upperCI = ci_upper,
                                           CI_95 = ci_95, P_Value = p_value, 
                                         variable = current_var, 
                                         stringsAsFactors = FALSE))
  }
  
  return(results)
}


extract_hr_coxph_multi <- function(models, model_names, variables) {
  if (length(models) != length(model_names) || length(models) != length(variables)) {
    stop("Length of models, model_names, and variables must be the same.")
  }
  
  results <- data.frame(Name = character(), HR = numeric(), lowerCI = numeric(),
                        upperCI = numeric(), CI_95 = character(), P_Value = numeric(), 
                        Variable = character(),  
                        stringsAsFactors = FALSE)
  
  for (i in seq_along(models)) {
    model <- models[[i]]
    name <- model_names[i]
    vars <- variables[[i]]
    
    # Extract coefficients and confidence intervals
    coefs <- summary(model)$coefficients
    conf_int <- summary(model)$conf.int
    
    for (var in vars) {
      if (!(var %in% rownames(coefs))) {
        warning(paste("Variable", var, "not in model", name))
        next
      }
      
      hr <- conf_int[var, "exp(coef)"]
      ci_lower <- conf_int[var, "lower .95"]
      ci_upper <- conf_int[var, "upper .95"]
      p_value <- coefs[var, "Pr(>|z|)"]
      
      # Create formatted 95% CI column
      ci_95 <- paste0("(", signif(ci_lower, 3), " - ", signif(ci_upper, 3), ")")
      
      # Append to results
      results <- rbind(results, data.frame(Name = name, HR = hr, lowerCI = ci_lower,
                                           upperCI = ci_upper, CI_95 = ci_95, 
                                           P_Value = p_value, Variable = var, 
                                           stringsAsFactors = FALSE))
    }
  }
  
  return(results)
}


convert_int_to_numeric <- function(data) {
  data %>%
    mutate(across(where(is.numeric), as.numeric))
}


####Connection and establish dataset####
cprdenvname <- "CPRD_diabetes_data"
yaml <- ""

cprd <- CPRDData$new(cprdEnv = cprdenvname, cprdConf = yaml)
analysis <- cprd$analysis("dh_hu")
hu <- hu %>% analysis$cached("hba1c_visits_dataset_cleaned_v13")
hu_local <- hu %>% collect() #bring local

#Figure 2 models.
#Format the table with the helper function
hu_local_formatted <- convert_int_to_numeric(hu_local)
hu_local_formatted <- data.frame(hu_local_formatted)


#perform some QC checks.
hu_local_formatted_round <- hu_local_formatted %>%
  mutate(
    count_lag1 = lag(count), 
    n_unique_hu_dates_lag = lag(n_unique_hu_dates)
  ) %>%
  filter(!is.na(count_lag1)) %>%
  mutate(count_cat_lag1 = case_when(
    count_lag1 == 0 ~ "0", 
    count_lag1 > 0 & count_lag1 < 3 ~ "1-2", 
    count_lag1 >= 3 & count_lag1 < 5 ~ "3-4", 
    count_lag1 >= 5 ~ "5+"
  ))

table(hu_local_formatted_round$count_cat_lag1)

hu_local_formatted_round <- hu_local_formatted_round %>% filter(start >= 2)


model_hu <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + n_unique_hu_dates_lag + ethnicity_bin + year_of_diagnosis + count_lag1, data = hu_local_formatted_round_short)
model_hu_unadj  <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + year_of_diagnosis + count_lag1 + cluster(patid), data = hu_local_formatted_round_short)

#Try adjusting for baseline Hb.
model_hu_unadj_hb <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + year_of_diagnosis + count_lag1 + hba1c_index, data = hu_local_formatted_round_short)
model_hu_adj_hb <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + n_unique_hu_dates_lag + ethnicity_bin + year_of_diagnosis + count_lag1 + hba1c_index + cluster(patid), 
                         data = hu_local_formatted_round_short)



summary(model_hu)
summary(model_hu_unadj)



#These models are run again later... this is more for the text in the paper to justify these models.
model_test <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_lag1, data = hu_local_formatted_round)
model_test_covs <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_lag1 + incident_comorbidity_score_timevar + preexisting_comorb_score, data = hu_local_formatted_round)

model_test_cat <- coxph(Surv(start, stop_ceil, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1, data = hu_local_formatted_round)
model_test_cat_covs <- coxph(Surv(start, stop_ceil, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_comorbidity_score_timevar + preexisting_comorb_score, data = hu_local_formatted_round)


#tidy models.
tidy(model_test, exponentiate = T, conf.int = T)
tidy(model_test_covs, exponentiate = T, conf.int = T)
tidy(model_test_cat_covs, exponentiate = T, conf.int = T)

#Further formatting, and create the incident and pre-existing depression tables.
nopredep_df <- hu_local_formatted_round %>% filter(pret2d_dep == 0) %>% 
  mutate(patid = as.factor(patid), 
         ethnicity_bin = as.factor(ethnicity_bin), 
         gender = as.factor(gender), 
         incident_depression_timevar = as.factor(incident_depression_timevar), 
         pret2d_dep = as.factor(pret2d_dep))

predep_df <- hu_local_formatted_round %>% 
  mutate(patid = as.factor(patid), 
         ethnicity_bin = as.factor(ethnicity_bin), 
         gender = as.factor(gender), 
         incident_depression_timevar = as.factor(incident_depression_timevar), 
         pret2d_dep = as.factor(pret2d_dep))


zst_test <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = predep_df[predep_df$start > 3,])
zph_res <- cox.zph(zst_test)

zst_test_inc <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start > 2,])
zph_res_inc <- cox.zph(zst_test_inc)


full_model_predep <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start >= 3, ])



#Stratify the model by 3 - 6yrs, 6 - 9yrs, and 9 - 12yrs
time_3 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start >= 3 & nopredep_df$start <= 5,])
time_6 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start > 5 & nopredep_df$start <= 8,])
time_9 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start > 8 & nopredep_df$start <= 11,])
time_12 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = nopredep_df[nopredep_df$start > 11 & nopredep_df$start <= 14,])


time_3_predep <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = predep_df[predep_df$start >= 3 & predep_df$start <= 5,])
time_6_predep <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = predep_df[predep_df$start > 5 & predep_df$start <= 8,])
time_9_predep <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = predep_df[predep_df$start > 8 & predep_df$start <= 11,])
time_12_predep <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid) + incident_comorbidity_score_timevar + preexisting_comorb_score, data = predep_df[predep_df$start > 11 & predep_df$start <= 14,])


get_hr_data <- function(model, time_label, variable) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(variable, term)) %>%
    mutate(time_period = time_label)
}

extract_hr_over_time <- function(models, time_labels, variable, var_label = NULL) {
  hr_data <- Map(function(model, label) {
    tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == variable) %>%
      mutate(time_period = label)
  }, models, time_labels) %>%
    bind_rows()
  
  hr_data$time_period <- factor(hr_data$time_period,
                                levels = time_labels,
                                ordered = TRUE)
  
  if (!is.null(var_label)) {
    hr_data$term <- var_label
  }
  
  return(hr_data)
}



time_labels <- c("2-5 years", "5-8 years", "8-11 years", "11-14 years")
models_incdep <- list(time_3, time_6, time_9, time_12)
models_predep <- list(time_3_predep, time_6_predep, time_9_predep, time_12_predep)

hr_incdep <- extract_hr_over_time(models_incdep, time_labels, "incident_depression_timevar1", "Incident Depression")
hr_predep <- extract_hr_over_time(models_predep, time_labels, "pret2d_dep1", "Pre-existing Depression")
hr_precom <- extract_hr_over_time(models_incdep, time_labels, "preexisting_comorb_score", "Pre-existing Comorbidity")
hr_inccom <- extract_hr_over_time(models_incdep, time_labels, "incident_comorbidity_score_timevar", "Incident Comorbidity")

plot_hr_variable <- function(hr_data, title_text) {
  ggplot(hr_data, aes(x = time_period, y = estimate)) +
    geom_line(group = 1, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "steelblue") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(
      title = title_text,
      x = "Time Period",
      y = "Hazard Ratio"
    ) +
    theme_minimal(base_size = 16)
}

plot_incdep <- plot_hr_variable(hr_incdep, "Hazard Ratio of Incident Depression Over Time")
plot_predep <- plot_hr_variable(hr_predep, "Hazard Ratio of Pre-existing Depression Over Time")
plot_precom <- plot_hr_variable(hr_precom, "Hazard Ratio of Pre-existing Comorbidity Score Over Time")
plot_inccom <- plot_hr_variable(hr_inccom, "Hazard Ratio of Incident Comorbidity Score Over Time")




# Ensure time_period is an ordered factor
hr_data <- bind_rows(
  get_hr_data(time_3, "2-5 years", "count_cat"),
  get_hr_data(time_6, "5-8 years", "count_cat"),
  get_hr_data(time_9, "8-11 years", "count_cat"),
  get_hr_data(time_12, "11-14 years", "count_cat")
)

hr_data$time_period <- factor(hr_data$time_period,
                              levels = c("2-5 years", "5-8 years", "8-11 years", "11-14 years"),
                              ordered = TRUE)

# Clean up variable names
hr_data <- hr_data %>%
  mutate(term = recode(term,
                       "count_cat_lag11-2" = "1–2",
                       "count_cat_lag13-4" = "3–4",
                       "count_cat_lag15+" = "5+"))

# Plot
ggplot(hr_data, aes(x = time_period, y = estimate, group = term, color = term)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Time Period",
    y = "Hazard Ratio",
    color = "HbA1c testing frequency"
  ) +
  theme_minimal(base_size = 16)



#Convert to table
zph_df <- as.data.frame(zph_res$table)
zph_df_inc <- as.data.frame(zph_res_inc$table)

zph_df_new <- zph_df %>%
  slice(1:7) %>%                              # Take rows 1 through 7
  bind_rows(zph_df_inc["incident_depression_timevar", ]) %>%  # Insert new row
  bind_rows(slice(zph_df, 8:n()))            # Then add the rest of the rows


zph_df_new$covariate <- rownames(zph_df)
rownames(zph_df_new) <- NULL
zph_df_new$variable <- c("Age at T2D diagnosis", "Gender", "Ethnicity", "IMD decile", 
                     "Year of T2D diagnosis", "HbA1c count category", "Pre-existing depression", "Incident depression", 
                     "Incident comorbidity score", "Pre-existing comorbidity score", "Global")



zph_df_new <- zph_df_new %>% select(variable, chisq, df, p)

#Write the table
write.table(zph_df_new, , row.names = F, quote = F, sep = ",")


#Get the plots
gg_list_main <- ggcoxzph(zph_res)
gg_list_inc  <- ggcoxzph(zph_res_inc)

zph_names_inc <- rownames(zph_res_inc$table)
incident_dep_index <- which(zph_names_inc == "incident_depression_timevar")

# Extract plot safely
incident_dep_plot <- ggcoxzph(zph_res_inc)[[incident_dep_index]]

# Original plots
gg_list_main <- ggcoxzph(zph_res)

# Insert incident_dep_plot into position 8 (after 7th plot)
gg_list_combined <- append(
  gg_list_main[1:7],
  list(incident_dep_plot),
  after = 7
)

# Add remaining plots
gg_list_combined <- c(
  gg_list_combined,
  gg_list_main[8:length(gg_list_main)]
)

for (i in seq_along(gg_list_combined)) {
  gg_list_combined[[i]] <- gg_list_combined[[i]] +
    ggtitle(names_to_use[i]) +
    ylab("Beta(t)")
}

plots <- gg_list_combined
n_plots <- length(plots)
ncol <- 3
n_filled <- ceiling(n_plots / ncol) * ncol

# Pad with empty grobs if needed
if (n_plots < n_filled) {
  n_pad <- n_filled - n_plots
  plots <- c(plots, rep(list(nullGrob()), n_pad))
}

# Now arrange and save
all_plots <- arrangeGrob(grobs = plots, ncol = ncol)


ggsave(
  filename = ,
  plot     = all_plots,
  width    = 11.69,
  height   = 8.27
)



#Raw
model_list_count <- list()
model_list_count[["count"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + cluster(patid), data = predep_df)

#postdep
model_list_count[["incdep_count"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + cluster(patid), data = nopredep_df)
model_list_count[["incdep_nocomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar+ cluster(patid), data = nopredep_df)
model_list_count[["incdep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score  + incident_depression_timevar + cluster(patid), data = nopredep_df)
model_list_count[["incdep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + incident_depression_timevar+ cluster(patid), data = nopredep_df)

#Predep
model_list_count[["predep_count"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + cluster(patid), data = predep_df)
model_list_count[["predep_nocomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid), data = predep_df)
model_list_count[["predep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + pret2d_dep + cluster(patid), data = predep_df)
model_list_count[["predep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + pret2d_dep + cluster(patid), data = predep_df)

#Check names of variables
vars_to_extract <- c("count_cat_lag11-2", "count_cat_lag13-4", "count_cat_lag15+")
variables <- rep(list(vars_to_extract), length(model_list_count))
model_names <- names(model_list_count)

#Get the models in a friendly format
figure_2_table <- extract_hr_coxph_multi(model_list_count, model_names, variables)


#Organise results
figure_2_table <- figure_2_table %>%
  mutate(
    Phase = str_extract(Name, "^(predep|incdep)"),
    Group = case_when(
      str_detect(Name, "count") ~ "count",
      str_detect(Name, "nocomorb") ~ "nocomorb",
      str_detect(Name, "icomorb") ~ "icomorb"
    ),
    Phase = factor(Phase, levels = c("predep", "incdep")),
    Group = factor(Group, levels = c("count", "nocomorb", "icomorb"))
  ) %>%
  arrange(Phase, Group)


#Add column for significance and change names to be more readable
figure_2_table <- figure_2_table %>% mutate(signif = if_else(P_Value < (0.05/8), "*", ""))
new_name_vector <- rep(c("Comorbidity unadjusted", "Comorbidity adjusted", "Comorbidity unadjusted", "Comorbidity adjusted"), each = 3)

figure_2_table$new_names <- new_name_vector
figure_2_table_names <- figure_2_table



#keep order
df <- figure_2_table_names %>%
  mutate(
    new_names = factor(new_names, levels = unique(new_names)),
    Variable = factor(Variable)
  )


n_lines <- length(levels(df$new_names))
hline_positions <- seq(1.5, n_lines - 0.5, by = 1)

# Dodge width
dodge_width <- 0.6
n_groups <- nlevels(df$Variable)

# Offsets centered around 0 (e.g., for 3 groups: -0.2, 0, 0.2)
offsets <- tibble(
  Variable = levels(df$Variable),
  offset = seq(
    from = -dodge_width / 2 + dodge_width / (2 * n_groups),
    to   =  dodge_width / 2 - dodge_width / (2 * n_groups),
    length.out = n_groups
  )
)

# Merge with df and compute y positions for asterisks
df <- df %>%
  left_join(offsets, by = "Variable") %>%
  mutate(
    y_numeric = as.numeric(new_names),
    y_pos_star = y_numeric + offset - 0.04
  )

n_lines <- length(levels(df$new_names))
hline_positions <- seq(1.5, n_lines - 0.5, by = 1)

df <- df %>% mutate(Phase = case_when(Phase == "predep" ~ "Pre-existing depression", 
                                      Phase == "incdep" ~ "Incident depression"), 
                    Phase = factor(Phase, levels = c("Pre-existing depression", "Incident depression")))


# Plot
ggplot(df, aes(x = HR, y = new_names, color = Variable)) +
  geom_errorbarh(
    aes(xmin = lowerCI, xmax = upperCI),
    height = 0.3,
    position = position_dodge(width = dodge_width),
    size = 0.7
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    shape = 15, size = 2.5
  ) +
  geom_text(
    data = df %>% filter(signif == "*"),
    aes(x = upperCI + 0.02, y = y_pos_star, label = "*"),
    color = "black",
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 1, linetype = "solid") +
  theme_bw() +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    color = "HbA1c measurements per year"
  ) +
  scale_color_manual(
    values = c(
      "count_cat_lag11-2" = "#F8766D",  # red
      "count_cat_lag13-4" = "#00BA38",  # green
      "count_cat_lag15+"  = "#619CFF"   # blue
    ),
    labels = c("1-2", "3-4", "5+")
  ) +
  scale_y_discrete(limits = levels(df$new_names)) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, hjust = 0)
  ) +
  geom_hline(yintercept = hline_positions, linetype = "dashed", color = "gray80", size = 0.3) +
  ggh4x::facet_wrap2(
    ~Phase,
    ncol = 1,
    strip.position = "top",
    axes = "x"
  )


#next, plot the list of variables that violate the ph assumption over time
#Using stratified models.
timestrat_1 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + pret2d_dep + cluster(patid), 
                     data = predep_df[predep_df$stop >= 2 & predep_df$stop <= 5,])

timestrat_2 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + pret2d_dep + cluster(patid), 
                     data = predep_df[predep_df$stop >= 6 & predep_df$stop <= 10,])

timestrat_3 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + pret2d_dep + cluster(patid), 
                     data = predep_df[predep_df$stop >= 10 & predep_df$stop <= 15,])

timestrat_4 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + incident_depression_timevar + cluster(patid), 
                     data = nopredep_df[nopredep_df$stop >= 2 & nopredep_df$stop <= 5,])

timestrat_5 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + incident_depression_timevar + cluster(patid), 
                     data = nopredep_df[nopredep_df$stop >= 6 & nopredep_df$stop <= 10,])

timestrat_6 <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + 
                       imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + 
                       incident_comorbidity_score_timevar + incident_depression_timevar + cluster(patid), 
                     data = nopredep_df[nopredep_df$stop >= 10 & nopredep_df$stop <= 15,])


hr1 <- tidy(timestrat_1, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "1–5")
hr2 <- tidy(timestrat_2, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "6–10")
hr3 <- tidy(timestrat_3, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "11–15")

hr4 <- tidy(timestrat_4, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "1–5")
hr5 <- tidy(timestrat_5, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "6–10")
hr6 <- tidy(timestrat_6, exponentiate = TRUE, conf.int = TRUE) %>% mutate(time_period = "11–15")

hr4 <- hr4[hr4$term == "incident_depression_timevar1",]
hr5 <- hr5[hr5$term == "incident_depression_timevar1",]
hr6 <- hr6[hr6$term == "incident_depression_timevar1",]

hr_combined <- bind_rows(hr1, hr2, hr3, hr4, hr5, hr6)

hr_combined$time_period <- factor(hr_combined$time_period, levels = unique(hr_combined$time_period))

hr_combined <- hr_combined %>% filter(!term %in% c("(Intercept)"))
hr_combined <- hr_combined %>% mutate(term = case_when(
                                      term == "count_cat_lag11-2" ~"HbMonFreq 1-2", 
                                      term == "count_cat_lag13-4" ~"HbMonFreq 3-4", 
                                      term == "count_cat_lag15+" ~ "HbMonFreq 5+", 
                                      term == "dm_diag_age" ~ "Age at T2D diagnosis", 
                                      term == "ethnicity_bin1" ~ "Ethnicity",
                                      term == "gender2" ~ "Gender", 
                                      term == "imd_decile" ~ "IMD decile", 
                                      term == "incident_comorbidity_score_timevar" ~ "IncComorbScore", 
                                      term == "preexisting_comorb_score" ~ "PreComorbScore", 
                                      term == "pret2d_dep1" ~ "PreDep", 
                                      term == "year_of_diagnosis" ~ "Year of diagnosis", 
                                      term == "incident_depression_timevar1" ~ "IncDep"))

ggplot(hr_combined, aes(x = time_period, y = estimate)) +
  geom_col(fill = "steelblue", width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ term, scales = "free", strip.position = "top", ncol = 3) +
  labs(
    x = "Years Since Diagnosis",
    y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 16) +  # Increase base font size
  theme(
    strip.text = element_text(size = 18, face = "bold"),   # Facet titles
    axis.title = element_text(size = 16),                  # Axis titles
    axis.text = element_text(size = 14),                   # Tick labels
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Main title
  )


#Generate the models and append them to a list.
model_list_predep <- list()

#Figure 2 part 1
model_list_predep[["predep_nocount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + pret2d_dep + cluster(patid), data = predep_df)
model_list_predep[["predep_nocomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep + cluster(patid), data = predep_df)
model_list_predep[["predep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + pret2d_dep + cluster(patid), data = predep_df)
model_list_predep[["predep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + pret2d_dep + cluster(patid), data = predep_df)

#Get the table ready
variables = c("pret2d_dep1", "pret2d_dep1", "pret2d_dep1")
figure_3a_table <- extract_hr_coxph(models = model_list_predep, model_names = names(model_list_predep), variable = variables)

#Get the data from the plots
figure_3a_plotready <- figure_3a_table %>% mutate(variable = str_replace_all(variable, "pret2d_dep1", "Pre-existing depression"), 
                                                variable = str_replace_all(variable, "incident_depression_timevar1", "Incident depression"), 
                                                Name = case_when(
                                                  str_detect(Name, "nocomorb") ~ "Comorbidity and HbA1c testing frequency unadjusted",
                                                  str_detect(Name, "nocount") ~ "HbA1c testing frequency adjusted",
                                                  str_detect(Name, "icomorb") ~ "Comorbidity and HbA1c testing frequency adjusted")) %>%
  mutate(star_label = ifelse(P_Value < 0.05, "*", "")) %>% mutate(Name = factor(Name, levels = c("Comorbidity and HbA1c testing frequency adjusted" ,
                                                                                              "HbA1c testing frequency adjusted", 
                                                                                              "Comorbidity and HbA1c testing frequency unadjusted")))





figure_3a <- ggplot(figure_3a_plotready, aes(x = Name, y = HR)) +
  # Point estimate
  geom_point(shape = 15, size = 2.5) +
  # Corrected error bars: use ymin/ymax for vertical error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Reference line at HR = 1
  geom_hline(yintercept = 1, linetype = "solid") +
  # Add star label above error bars
  geom_text(
    aes(label = star_label, y = upperCI * 1.02),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  coord_flip(clip = "off") + 
  theme_bw() +
  labs(
    x = NULL,
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  ) +
  ggtitle("Pre-existing depression") + ylim(NA, 1.5)
figure_3a



model_list_postdep <- list()
model_list_postdep[["incdep_nocount"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + incident_depression_timevar + cluster(patid), data = nopredep_df)
model_list_postdep[["incdep_nocomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar + cluster(patid), data = nopredep_df)
model_list_postdep[["incdep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score  + incident_depression_timevar + cluster(patid), data = nopredep_df)
model_list_postdep[["incdep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + incident_depression_timevar + cluster(patid), data = nopredep_df)
 

variables <- c("incident_depression_timevar1",  "incident_depression_timevar1", "incident_depression_timevar1")

figure_3b_table <- extract_hr_coxph(models = model_list_postdep, model_names = names(model_list_postdep), variable = variables)

figure_3b_plotready <- figure_3b_table %>% mutate(variable = str_replace_all(variable, "pret2d_dep1", "Pre-existing depression"), 
                                                variable = str_replace_all(variable, "incident_depression_timevar1", "Incident depression"), 
                                                Name = case_when(
                                                  str_detect(Name, "nocomorb") ~ "Comorbidity and HbA1c testing frequency unadjusted",
                                                  str_detect(Name, "nocount") ~ "HbA1c testing frequency adjusted",
                                                  str_detect(Name, "icomorb") ~ "Comorbidity and HbA1c testing frequency adjusted")) %>%
  mutate(star_label = ifelse(P_Value < 0.05, "*", "")) %>% mutate(Name = factor(Name, levels = c("Comorbidity and HbA1c testing frequency adjusted" ,
                                                                                                 "HbA1c testing frequency adjusted", 
                                                                                                 "Comorbidity and HbA1c testing frequency unadjusted")))



figure_3b <- ggplot(figure_3b_plotready, aes(x = Name, y = HR)) +
  # Point estimate
  geom_point(shape = 15, size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Dashed reference line at HR=1
  geom_hline(yintercept = 1, linetype = "solid") +
  # Add star to the right of error bars if p < 0.05
  geom_text(
    aes(label = star_label, y = upperCI * 1.02),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),  # X-axis text size
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18), 
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 14)
    # Y-axis text size
  ) + ggtitle("Incident depression")

figure_3b

figure_3 <- grid.arrange(figure_3a, figure_3b, ncol = 1)


#The variables to extract from each model. 

#Run the models with the function and pray.
figure_2_table <- extract_hr_coxph(models = model_list_postdep, model_names = names(model_list_postdep), variable = variables)

#Sort out the names for plotting
df <- figure_2_table %>%
  mutate(
    stacked_y = interaction(new_names, Variable, lex.order = TRUE),
    stacked_y = factor(stacked_y, levels = rev(levels(stacked_y)))  # flip for forest plot
  )

ggplot(df, aes(x = HR, y = stacked_y, xmin = lowerCI, xmax = upperCI, color = Variable)) +
  geom_pointrange(size = 0.7) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray50") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Forest Plot",
    x = "Hazard Ratio (95% CI)",
    y = "",
    color = "Variable"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# This is that plot.
figure_2_plotready <- figure_2_plotready %>%
  mutate(variable = factor(variable, 
                           levels = c("Pre-existing depression", 
                                      "Incident depression"))) %>%  # Ensure PED is above ID
  mutate(variable_name = paste(variable, Name, sep = " - ")) %>%
  mutate(
    variable_name = factor(variable_name, 
                           levels = c("Pre-existing depression - covariates + HMF + CS",
                                      "Pre-existing depression - covariates + HMF",
                                      "Pre-existing depress - covariates",
                                      "Incident depression - covariates + HMF + CS",
                                      "Incident depression - covariates + HMF",
                                      "Incident depression - covariates")))

write.table(figure_2_plotready, "", row.names = F, quote = F, sep = "\t")


#This is the big one
# These models have all been run before
# Covs + count
# Covs + count + predep
# Covs + count + PECS
# Covs + count + PECS + ICS
# Covs + count + postdep
# Covs + count + postdep + PECS
# Covs + count + postdep + PECS + ICS


model_list_count <- list()

predep_df_short <- predep_df[predep_df$start > 1, ]

#Just count
model_list_count[["count"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count, data = predep_df_short)

#Check for deviations from the proportional hazards assumptions.
#Re-run the model
sens_model <- model_list_count[["count"]]
sens_stats <- cox.zph(sens_model)
sens_plot <- plot(sens_stats, var = "count")



model_list_count[["count"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat, data = predep_df)



#postdep
model_list_count[["incdep_nocomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + incident_depression_timevar, data = nopredep_df)
model_list_count[["incdep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + preexisting_comorb_score  + incident_depression_timevar, data = nopredep_df)
model_list_count[["incdep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + preexisting_comorb_score + incident_comorbidity_score_timevar + incident_depression_timevar, data = nopredep_df)

#Predep
model_list_count[["predep_nocomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + pret2d_dep, data = predep_df)
model_list_count[["predep_pecomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + preexisting_comorb_score + pret2d_dep, data = predep_df)
model_list_count[["predep_icomorb"]] <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count + preexisting_comorb_score + incident_comorbidity_score_timevar + pret2d_dep, data = predep_df)


#Extract count
variables <- rep("count", 7)

figure_2_table <- extract_hr_coxph(models = model_list_count, model_names = names(model_list_count), variable = variables)

figure_2_plotready <- figure_2_table %>% mutate(variable = str_replace_all(variable, "pret2d_dep1", "Pre-existing depression"), 
                                                variable = str_replace_all(variable, "incident_depression_timevar1", "Incident depression"), 
                                                Name = case_when(
                                                  str_detect(Name, "count") ~ "Covariates",
                                                  str_detect(Name, "predep_nocomorb") ~ "Covariates + HMF + PED",
                                                  str_detect(Name, "predep_pecomorb") ~ "Covariates + HMF + PED + PECS", 
                                                  str_detect(Name, "predep_icomorb") ~ "Covariates + HMF + PECS + ICS", 
                                                  str_detect(Name, "incdep_nocomorb") ~ "Covariates +  HMF + ID", 
                                                  str_detect(Name, "incdep_pecomorb") ~ "Covariates + HMF + ID + PECS",
                                                  str_detect(Name, "incdep_icomorb") ~ "Covariates + HMF + ID + PECS + ICS")) %>%
                                                  mutate(star_label = ifelse(P_Value < 0.05, "*", ""))

figure_2_plotready$Name <- factor(figure_2_plotready$Name, levels = figure_2_plotready$Name)

figure_3 <- ggplot(figure_2_plotready, aes(x = Name, y = HR)) +
  # Point estimate
  geom_point(shape = 15, size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Dashed reference line at HR=1
  geom_hline(yintercept = 1, linetype = "solid") +
  # Add star to the right of error bars if p < 0.05
  geom_text(
    aes(label = star_label, y = upperCI * 1.05),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14),  # X-axis text size
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18), 
    legend.text = element_text(size = 14)
    # Y-axis text size
  )

figure_3




incdep <- read.xlsx("")
incdep$row <- 1:nrow(incdep)
incdep$row <- factor(incdep$row, levels = incdep$row)

figure_4 <- ggplot(incdep, aes(x = Name, y = IRR)) +
  # Point estimate
  geom_point(shape = 15, size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Dashed reference line at HR=1
  geom_hline(yintercept = 1, linetype = "solid") +
  facet_wrap(~grouping) + 
  # Add star to the right of error bars if p < 0.05
  geom_text(
    aes(label = star_label, y = upperCI * 1.05),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18),  # X-axis text size
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 18), 
    legend.text = element_text(size = 18), 
    axis.title = element_text(size = 20), 
    strip.text = element_text(size = 20)
    # Y-axis text size
  )
figure_4




#Actually run the depression severity variables. 
model_list_severity <- list()


#Create the predeponly
predeponly <- predep_df %>% filter(pret2d_dep == 1)

#Generate the models
model_list_sev <- list()


#Depression first within 10 years, and depression last within 10 years. 
predeponly$first_pret2d_dep_diff
predeponly$last_pret2d_dep_diff
predeponly$n_pret2d_dep_codes
predeponly$pret2d_dep_recurrent
predeponly$first_pret2d_dep_10yrs
predeponly$last_pret2d_dep_10yrs
predeponly$count_cat_lag1

#First depression
model_list_sev[["depsev_tsdd_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + first_pret2d_dep_diff + cluster(patid), data = predeponly)
model_list_sev[["depsev_tsdd_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + first_pret2d_dep_diff+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tsdd_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + first_pret2d_dep_diff+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tsdd_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + first_pret2d_dep_diff+ cluster(patid), data = predeponly)


#Last depression
model_list_sev[["depsev_tsld_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + last_pret2d_dep_diff+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tsld_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + last_pret2d_dep_diff+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tsld_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + last_pret2d_dep_diff+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tsld_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + last_pret2d_dep_diff+ cluster(patid), data = predeponly)

#total number of depression codes
model_list_sev[["depsev_tndc_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + n_pret2d_dep_codes+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tndc_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + n_pret2d_dep_codes+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tndc_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + n_pret2d_dep_codes+ cluster(patid), data = predeponly)
model_list_sev[["depsev_tndc_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + n_pret2d_dep_codes+ cluster(patid), data = predeponly)

#Recurrent depression
model_list_sev[["depsev_rd_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + pret2d_dep_recurrent+ cluster(patid), data = predeponly)
model_list_sev[["depsev_rd_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + pret2d_dep_recurrent+ cluster(patid), data = predeponly)
model_list_sev[["depsev_rd_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + pret2d_dep_recurrent+ cluster(patid), data = predeponly)
model_list_sev[["depsev_rd_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + pret2d_dep_recurrent+ cluster(patid), data = predeponly)

model_list_sev[["depsev_firstdep_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + first_pret2d_dep_10yrs+ cluster(patid), data = predeponly)
model_list_sev[["depsev_firstdep_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + first_pret2d_dep_10yrs + cluster(patid), data = predeponly)
model_list_sev[["depsev_firstdep_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + first_pret2d_dep_10yrs+ cluster(patid), data = predeponly)
model_list_sev[["depsev_firstdep_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + first_pret2d_dep_10yrs + cluster(patid), data = predeponly)

model_list_sev[["depsev_lastdep_covs"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + last_pret2d_dep_10yrs+ cluster(patid), data = predeponly)
model_list_sev[["depsev_lastdep_covscount"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + last_pret2d_dep_10yrs + cluster(patid), data = predeponly)
model_list_sev[["depsev_lastdep_precomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + last_pret2d_dep_10yrs+ cluster(patid), data = predeponly)
model_list_sev[["depsev_lastdep_allcomorb"]]<- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + imd_decile + year_of_diagnosis + count_cat_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar + last_pret2d_dep_10yrs + cluster(patid), data = predeponly)


#Two figures:
variables <- c(rep("first_pret2d_dep_diff", 4), rep("last_pret2d_dep_diff", 4), rep("n_pret2d_dep_codes", 4), rep("pret2d_dep_recurrent", 4), rep("first_pret2d_dep_10yrs", 4), rep("last_pret2d_dep_10yrs", 4))
depsev_figure_table <- extract_hr_coxph(models = model_list_sev, model_names = names(model_list_sev), variable = variables)

#Create new naming values
depsev_figure_table$covs <- rep(c("Covs", "Covs + HbMonFreq", "Covs + HbMonFreq \n+ PreComorbScore", "Covs + HbMonFreq \n+ PreComorbScore \n+ IncComorbScore"), 6)
depsev_figure_table$variable <- rep(c("Time since first depression diagnosis", "Time since last depression diagnosis", "Number of pre-existing depression codes", "Recurrent depression", "Depression diagnosed within 10 years prior to T2D", "Last depression code within 10 years prior to T2D"), each = 4)
depsev_figure_table$star_label <- ifelse(depsev_figure_table$P_Value < 0.05/8, "*", "")
depsev_figure_table$marker <- c(rep(0, 12), rep(1, 12))
depsev_figure_cont <- depsev_figure_table %>% filter(marker == 0)
depsev_figure_di <- depsev_figure_table %>% filter(marker == 1)


supfig_4_1 <- ggplot(depsev_figure_cont, aes(x = covs, y = HR)) +
  # Point estimate
  geom_point(shape = 15, size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Dashed reference line at HR=1
  geom_hline(yintercept = 1, linetype = "solid") +
  facet_wrap(~variable) + 
  # Add star to the right of error bars if p < 0.05
  geom_text(
    aes(label = star_label, y = upperCI * 1.05),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18),  # X-axis text size
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 18), 
    legend.text = element_text(size = 18), 
    axis.title = element_text(size = 20), 
    strip.text = element_text(size = 18)
    # Y-axis text size
  )
supfig_4_1

supfig_4_2 <- ggplot(depsev_figure_di, aes(x = covs, y = HR)) +
  # Point estimate
  geom_point(shape = 15, size = 3) +
  # Error bars
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2) +
  # Dashed reference line at HR=1
  geom_hline(yintercept = 1, linetype = "solid") +
  facet_wrap(~variable) + 
  # Add star to the right of error bars if p < 0.05
  geom_text(
    aes(label = star_label, y = upperCI * 1.05),
    vjust = 0.75,
    size = 6,
    color = "black",
    na.rm = TRUE
  ) +
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Hazard Ratio"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18),  # X-axis text size
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 18), 
    legend.text = element_text(size = 18), 
    axis.title = element_text(size = 20), 
    strip.text = element_text(size = 18)
    # Y-axis text size
  )

supfig_4_2



#First, generate the time-varying plot for 
df <- hu_local_formatted_round_short

# Create a depression status variable
df <- df %>%
  mutate(depression_status = case_when(
    pret2d_dep == 1 ~ "Pre-existing depression",
    postt2d_dep == 1 ~ "Incident depression",
    TRUE ~ "No depression"
  ))

# Summarise mean & 95% CI of n_unique_hu_dates by year & depression group
summary_df <- df %>%
  filter(hba1c_date_yrs <= 10) %>%   # limit to first 10 years
  group_by(hba1c_date_yrs, depression_status) %>%
  summarise(
    mean_val = mean(n_unique_hu_dates, na.rm = TRUE),
    se = sd(n_unique_hu_dates, na.rm = TRUE) / sqrt(n()),
    n = n()
  ) %>%
  mutate(
    ci_lower = mean_val - 1.96 * se,
    ci_upper = mean_val + 1.96 * se
  )

# Plot
ggplot(summary_df, aes(x = hba1c_date_yrs, y = mean_val, color = depression_status)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    title = "Mean number of unique HU dates (first 10 years)",
    x = "Years since T2D diagnosis",
    y = "Mean n_unique_hu_dates",
    color = "Depression Status"
  ) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


df <- hu_local_formatted_round


df <- df %>%
  mutate(depression_status = case_when(
    pret2d_dep == 1 ~ "Pre-existing depression",
    postt2d_dep == 1 ~ "Incident depression",
    TRUE ~ "No depression"
  ), 
  depression_status = factor(
    depression_status,
    levels = c("No depression", "Pre-existing depression", "Incident depression")
  )
)
# Bin by STOP -> whole years (completed years), cap to 1..10
df_binned <- df %>%
  mutate(
    year_bin = pmin(10, pmax(1, floor(stop)))  # use round()/ceiling() if desired
  ) %>%
  filter(year_bin <= 10, !is.na(n_unique_hu_dates))

# Ensure one record per person-year (if duplicates exist)
df_binned <- df_binned %>%
  distinct(patid, year_bin, depression_status, .keep_all = TRUE)

# Summarise across distinct patients in each year × group
summary_df <- df_binned %>%
  group_by(year_bin, depression_status) %>%
  summarise(
    mean_val = mean(n_unique_hu_dates, na.rm = TRUE),
    sd_val   = sd(n_unique_hu_dates, na.rm = TRUE),
    n_pat    = n_distinct(patid),
    se       = sd_val / sqrt(n_pat),
    ci_lower = mean_val - 1.96 * se,
    ci_upper = mean_val + 1.96 * se,
    .groups = "drop"
  )

# Plot (years 1–10 on x-axis), larger text
ggplot(summary_df, aes(x = year_bin, y = mean_val, color = depression_status)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    x = "Years since T2D diagnosis",
    y = "Mean unique healthcare dates",
    color = "Depression Status"
  ) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16)  # make all text larger
  )


inc_dx <- df %>%
  filter(incident_depression_timevar == 1) %>%
  group_by(patid) %>%
  summarise(first_inc_stop = min(stop, na.rm = TRUE), .groups = "drop")

patients <- df %>%
  distinct(patid) %>%
  left_join(inc_dx, by = "patid") %>%
  mutate(
    has_incident_dep = !is.na(first_inc_stop),
    # facet bin for those with incident depression (0–2, 2–4, 4–6, 6–8, 8–10)
    facet_bin = cut(
      first_inc_stop,
      breaks = c(0, 2, 4, 6, 8, 10),
      include.lowest = TRUE, right = TRUE,
      labels = c("0–2 years", "2–4 years", "4–6 years", "6–8 years", "8–10 years")
    )
  )

dat <- df %>%
  inner_join(patients, by = "patid") %>%
  mutate(
    # group for plotting (time-invariant membership)
    grp = if_else(has_incident_dep, "Incident depression", "No depression"),
    grp = factor(grp, levels = c("No depression", "Incident depression")),
    # whole-year bin from STOP; cap to 1..10
    year_bin = pmin(10, pmax(1, floor(stop)))
  ) %>%
  filter(year_bin <= 10, !is.na(n_unique_hu_dates))


no_dep <- dat %>%
  filter(grp == "No depression") %>%
  mutate(facet_bin = factor("0–2 years",
                            levels = c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years"))) %>%
  bind_rows(
    dat %>% filter(grp == "No depression") %>% mutate(facet_bin = factor("2–4 years",
                                                                         levels = c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years"))),
    dat %>% filter(grp == "No depression") %>% mutate(facet_bin = factor("4–6 years",
                                                                         levels = c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years"))),
    dat %>% filter(grp == "No depression") %>% mutate(facet_bin = factor("6–8 years",
                                                                         levels = c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years"))),
    dat %>% filter(grp == "No depression") %>% mutate(facet_bin = factor("8–10 years",
                                                                         levels = c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years")))
  )

dat_faceted <- bind_rows(
  dat %>% filter(grp == "Incident depression", !is.na(facet_bin)),
  no_dep
)

summ <- dat_faceted %>%
  group_by(facet_bin, year_bin, grp) %>%
  summarise(
    mean_val = mean(n_unique_hu_dates, na.rm = TRUE),
    sd_val   = sd(n_unique_hu_dates, na.rm = TRUE),
    n_pat    = n_distinct(patid),
    se       = sd_val / sqrt(n_pat),
    lo       = mean_val - 1.96 * se,
    hi       = mean_val + 1.96 * se,
    .groups = "drop"
  )

facet_counts <- patients %>%
  filter(has_incident_dep, !is.na(facet_bin)) %>%
  count(facet_bin, name = "n_inc") %>%
  mutate(label = paste0(facet_bin, " (n = ", scales::comma(n_inc), ")"))

vline_data <- data.frame(
  facet_bin = factor(
    rep(c("0–2 years","2–4 years","4–6 years","6–8 years","8–10 years"), each = 2),
    levels = levels(dat_faceted$facet_bin)
  ),
  x = c(0,2, 2,4, 4,6, 6,8, 8,10)
)

lab_map <- setNames(facet_counts$label, facet_counts$facet_bin)


ggplot(summ, aes(x = year_bin, y = mean_val, color = grp, group = grp)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_vline(data = vline_data, aes(xintercept = x), linetype = "dashed", inherit.aes = FALSE) +
  facet_wrap(~ facet_bin, ncol = 3, labeller = labeller(facet_bin = lab_map)) +
  scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
  scale_color_manual(values = c("No depression" = "black", "Incident depression" = "red")) +
  labs(
    x = "Years since T2D diagnosis",
    y = "Mean number of unique healthcare dates",
    color = "Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.position = "bottom"
  )


model_ttdep <- coxph(Surv(start, stop, incident_depression_timevar) ~ dm_diag_age + gender + n_unique_hu_dates_lag + ethnicity_bin + year_of_diagnosis + count_lag1 + preexisting_comorb_score + incident_comorbidity_score_timevar, data = hu_local_formatted_round %>% filter(pret2d_dep != 1))
mat <- model.matrix(model_ttdep)
X <- mat[, colnames(mat) != "(Intercept)", drop = FALSE]

model_ttdep_cat <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + n_unique_hu_dates_lag + ethnicity_bin + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar, data = hu_local_formatted_round %>% filter(pret2d_dep != 1))
model_ttdep_cat_unadj <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + year_of_diagnosis + count_cat_lag1 + incident_depression_timevar, data = hu_local_formatted_round_short %>% filter(pret2d_dep != 1))

model_ttdep_cat_pre <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + n_unique_hu_dates_lag + ethnicity_bin + year_of_diagnosis + count_cat_lag1 + pret2d_dep, data = hu_local_formatted_round)
model_ttdep_cat_unadj_pre <- coxph(Surv(start, stop, event) ~ dm_diag_age + gender + ethnicity_bin + year_of_diagnosis + count_cat_lag1 + pret2d_dep, data = hu_local_formatted_round)

tidy_unadj <- tidy(model_ttdep_cat_unadj, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(str_starts(term, "count_cat_lag1")) %>%
  mutate(
    level = str_remove(term, "^count_cat_lag1"),
    level = recode(level, "1-2" = "1–2", "3-4" = "3–4", "5+" = "5+"),
    model = "Unadjusted", 
    depression = "Incident depression"
  )

tidy_adj <- tidy(model_ttdep_cat, conf.int = TRUE, exponentiate = TRUE) 

tidy_unadj_pre <- tidy(model_ttdep_cat_unadj_pre, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(str_starts(term, "count_cat_lag1")) %>%
  mutate(
    level = str_remove(term, "^count_cat_lag1"),
    level = recode(level, "1-2" = "1–2", "3-4" = "3–4", "5+" = "5+"),
    model = "Unadjusted", 
    depression = "Pre-existing depression"
  )

tidy_adj_pre <- tidy(model_ttdep_cat_pre, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(str_starts(term, "count_cat_lag1")) %>%
  mutate(
    level = str_remove(term, "^count_cat_lag1"),
    level = recode(level, "1-2" = "1–2", "3-4" = "3–4", "5+" = "5+"),
    model = "Adjusted", 
    depression = "Pre-existing depression"
  )




df <- bind_rows(tidy_unadj, tidy_adj, tidy_unadj_pre, tidy_adj_pre) %>%
  mutate(level = factor(level, levels = c("1–2", "3–4", "5+")),
         model = factor(model, levels = c("Unadjusted", "Adjusted")), 
         depression = factor(depression, levels = c("Pre-existing depression", "Incident depression")))



# x positions between categories (don't draw at the very ends)
sep_x <- seq_len(nlevels(df$level) - 1) + 0.5

ggplot(
  df,
  aes(
    x    = level,
    y    = estimate,
    ymin = conf.low,
    ymax = conf.high,
    color = model
  )
) +
  geom_hline(yintercept = 1, linetype = "solid") +
  
  geom_errorbar(
    position = position_dodge(width = dodge_width),
    width = 0.3, linewidth = 0.7
  ) +
  geom_point(
    position = position_dodge(width = dodge_width),
    shape = 15, size = 2.5
  ) +
  
  # full-width dashed separators (vertical before flip -> horizontal after flip)
  geom_vline(
    xintercept = sep_x,
    linetype = "dashed",
    color = "gray70",
    linewidth = 0.3
  ) +
  
  coord_flip() +
  scale_y_log10(
    name   = "Hazard Ratio (95% CI)",
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2),
    limits = c(y_low, y_high)
  ) +
  scale_color_manual(values = c("Unadjusted" = "#E41A1C", "Adjusted" = "#377EB8")) +
  labs(x = "HbA1c testing frequency category", color = "Model") +
  facet_wrap(~ depression, ncol = 1, strip.position = "top") +
  theme_bw(base_size = 16) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 16, face = "plain", hjust = 0.5),
    strip.placement    = "outside"
  )
