# Temperature curves models - Non linear models using rTPC
# Fitting quadratic rTPC curve for temperature curves
# Adapted from rTPC vignette pipeline
# 26/01/2024
# Jiri Subrt

# Need to separate my data to four independent data frames
# np_control, np_treatment, dr_control, dr_treatment

# Pipeline
# 1. Model selection:
# I am using models from rTPC that are able to fit data points that are negative
# ar both before and after Topt. Then I compare their AIC values and choose the best one.
# The model compared are:
# joehnk_2008(), kamykowski 1985(), lactin2_1995()
# quadratic_2008(), thomas_2021, thomas 2017

# 2. Fit best model to my data and extract predictions
# 3. Bootstrap model uncertainty and extract paramaters with their CIs
# 4. Create final figure with model predictions and comparison of the parameters 

# Libraries
library(nls.multstart)
library(rTPC)
library(broom)
library(tidyverse)
library(ggrepel)
library(boot)
library(car)
library(patchwork)

# Functions
get_value_at_temperature <- function(model, target_temperature = 5) {
  
  # Capture environment from model - contains data
  x <- model$m$getEnv()
  
  # Get the name of the temperature column
  formula <- stats::as.formula(model$m$formula())
  param_ind <- all.vars(formula[[3]])[!all.vars(formula[[3]]) %in% names(model$m$getPars())]
  
  # Extract the temperature values
  vals <- x[[param_ind]]
  
  # Create a data frame with the target temperature
  newdata <- data.frame(x = target_temperature, stringsAsFactors = FALSE)
  
  # Rename to be the correct column name
  names(newdata) <- param_ind
  
  # Predict at the specified temperature
  newdata$preds <- stats::predict(model, newdata = newdata)
  
  # Extract the value at the specified temperature
  value_at_temperature = newdata$preds
  
  return(value_at_temperature)
}

calc_temp_steps <- function(model){
  t <- data.frame(zero = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 0), error = function(err) NA)),
                  five = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 5), error = function(err) NA)),
                  ten = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 10), error = function(err) NA)),
                  fifteen = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 15), error = function(err) NA)),
                  twenty = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 20), error = function(err) NA)),
                  twentyfive = suppressWarnings(tryCatch(get_value_at_temperature(model, target_temperature = 25), error = function(err) NA)),
                  stringsAsFactors = FALSE)
  
  return(t)
}


# MODEL SELECTION ----
# Load in the DR NP GP dataset with calculated max NP and max DR values
np_gr_dp_raw_long <- read.csv("data/temperature_response/np_dr_gp_raw_long_zkouska.csv")

# np_control model selection ----
np_control_data <-  np_gr_dp_raw_long %>% 
  filter(type == "NP") %>% 
  filter(treatment_type == "control") %>% 
  rename(temp = "temperature") %>% 
  rename(rate = "np_DW") %>% 
  rename(curve_id = "sample") %>% 
  select(-X)

# 1. Model selection
# fit five chosen models that were suitable from rTPC
np_control_fits <- nest(np_control_data, data = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                                       data = np_control_data,
                                                       iter = c(4,4,4),
                                                       start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008') - 10,
                                                       start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008') + 10,
                                                       lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                                                       upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                                                       supp_errors = 'Y',
                                                       convergence_count = FALSE)),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = np_control_data,
                                               iter = c(3,3,3,3,3),
                                               start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                           data = np_control_data,
                                           iter = c(3,3,3,3),
                                           start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'lactin2_1995') - 10,
                                           start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'lactin2_1995') + 10,
                                           lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'lactin2_1995'),
                                           upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'lactin2_1995'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         thomas17 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a, b, c, d, e),
                                             data = np_control_data,
                                             iter = c(3,3,3,3,3),
                                             start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2017') - 10,
                                             start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2017') + 10,
                                             lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2017'),
                                             upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2017'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         thomas12 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a, b, c, tref),
                                             data = np_control_data,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2012') - 1,
                                             start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2012') + 2,
                                             lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2012'),
                                             upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'thomas_2012'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)))

# stack models
np_control_stack <- select(np_control_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', quadratic:thomas12)

# get predictions using augment
np_control_new <- tibble(temp = seq(min(np_control_data$temp), max(np_control_data$temp), length.out = 100))
np_control_preds <- np_control_stack %>%
  mutate(., preds = map(fit, augment, newdata = np_control_new)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
np_control_labs <- filter(np_control_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(np_control_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', np_control_labs) +
  geom_point(aes(temp, rate), np_control_data) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'rate') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

np_control_aic <- np_control_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

np_control_aic # thomas12 best here, although it does not capture the curve, choosing quadratic

# filter for best model
np_best_model = filter(np_control_aic, AICc == min(AICc)) %>% pull(model_name)

# np_treatment model selection ----
np_treatment_data <-  np_gr_dp_raw_long %>% 
  filter(type == "NP") %>% 
  filter(treatment_type == "treatment") %>% 
  rename(temp = "temperature") %>% 
  rename(rate = "np_DW") %>% 
  rename(curve_id = "sample") %>% 
  select(-X)

# show the data
ggplot(np_treatment_data, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# 1. Model selection
# fit five chosen models that were suitable from rTPC
np_treatment_fits <- nest(np_treatment_data, data = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = np_treatment_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = np_treatment_data,
                                               iter = c(3,3,3,3,3),
                                               start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                           data = np_treatment_data,
                                           iter = c(3,3,3,3),
                                           start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'lactin2_1995') - 10,
                                           start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'lactin2_1995') + 10,
                                           lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'lactin2_1995'),
                                           upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'lactin2_1995'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         thomas17 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a, b, c, d, e),
                                             data = np_treatment_data,
                                             iter = c(3,3,3,3,3),
                                             start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2017') - 10,
                                             start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2017') + 10,
                                             lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2017'),
                                             upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2017'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         thomas12 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a, b, c, tref),
                                             data = np_treatment_data,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2012') - 1,
                                             start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2012') + 2,
                                             lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2012'),
                                             upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'thomas_2012'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)))

# stack models
np_treatment_stack <- select(np_treatment_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', quadratic:thomas12)

# get predictions using augment
np_treatment_new <- tibble(temp = seq(min(np_treatment_data$temp), max(np_treatment_data$temp), length.out = 100))
np_treatment_preds <- np_treatment_stack %>%
  mutate(., preds = map(fit, augment, newdata = np_treatment_new)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
np_treatment_labs <- filter(np_treatment_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(np_treatment_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', np_treatment_labs) +
  geom_point(aes(temp, rate), np_treatment_data) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

np_treatment_aic <- np_treatment_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

np_treatment_aic # quadratic here

# filter for best model
np_best_model = filter(np_treatment_aic, AICc == min(AICc)) %>% pull(model_name)

# dr_control model selection ----
dr_control_data <-  np_gr_dp_raw_long %>% 
  filter(type == "DR") %>% 
  filter(treatment_type == "control") %>% 
  rename(temp = "temperature") %>% 
  rename(rate = "np_DW") %>% 
  rename(curve_id = "sample") %>% 
  select(-X) %>% 
  mutate(rate = abs(rate))

# show the data
ggplot(dr_control_data, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# 1. Model selection
# fit five chosen models that were suitable from rTPC
dr_control_fits <- nest(dr_control_data, data = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = dr_control_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = dr_control_data,
                                               iter = c(3,3,3,3,3),
                                               start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                           data = dr_control_data,
                                           iter = c(3,3,3,3),
                                           start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'lactin2_1995') - 10,
                                           start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'lactin2_1995') + 10,
                                           lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'lactin2_1995'),
                                           upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'lactin2_1995'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         thomas17 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a, b, c, d, e),
                                             data = dr_control_data,
                                             iter = c(3,3,3,3,3),
                                             start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2017') - 10,
                                             start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2017') + 10,
                                             lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2017'),
                                             upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2017'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         thomas12 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a, b, c, tref),
                                             data = dr_control_data,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2012') - 1,
                                             start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2012') + 2,
                                             lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2012'),
                                             upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'thomas_2012'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)))

# stack models
dr_control_stack <- select(dr_control_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', quadratic:thomas12)

# get predictions using augment
dr_control_new <- tibble(temp = seq(min(dr_control_data$temp), max(dr_control_data$temp), length.out = 100))
dr_control_preds <- dr_control_stack %>%
  mutate(., preds = map(fit, augment, newdata = dr_control_new)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
dr_control_labs <- filter(dr_control_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(dr_control_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', dr_control_labs) +
  geom_point(aes(temp, rate), dr_control_data) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

dr_control_aic <- dr_control_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

dr_control_aic # quadratic here

# filter for best model
dr_best_model = filter(dr_control_aic, AICc == min(AICc)) %>% pull(model_name)

# dr_treatment model selection ----
dr_treatment_data <-  np_gr_dp_raw_long %>% 
  filter(type == "DR") %>% 
  filter(treatment_type == "treatment") %>% 
  rename(temp = "temperature") %>% 
  rename(rate = "np_DW") %>% 
  rename(curve_id = "sample") %>% 
  select(-X) %>% 
  mutate(rate = abs(rate))

# show the data
ggplot(dr_treatment_data, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# 1. Model selection
# fit five chosen models that were suitable from rTPC
dr_treatment_fits <- nest(dr_treatment_data, data = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = dr_treatment_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = dr_treatment_data,
                                               iter = c(3,3,3,3,3),
                                               start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                           data = dr_treatment_data,
                                           iter = c(3,3,3,3),
                                           start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'lactin2_1995') - 10,
                                           start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'lactin2_1995') + 10,
                                           lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'lactin2_1995'),
                                           upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'lactin2_1995'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         thomas17 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a, b, c, d, e),
                                             data = dr_treatment_data,
                                             iter = c(3,3,3,3,3),
                                             start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2017') - 10,
                                             start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2017') + 10,
                                             lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2017'),
                                             upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2017'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         thomas12 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a, b, c, tref),
                                             data = dr_treatment_data,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2012') - 1,
                                             start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2012') + 2,
                                             lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2012'),
                                             upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'thomas_2012'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)))

# stack models
dr_treatment_stack <- select(dr_treatment_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', quadratic:thomas12)

# get predictions using augment
dr_treatment_new <- tibble(temp = seq(min(dr_treatment_data$temp), max(dr_treatment_data$temp), length.out = 100))
dr_treatment_preds <- dr_treatment_stack %>%
  mutate(., preds = map(fit, augment, newdata = dr_treatment_new)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
dr_treatment_labs <- filter(dr_treatment_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(dr_treatment_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', dr_treatment_labs) +
  geom_point(aes(temp, rate), dr_treatment_data) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

dr_treatment_aic <- dr_treatment_stack %>%
  mutate(., info = map(fit, glance),
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, sigma, AIC, AICc, BIC, df.residual)

dr_treatment_aic # quadratic here

# filter for best model
dr_best_model = filter(dr_treatment_aic, AICc == min(AICc)) %>% pull(model_name)

# BOOTSTRAPPING 
# np_control bootstrapping ----
np_control_fit <- nest(np_control_data, data  = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = np_control_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))


# unnest predictions
np_control_fit_preds <- select(np_control_fit, preds, curve_id, treatment_type) %>%
  unnest(preds)

# plot data and predictions
ggplot() +
  geom_line(aes(temp, .fitted), np_control_fit_preds, col = 'blue') +
  geom_point(aes(temp, rate), np_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'NP',
       title = 'NP across temperatures')

# refit using nlsLM
fit_nlsLM_np_control <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                               data = np_control_data,
                               start = coef(np_control_fit$quadratic[[1]]),
                               lower = get_lower_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(np_control_data$temp, np_control_data$rate, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(np_control_data)))

# bootstrap using case resampling
boot1_np_control <- Boot(fit_nlsLM_np_control, method = 'case')

hist(boot1_np_control, layout = c(2,2))

# create predictions of each bootstrapped model
boot1_preds_np_control <- boot1_np_control$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(np_control_data$temp), max(np_control_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_np_control <- group_by(boot1_preds_np_control, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>% 
  mutate(treatment_type = "control")
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), np_control_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_np_control, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), np_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), np_control_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot1_preds_np_control, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), np_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# bootstrap using residual resampling
boot_resid_np_control <- Boot(fit_nlsLM_np_control, method = 'residual')

# predict over new data
boot_resid_np_control_preds <- boot_resid_np_control$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(np_control_data$temp), max(np_control_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot_resid_np_control_conf_preds <- group_by(boot_resid_np_control_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), np_control_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_resid_np_control_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), np_control_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), np_control_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot_resid_np_control_preds, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), np_control_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

p1 + p2

# np treatment bootstrapping ----
np_treatment_fit <- nest(np_treatment_data, data  = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = np_treatment_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 1000))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))


# unnest predictions
np_treatment_fit_preds <- select(np_treatment_fit, preds, curve_id, treatment_type) %>%
  unnest(preds)

# plot data and predictions
ggplot() +
  geom_line(aes(temp, .fitted), np_treatment_fit_preds, col = 'blue') +
  geom_point(aes(temp, rate), np_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'NP',
       title = 'NP across temperatures')

# refit using nlsLM
fit_nlsLM_np_treatment <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                                          data = np_treatment_data,
                                          start = coef(np_treatment_fit$quadratic[[1]]),
                                          lower = get_lower_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                          upper = get_upper_lims(np_treatment_data$temp, np_treatment_data$rate, model_name = 'quadratic_2008'),
                                          weights = rep(1, times = nrow(np_treatment_data)))

# bootstrap using case resampling
boot1_np_treatment <- Boot(fit_nlsLM_np_treatment, method = 'case')

hist(boot1_np_treatment, layout = c(2,2))

# create predictions of each bootstrapped model
boot1_preds_np_treatment <- boot1_np_treatment$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(np_treatment_data$temp), max(np_treatment_data$temp), length.out = 1000))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_np_treatment <- group_by(boot1_preds_np_treatment, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  mutate(treatment_type = "treatment") %>% 
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), np_treatment_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_np_treatment, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), np_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), np_treatment_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot1_preds_np_treatment, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), np_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# bootstrap using residual resampling
boot_resid_np_treatment <- Boot(fit_nlsLM_np_treatment, method = 'residual')

# predict over new data
boot_resid_np_treatment_preds <- boot_resid_np_treatment$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(np_treatment_data$temp), max(np_treatment_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot_resid_np_treatment_conf_preds <- group_by(boot_resid_np_treatment_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), np_treatment_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_resid_np_treatment_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), np_treatment_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), np_treatment_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot_resid_np_treatment_preds, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), np_treatment_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

p1 + p2

# dr control bootstrapping ----
dr_treatment_fit <- nest(dr_treatment_data, data  = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = dr_treatment_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))


# unnest predictions
dr_treatment_fit_preds <- select(dr_treatment_fit, preds, curve_id) %>%
  mutate(treatment_type = "treatment") %>% 
  unnest(preds)

# plot data and predictions
ggplot() +
  geom_line(aes(temp, .fitted), dr_treatment_fit_preds, col = 'blue') +
  geom_point(aes(temp, rate), dr_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'dr',
       title = 'dr across temperatures')

# refit using nlsLM
fit_nlsLM_dr_treatment <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                                            data = dr_treatment_data,
                                            start = coef(dr_treatment_fit$quadratic[[1]]),
                                            lower = get_lower_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                            upper = get_upper_lims(dr_treatment_data$temp, dr_treatment_data$rate, model_name = 'quadratic_2008'),
                                            weights = rep(1, times = nrow(dr_treatment_data)))

# bootstrap using case resampling
boot1_dr_treatment <- Boot(fit_nlsLM_dr_treatment, method = 'case')

hist(boot1_dr_treatment, layout = c(2,2))

# create predictions of each bootstrapped model
boot1_preds_dr_treatment <- boot1_dr_treatment$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(dr_treatment_data$temp), max(dr_treatment_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_dr_treatment <- group_by(boot1_preds_dr_treatment, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  mutate(treatment_type = "treatment") %>% 
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_treatment_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_dr_treatment, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), dr_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_treatment_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot1_preds_dr_treatment, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), dr_treatment_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# bootstrap using residual resampling
boot_resid_dr_treatment <- Boot(fit_nlsLM_dr_treatment, method = 'residual')

# predict over new data
boot_resid_dr_treatment_preds <- boot_resid_dr_treatment$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(dr_treatment_data$temp), max(dr_treatment_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot_resid_dr_treatment_conf_preds <- group_by(boot_resid_dr_treatment_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_treatment_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_resid_dr_treatment_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), dr_treatment_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_treatment_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot_resid_dr_treatment_preds, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), dr_treatment_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

p1 + p2

# dr control bootstrapping ----
dr_control_fit <- nest(dr_control_data, data  = c(temp, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = dr_control_data,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))


# unnest predictions
dr_control_fit_preds <- select(dr_control_fit, preds, curve_id) %>%
  mutate(treatment_type = "control") %>% 
  unnest(preds)

# plot data and predictions
ggplot() +
  geom_line(aes(temp, .fitted), dr_control_fit_preds, col = 'blue') +
  geom_point(aes(temp, rate), dr_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'dr',
       title = 'dr across temperatures')

# refit using nlsLM
fit_nlsLM_dr_control <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                                            data = dr_control_data,
                                            start = coef(dr_control_fit$quadratic[[1]]),
                                            lower = get_lower_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                            upper = get_upper_lims(dr_control_data$temp, dr_control_data$rate, model_name = 'quadratic_2008'),
                                            weights = rep(1, times = nrow(dr_control_data)))

# bootstrap using case resampling
boot1_dr_control <- Boot(fit_nlsLM_dr_control, method = 'case')

hist(boot1_dr_control, layout = c(2,2))

# create predictions of each bootstrapped model
boot1_preds_dr_control <- boot1_dr_control$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(dr_control_data$temp), max(dr_control_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_dr_control <- group_by(boot1_preds_dr_control, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  mutate(treatment_type = "control") %>% 
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_control_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_dr_control, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), dr_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_control_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot1_preds_dr_control, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), dr_control_data, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# bootstrap using residual resampling
boot_resid_dr_control <- Boot(fit_nlsLM_dr_control, method = 'residual')

# predict over new data
boot_resid_dr_control_preds <- boot_resid_dr_control$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(dr_control_data$temp), max(dr_control_data$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = temp, a, b, c))

# calculate bootstrapped confidence intervals
boot_resid_dr_control_conf_preds <- group_by(boot_resid_dr_control_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_control_fit_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot_resid_dr_control_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), dr_control_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), dr_control_fit_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot_resid_dr_control_preds, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), dr_control_data, size = 2) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

p1 + p2

# CALCULATING CIs OF PARAMETERS using different methods ----
# Chose residual method afterall
# np_control_calculation ----
# get parameters of fitted model
param_np_control <- broom::tidy(fit_nlsLM_np_control) %>%
  select(param = term, estimate)

# calculate confidence intervals of models
ci_np_control1 <- nlstools::confint2(fit_nlsLM_np_control, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

ci_np_control2 <- confint(fit_nlsLM_np_control) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')

# CIs from case resampling
ci_np_control3 <- confint(boot1_np_control, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# CIs from residual resampling
ci_np_control4 <- Boot(fit_nlsLM_np_control, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_np_control1, ci_np_control2, ci_np_control3, ci_np_control4) %>%
  left_join(., param_np_control)

# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'dashed lines are CI of profiling method')

extra_params_np_control <- calc_params(fit_nlsLM_np_control) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

extra_temp_params_np_control <- calc_temp_steps(fit_nlsLM_np_control) %>% 
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params_np_control <- Boot(fit_nlsLM_np_control, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_np_control)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_temp_params_np_control <- Boot(fit_nlsLM_np_control, f = function(x){unlist(calc_temp_steps(x))}, labels = names(calc_temp_steps(fit_nlsLM_np_control)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_np_control <- left_join(ci_extra_params_np_control, extra_params_np_control) %>% 
  mutate(treatment_type = "control")

ci_extra_temp_params_np_control <- left_join(ci_extra_temp_params_np_control, extra_temp_params_np_control) %>% 
  mutate(treatment_type = "control")

ggplot(ci_extra_temp_params_np_control, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'using case resampling')

# np_treatment_calculation ----
# get parameters of fitted model
param_np_treatment <- broom::tidy(fit_nlsLM_np_treatment) %>%
  select(param = term, estimate)

# calculate confidence intervals of models
ci_np_treatment1 <- nlstools::confint2(fit_nlsLM_np_treatment, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

ci_np_treatment2 <- confint(fit_nlsLM_np_treatment) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')

# CIs from case resampling
ci_np_treatment3 <- confint(boot1_np_treatment, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# CIs from residual resampling
ci_np_treatment4 <- Boot(fit_nlsLM_np_treatment, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_np_treatment1, ci_np_treatment2, ci_np_treatment3, ci_np_treatment4) %>%
  left_join(., param_np_treatment)

# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'dashed lines are CI of profiling method')

# Normal params
extra_params_np_treatment <- calc_params(fit_nlsLM_np_treatment) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

# temp params
extra_temp_params_np_treatment <- calc_temp_steps(fit_nlsLM_np_control) %>% 
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params_np_treatment <- Boot(fit_nlsLM_np_treatment, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM_np_treatment)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_temp_params_np_treatment <- Boot(fit_nlsLM_np_control, f = function(x){unlist(calc_temp_steps(x))}, labels = names(calc_temp_steps(fit_nlsLM_np_control)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

ci_extra_params_np_treatment <- left_join(ci_extra_params_np_treatment, extra_params_np_treatment) %>% 
  mutate(treatment_type = "treatment")

ci_extra_temp_params_np_treatment <- left_join(ci_extra_temp_params_np_treatment, extra_temp_params_np_treatment) %>% 
  mutate(treatment_type = "treatment")

ggplot(ci_extra_params_np_treatment, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'using case resampling')

# dr_treatment_calculation ----
param_dr_treatment <- broom::tidy(fit_nlsLM_dr_treatment) %>%
  select(param = term, estimate)

# calculate confidence intervals of models
ci_dr_treatment1 <- nlstools::confint2(fit_nlsLM_dr_treatment, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

ci_dr_treatment2 <- confint(fit_nlsLM_dr_treatment) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')

# CIs from case resampling
ci_dr_treatment3 <- confint(boot1_dr_treatment, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# CIs from residual resampling
ci_dr_treatment4 <- Boot(fit_nlsLM_dr_treatment, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_dr_treatment1, ci_dr_treatment2, ci_dr_treatment3, ci_dr_treatment4) %>%
  left_join(., param_dr_treatment)

# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'dashed lines are CI of profiling method')

extra_params_dr_treatment <- calc_temp_steps(fit_nlsLM_dr_treatment) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate') 

ci_extra_params_dr_treatment <- Boot(fit_nlsLM_dr_treatment, f = function(x){unlist(calc_temp_steps(x))}, labels = names(calc_temp_steps(fit_nlsLM_dr_treatment)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap') 

ci_extra_params_dr_treatment_connected <- left_join(ci_extra_params_dr_treatment, extra_params_dr_treatment) %>% 
  mutate(treatment_type = "treatment")

ggplot(ci_extra_params_dr_treatment_connected, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'using case resampling')

# dr_control_calculation ----
param_dr_control <- broom::tidy(fit_nlsLM_dr_control) %>%
  select(param = term, estimate)

# calculate confidence intervals of models
ci_dr_control1 <- nlstools::confint2(fit_nlsLM_dr_control, method = 'asymptotic') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'asymptotic')

ci_dr_control2 <- confint(fit_nlsLM_dr_control) %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'profile')

# CIs from case resampling
ci_dr_control3 <- confint(boot1_dr_control, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap')

# CIs from residual resampling
ci_dr_control4 <- Boot(fit_nlsLM_dr_control, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_bact <- bind_rows(ci_dr_control1, ci_dr_control2, ci_dr_control3, ci_dr_control4) %>%
  left_join(., param_dr_control)

# plot
ggplot(ci_bact, aes(forcats::fct_relevel(method, c('profile', 'asymptotic')), estimate, col = method)) +
  geom_hline(aes(yintercept = conf_lower), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_hline(aes(yintercept = conf_upper), linetype = 2, filter(ci_bact, method == 'profile')) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('', labels = function(x) stringr::str_wrap(x, width = 10)) +
  labs(title = 'Calculation of confidence intervals for model parameters',
       subtitle = 'dashed lines are CI of profiling method')

extra_params_dr_control <- calc_temp_steps(fit_nlsLM_dr_control) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params_dr_control <- Boot(fit_nlsLM_dr_control, f = function(x){unlist(calc_temp_steps(x))}, labels = names(calc_temp_steps(fit_nlsLM_dr_control)), R = 200, method = 'case') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'case bootstrap') 

ci_extra_params_dr_control_connected <- left_join(ci_extra_params_dr_control, extra_params_dr_control) %>% 
  mutate(treatment_type = "control")

ggplot(ci_extra_params_dr_control_connected, aes(param, estimate)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw() +
  facet_wrap(~param, scales = 'free') +
  scale_x_discrete('') +
  labs(title = 'Calculation of confidence intervals for extra parameters',
       subtitle = 'using case resampling')

# FIGURES ----
# connect parameters from np_control and np_treatment and plot them together
ci_extra_temp_params_np_treatment <- ci_extra_temp_params_np_treatment %>% 
  mutate(temperature = case_when(param == "zero" ~ 0,
                                 param == "five" ~ 5,
                                 param == "ten" ~ 10,
                                 param == "fifteen" ~ 15,
                                 param == "twenty" ~ 20,
                                 param == "twentyfive" ~ 25))

ci_extra_temp_params_np_control <-  ci_extra_temp_params_np_control %>% 
  mutate(temperature = case_when(param == "zero" ~ 0,
                                 param == "five" ~ 5,
                                 param == "ten" ~ 10,
                                 param == "fifteen" ~ 15,
                                 param == "twenty" ~ 20,
                                 param == "twentyfive" ~ 25))

extra_params_np <- rbind(ci_extra_temp_params_np_treatment, ci_extra_temp_params_np_control)

np_params <- rbind(ci_extra_params_np_control,ci_extra_params_np_treatment) %>% 
  mutate(treatment_type = case_when(grepl("control", treatment_type) ~ "Control",
                                    grepl("treatment", treatment_type) ~ "Treatment"))

(plot_topt <- filter(np_params, param %in% c('topt')) %>%
  ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
  geom_point(size = 3) + 
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('topt'))) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
    ylab(expression(T[opt]~"(ºC)")) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('(c)') +
  xlab(''))

(plot_rmax <- filter(np_params, param %in% c('rmax')) %>%
    ggplot(., aes(treatment_type, estimate/1000, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower/1000, ymax = conf_upper/1000, col = treatment_type), data = filter(np_params, param %in% c('rmax'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab(expression(NP[max]~(mu*mol~CO[2]~g^-1~s^-1))) +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_breadth <- filter(np_params, param %in% c('breadth')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('breadth'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab(expression(T[breadth]~"(ºC)")) +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(d)') +
    xlab(''))

(plot_ct_min <- filter(np_params, param %in% c('ctmin')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('ctmin'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab(expression(CT[min]~"(ºC)")) +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(e)') +
    xlab(''))

(plot_ct_max <- filter(np_params, param %in% c('ctmax')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('ctmax'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab(expression(CT[max]~"(ºC)")) +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(f)') +
    xlab(''))

(plot_e <- filter(np_params, param %in% c('e')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('e'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab('e') +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_eh <- filter(np_params, param %in% c('eh')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('eh'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab('eh') +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_q10 <- filter(np_params, param %in% c('q10')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('q10'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab('q10') +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_thermal_safety_margin <- filter(np_params, param %in% c('thermal_safety_margin')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('thermal_safety_margin'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab('thermal_safety_margin') +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_skewness <- filter(np_params, param %in% c('skewness')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('skewness'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab('skewness') +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(b)') +
    xlab(''))

(plot_thermal_tolerance <- filter(np_params, param %in% c('thermal_tolerance')) %>%
    ggplot(., aes(treatment_type, estimate, col = treatment_type)) + 
    geom_point(size = 3) + 
    geom_linerange(aes(ymin = conf_lower, ymax = conf_upper, col = treatment_type), data = filter(np_params, param %in% c('thermal_tolerance'))) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    ylab(expression(T[tol]~"(ºC)")) +
    scale_color_manual(values = c('dark blue', 'red')) +
    ggtitle('(g)') +
    xlab(''))

combined_figure <- patchwork::wrap_plots(plot_thermal_safety_margin,plot_thermal_tolerance, plot_skewness, plot_q10,plot_eh,plot_e,plot_ct_max,plot_ct_min, plot_breadth, plot_rmax, plot_topt)

# Combined figure for NP rates
combined_pred_np <- rbind(boot1_conf_preds_np_control, boot1_conf_preds_np_treatment)
combined_preds_fit_np <- rbind(np_control_fit_preds, np_treatment_fit_preds)
np_c <-  np_gr_dp_raw_long %>% 
  filter(type == "NP") 

(plot_bootstrap <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +  # Adding horizontal line at y = 0
    geom_point(aes(x = temperature, y = np_DW/1000, col = treatment_type), data = np_c, size = 3, alpha = 0.5) +
    geom_line(aes(x = temp, y = .fitted/1000, col = treatment_type, group = interaction(treatment_type)), data = combined_preds_fit_np) +
    geom_ribbon(aes(temp, ymin = conf_lower/1000, ymax = conf_upper/1000, fill = treatment_type), data = combined_pred_np, alpha = 0.2) +
    theme_bw(base_size = 14) +
    ylim(-1, 7) +
    scale_fill_manual(values = c('dark blue', 'red')) +
    scale_color_manual(values = c('dark blue', 'red')) +
    xlab('Temperature (ºC)') +
    ylab(label = expression(paste(
      "Net Photosynthesis µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")"))) + 
    ggtitle('(a)') +
    theme_bw(base_size = 15) +
    theme(legend.position = 'none'))
#geom_vline(aes(xintercept = 29.2), linetype = 2, col = 'black'))

np_figure <- plot_bootstrap + {
  plot_rmax + plot_topt + plot_breadth + plot_ct_min + plot_ct_max + plot_thermal_tolerance + plot_layout(ncol = 6)
} +
  plot_layout(nrow = 2, heights = c(0.7, 0.3))

## FIGURES DR ----
combined_pred_dr <- rbind(boot1_conf_preds_dr_control, boot1_conf_preds_dr_treatment)
combined_preds_fit_dr <- rbind(dr_control_fit_preds, dr_treatment_fit_preds)
dr_data <-  np_gr_dp_raw_long %>% 
  filter(type == "DR") %>% 
  mutate(np_DW = abs(np_DW))


ci_extra_params_dr_treatment_connected <- ci_extra_params_dr_treatment_connected %>% 
  mutate(temperature = case_when(param == "zero" ~ 0,
                                 param == "five" ~ 5,
                                 param == "ten" ~ 10,
                                 param == "fifteen" ~ 15,
                                 param == "twenty" ~ 20,
                                 param == "twentyfive" ~ 25))

ci_extra_params_dr_control_connected <-  ci_extra_params_dr_control_connected %>% 
  mutate(temperature = case_when(param == "zero" ~ 0,
                                 param == "five" ~ 5,
                                 param == "ten" ~ 10,
                                 param == "fifteen" ~ 15,
                                 param == "twenty" ~ 20,
                                 param == "twentyfive" ~ 25))

extra_params_dr <- rbind(ci_extra_params_dr_treatment_connected, ci_extra_params_dr_control_connected)

(plot_bootstrap_respiration <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.2) +  # Adding horizontal line at y = 0
    geom_point(aes(x = temperature, y = np_DW/1000, col = treatment_type, group = interaction(treatment_type)), data = dr_data, size = 3, alpha = 0.5) +
    geom_line(aes(temp, .fitted/1000, col = treatment_type, group = interaction(treatment_type)), data = combined_preds_fit_dr) +
    geom_ribbon(aes(temp, ymin = conf_lower/1000, ymax = conf_upper/1000, fill = treatment_type), data = combined_pred_dr, alpha = 0.2) +
    theme_bw(base_size = 14) +
    ylim(-1, 7) +
    scale_fill_manual(values = c('dark blue', 'red')) +
    scale_color_manual(values = c('dark blue', 'red')) +
    xlab('Temperature (ºC)') +
    ylab(label = expression(paste(
      "Dark Respiration (µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")"))) +
    theme_bw(base_size = 15) +
    ggtitle('(h)') +
    #ggtitle('(a) Thermal performance of respiration of Lichens in and out the OTCs') +
    theme(legend.position = 'none'))

### Both NP and DR together
figure_4 <- np_figure + {
  plot_bootstrap_respiration + plot_layout(ncol = 1)
} +
  plot_layout(nrow = 3, heights = c(20, 10))

ggsave("outcomes/temperature_curves/figure_4.png", plot = figure_4, width = 11, height = 11, dpi = 300)



np_figure <- plot_bootstrap + {
  plot_rmax + plot_topt + plot_breadth + plot_ct_min + plot_ct_max + plot_thermal_tolerance + plot_layout(ncol = 6)
} +
  plot_layout(nrow = 2, heights = c(0.7, 0.3))

dr_figure <- plot_bootstrap_respiration + {
  plot_rmax + plot_topt + plot_breadth + plot_ct_min + plot_ct_max + plot_thermal_tolerance + plot_layout(ncol = 6)
} +
  plot_layout(nrow = 2, heights = c(0.7, 0.3))

