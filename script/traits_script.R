# Lichen Functional Traits Analysis

## Algal layer analysis ## ----
## Jiri Subrt ##
## 25th of September 2023 ##

# Library
library(tidyverse)
library(rstatix)
library(plotrix)
library(ggpubr)
library(car) 
library(patchwork)

# Colours for figures
colours <- c("dark blue", "red")

# Specific algal area ----
# Load data
algal_layer <- read_csv("data/traits/algal_layer.csv")

# Select columns I want and summarise the sum of algal layer by treatment type
algal_layer_sum <- algal_layer %>% 
  select(-mean, -min,-max, -number) %>% 
  group_by(treatment_type, sample) %>% 
  filter(type == "algae") %>% 
  summarise(algal_sum = sum(area)) %>% 
  mutate(treatment_type = case_when(grepl("control", treatment_type) ~ "Control",
                                    grepl("treatment", treatment_type) ~ "Treatment")) %>% 
  ungroup()

# Select only whole lichen surface area
lichen_whole <- algal_layer %>% 
  select(-mean, -min,-max, -number) %>% 
  filter(type == "whole_lichen")

# Merge algal layer and whole surface area
algal_layer_summary <- merge(algal_layer_sum, lichen_whole, by = "sample") %>% 
  rename(treatment_type = "treatment_type.x")

# Calculate specific algal area
algal_area <- algal_layer_summary %>% 
  select(-treatment_type.y, -type) %>% 
  group_by(treatment_type, sample) %>% 
  summarise(specific_algal_area = algal_sum/area*100) %>% 
  ungroup()

# Stats
# Compute Shapiro wilk test by groups - normal distribution for both groups
algal_area %>%
  group_by(treatment_type) %>%
  shapiro_test(specific_algal_area)

# Draw a qq plot by group - Looks good
ggqqplot(algal_area, x = "specific_algal_area", facet.by = "treatment_type")

# equality of variances - Looks good
leveneTest(specific_algal_area ~ treatment_type, algal_area)

# t-test - Significant differences
saa_t_test <- algal_area %>%
  t_test(specific_algal_area ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# Get summary stats
average_saa_per_treatment <- algal_area %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(specific_algal_area),
            SE = std.error(specific_algal_area))

# Plot
algal_boxplot_mean <- ggplot(algal_area, aes(x = treatment_type, y = specific_algal_area, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(a)") +
  ylab("Specific algal area (%)") +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

### chlorophyll concentration ### ----
chlorophyll <- read_csv("data/traits/chlorophyll.csv")

# Stats
# Identify outliers 
chl_outliers <- chlorophyll %>% 
  group_by(treatment_type) %>% 
  identify_outliers("chl_dw") # Will remove those

chlorophyll_no_outlier <- chlorophyll %>% 
  anti_join(chl_outliers) 

# Stats
# Compute Shapiro wilk test by groups - normal distribution for both groups
chlorophyll_no_outlier %>%
  group_by(treatment_type) %>%
  shapiro_test(chl_dw)

# Draw a qq plot by group - Looks good
ggqqplot(chlorophyll_no_outlier, x = "chl_dw", facet.by = "treatment_type")

# equality of variances - Looks good
chlorophyll_no_outlier %>% levene_test(chl_dw ~ treatment_type)

# t-test - Significant differences!
chl_t_test <- chlorophyll_no_outlier %>%
  t_test(chl_dw ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# Get summary stats
average_chlorophyll_per_treatment <- chlorophyll_no_outlier %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(chl_dw),
            SE = std.error(chl_dw))

# Plots
chlorophyll$treatment_type <- as.factor(chlorophyll$treatment_type)

chlorophyll_boxplot_mean <- ggplot(chlorophyll_no_outlier, aes(x = treatment_type, y = chl_dw, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(b)") +
  ylab(label = expression(paste("Chlorophyll content "("mg g"^"−1")))) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

## Water holding capacity ----
whc <- read_csv("data/traits/STM_WHC.csv")

# Summary stats
whc %>%
  group_by(treatment_type) %>%
  get_summary_stats(WHC, type = "mean_sd")

# Identify outliers - Will remove those
whc_outliers <- whc %>% 
  group_by(treatment_type) %>% 
  identify_outliers("WHC")

whc_no_outlier <- whc %>% 
  anti_join(whc_outliers) %>% 
  mutate(treatment_type = case_when(grepl("control", treatment_type) ~ "Control",
                                    grepl("treatment", treatment_type) ~ "Treatment")) 

# Stats
# Compute Shapiro wilk test by groups - normal distribution for both groups
whc_no_outlier %>%
  group_by(treatment_type) %>%
  shapiro_test(WHC)

# Draw a qq plot by group - Looks good
ggqqplot(whc, x = "WHC", facet.by = "treatment_type")

# equality of variances - Looks good
whc_no_outlier %>% levene_test(WHC ~ treatment_type)

# t-test - not sig
whc_test <- whc_no_outlier %>%
  t_test(WHC ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# Get summary stats
average_chlorophyll_per_treatment <- whc_no_outlier %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(WHC),
            SE = std.error(WHC))

# Plots
whc_boxplot_mean <- ggplot(whc_no_outlier, aes(x = treatment_type, y = WHC, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(c)") +
  ylab(label = expression(paste("Water holding capacity "("g cm"^"−2")))) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))


# Specific Thallus Mass ----
stm <- read_csv("data/traits/STM_WHC.csv") %>% 
  mutate(treatment_type = case_when(grepl("control", treatment_type) ~ "Control",
                                    grepl("treatment", treatment_type) ~ "Treatment")) 

# Stats
# Compute Shapiro wilk test by groups - normal distribution for both groups
stm %>%
  group_by(treatment_type) %>%
  shapiro_test(STM) 

# t-test - not different
stm_test <- stm %>%
  t_test(STM ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# Get summary stats
average_stm_per_treatment <- stm %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(WHC),
            SE = std.error(WHC))

# plots
stm_boxplot_mean <- ggplot(stm, aes(x = treatment_type, y = STM, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(d)") +
  ylab(label = expression(paste("Specific thallus mass "("g cm"^"−2")))) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))


# Summary figure
figure_6 <- {
  algal_boxplot_mean + chlorophyll_boxplot_mean + whc_boxplot_mean + stm_boxplot_mean + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 2, heights = c(1, 1))

ggsave("outcomes/traits/figure_6.png", plot = figure_6, width = 10, height = 10, dpi = 300)


