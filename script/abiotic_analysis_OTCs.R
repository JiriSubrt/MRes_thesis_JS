### Abiotic measurements analysis from TOMST loggers and iButtons ###
### Jiri Subrt
### 31/08/23

# Libraries
library(tidyverse)
library(stringr)
library(plotrix)
library(lubridate)
library(ggblend)
library(rstatix)
library(ggpubr)
library(pwr)
library(rstatix)
library(coin)
library(patchwork)
library(nortest)
library(tools)
library(data.table)
library(plantecophys)

### TOMST loggers data ### ----

# Read in TMS data function (Adapted from Elise Gallois)
read_tms4 <- function(file) {      
  
  # Extract serial number from filename
  serial <- file 
  print(file)
  
  # Read the data file
  data <- read_delim(file, delim = ";",
                     col_names = F, 
                     locale=locale(decimal_mark = ",")) 
  
  
  # Check file has contents. Empty files due to bad data download only have "File is empty" as text. 
  if (ncol(data) > 1) {
    # Create vector of column names
    vars <- c("Index", "Datetime_UTC", "TimeZone", "soil_temp", "surface_temp", "top_temp", "moisture_soil", "shake",
              "errFlag", "empty")
    
    # Format data for output
    names(data) <- vars
    
    data_with_ID <- data  %>% 
      mutate(SerialID = serial) %>% 
      select(SerialID, everything()) %>% 
      mutate(Datetime_UTC = lubridate::parse_date_time(Datetime_UTC,orders = c("%Y.%m.%d %H:%M")))
    
  } else {
    print("empty file")
    data_with_ID <- NULL
  }
  
  
  return(data_with_ID)
}

## Read-in data files
tomst <- "data/microclimate/TOMST"
files <- list.files(path = tomst, pattern = "^data_*", full.names = T)
OTC_data <- map_dfr(files, read_tms4)

# Filter correct date, correct sensor, add control and treatment groups
tomst_otc_temp <-  OTC_data %>% 
  pivot_longer(cols = 5:8,
               names_to = "Variable",
               values_to = "Value") %>% 
  filter(format(Datetime_UTC, "%m-%d") != "08-30") %>% 
  filter(format(Datetime_UTC, "%m-%d") != "06-28") %>% 
  filter(Variable == "surface_temp") %>%  
  mutate(treatment_type = case_when(grepl("C", SerialID) ~ "Control",
                                    grepl("T", SerialID) ~ "Treatment"))

# Need to change hour to 1 hour plus, TOMSTS are for example 17:00, but that is UK time
tomst_otc$Datetime_UTC <- tomst_otc$Datetime_UTC + hours(1)

# TOMST TEMPERATURE summary stats
# create min max and mean values before doing summary stats
OTC_temp_stats <- tomst_otc_temp %>%
  mutate(date = as.Date(Datetime_UTC)) %>%
  group_by(date, treatment_type) %>%
  summarise(
    max_temperature = max(Value),
    min_temperature = min(Value),
    mean_temperature = mean(Value)
  )

# Summary stats
# Summary mean temp points from RAW points
# Control -  9.36 ± 0.0236
# Treatment - 10.6 ± 0.0292
tomst_otc_temp %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(Value),
            SE = std.error(Value))

# Summary max temperatures
#  Control n = 62 | mean daily max temp and SE = 12.8 0.442
#  Treatment  n = 62  | mean daily max temp and SE =  15.5 0.582
OTC_temp_stats %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(max_temperature),
            SE = std.error(max_temperature))

# Summary min temperatures
# Control  n = 62 | mean daily min temp and SE =  6.55  0.252
# Treatment n = 62 | mean daily min temp and SE =  7.07  0.288
OTC_temp_stats %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(min_temperature),
            SE = std.error(min_temperature))


# Normality - raw data - too many points - using Anderson-Darling normality test
normality_raw_temp_otc_control <- tomst_otc_temp %>% 
  filter(treatment_type == "Control") 

ad.test(normality_raw_temp_otc_control$Value) # Not normal

normality_raw_temp_otc_treatment <- tomst_otc_temp %>% 
  filter(treatment_type == "Treatment") 

ad.test(normality_raw_temp_otc_treatment$Value) # Not normal

# Normality - max temp - assumption not met
OTC_temp_stats %>%
  group_by(treatment_type) %>%
  shapiro_test(max_temperature)

# Normality - min temp - assumption not met
OTC_temp_stats %>%
  group_by(treatment_type) %>%
  shapiro_test(min_temperature)

# Draw a qq plot by group - does not look good
ggqqplot(tomst_otc_temp, x = "Value", facet.by = "treatment_type")
ggqqplot(OTC_temp_stats, x = "max_temperature", facet.by = "treatment_type")
ggqqplot(OTC_temp_stats, x = "min_temperature", facet.by = "treatment_type")

# equality of variances - assumption not met
tomst_otc_temp %>% levene_test(Value ~ treatment_type)

# Raw data test - Significant - W = 131962267, p-value < 2.2e-16
result_OTC_temp_raw <- wilcox.test(Value ~ treatment_type, data = tomst_otc_temp)

# Wilcoxon test - Significant - W = 1282, p-value = 0.001393
result_OTC_temp_daily_max <- wilcox.test(max_temperature ~ treatment_type, data = OTC_temp_stats)

# Wilcoxon test - Non-sig - W = 1658.5, p-value = 0.1886
result_OTC_temp_daily_min <- wilcox.test(min_temperature ~ treatment_type, data = OTC_temp_stats)

# TOMST temperature plots ----
my_colors <- c("Control" = "dark blue", "Treatment" = "red")

# Linegraphs
(tomst_mean_temp_line <- ggplot(OTC_temp_stats, aes(x = date, y = mean_temperature, group = treatment_type)) +
    geom_point(data = tomst_otc_temp, aes(x = as.Date(Datetime_UTC), y = Value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    # Move legend to the top right and adjust legend labels
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Mean Daily Temperature (˚C)"))

(tomst_max_temp_line <- ggplot(OTC_temp_stats, aes(x = date, y = max_temperature, group = treatment_type)) +
    geom_point(data = tomst_otc_temp, aes(x = as.Date(Datetime_UTC), y = Value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    # Move legend to the top right and adjust legend labels
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Average Max Daily Temperature (˚C)"))

(tomst_min_temp_line <- ggplot(OTC_temp_stats, aes(x = date, y = min_temperature, group = treatment_type)) +
    geom_point(data = tomst_otc_temp, aes(x = as.Date(Datetime_UTC), y = Value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    # Move legend to the top right and adjust legend labels
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size

    xlab("Date") +
    ylab("Average Max Daily Temperature (˚C)"))

# Mean boxplots
tomst_mean_temp_boxplot <- ggplot(OTC_temp_stats, aes(x = treatment_type, y = mean_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(3, 25)) +
  
  # Other styling options
  ylab("") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

tomst_max_temp_boxplot <- ggplot(OTC_temp_stats, aes(x = treatment_type, y = max_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(3, 25)) +
  
  # Other styling options
  ggtitle("") +
  ylab("") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

tomst_min_temp_boxplot <- ggplot(OTC_temp_stats, aes(x = treatment_type, y = min_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(3, 25)) +
  
  # Other styling options
  ggtitle("") +
  ylab("") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Combined TOMST temp plot - Supplementary Information S5
tomst_temp_plot <- {
  tomst_mean_temp_line + tomst_mean_temp_boxplot + tomst_max_temp_line + tomst_max_temp_boxplot + tomst_min_temp_line + tomst_min_temp_boxplot + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))


### iButtons data ### ----
# Extraction of useful iButton data. Code adapted from Eleanor Walker (2017)
`%notin%` <- function(x,y) !(x %in% y)

#Read Svalbard iButton Data. Create list of temperature data files
filenames <- list.files(path="data/microclimate/iButtons_OTCs") 
# read in all data files into large list
tempdata <- lapply(filenames, function(x) read.csv(paste("data/microclimate/iButtons_OTCs/",x,sep=""), sep=",", header=F, skip=20)) 

for(i in 1:length(filenames)) {
  #colnames(tempdata[[i]]) <- tempdata[[i]][-1,]; # put first line as header
  tempdata[[i]] <- tempdata[[i]][,-2]; # delete dispensable second column
  tempdata[[i]][,2] <- as.numeric(as.character(tempdata[[i]][,2]))
  names(tempdata[[i]])[c(1:2)]<-c("Date/Time","Value")
  tempdata[[i]][,1] <- as.POSIXct(tempdata[[i]][,1],format="%d/%m/%Y %H:%M") # convert date column into date format
  #days <- tempdata[[i]][,1]
  #tempdata[[i]] <- aggregate(tempdata[[i]][,2], by = list(days), mean); # compute mean temperatures per day
  
  fileinfo <- strsplit(file_path_sans_ext(filenames[i]),"_")[[1]]; # extract plot infos from file name
  
  for (j in 1:3) {
    tempdata[[i]][,2+j] <- fileinfo[j]
  } # add columns for plot infos
  colnames(tempdata[[i]])[c(1:5)] <- c("date","value","plot","position","type")
}

# Raw Data conversion into usablw dataset
microclimate <- do.call("rbind", tempdata)
microclimate$date<-microclimate$date + years(2000) # Correct years
microclimate$day<-cut(microclimate$date, "day") # Extract days
microclimate$time<-format(microclimate$date, "%H:%M:%S") # Extract days
microclimate <- microclimate %>% 
  mutate(treatment_type = case_when(grepl("c", plot) ~ "Control",
                                    grepl("t", plot) ~ "Treatment")) %>% 
  filter(as.Date(day) <= as.Date("2023-08-29")) # Cut days thet were not valid

microclimate$datetime <- as.POSIXct(paste(microclimate$day, microclimate$time), format = "%Y-%m-%d %H:%M:%S")

# Filter temperature and humidity and calculate VPD
temperature <- microclimate %>% 
  filter(type == "temp") %>% 
  filter(position == "down") 
  

humidity <- microclimate %>% 
  filter(type == "hum") %>% 
  filter(position == "down")


# VPD needs a wide format
microclimate_wide <- microclimate %>%
  pivot_wider(
    id_cols = c(date, plot, position, day, time, treatment_type),
    names_from = type,
    values_from = value
  )

vpd <- microclimate_wide %>% 
  mutate(VPD = RHtoVPD(microclimate_wide$hum, microclimate_wide$temp, Pa = 101)) %>% 
  filter(position == "down")


# iButtons temperature ----
# Only did plots, as stats were done with TOMST data for temperature
# Stats
daily_temp <- temperature %>%
  filter(as.Date(day) <= as.Date("2023-08-29")) %>% 
  group_by(day,treatment_type,position) %>%
  summarise(
    max_temperature = max(value),
    min_temperature = min(value),
    mean_temperature = mean(value)
  )

# Plots 
(ibutt_mean_temp <- ggplot(daily_temp, aes(x = as.Date(day), y = mean_temperature, group = treatment_type)) +
    geom_point(data = temperature, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Mean Daily Temperature (˚C)"))


(ibutt_max_temp <- ggplot(daily_temp, aes(x = as.Date(day), y = max_temperature, group = treatment_type)) +
    geom_point(data = temperature, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size

    xlab("Date") +
    ylab("Average Max Daily Temperature (˚C)"))

(ibutt_min_temp <- ggplot(daily_temp, aes(x = as.Date(day), y = min_temperature, group = treatment_type)) +
    geom_point(data = temperature, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size

    xlab("Date") +
    ylab("Average Min Daily Temperature (˚C)"))

# Boxplot with jittered points
ibutt_temp_boxplot_mean <- ggplot(daily_temp, aes(x = treatment_type, y = mean_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 40)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_temp_boxplot_max <- ggplot(daily_temp, aes(x = treatment_type, y = max_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 40)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_temp_boxplot_min <- ggplot(daily_temp, aes(x = treatment_type, y = min_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 40)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Combined iButtons temp plot - Supplementary Information S6
ibuttons_temp_plot <- {
  ibutt_mean_temp + ibutt_temp_boxplot_mean + ibutt_max_temp + ibutt_temp_boxplot_max + ibutt_min_temp + ibutt_temp_boxplot_min + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))

ggsave("outcomes/microclimate/figure_S6.png", plot = ibuttons_temp_plot, width = 13, height = 13, dpi = 300)

# iButtons humidity ----
# Summary stats
daily_hum <- humidity %>%
  filter(position == "down") %>% 
  filter(as.Date(day) <= as.Date("2023-08-29")) %>% 
  group_by(day,treatment_type,position) %>%
  summarise(
    max_humidity = max(value),
    min_humidity = min(value),
    mean_humidity = mean(value)
  )

# Summary mean hum points from RAW points
# Control -   87.9 ± 0.123
# Treatment -  85.7  ± 0.138
humidity %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(value),
            SE = std.error(value))

# Summary max hum
#  Control 98.0 ± 0.665
#  Treatment   98.0 ± 0.536
daily_hum %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(max_humidity),
            SE = std.error(max_humidity))

# Summary min hum
# Control   64.7 ± 3.98
# Treatment 60.3 ± 3.29
daily_hum %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(min_humidity),
            SE = std.error(min_humidity))


# Normality - raw data - too many points - using Anderson-Darling normality test
normality_raw_hum_otc_control <- humidity %>% 
  filter(treatment_type == "Control") 

ad.test(normality_raw_hum_otc_control$value) # Not normal

normality_raw_hum_otc_treatment <- humidity %>% 
  filter(treatment_type == "Treatment") 

ad.test(normality_raw_hum_otc_treatment$value) # Not normal

# Normality - max temp - assumption not met
daily_hum %>%
  group_by(treatment_type) %>%
  shapiro_test(max_humidity)

# Normality - min temp - assumption not met
daily_hum %>%
  group_by(treatment_type) %>%
  shapiro_test(min_humidity)

# Draw a qq plot by group - does not look good
ggqqplot(humidity, x = "value", facet.by = "treatment_type")
ggqqplot(daily_hum, x = "max_humidity", facet.by = "treatment_type")
ggqqplot(daily_hum, x = "min_humidity", facet.by = "treatment_type")

# equality of variances - assumption not met
humidity %>% levene_test(value ~ treatment_type)

# Raw data test - Significant - W = 72887954, p-value < 2.2e-16
result_OTC_hum_raw <- wilcox.test(value ~ treatment_type, data = humidity)

# Wilcoxon test - Significant - W = 911.5, p-value = 0.4392
result_OTC_hum_daily_max <- wilcox.test(max_humidity ~ treatment_type, data = daily_hum)

# Wilcoxon test - Non-sig - W = 977, p-value = 0.2072
result_OTC_hum_daily_min <- wilcox.test(min_humidity ~ treatment_type, data = daily_hum)

# Plots 
(ibutt_mean_hum <- ggplot(daily_hum, aes(x = as.Date(day), y = mean_humidity, group = treatment_type)) +
    geom_point(data = humidity, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Mean Daily RH (%)"))


(ibutt_max_hum <- ggplot(daily_hum, aes(x = as.Date(day), y = max_humidity, group = treatment_type)) +
    geom_point(data = humidity, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Average Max Daily RH (%)"))

(ibutt_min_hum <- ggplot(daily_hum, aes(x = as.Date(day), y = min_humidity, group = treatment_type)) +
    geom_point(data = humidity, aes(x = as.Date(datetime), y = value, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Average Min Daily RH (%)"))

# Boxplot with jittered points
ibutt_hum_boxplot_mean <- ggplot(daily_hum, aes(x = treatment_type, y = mean_humidity, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_hum_boxplot_max <- ggplot(daily_hum, aes(x = treatment_type, y = max_humidity, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_hum_boxplot_min <- ggplot(daily_hum, aes(x = treatment_type, y = min_humidity, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Combined iButtons temp plot - Supplementary Information S7
ibuttons_hum_plot <- {
  ibutt_mean_hum + ibutt_hum_boxplot_mean + ibutt_max_hum + ibutt_hum_boxplot_max + ibutt_min_hum + ibutt_hum_boxplot_min + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))

ggsave("outcomes/microclimate/figure_S7.png", plot = ibuttons_hum_plot, width = 13, height = 13, dpi = 300)


# iButtons VPD ----
# Summary stats
daily_vpd <- vpd %>%
  filter(position == "down") %>% 
  filter(as.Date(day) <= as.Date("2023-08-29")) %>% 
  group_by(day,treatment_type,position) %>%
  summarise(
    max_vpd = max(VPD),
    min_vpd = min(VPD),
    mean_vpd = mean(VPD)
  )

# Summary mean VPD points from RAW points
# Control -   0.167 ± 0.00199
# Treatment -  0.252 ± 0.00371
vpd %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(VPD),
            SE = std.error(VPD))

# Summary max VPD
#  Control  0.639 ± 0.0904
#  Treatment  1.09 ± 0.165 
daily_vpd %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(max_vpd),
            SE = std.error(max_vpd))

# Summary min VPD
# Control  0.639  ± 0.0904
# Treatment  1.09  ± 0.165 
daily_vpd %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(min_vpd),
            SE = std.error(min_vpd))


# Normality - raw data - too many points - using Anderson-Darling normality test
normality_raw_vpd_otc_control <- vpd %>% 
  filter(treatment_type == "Control") 

ad.test(normality_raw_vpd_otc_control$VPD) # Not normal

normality_raw_vpd_otc_treatment <- vpd %>% 
  filter(treatment_type == "Treatment") 

ad.test(normality_raw_vpd_otc_treatment$VPD) # Not normal

# Normality - max temp - assumption not met
daily_vpd %>%
  group_by(treatment_type) %>%
  shapiro_test(max_vpd)

# Normality - min temp - assumption not met
daily_vpd %>%
  group_by(treatment_type) %>%
  shapiro_test(min_vpd)

# Draw a qq plot by group - does not look good
ggqqplot(vpd, x = "VPD", facet.by = "treatment_type")
ggqqplot(daily_vpd, x = "max_vpd", facet.by = "treatment_type")
ggqqplot(daily_vpd, x = "min_vpd", facet.by = "treatment_type")

# equality of variances - assumption not met
vpd %>% levene_test(VPD ~ treatment_type)

# Raw data test - Significant - W = 60388064, p-value < 2.2e-16
result_OTC_vpd_raw <- wilcox.test(VPD ~ treatment_type, data = vpd)

# Wilcoxon test - non-Significant - W = 640, p-value = 0.06362
result_OTC_vpd_daily_max <- wilcox.test(max_vpd ~ treatment_type, data = daily_vpd)

# Wilcoxon test - Non-sig - W = 766.5, p-value = 0.42
result_OTC_vpd_daily_min <- wilcox.test(min_vpd ~ treatment_type, data = daily_vpd)

# Plots 
(ibutt_mean_vpd <- ggplot(daily_vpd, aes(x = as.Date(day), y = mean_vpd, group = treatment_type)) +
    geom_point(data = vpd, aes(x = as.Date(date), y = VPD, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Mean Daily VPD (kPa)"))


(ibutt_max_vpd <- ggplot(daily_vpd, aes(x = as.Date(day), y = max_vpd, group = treatment_type)) +
    geom_point(data = vpd, aes(x = as.Date(date), y = VPD, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Average Max Daily VPD (kPa)"))

(ibutt_min_vpd <- ggplot(daily_vpd, aes(x = as.Date(day), y = min_vpd, group = treatment_type)) +
    geom_point(data = vpd, aes(x = as.Date(date), y = VPD, color = treatment_type), size = 1, alpha = 0.05, colour = "grey") +
    geom_line(aes(color = treatment_type)) +
    scale_color_manual(values = my_colors) +  # Manually set colors
    theme_bw(base_size = 18) +
    
    theme(axis.title.x = element_blank(),
          legend.position = "none",  # Adjust the values as needed
          legend.title = element_blank(),    # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    
    xlab("Date") +
    ylab("Average Min Daily VPD (kPa)"))

# Boxplot with jittered points
ibutt_vpd_boxplot_mean <- ggplot(daily_vpd, aes(x = treatment_type, y = mean_vpd, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 4.7)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_vpd_boxplot_max <- ggplot(daily_vpd, aes(x = treatment_type, y = max_vpd, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 4.7)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_vpd_boxplot_min <- ggplot(daily_vpd, aes(x = treatment_type, y = min_vpd, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 4.7)) +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Combined iButtons temp plot - Supplementary Information S8
ibuttons_vpd_plot <- {
  ibutt_mean_vpd + ibutt_vpd_boxplot_mean + ibutt_max_vpd + ibutt_vpd_boxplot_max + ibutt_min_vpd + ibutt_vpd_boxplot_min + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))

ggsave("outcomes/microclimate/figure_S8.png", plot = ibuttons_vpd_plot, width = 13, height = 13, dpi = 300)

### VPD COUNTS DAYS 0.3 kPa
# From 41 measured days, there were 15 days inside the OTCs when the mean VPD was higher than 0.3 kPa, 
# compared to 8 days in control plots.
counts_below_0.3 <- table(daily_vpd$treatment_type, daily_vpd$mean_vpd < 0.3)


#### Raw data figures ####
# Figures for raw data TOMST temperature
tomst_raw_temp_boxplot <- ggplot(tomst_otc_temp, aes(x = treatment_type, y = Value, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  #geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(3, 25)) +
  
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Histogram raw data temperature tomst
hist_temp_otc <- ggplot(tomst_otc_temp, aes(x=Value, fill=treatment_type)) +
  geom_histogram(binwidth=0.5, position="dodge", alpha=0.7) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  
  # Other styling options
  ylab("Frequency") +
  xlab("Temperature (˚C)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Raw data for ibuttons humidity
raw_hum_boxplot <- ggplot(humidity, aes(x = treatment_type, y = value, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  #geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Histogram raw data humidity
hist_hum <- ggplot(humidity, aes(x=value, fill=treatment_type)) +
  geom_histogram(binwidth=2.5, position="dodge", alpha=0.7) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  
  # Other styling options
  ylab("Frequency") +
  xlab("Relative humidity (%)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Raw data for ibuttons VPD
raw_vpd_boxplot <- ggplot(vpd, aes(x = treatment_type, y = VPD, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  #geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Histogram raw data humidity
hist_vpd <- ggplot(vpd, aes(x=VPD, fill=treatment_type)) +
  geom_histogram(binwidth=0.2, position="dodge", alpha=0.7) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  
  # Other styling options
  ylab("Frequency") +
  xlab("VPD (kPa)") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

raw_microclimate_plot <- {
  hist_temp_otc + tomst_raw_temp_boxplot + hist_hum + raw_hum_boxplot + hist_vpd + raw_vpd_boxplot + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))

ggsave("outcomes/microclimate/figure_S4.png", plot = raw_microclimate_plot, width = 13, height = 13, dpi = 300)


#### LIGHT OTCS ----
### Part 1: Light measurements ### ----
light_OTC <- read.csv("data/light_OTC.csv") %>% 
  mutate(in_out = case_when(grepl("out", in_out) ~ "Control",
                            grepl("in", in_out) ~ "Treatment"))

light_OTC$plot_n <- as.factor(light_OTC$plot_n)
light_OTC$in_out <- as.factor(light_OTC$in_out)

# light means
light_OTC_means <- light_OTC %>% 
  group_by(plot_n, in_out) %>% 
  summarise(mean_light = mean(value)) 

average_per_OTC <- light_OTC %>% 
  group_by(in_out) %>% 
  summarise(mean_value = mean(value),
            SE = std.error(value))

# Light stats
shapiro.test(light_OTC$value)   # is normally distributed, p > 0.05
t.test(value ~ in_out, data = light_OTC) # t = 4.5459, df = 47.718, p-value = 3.745e-05

# plot
light_OTC_boxplot_mean <- ggplot(light_OTC_means, aes(x = in_out, y = mean_light, fill = in_out)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  ggtitle("(d)") +
  
  # Other styling options
  #ggtitle("Mean points = Non sig") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  ylab(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")"))  # Using expression() to include superscripts


#### FIGURE OTCS MAIN MANUSCRIPT FIGURE 2 ----
# Maximum daily temp tomst + minimum RH ibut + max daily vpd ibut + light levels
tomst_max_temp_box_main_text <- ggplot(OTC_temp_stats, aes(x = treatment_type, y = max_temperature, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(3, 25)) +
  
  # Other styling options
  ylab("Average maximum daily temperature (˚C)") +
  ggtitle("(a)") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_hum_boxplot_min_main_text <- ggplot(daily_hum, aes(x = treatment_type, y = min_humidity, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 100)) +
  ylab("Average minimum daily RH (%)") +
  ggtitle("(b)") +
  # Other styling options
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

ibutt_vpd_boxplot_max_main_text <- ggplot(daily_vpd, aes(x = treatment_type, y = max_vpd, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4,outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 4.7)) +
  # Other styling options
  ylab("Average maximum daily VPD (kPa)") +
  ggtitle("(c)") +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))


figure_2 <- {
  tomst_max_temp_box_main_text + ibutt_hum_boxplot_min_main_text + ibutt_vpd_boxplot_max_main_text + light_OTC_boxplot_mean + plot_layout(ncol = 2)
} +
  plot_layout(nrow = 3, heights = c(10, 10))

ggsave("outcomes/microclimate/figure_2.png", plot = figure_2, width = 10, height = 13.5, dpi = 300)


