# CUE script
# Jiri Subrt
# 28/07/2023

library(tidyverse)

# Getting max NP values, and corresponsing DR values
cuedata <- read.csv("data/cue/np_dr_gp_corresponding.csv")

# Set values that are negative to 0, because CUE cannot be zero
cue_capped <- cuedata %>% 
  mutate(np_dw = np_dw/1000) %>% 
  pivot_wider(names_from = measurement, values_from = np_dw) %>% 
  group_by(sample, temperature, treatment_type) %>% 
  mutate(CUE = NP / GP * 100,
         CUE = ifelse(CUE < 0, 0, CUE)) %>%  # Set negative values of CUE to 0
  dplyr::select(!c(GP, DR, NP)) 

cue_capped$temperature <- as.factor(cue_capped$temperature)
cue_capped$treatment_type <- as.factor(cue_capped$treatment_type)

colours <- c("dark blue", "red")

figure_5 <- ggplot(cue_capped, aes(x = temperature, y = CUE, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  #ggtitle("(CUE)") +
  ylab("Carbon use efficiency (%)") +
  xlab("Temperature (˚C)") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

ggsave("outcomes/cue/figure_5.png", plot = figure_5, width = 10, height = 5, dpi = 300)


