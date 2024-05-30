# Raw desiccation curves from gas-exchange
# Jiri Subrt
# 24/07/2023

library(tidyverse)

# Load data
wc <- read.csv("data/lcp_wc/correct_wc_cleaned.csv")

wc$temperature <- as.factor(wc$temperature)
wc$sample <- as.factor(wc$sample)
wc$light_dark <- as.factor(wc$light_dark)
wc$CO2_PER_WEIGHT <- as.numeric(wc$CO2_PER_WEIGHT)

# converting relative time to total minutes
convert_to_minutes <- function(time_str) {
  parts <- strsplit(time_str, ":")[[1]]
  hours <- as.numeric(parts[1])
  minutes <- as.numeric(parts[2])
  total_minutes <- hours * 60 + minutes
  return(total_minutes)
}

wc <- wc %>%
  mutate(total.minutes = sapply(relative.minutes, convert_to_minutes))

# renaming t3 t4 t5 t6 to t1 t2 t3 t4
wc_new <- wc %>% 
  mutate(sample = case_when(sample == "T3" ~ "T1",
                            sample == "T4" ~ "T2",
                            sample == "T5" ~ "T3",
                            sample == "T6" ~ "T4",
                            TRUE ~ sample)) %>% 
  mutate(CO2_PER_WEIGHT = CO2_PER_WEIGHT/1000)


wc_25 <- ggplot(filter(wc_new, temperature=="25"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) + 
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab("") + 
  theme(legend.position='none') +
  ggtitle("25˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )

wc_20 <- ggplot(filter(wc_new, temperature=="20"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) +
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab("") +  
  theme(legend.position='none') +
  ggtitle("20˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )

wc_15 <- ggplot(filter(wc_new, temperature=="15"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) +
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab("") + 
  theme(legend.position='none') +
  ggtitle("15˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )

wc_10 <- ggplot(filter(wc_new, temperature=="10"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) +
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab("") + 
  theme(legend.position='none') +
  ggtitle("10˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )

wc_5 <- ggplot(filter(wc_new, temperature=="5"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) +
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab("") + 
  theme(legend.position='none') +
  ggtitle("5˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )

wc_0 <- ggplot(filter(wc_new, temperature=="0"), aes(x=total.minutes, y=CO2_PER_WEIGHT, color = light_dark)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample, ncol = 1) +
  #scale_x_reverse() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_y_continuous(limits = c(-6.6,3.7)) +
  scale_x_continuous(limits = c(0,400)) +
  scale_color_manual(values = c('brown', 'light green')) +
  xlab("") +
  ylab(label = expression(paste(
    "Net Photosynthesis µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")"))) + 
  theme(legend.position='none') +
  ggtitle("0˚C") +
  theme(
    legend.position = 'none',
    strip.text = element_blank()  # Remove the facet labels
  )


figure_s10 <- wc_0 + wc_5 + wc_10 + wc_15 + wc_20 + wc_25  +
  plot_layout(ncol  = 6)

ggsave("outcomes/dessication_curves/figure_s10.png", plot = figure_s10, width = 11, height = 6, dpi = 300)
