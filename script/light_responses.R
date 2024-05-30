##%######################################################%##
#                                                          #
####        Light Response Curve analyses # Jiri        ####
####              Subrt # 6th of July 2023              ####
#                                                          #
##%######################################################%##

# Libraries
library(tidyverse)
library(onls)
library(viridis)
library(ggpubr)
library(ggalt)
library(grid)
library(stats)
library(plotrix)
library(ggnewscale)
library(broom)
library(nls.multstart)
library(cowplot)

# Function for the LRC fitting - from Tomeo 2017
diagnostic_AQ_plot <- function(curve_data, fit_data, Photo, PARi, group_id,
                               save_to_pdf = FALSE, save_path, file_name){
  if(save_to_pdf == FALSE){ 
    par(mar = c(3, 3, 1, 1), oma = c(1, 1, 1, 1))
    for(i in seq_along(1:length(unique(curve_data[[group_id]])))){
      single_curve <- 
        curve_data[curve_data[[group_id]] == 
                     unique(curve_data[[group_id]])[i],]
      plot(
        single_curve[[Photo]] ~ single_curve[[PARi]] ,
        xlim = c(-2, max(curve_data[[PARi]])), 
        ylim = c(min(curve_data[[Photo]]) - 2,
                 max(curve_data[[Photo]]) + 2),
        pch = 3,
        cex = 2,
        xlab = "",
        ylab = "",
        main = paste("Data from curve ",
                     as.character(
                       unique(single_curve[[group_id]])))
      )
      mtext(expression("Photo (µmol "*CO[2]*" "*m^-2*" "*s^-1*")"),
            line = 2.4, side = 2)
      mtext(expression("PARi (µmol photons "*m^-2*" "*s^-1*")"),
            line = 2.4, side = 1)
      par(new = TRUE)
      curve(((
        fit_data$Phi[i] * PARi + fit_data$Asat[i] - 
          sqrt((fit_data$Phi[i] * PARi + fit_data$Asat[i])^2 - 4 *
                 fit_data$Phi[i] * fit_data$theta[i] * PARi *
                 fit_data$Asat[i])
      ) / (2*fit_data$theta[i]) - fit_data$Rd[i]),
      from = 0, to = 1600, 
      xname = "PARi",
      xlab = "", ylab = "", 
      xlim = c(-2, max(curve_data[[PARi]])), 
      ylim = c(min(curve_data[[Photo]]) - 2,
               max(curve_data[[Photo]]) + 2),
      axes = FALSE,
      col = "red",
      lwd = 2
      )
    }} else{
      if(dir.exists(save_path)){
        pdf(paste0(save_path, file_name, ".pdf"))
        par(mar = c(3, 3, 1, 1), oma = c(1, 1, 1, 1))
        for(i in seq_along(1:length(unique(curve_data[[group_id]])))){
          single_curve <- 
            curve_data[curve_data[[group_id]] == 
                         unique(curve_data[[group_id]])[i],]
          plot(
            single_curve[[Photo]] ~ single_curve[[PARi]] ,
            xlim = c(-2, max(curve_data[[PARi]])), 
            ylim = c(min(curve_data[[Photo]]) - 2,
                     max(curve_data[[Photo]]) + 2),
            pch = 3,
            cex = 2,
            xlab = "",
            ylab = "",
            main = paste("Data from curve ",
                         as.character(
                           unique(single_curve[[group_id]])))
          )
          mtext(expression("Photo (µmol "*CO[2]*" "*m^-2*" "*s^-1*")"),
                line = 2.4, side = 2)
          mtext(expression("PARi (µmol photons "*m^-2*" "*s^-1*")"),
                line = 2.4, side = 1)
          par(new = TRUE)
          curve(((
            fit_data$Phi[i] * PARi + fit_data$Asat[i] - 
              sqrt((fit_data$Phi[i] * PARi + fit_data$Asat[i])^2 - 4 *
                     fit_data$Phi[i] * fit_data$theta[i] * PARi *
                     fit_data$Asat[i])
          ) / (2*fit_data$theta[i]) - fit_data$Rd[i]),
          from = 0, to = 1600, 
          xname = "PARi",
          xlab = "", ylab = "", 
          xlim = c(-2, max(curve_data[[PARi]])), 
          ylim = c(min(curve_data[[Photo]]) - 2,
                   max(curve_data[[Photo]]) + 2),
          axes = FALSE,
          col = "red",
          lwd = 2
          )
        }
        dev.off()
      } else {
        return(
          "Warning: the file path provided to save_path does not exist"
        )}
      
    }
}

# Data import
lrc <- read.csv("data/lcp_wc/lrc_data.csv")

# Extract important columns from the data: Light,	NP_DW,	sample, treatment_type
lrc_short <- lrc %>% 
  select(light, CO2_PER_WEIGHT, sample, temperature)

# Add treatment column based on conditions
lrc_short$treatment <- ifelse(grepl("T", lrc_short$sample), "treatment",
                              ifelse(grepl("C", lrc_short$sample), "control", NA))

# Change samples for just numbers
lrc_short <- lrc_short %>% 
  rename(NP_DW = CO2_PER_WEIGHT) %>% 
  mutate(sample_num = case_when(sample == "C1" ~ "C1",
                                sample == "C2" ~ "C2",
                                sample == "C3" ~ "C3",
                                sample == "C4" ~ "C4",
                                sample == "T3" ~ "T3",
                                sample == "T4" ~ "T4",
                                sample == "T5" ~ "T5",
                                sample == "T6" ~ "T6")) %>% 
  mutate(sample = as.factor(sample)) %>% 
  mutate(sample_num = as.character(sample_num))


# split into different temperatures
temp_0 <- lrc_short %>% 
  filter(temperature == 0)

temp_5 <- lrc_short %>% 
  filter(temperature == 5) 

temp_10 <- lrc_short %>% 
  filter(temperature == 10)

temp_15 <- lrc_short %>% 
  filter(temperature == 15)

#temp_17 <- lrc_short %>% 
#  filter(temperature == 17.5)

temp_20 <- lrc_short %>% 
  filter(temperature == 20)

temp_25 <- lrc_short %>% 
  filter(temperature == 25)

# creating multiple photosynthetic light response (A-Q) curves 
fit_AQ_curve <- function(df, group_id, Photo, PARi, fit_type = "onls"){
  AQ_curve_fits <- data.frame(ID = character(),
                              Asat = numeric(),
                              Phi = numeric(),
                              Rd = numeric(),
                              theta = numeric(),
                              resid_SSs = numeric(),
                              LCP = numeric(),
                              Q_sat_75 = numeric(),
                              Q_sat_85 = numeric(),  
                              stringsAsFactors = FALSE
  )
  if(fit_type == "onls"){
    if(require("onls")){
      print("onls is loaded correctly")
    } else {
      print("trying to install onls")
      install.packages("onls")
      if(require("onls")){
        print("onls installed and loaded")
      } else {
        stop("could not install onls")
      }
    }
    library("onls")      
    for(i in seq_along(unique(df[[group_id]]))){
      tryCatch({
        AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
        # Subset by group_ID iteratively:
        single_curve1 <- df[df[[group_id]] == unique(df[[group_id]])[i],]
        single_curve1$assim <- single_curve1[[Photo]]
        single_curve1$PAR <- single_curve1[[PARi]]
        single_curve = single_curve1[order(single_curve1$PAR),]
        phi.as.slope <- with(single_curve,
                             as.numeric(coef(lm(
                               assim[1:5] ~ PAR[1:5]))[2]))
        # Fit the curve:
        temp.fit <- with(single_curve, # use the subset of a single curve
                         onls(assim ~ ((Phi * PAR + Asat - 
                                          sqrt((Phi * PAR + Asat)^2 - 
                                                 4 * Phi * theta * 
                                                 Asat * PAR ))
                         )/(2*theta) - Rd,
                         start=list(
                           Asat = (max(assim)),
                           Phi = phi.as.slope,
                           Rd = -min(assim),
                           theta = 0.5),
                         control = list(maxiter = 50)#,
                         #algorithm = "port"
                         )
        )
        AQ_curve_fits[i, 2] <- as.numeric(coef(temp.fit)[1]) # asat 
        AQ_curve_fits[i, 3] <- as.numeric(coef(temp.fit)[2]) # Phi
        AQ_curve_fits[i, 4] <- as.numeric(coef(temp.fit)[3]) # Rd
        AQ_curve_fits[i, 5] <- as.numeric(coef(temp.fit)[4]) # theta
        AQ_curve_fits[i, 6] <- sum(resid(temp.fit)^2)
        AQ_curve_fits[i, 7] <- (as.numeric(coef(temp.fit)[3]) *(
          as.numeric(coef(temp.fit)[3]) * as.numeric(coef(temp.fit)[4]) - 
            as.numeric(coef(temp.fit)[1]))
        ) / (as.numeric(coef(temp.fit)[2]) * (
          as.numeric(coef(temp.fit)[3]) - as.numeric(coef(temp.fit)[1])
        ))
        AQ_curve_fits[i, 8] <- (
          (as.numeric(coef(temp.fit)[1]) * 0.75 + 
             (as.numeric(coef(temp.fit)[3]))) * (
               as.numeric(coef(temp.fit)[1]) * 0.75 *
                 as.numeric(coef(temp.fit)[4]) +
                 as.numeric(coef(temp.fit)[3]) *
                 as.numeric(coef(temp.fit)[4]) -
                 as.numeric(coef(temp.fit)[1])
             )) / (
               as.numeric(coef(temp.fit)[2])* (
                 as.numeric(coef(temp.fit)[1]) * 0.75 +
                   as.numeric(coef(temp.fit)[3]) -
                   as.numeric(coef(temp.fit)[1])
               ))
        
        AQ_curve_fits[i, 9] <- (
          (as.numeric(coef(temp.fit)[1]) * 0.85 + 
             (as.numeric(coef(temp.fit)[3]))) * (
               as.numeric(coef(temp.fit)[1]) * 0.85 *
                 as.numeric(coef(temp.fit)[4]) +
                 as.numeric(coef(temp.fit)[3]) *
                 as.numeric(coef(temp.fit)[4]) -
                 as.numeric(coef(temp.fit)[1])
             )) / (
               as.numeric(coef(temp.fit)[2])* (
                 as.numeric(coef(temp.fit)[1]) * 0.85 +
                   as.numeric(coef(temp.fit)[3]) -
                   as.numeric(coef(temp.fit)[1])
               ))
      }, error = function(E){cat("Error: ", conditionMessage(E), "\n")})
    }
    return(AQ_curve_fits)
  } else{
    if(fit_type == "nls"){
      for(i in seq_along(unique(df[[group_id]]))){
        tryCatch({
          AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
          # Subset by group_ID iteratively:
          single_curve1 <- df[df[[group_id]] == unique(df[[group_id]])[i],]
          single_curve1$assim <- single_curve1[[Photo]]
          single_curve1$PAR <- single_curve1[[PARi]]
          single_curve = single_curve1[order(single_curve1$PAR),]
          phi.as.slope <- with(single_curve,
                               as.numeric(coef(lm(
                                 assim[1:5] ~ PAR[1:5]))[2]))
          # Fit the curve:
          temp.fit <- with(single_curve, 
                           nls(assim ~ ((Phi * PAR + Asat - 
                                           sqrt((Phi * PAR + Asat)^2 - 
                                                  4 * Phi * theta * 
                                                  Asat * PAR ))
                           )/(2*theta) - Rd,
                           start=list(
                             Asat = (max(assim)),
                             Phi = phi.as.slope,
                             Rd = -min(assim),
                             theta = 0.5),
                           control = list(maxiter = 50),
                           algorithm = "port")
          )
          AQ_curve_fits[i, 2] <- as.numeric(coef(temp.fit)[1]) # asat 
          AQ_curve_fits[i, 3] <- as.numeric(coef(temp.fit)[2]) # Phi
          AQ_curve_fits[i, 4] <- as.numeric(coef(temp.fit)[3]) # Rd
          AQ_curve_fits[i, 5] <- as.numeric(coef(temp.fit)[4]) # theta
          AQ_curve_fits[i, 6] <- sum(resid(temp.fit)^2)
          AQ_curve_fits[i, 7] <- (as.numeric(coef(temp.fit)[3]) *(
            as.numeric(coef(temp.fit)[3]) * 
              as.numeric(coef(temp.fit)[4]) - 
              as.numeric(coef(temp.fit)[1]))
          ) / (as.numeric(coef(temp.fit)[2]) * (
            as.numeric(coef(temp.fit)[3]) - 
              as.numeric(coef(temp.fit)[1])
          ))
          AQ_curve_fits[i, 8] <- (
            (as.numeric(coef(temp.fit)[1]) * 0.75 + 
               (as.numeric(coef(temp.fit)[3]))) * (
                 as.numeric(coef(temp.fit)[1]) * 0.75 *
                   as.numeric(coef(temp.fit)[4]) +
                   as.numeric(coef(temp.fit)[3]) *
                   as.numeric(coef(temp.fit)[4]) -
                   as.numeric(coef(temp.fit)[1])
               )) / (
                 as.numeric(coef(temp.fit)[2])* (
                   as.numeric(coef(temp.fit)[1]) * 0.75 +
                     as.numeric(coef(temp.fit)[3]) -
                     as.numeric(coef(temp.fit)[1])
                 ))
          AQ_curve_fits[i, 9] <- (
            (as.numeric(coef(temp.fit)[1]) * 0.85 + 
               (as.numeric(coef(temp.fit)[3]))) * (
                 as.numeric(coef(temp.fit)[1]) * 0.85 *
                   as.numeric(coef(temp.fit)[4]) +
                   as.numeric(coef(temp.fit)[3]) *
                   as.numeric(coef(temp.fit)[4]) -
                   as.numeric(coef(temp.fit)[1])
               )) / (
                 as.numeric(coef(temp.fit)[2])* (
                   as.numeric(coef(temp.fit)[1]) * 0.85 +
                     as.numeric(coef(temp.fit)[3]) -
                     as.numeric(coef(temp.fit)[1])
                 ))
        }, error = function(E){
          cat("Error: ", conditionMessage(E), "\n")})
      }
      return(AQ_curve_fits)      
    } else{print("ERROR: 'fit_type' specified incorrectly.")}
  }
}

# temp_0 ----
temp_0 <- fit_AQ_curve(df = filter(temp_0),
                       Photo = "NP_DW",
                       PARi = "light",
                       group_id = "sample",
                       fit_type = "onls")

temp_0 <- temp_0 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 0) 
  

# 0˚C C1
curve_c1_0 <- function(PAR){
  (7.550395e-03 * PAR + 1.427091e+00 - 
     sqrt((7.550395e-03 * PAR + 1.427091e+00)^2 - 4 *
            7.550395e-03 * 7.107184e-01 * PAR *
            1.427091e+00)
  ) / (2*7.107184e-01) - 0.7741790
}

par_c1_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c1_0(PAR = 0:1200),
                       sample = "C1")

# 0˚C C2
curve_c2_0 <- function(PAR){
  (3.035465e-02 * PAR + 1.982598e+00 - 
     sqrt((3.035465e-02 * PAR + 1.982598e+00)^2 - 4 *
            3.035465e-02 * -7.189349e+00 * PAR *
            1.982598e+00)
  ) / (2*-7.189349e+00) - 0.6450219
}

par_c2_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c2_0(PAR = 0:1200),
                       sample = "C2")

# 0˚C C3
curve_c3_0 <- function(PAR){
  (1.805747e-02 * PAR + 1.582488e+00 - 
     sqrt((1.805747e-02 * PAR + 1.582488e+00)^2 - 4 *
            1.805747e-02 * 8.161075e-04 * PAR *
            1.582488e+00)
  ) / (2*8.161075e-04) - 0.9637622
}

par_c3_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c3_0(PAR = 0:1200),
                       sample = "C3")

# 0˚C C4
curve_c4_0 <- function(PAR){
  (2.139457e-02 * PAR + 2.305223e+00  - 
     sqrt((2.139457e-02 * PAR + 2.305223e+00)^2 - 4 *
            2.139457e-02* 6.311851e-01 * PAR *
            2.305223e+00)
  ) / (2*6.311851e-01) - 1.3447847
}

par_c4_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c4_0(PAR = 0:1200),
                       sample = "C4")

# 0˚C T3
curve_t3_0 <- function(PAR){
  (5.768933e-05 * PAR + -9.137169e-03  - 
     sqrt((5.768933e-05 * PAR + -9.137169e-03)^2 - 4 *
            5.768933e-05* 6.094533e-03 * PAR *
            -9.137169e-03)
  ) / (2*6.094533e-03) - 0.5333318
}

par_t3_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t3_0(PAR = 0:1200),
                       sample = "T3")

# 0˚C T4
curve_t4_0 <- function(PAR){
  (1.863563e-11 * PAR + -2.463994e-09  - 
     sqrt((1.863563e-11 * PAR + -2.463994e-09)^2 - 4 *
            1.863563e-11* 3.504169e-09 * PAR *
            -2.463994e-09)
  ) / (2*3.504169e-09) - 0.4441560
}

par_t4_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t4_0(PAR = 0:1200),
                       sample = "T4")

# 0˚C T5
curve_t5_0 <- function(PAR){
  (8.980728e-03* PAR + 8.760918e-01  - 
     sqrt((8.980728e-03 * PAR + 8.760918e-01)^2 - 4 *
            8.980728e-03* 9.954384e-01 * PAR *
            8.760918e-01)
  ) / (2*9.954384e-01) - 1.0113035
}

par_t5_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t5_0(PAR = 0:1200),
                       sample = "T5")

# 0˚C T6
curve_t6_0 <- function(PAR){
  (1.228884e-12 * PAR + -1.402095e-10  - 
     sqrt((1.228884e-12 * PAR + -1.402095e-10)^2 - 4 *
            1.228884e-12* 3.128248e-10 * PAR *
            -1.402095e-10)
  ) / (2*3.128248e-10) - 0.3700296
}

par_t6_0 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t6_0(PAR = 0:1200),
                       sample = "T6")

# all 0 together
par.0 <- rbind(par_c1_0, par_c2_0, par_c3_0, par_c4_0, par_t3_0, par_t4_0, par_t5_0, par_t6_0)
par.0 <- par.0 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 0)

lrc_short_0 <- lrc_short %>% 
  filter(temperature == 0)

# horizontal plot 
(light_curves_0 <- ggplot(par.0, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_0, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('0˚C') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(-7, 3))

# temp_5 curve ---- 
temp_5 <- fit_AQ_curve(df = temp_5,
                       Photo = "NP_DW",
                       PARi = "light",
                       group_id = "sample_num",
                       fit_type = "onls")

try <- fit_AQ_curve(df = temp_5, group_id = "sample_num", Photo = "NP_DW", PARi = "light", fit_type = "onls")

temp_5 <- temp_5 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 5)

# 5˚C C1 - DID NOT WORK HAVE NAs


# 5˚C C2 
curve_c2_5 <- function(PAR){
  (0.04933215 * PAR + 1.994893 - 
     sqrt((0.04933215 * PAR + 1.994893)^2 - 4 *
            0.04933215 * -0.7072987 * PAR *
            1.994893)
  ) / (2*-0.7072987) - 0.9447595
}


par_c2_5 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_5(PAR = 0:1200),
                        sample = "C2")

# 5˚C C3
curve_c3_5 <- function(PAR){
  (0.02523834 * PAR + 2.333988 - 
     sqrt((0.02523834 * PAR + 2.333988)^2 - 4 *
            0.02523834 * 0.5574019 * PAR *
            2.333988)
  ) / (2*0.5574019) - 1.2820530
}


par_c3_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c3_5(PAR = 0:1200),
                       sample = "C3")

# 5˚C C4
curve_c4_5 <- function(PAR){
  (0.03154612 * PAR + 3.393663 - 
     sqrt((0.03154612 * PAR + 3.393663)^2 - 4 *
            0.03154612 * 0.5598045 * PAR *
            3.393663)
  ) / (2*0.5598045) - 2.0263203
}


par_c4_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_c4_5(PAR = 0:1200),
                       sample = "C4")

# 5˚C T5
curve_t5_5 <- function(PAR){
  (0.44567584 * PAR + 2.663236 - 
     sqrt((0.44567584 * PAR + 2.663236)^2 - 4 *
            0.44567584 * -5.8581354 * PAR *
            2.663236)
  ) / (2*-5.8581354) - 2.2142130
}


par_t5_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t5_5(PAR = 0:1200),
                       sample = "T5")

# 5˚C T6
curve_t6_5 <- function(PAR){
  (0.02932903 * PAR + 2.280328 - 
     sqrt((0.02932903 * PAR + 2.280328)^2 - 4 *
            0.02932903 * 0.9734507 * PAR *
            2.280328)
  ) / (2*0.9734507) - 1.9053821
}

par_t6_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t6_5(PAR = 0:1200),
                       sample = "T6")

# 5˚C T3
curve_t3_5 <- function(PAR){
  (0.01960689 * PAR + 2.942371 - 
     sqrt((0.01960689 * PAR + 2.942371)^2 - 4 *
            0.01960689 * 0.6800631 * PAR *
            2.942371)
  ) / (2*0.6800631) - 2.3020977
}


par_t3_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t3_5(PAR = 0:1200),
                       sample = "T3")

# 5˚C T4
curve_t4_5 <- function(PAR){
  (0.02234915 * PAR + 1.748811 - 
     sqrt((0.02234915 * PAR + 1.748811)^2 - 4 *
            0.02234915 * 0.9221760 * PAR *
            1.748811)
  ) / (2*0.9221760) - 1.5454751
}


par_t4_5 <- data.frame(Lcuv = 0:1200,
                       curve = curve_t4_5(PAR = 0:1200),
                       sample = "T4")

# all 5 together 
par.5 <- rbind(par_c2_5, par_c3_5, par_c4_5, par_t3_5, par_t4_5, par_t5_5, par_t6_5)
par.5 <- par.5 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 5)

lrc_short_5 <- lrc_short %>% 
  filter(temperature == 5)

# horizontal plot 
(light_curves_5 <- ggplot(par.5, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_5, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('5˚C') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(-7, 3))


# temp_10 curve ----
temp_10 <- fit_AQ_curve(df = temp_10,
                       Photo = "NP_DW",
                       PARi = "light",
                       group_id = "sample",
                       fit_type = "onls")

temp_10 <- temp_10 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 10)

# 10˚C C1
curve_c1_10 <- function(PAR){
  (0.04027955 * PAR + 3.463369 - 
     sqrt((0.04027955 * PAR + 3.463369)^2 - 4 *
            0.04027955 * 0.4215458 * PAR *
            3.463369)
  ) / (2*0.4215458) - 1.674355
}

par_c1_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c1_10(PAR = 0:1200),
                        sample = "C1")

# 10˚C C2
curve_c2_10 <- function(PAR){
  (0.05281387 * PAR + 4.104246 - 
     sqrt((0.05281387 * PAR + 4.104246)^2 - 4 *
            0.05281387 * -0.1461889 * PAR *
            4.104246)
  ) / (2*-0.1461889) - 1.504160
}

par_c2_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_10(PAR = 0:1200),
                        sample = "C2")

# 10˚C T3
curve_t3_10 <- function(PAR){
  (0.05101219 * PAR + 3.927552 - 
     sqrt((0.05101219 * PAR + 3.927552)^2 - 4 *
            0.05101219 * 0.3174288 * PAR *
            3.927552)
  ) / (2*0.3174288) - 2.793004
}

par_t3_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t3_10(PAR = 0:1200),
                        sample = "T3")

# 10˚C T4
curve_t4_10 <- function(PAR){
  (0.03654339 * PAR + 2.211857 - 
     sqrt((0.03654339 * PAR + 2.211857)^2 - 4 *
            0.03654339 * 0.7278627 * PAR *
            2.211857)
  ) / (2*0.7278627) - 2.365316
}

par_t4_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t4_10(PAR = 0:1200),
                        sample = "T4")

# 10˚C T5
curve_t5_10 <- function(PAR){
  (0.04912724 * PAR + 2.706180 - 
     sqrt((0.04912724 * PAR + 2.706180)^2 - 4 *
            0.04912724 * -0.2241861 * PAR *
            2.706180)
  ) / (2*-0.2241861) - 2.640106
}

par_t5_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t5_10(PAR = 0:1200),
                        sample = "T5")

# 10˚C T6
curve_t6_10 <- function(PAR){
  (0.08139189 * PAR + 3.488607 - 
     sqrt((0.08139189 * PAR + 3.488607)^2 - 4 *
            0.08139189 * 0.3626814 * PAR *
            3.488607)
  ) / (2*0.3626814) - 2.644533
}

par_t6_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t6_10(PAR = 0:1200),
                        sample = "T6")

# 10˚C C3
curve_c3_10 <- function(PAR){
  (0.08526146 * PAR + 4.837332 - 
     sqrt((0.08526146 * PAR + 4.837332)^2 - 4 *
            0.08526146 * -0.8314027 * PAR *
            4.837332)
  ) / (2*-0.8314027) - 2.671903
}

par_c3_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c3_10(PAR = 0:1200),
                        sample = "C3")

# 10˚C C4
curve_c4_10 <- function(PAR){
  (0.07711163 * PAR + 5.122920 - 
     sqrt((0.07711163 * PAR + 5.122920)^2 - 4 *
            0.07711163 * -2.6228638 * PAR *
            5.122920)
  ) / (2*-2.6228638) - 2.329931
}

par_c4_10 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c4_10(PAR = 0:1200),
                        sample = "C4")

# all 10 together
par.10 <- rbind(par_c1_10, par_c2_10, par_c3_10, par_c4_10, par_t3_10, par_t4_10, par_t5_10, par_t6_10)
par.10 <- par.10 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 10)

lrc_short_10 <- lrc_short %>% 
  filter(temperature == 10)

# horizontal plot 
(light_curves_10 <- ggplot(par.10, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_10, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('10˚C') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(-7, 3))

# temp_15 curve ----
temp_15 <- fit_AQ_curve(df = temp_15,
                        Photo = "NP_DW",
                        PARi = "light",
                        group_id = "sample",
                        fit_type = "onls")

temp_15 <- temp_15 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 15)

# 15˚C C1 
curve_c1_15 <- function(PAR){
  (0.03561512 * PAR + 4.333145 - 
     sqrt((0.03561512 * PAR + 4.333145)^2 - 4 *
            0.03561512 * 0.59149471 * PAR *
            4.333145)
  ) / (2*0.59149471) - 2.041315
}

par_c1_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c1_15(PAR = 0:1200),
                        sample = "C1")

# 15˚C C2
curve_c2_15 <- function(PAR){
  (0.05355430 * PAR + 5.759149 - 
     sqrt((0.05355430 * PAR + 5.759149)^2 - 4 *
            0.05355430 * -0.19151041 * PAR *
            5.759149)
  ) / (2*-0.19151041) - 2.506958
}

par_c2_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_15(PAR = 0:1200),
                        sample = "C2")

# 15˚C C3
curve_c3_15 <- function(PAR){
  (0.06193610 * PAR + 5.460485 - 
     sqrt((0.06193610 * PAR + 5.460485)^2 - 4 *
            0.06193610 * -0.23420521 * PAR *
            5.460485)
  ) / (2*-0.23420521) - 3.580046
}

par_c3_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c3_15(PAR = 0:1200),
                        sample = "C3")

# 15˚C C4
curve_c4_15 <- function(PAR){
  (0.04917961 * PAR + 6.236724 - 
     sqrt((0.04917961 * PAR + 6.236724)^2 - 4 *
            0.04917961 * -0.90271316 * PAR *
            6.236724)
  ) / (2*-0.90271316) - 2.827840
}

par_c4_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c4_15(PAR = 0:1200),
                        sample = "C4")

# 15˚C T5
curve_t5_15 <- function(PAR){
  (0.09288995 * PAR + 4.647854 - 
     sqrt((0.09288995 * PAR + 4.647854)^2 - 4 *
            0.09288995 * -2.10215239 * PAR *
            4.647854)
  ) / (2*-2.10215239) - 3.747157
}

par_t5_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t5_15(PAR = 0:1200),
                        sample = "T5")

# 15˚C T6
curve_t6_15 <- function(PAR){
  (0.04950741 * PAR + 3.743937 - 
     sqrt((0.04950741 * PAR + 3.743937)^2 - 4 *
            0.04950741 *  0.79842629 * PAR *
            3.743937)
  ) / (2*0.79842629) - 3.347763 
}

par_t6_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t6_15(PAR = 0:1200),
                        sample = "T6")

# 15˚C T3
curve_t3_15 <- function(PAR){
  (0.09411343 * PAR + 6.261386 - 
     sqrt((0.09411343 * PAR + 6.261386)^2 - 4 *
            0.09411343 *  0.47057785 * PAR *
            6.261386)
  ) / (2*0.47057785) - 4.520277 
}

par_t3_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t3_15(PAR = 0:1200),
                        sample = "T3")

# 15˚C T4
curve_t4_15 <- function(PAR){
  (0.06864855 * PAR + 3.696427 - 
     sqrt((0.06864855 * PAR + 3.696427)^2 - 4 *
            0.06864855 *  -0.04518261 * PAR *
            3.696427)
  ) / (2*-0.04518261) - 3.116718 
}

par_t4_15 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t4_15(PAR = 0:1200),
                        sample = "T4")

# all 15 together
par.15 <- rbind(par_c1_15, par_c2_15, par_c3_15, par_c4_15, par_t3_15, par_t4_15, par_t5_15, par_t6_15)
par.15 <- par.15 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 15)

lrc_short_15 <- lrc_short %>% 
  filter(temperature == 15)

# horizontal plot 
(light_curves_15 <- ggplot(par.15, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_15, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('15˚C') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(-7, 3))

# temp_17 - NOT USING FOR MY ANALYSIS ----

temp_17 <- fit_AQ_curve(df = temp_17,
                        Photo = "NP_DW",
                        PARi = "light",
                        group_id = "sample",
                        fit_type = "onls")

temp_17 <- temp_17 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 17.5)

# 17˚C C1
curve_c1_17 <- function(PAR){
  (0.04621577 * PAR + 4.400506 - 
     sqrt((0.04621577 * PAR + 4.400506)^2 - 4 *
            0.04621577 *  0.22949616 * PAR *
            4.400506)
  ) / (2*0.22949616) - 2.530811 
}

par_c1_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c1_17(PAR = 0:1200),
                        sample = "C1")

# 17˚C C2
curve_c2_17 <- function(PAR){
  (0.04645883 * PAR + 6.083790 - 
     sqrt((0.04645883 * PAR + 6.083790)^2 - 4 *
            0.04645883 *  0.05663749 * PAR *
            6.083790)
  ) / (2*0.05663749) - 2.627307 
}

par_c2_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_17(PAR = 0:1200),
                        sample = "C2")

# 17˚C C3
curve_c3_17 <- function(PAR){
  (0.09297241 * PAR + 7.588515 - 
     sqrt((0.09297241 * PAR + 7.588515)^2 - 4 *
            0.09297241 *  -0.90796841 * PAR *
            7.588515)
  ) / (2*-0.90796841) - 5.199819 
}

par_c3_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c3_17(PAR = 0:1200),
                        sample = "C3")


# 17˚C C4
curve_c4_17 <- function(PAR){
  (0.05430653 * PAR + 6.945452 - 
     sqrt((0.05430653 * PAR + 6.945452)^2 - 4 *
            0.05430653 *  -0.33805715 * PAR *
            6.945452)
  ) / (2*-0.33805715) - 4.354104 
}

par_c4_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c4_17(PAR = 0:1200),
                        sample = "C4")

# 17˚C T5
curve_t5_17 <- function(PAR){
  (0.09261711 * PAR + 5.164643 - 
     sqrt((0.09261711 * PAR + 5.164643)^2 - 4 *
            0.09261711 *  -1.19066144  * PAR *
            5.164643)
  ) / (2*-1.19066144 ) - 5.276009 
}

par_t5_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t5_17(PAR = 0:1200),
                        sample = "T5")

# 17˚C T6
curve_t6_17 <- function(PAR){
  (0.05540621 * PAR + 5.225643 - 
     sqrt((0.05540621 * PAR + 5.225643)^2 - 4 *
            0.05540621 *  0.33848416  * PAR *
            5.225643)
  ) / (2*0.33848416) - 5.556099 
}

par_t6_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t6_17(PAR = 0:1200),
                        sample = "T6")

# 17˚C T3
curve_t3_17 <- function(PAR){
  (0.07174677 * PAR + 6.261555 - 
     sqrt((0.07174677 * PAR + 6.261555)^2 - 4 *
            0.07174677 *  0.46047444  * PAR *
            6.261555)
  ) / (2*0.46047444) - 4.065075 
}

par_t3_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t3_17(PAR = 0:1200),
                        sample = "T3")

# 17˚C T4
curve_t4_17 <- function(PAR){
  (0.07169892 * PAR + 4.762740 - 
     sqrt((0.07169892 * PAR + 4.762740)^2 - 4 *
            0.07169892 *  -0.15412619  * PAR *
            4.762740)
  ) / (2*-0.15412619) - 3.501585 
}

par_t4_17 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t4_17(PAR = 0:1200),
                        sample = "T4")

# all 17 together
par.17 <- rbind(par_c1_17, par_c2_17, par_c3_17, par_c4_17, par_t3_17, par_t4_17, par_t5_17, par_t6_17)
par.17 <- par.17 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment"))

lrc_short_17 <- lrc_short %>% 
  filter(temperature == 17.5)

# horizontal plot 
(light_curves_17 <- ggplot(par.17, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_17, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw() +
  theme(axis.title.x = 
          element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = 
          element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)), 
        panel.grid = element_blank(),
        legend.key.width = unit(1,"cm"),
        legend.position = "right") +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  scale_y_continuous(limits = c(-7, 3))


# temp_20 ----
temp_20 <- fit_AQ_curve(df = temp_20,
                        Photo = "NP_DW",
                        PARi = "light",
                        group_id = "sample",
                        fit_type = "onls")

temp_20 <- temp_20 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 20)

# 20˚C C1
curve_c1_20 <- function(PAR){
  (0.04264801 * PAR + 4.441757 - 
     sqrt((0.04264801 * PAR + 4.441757)^2 - 4 *
            0.04264801 *  0.2475178  * PAR *
            4.441757)
  ) / (2*0.2475178) - 2.548217 
}

par_c1_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c1_20(PAR = 0:1200),
                        sample = "C1")

# 20˚C C2
curve_c2_20 <- function(PAR){
  (0.04453881 * PAR + 5.531482 - 
     sqrt((0.04453881 * PAR + 5.531482)^2 - 4 *
            0.04453881 *  0.2351505  * PAR *
            5.531482)
  ) / (2*0.2351505) - 2.597481 
}

par_c2_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_20(PAR = 0:1200),
                        sample = "C2")

# 20˚C C3
curve_c3_20 <- function(PAR){
  (0.08739557 * PAR + 6.893084 - 
     sqrt((0.08739557 * PAR + 6.893084)^2 - 4 *
            0.08739557 *  -2.6543509  * PAR *
            6.893084)
  ) / (2*-2.6543509) - 4.676986 
}

par_c3_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c3_20(PAR = 0:1200),
                        sample = "C3")

# 20˚C C4
curve_c4_20 <- function(PAR){
  (0.05264463 * PAR + 6.508801 - 
     sqrt((0.05264463 * PAR + 6.508801)^2 - 4 *
            0.05264463 *  -0.1482757  * PAR *
            6.508801)
  ) / (2*-0.1482757) - 4.565908 
}

par_c4_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c4_20(PAR = 0:1200),
                        sample = "C4")

# 20˚C T5
curve_t5_20 <- function(PAR){
  (0.04870653 * PAR + 5.533976 - 
     sqrt((0.04870653 * PAR + 5.533976)^2 - 4 *
            0.04870653 *  0.3811350  * PAR *
            5.533976)
  ) / (2*0.3811350) - 5.675767 
}

par_t5_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t5_20(PAR = 0:1200),
                        sample = "T5")

# 20˚C T6 - GIVES NA - NEED TO FIND THE ISSUE
curve_t6_20 <- function(PAR){
  (0.05540621 * PAR + 5.225643 - 
     sqrt((0.05540621 * PAR + 5.225643)^2 - 4 *
            0.05540621 *  0.33848416  * PAR *
            5.225643)
  ) / (2*0.33848416) - 5.556099 
}

par_t6_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t6_20(PAR = 0:1200),
                        sample = "T6")

# 20˚C T3
curve_t3_20 <- function(PAR){
  (0.11463993 * PAR + 6.924705 - 
     sqrt((0.11463993 * PAR + 6.924705)^2 - 4 *
            0.11463993 *  -0.1442066  * PAR *
            6.924705)
  ) / (2*-0.1442066) - 5.770474 
}

par_t3_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t3_20(PAR = 0:1200),
                        sample = "T3")

# 20˚C T4
curve_t4_20 <- function(PAR){
  (0.04656819 * PAR + 4.257034 - 
     sqrt((0.04656819 * PAR + 4.257034)^2 - 4 *
            0.04656819 *  0.3763891  * PAR *
            4.257034)
  ) / (2*0.3763891) - 4.164198 
}

par_t4_20 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t4_20(PAR = 0:1200),
                        sample = "T4")

# all 20 together
par.20 <- rbind(par_c1_20, par_c2_20, par_c3_20, par_c4_20, par_t3_20, par_t5_20, par_t6_20)
par.20 <- par.20 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 20)

lrc_short_20 <- lrc_short %>% 
  filter(temperature == 20)

# horizontal plot 
(light_curves_20 <- ggplot(par.20, aes(x = Lcuv, y = curve, group = sample))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey30") +
    geom_point(data = lrc_short_20, aes(x = light, y = NP_DW, color = treatment,
                                        shape = treatment, fill = treatment), 
               alpha = 0.5, size = 4) +
    geom_line(aes(color = treatment, linetype = treatment),
              size = 0.8)) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('20˚C') +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(-7, 3))

# temp_25 ----

temp_25 <- fit_AQ_curve(df = temp_25,
                        Photo = "NP_DW",
                        PARi = "light",
                        group_id = "sample",
                        fit_type = "onls")

temp_25 <- temp_25 %>% 
  mutate(sample = case_when(ID == "1" ~ "C1",
                            ID == "2" ~ "C2",
                            ID == "3" ~ "C3",
                            ID == "4" ~ "C4",
                            ID == "5" ~ "T3",
                            ID == "6" ~ "T4",
                            ID == "7" ~ "T5",
                            ID == "8" ~ "T6")) %>% 
  mutate(temperature = 25)

# 25˚C C1
curve_c1_25 <- function(PAR){
  (0.06130804 * PAR + 4.219164 - 
     sqrt((0.06130804 * PAR + 4.219164)^2 - 4 *
            0.06130804 *  0.3840358  * PAR *
            4.219164)
  ) / (2*0.3840358) - 3.896985 
}

par_c1_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c1_25(PAR = 0:1200),
                        sample = "C1")

# 25˚C C2
curve_c2_25 <- function(PAR){
  (0.07547782 * PAR + 5.171124 - 
     sqrt((0.07547782 * PAR + 5.171124)^2 - 4 *
            0.07547782 *  -0.9175394  * PAR *
            5.171124)
  ) / (2*-0.9175394) - 3.924926 
}

par_c2_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c2_25(PAR = 0:1200),
                        sample = "C2")

# 25˚C C3
curve_c3_25 <- function(PAR){
  (0.08554596 * PAR + 7.231979 - 
     sqrt((0.08554596 * PAR + 7.231979)^2 - 4 *
            0.08554596 *  -1.4685783  * PAR *
            7.231979)
  ) / (2*-1.4685783) - 6.062066 
}

par_c3_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c3_25(PAR = 0:1200),
                        sample = "C3")

# 25˚C C4
curve_c4_25 <- function(PAR){
  (0.04480221 * PAR + 7.253526 - 
     sqrt((0.04480221 * PAR + 7.253526)^2 - 4 *
            0.04480221 *  -0.1268348  * PAR *
            7.253526)
  ) / (2*-0.1268348) - 5.742690 
}

par_c4_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_c4_25(PAR = 0:1200),
                        sample = "C4")

# 25˚C T5
curve_t5_25 <- function(PAR){
  (64.98305470 * PAR + -7695.054425 - 
     sqrt((64.98305470 * PAR + -7695.054425)^2 - 4 *
            64.98305470 *  0.9992626  * PAR *
            -7695.054425)
  ) / (2*0.9992626) - -7694.035508 
}

par_t5_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t5_25(PAR = 0:1200),
                        sample = "T5")

# 25˚C T6
curve_t6_25 <- function(PAR){
  (28.97634323 * PAR + -2403.518587 - 
     sqrt((28.97634323 * PAR + -2403.518587)^2 - 4 *
            28.97634323 *  0.9976309  * PAR *
            -2403.518587)
  ) / (2*0.9976309) - -2402.887763 
}

par_t6_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t6_25(PAR = 0:1200),
                        sample = "T6")

# 25˚C T3
curve_t3_25 <- function(PAR){
  (0.16153343 * PAR + 8.200713 - 
     sqrt((0.16153343 * PAR + 8.200713)^2 - 4 *
            0.16153343 *  -0.7207734  * PAR *
            8.200713)
  ) / (2*-0.7207734) - 6.896605 
}

par_t3_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t3_25(PAR = 0:1200),
                        sample = "T3")

# 25˚C T4
curve_t4_25 <- function(PAR){
  (0.06541427 * PAR + 5.308415 - 
     sqrt((0.06541427 * PAR + 5.308415)^2 - 4 *
            0.06541427 *  -0.2997131  * PAR *
            5.308415)
  ) / (2*-0.2997131) - 4.897868 
}

par_t4_25 <- data.frame(Lcuv = 0:1200,
                        curve = curve_t4_25(PAR = 0:1200),
                        sample = "T4")

# all 25 together
par.25 <- rbind(par_c1_25, par_c2_25, par_c3_25, par_c4_25, par_t3_25, par_t4_25, par_t5_25, par_t6_25)
par.25 <- par.25 %>% 
  mutate(treatment = ifelse(grepl("C", sample), "control", "treatment")) %>% 
  mutate(temperature = 25)

lrc_short_25 <- lrc_short %>% 
  filter(temperature == 25)

# horizontal plots
plot_lrc_0 <- ggplot(par.0, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_0, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('0˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))

      
plot_lrc_5 <- ggplot(par.5, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_5, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('5˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))


plot_lrc_10 <- ggplot(par.10, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_10, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('10˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))


plot_lrc_15 <- ggplot(par.15, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_15, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('15˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))


plot_lrc_20 <- ggplot(par.20, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_20, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('20˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))

plot_lrc_25 <- ggplot(par.25, aes(x = Lcuv, y = curve, group = sample))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1.5, alpha = 0.7) +
  geom_point(data = lrc_short_25, aes(x = light, y = NP_DW, col = treatment), size = 3, alpha = 0.5) +
  geom_line(aes(color = treatment, linetype = treatment),
            size = 0.8) +
  # facet_wrap(~treatment) +
  ylab(label = expression(NP~(µmol~CO2~m^-2~s^-1))) +  
  xlab(label = expression(paste(
    "PPFD ", "(µmol photons ", "m"^-2, " s"^-1, ")"))) +
  theme_bw(base_size = 15) +
  scale_fill_manual(values = c('dark blue', 'red')) +
  scale_color_manual(values = c('dark blue', 'red')) +
  ggtitle('25˚C') +
  theme(legend.position = 'none', axis.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_y_continuous(limits = c(-7, 3)) +
  scale_x_continuous(breaks = c(0,1200), 
                     labels = c(0,1200))

# ALL LRCs together plot ----
LRC_plot <- plot_grid(plot_lrc_0, plot_lrc_5, plot_lrc_10, plot_lrc_15, plot_lrc_20, plot_lrc_25,
                      nrow = 1, align = "h")

LRC_plot <- annotate_figure(LRC_plot, left = textGrob(expression(paste(
  "NP (µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")")), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
  bottom = textGrob(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")"), gp = gpar(cex = 1.3)))

# Adding common legend at the bottom
LRC_plot <- LRC_plot + theme(legend.position = "bottom") +
  plot_annotation(title = "(a)")

ggsave("outcomes/figure_3a.png", plot = LRC_plot, width = 12, height = 5, dpi = 300)


### Extracted parameters from LRCs ----

# LIGHT COMPENSATION POINTS ----

# Option 1: Extract from code
#lrc_merged <- rbind(temp_0, temp_5, temp_10, temp_15, temp_20, temp_25) 

#lrc_merged <- lrc_merged %>% 
#  select(LCP, sample, temperature) %>% 
#  mutate(treatment_type = case_when(grepl("C", sample) ~ "control",
#                                    grepl("T", sample) ~ "treatment")) %>% 
#  na.omit()

# Option 2: Collide in Excel and use the csv
lcp_merged <- read.csv('data/lcp_wc/light_compensation_points.csv')

# See how much crossed zero
lcp_merged_summary <- lcp_merged %>% 
  mutate(treatment = as.factor(treatment),
         curve = as.factor(curve),
         temp = as_factor(temp),
         sample = as.factor(sample)) %>% 
  group_by(temp,treatment) %>%
  summarise(
    yes_count = sum(curve == "yes")
  )

ggplot(data = lcp_valid, aes(x = temp, y = lcp_valid$function., color = sample)) +
  geom_point(size = 5) +
  labs(x = "Temperature", y = "Light Compensation Point") +
  facet_wrap(~sample)

# How many of which temps and treatments crossed the zero line?
# 0 - All controls, none treatments
# 5 - All controls, all treatments
# 10 - All controls, 2/4 treatments (T6 and T3)
# 15 - All controls, all treatments
# 20 - All controls, 1/4 treatments (T6, T3, T4)
# 25 - All controls, 1/4 treatments

# testing t-test LCPs at 5˚C and at 15˚C because only these had all positive LCPs, could not do stats on other temp steps
lcp_only_15 <- lcp_merged %>% 
  filter(temp == 15)

lcp_only_5 <- lcp_merged %>% 
  filter(temp == 5)

t.test(function. ~ treatment,data = lcp_only_5)

# testing t-test for all that crossed 
lcp_only_yes <- lcp_merged %>% 
  filter(curve == "yes") %>% 
  mutate(treatment = case_when(grepl("control", treatment) ~ "Control",
                                    grepl("treatment", treatment) ~ "Treatment"))

# t-test
# mean in group Control mean in group Treatment 
#            222.4031                231.7662 
t.test(function. ~ treatment,data = lcp_only_yes)

average_per_treatment_lcp <- lcp_only_yes %>% 
  group_by(treatment) %>% 
  summarise(mean_value = mean(function.),
            SE = std.error(function.))

# Plots for LCPs
colours <- c("dark blue", "red")

# All temps together
LCP_boxplot_mean <- ggplot(lcp_only_yes, aes(x = treatment, y = function., fill = treatment)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  #facet_wrap(~temperature) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(b)") +
  ylab(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Indivisual temp steps
LCP_boxplot_split <- ggplot(lcp_only_yes, aes(x = treatment, y = function., fill = treatment)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  facet_wrap(~temp, ncol = 6) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ylab(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Light Saturation Points - extracted from modelled values ----
# careful with the negative LSPs in treatment ones, the percentage needs to be calculated from positive value and then converted back to negative
par.all <- rbind(par.0, par.5, par.10, par.15, par.20, par.25)

# splitting calculations for 90% to positive and negative
par_sum <- par.all %>%
  group_by(temperature, sample) %>%
  summarize(maxNP = max(curve)) %>%
  mutate(
    maxNP_modified = ifelse(maxNP >= 0, maxNP, -maxNP), # Replace negative values with their positive counterparts
    LSPx90 = ifelse(maxNP >= 0, maxNP_modified / 100 * 90, maxNP_modified * 100 / 90), # Perform different calculations for positive and negative values
    LSPx90_final = ifelse(maxNP < 0, -LSPx90, LSPx90) # Convert back to negative if original value was negative
  ) %>% 
  select(-maxNP_modified, -LSPx90)

# Need to get an intersection of the modelled 90% maxNP with light intensity
# Manuslly using values from par_sum
# 0˚C LSP 
c1_0_intersection <- 0.52530101
c2_0_intersection <-0.74583116
c3_0_intersection <- 0.45998857
c4_0_intersection <- 0.79416922
t3_0_intersection <- -0.60427459
t4_0_intersection <- -0.49350667
t5_0_intersection <- -0.15062781
t6_0_intersection <- -0.41114400

# Interpolate the light value corresponding to the 90th percentile of maxNP
c1_0_lsp <- approx(par_c1_0$curve, par_c1_0$Lcuv, xout = c1_0_intersection)$y
c2_0_lsp <- approx(par_c2_0$curve, par_c2_0$Lcuv, xout = c2_0_intersection)$y
c3_0_lsp <- approx(par_c3_0$curve, par_c3_0$Lcuv, xout = c3_0_intersection)$y
c4_0_lsp <- approx(par_c4_0$curve, par_c4_0$Lcuv, xout = c4_0_intersection)$y
t3_0_lsp <- approx(par_t3_0$curve, par_t3_0$Lcuv, xout = t3_0_intersection)$y
t4_0_lsp <- approx(par_t4_0$curve, par_t4_0$Lcuv, xout = t4_0_intersection)$y
t5_0_lsp <- approx(par_t5_0$curve, par_t5_0$Lcuv, xout = t5_0_intersection)$y
t6_0_lsp <- approx(par_t6_0$curve, par_t6_0$Lcuv, xout = t6_0_intersection)$y

# joining 0˚C LSP
lsp_0 <- par_sum %>% 
  filter(temperature == 0)
lsp_0$LSPcuv <- c(679.1011, 927.7867, 785.4524, 574.324, 173.2026, 122.9398, 120.249, 103.6289)
lsp_0$treatment_type <- c("control", "control", "control", "control",
                              "treatment", "treatment", "treatment", "treatment")

# Normal distribution
lsp_0 %>%
 group_by(treatment_type) %>%
 shapiro_test(LSPcuv)

# t-test - sig - 0.000212
lsp_0_test <- lsp_0 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# 5˚C LSP 
# Need to get an intersection of the modelled 90% maxNP with light intensity
c2_5_intersection <-0.84948519
c3_5_intersection <- 0.87456091
c4_5_intersection <- 1.10901668
t3_5_intersection <- 0.46570751
t4_5_intersection <- 0.17455251
t5_5_intersection <- 0.32706850
t6_5_intersection <- 0.33369073

# Interpolate the light value corresponding to the 95th percentile of maxNP
#c1_5_lsp <- approx(par.all$curve, par.all$Lcuv, xout = c1_5_intersection)$y
c2_5_lsp <- approx(par_c2_5$curve, par_c2_5$Lcuv, xout = c2_5_intersection)$y
c3_5_lsp <- approx(par_c3_5$curve, par_c3_5$Lcuv, xout = c3_5_intersection)$y
c4_5_lsp <- approx(par_c4_5$curve, par_c4_5$Lcuv, xout = c4_5_intersection)$y
t3_5_lsp <- approx(par_t3_5$curve, par_t3_5$Lcuv, xout = t3_5_intersection)$y
t4_5_lsp <- approx(par_t4_5$curve, par_t4_5$Lcuv, xout = t4_5_intersection)$y
t5_5_lsp <- approx(par_t5_5$curve, par_t5_5$Lcuv, xout = t5_5_intersection)$y
t6_5_lsp <- approx(par_t6_5$curve, par_t6_5$Lcuv, xout = t6_5_intersection)$y

# joining 5˚C LSP
lsp_5 <- par_sum %>% 
  filter(temperature == 5)
lsp_5$LSPcuv <- c(591.645, 545.2867, 630.396, 857.255, 434.8782, 820.5843, 186.3504)
lsp_5$treatment_type <- c("control", "control", "control",
                          "treatment", "treatment", "treatment", "treatment")   

# Normal distribution
lsp_5 %>%
  group_by(treatment_type) %>%
  shapiro_test(LSPcuv)

# t-test - not sig - 0.943
lsp_5_test <- lsp_5 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

t.test(LSPcuv ~ treatment_type, data = lsp_5)

# 10˚C LSP
c1_10_intersection <- 1.48251075
c2_10_intersection <-2.08692053
c3_10_intersection <- 1.61316647
c4_10_intersection <- 1.81511485
t3_10_intersection <- 0.86996252
t4_10_intersection <- -0.20501959
t5_10_intersection <- -0.08510308
t6_10_intersection <- 0.68890856

# Interpolate the light value corresponding to the 910th percentile of maxNP
c1_10_lsp <- approx(par_c1_10$curve, par_c1_10$Lcuv, xout = c1_10_intersection)$y
c2_10_lsp <- approx(par_c2_10$curve, par_c2_10$Lcuv, xout = c2_10_intersection)$y
c3_10_lsp <- approx(par_c3_10$curve, par_c3_10$Lcuv, xout = c3_10_intersection)$y
c4_10_lsp <- approx(par_c4_10$curve, par_c4_10$Lcuv, xout = c4_10_intersection)$y
t3_10_lsp <- approx(par_t3_10$curve, par_t3_10$Lcuv, xout = t3_10_intersection)$y
t4_10_lsp <- approx(par_t4_10$curve, par_t4_10$Lcuv, xout = t4_10_intersection)$y
t5_10_lsp <- approx(par_t5_10$curve, par_t5_10$Lcuv, xout = t5_10_intersection)$y
t6_10_lsp <- approx(par_t6_10$curve, par_t6_10$Lcuv, xout = t6_10_intersection)$y

# joining 10˚C LSP
lsp_10 <- par_sum %>% 
  filter(temperature == 10)
lsp_10$LSPcuv <- c(545.3148, 613.3777, 764.4277, 879.2364, 750.3435, 733.1598, 1128.032, 601.7018)
lsp_10$treatment_type <- c("control", "control", "control", "control", 
                          "treatment", "treatment", "treatment", "treatment")   

# Normal distribution
lsp_10 %>%
  group_by(treatment_type) %>%
  shapiro_test(LSPcuv)

# t-test - not sig - 0.478
lsp_10_test <- lsp_10 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# 15˚C LSP 
c1_15_intersection <- 1.89849411
c2_15_intersection <-2.43383255
c3_15_intersection <- 1.28965359
c4_15_intersection <- 2.18726809
t3_15_intersection <- 1.40225466
t4_15_intersection <- 0.37298058
t5_15_intersection <- 0.36257114
t6_15_intersection <- 0.31211056

# Interpolate the light value corresponding to the 915th percentile of maxNP
c1_15_lsp <- approx(par_c1_15$curve, par_c1_15$Lcuv, xout = c1_15_intersection)$y
c2_15_lsp <- approx(par_c2_15$curve, par_c2_15$Lcuv, xout = c2_15_intersection)$y
c3_15_lsp <- approx(par_c2_15$curve, par_c2_15$Lcuv, xout = c3_15_intersection)$y
c4_15_lsp <- approx(par_c2_15$curve, par_c2_15$Lcuv, xout = c4_15_intersection)$y
t3_15_lsp <- approx(par_t3_15$curve, par_t3_15$Lcuv, xout = t3_15_intersection)$y
t4_15_lsp <- approx(par_t4_15$curve, par_t4_15$Lcuv, xout = t4_15_intersection)$y
t5_15_lsp <- approx(par_t5_15$curve, par_t5_15$Lcuv, xout = t5_15_intersection)$y
t6_15_lsp <- approx(par_t6_15$curve, par_t6_15$Lcuv, xout = t6_15_intersection)$y

# joining 15˚C LSP
lsp_15 <- par_sum %>% 
  filter(temperature == 15)
lsp_15$LSPcuv <- c(563.2592, 755.9286, 234.3031, 548.0303, 645.2373, 947.7188, 1092.425, 722.693)
lsp_15$treatment_type <- c("control", "control", "control", "control", 
                           "treatment", "treatment", "treatment", "treatment")   

# Normal distribution
lsp_15 %>%
  group_by(treatment_type) %>%
  shapiro_test(LSPcuv)

# t-test - not sig - 0.0708
lsp_15_test <- lsp_15 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()


# 20˚C LSP 
c1_20_intersection <- 1.45439205
c2_20_intersection <-2.26764913
c3_20_intersection <- 0.91813888
c4_20_intersection <- 1.13637504
t3_20_intersection <- 0.70159642
t5_20_intersection <- -0.50920053
t6_20_intersection <- -0.66119455

# Interpolate the light value corresponding to the 920th percentile of maxNP
c1_20_lsp <- approx(par_c1_20$curve, par_c1_20$Lcuv, xout = c1_20_intersection)$y
c2_20_lsp <- approx(par_c2_20$curve, par_c2_20$Lcuv, xout = c2_20_intersection)$y
c3_20_lsp <- approx(par_c3_20$curve, par_c3_20$Lcuv, xout = c3_20_intersection)$y
c4_20_lsp <- approx(par_c4_20$curve, par_c4_20$Lcuv, xout = c4_20_intersection)$y
t3_20_lsp <- approx(par_t3_20$curve, par_t3_20$Lcuv, xout = t3_20_intersection)$y
t5_20_lsp <- approx(par_t5_20$curve, par_t5_20$Lcuv, xout = t5_20_intersection)$y
t6_20_lsp <- approx(par_t6_20$curve, par_t6_20$Lcuv, xout = t6_20_intersection)$y

# joining 20˚C LSP
lsp_20 <- par_sum %>% 
  filter(temperature == 20)
lsp_20$LSPcuv <- c(737.5368, 719.2237, 1072.529, 987.6947, 980.1054, 1029.205, 953.2849)
lsp_20$treatment_type <- c("control", "control", "control", "control", 
                          "treatment", "treatment", "treatment")   

# Normal distribution
lsp_20 %>%
  group_by(treatment_type) %>%
  shapiro_test(LSPcuv)

# t-test - not sig - 0.357
lsp_20_test <- lsp_20 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# 25˚C LSP 
c1_25_intersection <- 0.15768390
c2_25_intersection <-0.68116750
c3_25_intersection <- 0.15678018
c4_25_intersection <- 0.50923718
t3_25_intersection <- 0.68579935
t4_25_intersection <- -0.01233609
t5_25_intersection <- -1.69848178
t6_25_intersection <- -1.11009341

# Interpolate the light value corresponding to the 925th percentile of maxNP
c1_25_lsp <- approx(par_c1_25$curve, par_c1_25$Lcuv, xout = c1_25_intersection)$y
c2_25_lsp <- approx(par_c2_25$curve, par_c2_25$Lcuv, xout = c2_25_intersection)$y
c3_25_lsp <- approx(par_c3_25$curve, par_c3_25$Lcuv, xout = c3_25_intersection)$y
c4_25_lsp <- approx(par_c4_25$curve, par_c4_25$Lcuv, xout = c4_25_intersection)$y
t3_25_lsp <- approx(par_t3_25$curve, par_t3_25$Lcuv, xout = t3_25_intersection)$y
t4_25_lsp <- approx(par_t4_25$curve, par_t4_25$Lcuv, xout = t4_25_intersection)$y
t5_25_lsp <- approx(par_t5_25$curve, par_t5_25$Lcuv, xout = t5_25_intersection)$y
t6_25_lsp <- approx(par_t6_25$curve, par_t6_25$Lcuv, xout = t6_25_intersection)$y

# joining 25˚C LSP 
lsp_25 <- par_sum %>% 
  filter(temperature == 25)
lsp_25$LSPcuv <- c(1070.28, 1014.96, 1174.237, 1121.055, 1037.475, 1196.132, 870.5187, 902.9255)
lsp_25$treatment_type <- c("control", "control", "control", "control", 
                           "treatment", "treatment", "treatment", "treatment")   

# Normal distribution
lsp_25 %>%
  group_by(treatment_type) %>%
  shapiro_test(LSPcuv)

# t-test - not sig - 0.296
lsp_25_test <- lsp_25 %>%
  t_test(LSPcuv ~ treatment_type, var.equal = TRUE) %>%
  add_significance()

# All LSPs together ----
all_lsp <- rbind(lsp_0,lsp_5,lsp_10,lsp_15,lsp_20,lsp_25) %>% 
  mutate(treatment_type = case_when(grepl("control", treatment_type) ~ "Control",
                                    grepl("treatment", treatment_type) ~ "Treatment"))

t.test(LSPcuv ~ treatment_type, data = all_lsp)  # p = 0.5406 (not sig.)

# Control              762. +- 49.1
#2 Treatment            713. +- 74.0
average_lsp_per_treatment <- all_lsp %>% 
  group_by(treatment_type) %>% 
  summarise(mean_value = mean(LSPcuv),
            SE = std.error(LSPcuv))

# Plots
# All together
LSP_boxplot_mean <- ggplot(all_lsp, aes(x = treatment_type, y = LSPcuv, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  #facet_wrap(~temperature) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(c)") +
  ylab(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(), legend.text = element_text(size = 12))

# Separate temps
LSP_boxplot_split <- ggplot(all_lsp, aes(x = treatment_type, y = LSPcuv, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  facet_wrap(~temperature, ncol = 6) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ylab(expression("PPFD (μmol photons m"^"-2"*" s"^"-1"~")")) +
  theme_bw(base_size = 15) +
  theme(axis.title.x = element_blank(), legend.position = "none", legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# MEAN MAXIMAL RATES OF PS from Light curves ----
max_rates_lrc <- par.all %>%
  group_by(temperature, sample, treatment) %>% 
  summarise(max_np_curve = max(curve))

max_rates_lrc %>%
  group_by(treatment) %>%
  shapiro_test(max_np_curve)

t.test(max_np_curve ~ treatment, data = max_rates_lrc)   

# control       1.33  +- 0.158
# treatment     0.0911 +- 0.144
average_per_treatment_max_np <- max_rates_lrc %>% 
  group_by(treatment) %>% 
  summarise(mean_value = mean(max_np_curve),
            SE = std.error(max_np_curve))

maxNP_boxplot_mean <- ggplot(all_lsp, aes(x = treatment_type, y = maxNP, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  #facet_wrap(~temperature) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ggtitle("(d)") +
  ylab(label = expression(paste(
    "NP (µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")"))) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "none", axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))

maxNP_boxplot_split <- ggplot(all_lsp, aes(x = treatment_type, y = maxNP, fill = treatment_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.4, outlier.shape = NA) +
  facet_wrap(~temperature, ncol = 6) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 3) +
  
  # Set colors manually
  scale_fill_manual(values = colours) +
  
  # Other styling options
  ylab(label = expression(paste(
    "NP (µmol ", "CO"[2], " g"^"−1", " s"^"−1", ")"))) + 
  theme_bw(base_size = 15) +
  theme(legend.position = "none", axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 12))


#  extraxted parameters plot - Figure 3bcd
LRC_extract_plot <- ggarrange(LCP_boxplot_mean, LSP_boxplot_mean, maxNP_boxplot_mean, 
                      common.legend = TRUE, legend = "none")

LRC_extract_plot <- {
  LCP_boxplot_mean + LSP_boxplot_mean + maxNP_boxplot_mean + plot_layout(ncol = 3)
} +
  plot_layout(nrow = 1, heights = c(5, 1))

ggsave("outcomes/figure_3bcd.png", plot = LRC_extract_plot, width = 12, height = 5, dpi = 300)


# Split extracted parameters, split by temperature steps - Figure S9
split_temps <- {
  LCP_boxplot_split + LSP_boxplot_split + maxNP_boxplot_split + plot_layout(ncol = 1)
} 

ggsave("outcomes/figure_S9.png", plot = split_temps, width = 10, height = 10, dpi = 300)




