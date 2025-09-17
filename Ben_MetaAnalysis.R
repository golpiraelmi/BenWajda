library(tidyverse)
library(meta)
library(metafor)
library(readxl)
library(dplyr)

setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda/')

# Get sheet names excluding first two
all_sheets <- excel_sheets('data/GLP_Reformat.xlsx')
sheet_names <- setdiff(all_sheets[-c(1,2)], "LOS")

for(sheet_name in sheet_names){
  print(sheet_name)
  
  data <- read_excel('data/GLP_Reformat.xlsx', sheet = sheet_name)
  data <- na.omit(data)
  
  # Run fixed effect first to get I2 and Tau2
  m.fixed <- metabin(
    event.e = n_events_GLP_RA, 
    n.e = n_total_GLP_RA,
    event.c = n_events_non_GLP_RA,
    n.c = n_total_non_GLP_RA,
    studlab = StudyID,
    data = data,
    sm = "OR",
    method = "MH",
    MH.exact = FALSE,
    common = TRUE,
    random = FALSE,
    method.random.ci = TRUE,
    subgroup = Surgery,
    title = sheet_name,
    incr = 0.5,            # continuity correction
    allstudies = TRUE
  )
  
  # Decide model type
  I2 <- m.fixed$I2
  Tau2 <- m.fixed$tau^2
  
  if(I2 > .50 ){
    cat("Switching to random effects for sheet:", sheet_name, "(I2 =", I2, "Tau2 =", Tau2, ")\n")
    # Create a new random effects model
    m.bin <- metabin(
      event.e = n_events_GLP_RA, 
      n.e = n_total_GLP_RA,
      event.c = n_events_non_GLP_RA,
      n.c = n_total_non_GLP_RA,
      studlab = StudyID,
      data = data,
      sm = "OR",
      method = "MH",
      MH.exact = FALSE,
      common = FALSE,     # random effect
      random = TRUE,
      method.tau = "DL",
      method.ci = "knha" ,      # Knapp-Hartung CI adjustment
      method.random.ci = TRUE,
      subgroup  = Surgery,
      title = sheet_name,
      incr = 0.5,            # continuity correction
      allstudies = TRUE
    )
  } else {
    cat("Using fixed effect for sheet:", sheet_name, "(I2 =", I2, "Tau2 =", Tau2, ")\n")
    m.bin <- m.fixed
  }
  
  # # Print summary
  # print(summary(m.bin))

  
  # Save forest plot
  png(file = paste0("Results/", sheet_name, "_Forest.png"), 
      width = 2800, height = 2200, res = 300)
  forest(
    m.bin,
    label.e = 'GLP RA',
    label.c = 'Non GLP RA',
    fontsize = 10,
    col.common='black',
    col.square ='black',
    col.diamond.common ='#F5EB27',
    col.diamond.random ='blue',
    layout = "RevMan5"
  )
  dev.off()
  
  cat("Finished sheet:", sheet_name, "\n")
}

######################################## Publication bias ######################################## 
# Create overall model without subgroups
for(sheet_name in sheet_names){
  print(sheet_name)
  print('======================================')
  
  data <- read_excel('data/GLP_Reformat.xlsx', sheet = sheet_name)
  data <- na.omit(data)
  
  # Fit model without subgroups (for bias tests)
  m.overall <- metabin(
    event.e = n_events_GLP_RA, 
    n.e = n_total_GLP_RA,
    event.c = n_events_non_GLP_RA,
    n.c = n_total_non_GLP_RA,
    studlab = StudyID,
    data = data,
    sm = "OR",
    method = "MH",
    MH.exact = FALSE,
    common = FALSE,
    random = TRUE,
    method.tau = "DL",
    method.ci = "knha",
    incr = 0.5,
    allstudies = TRUE
  )
  
  # Publication bias (Egger’s test)
  print ("Egger's Test")
  print('======================================')
  if (nrow(data) >= 7) {
    print(metabias(m.overall, k.min = 7, method.bias = "linreg"))
    
    # Trim-and-Fill
    tf <- trimfill(m.overall)
    print('======================================')
    print ("Trim-and-Fill")
    print('======================================')
    print(tf)
    
    # Funnel plot
    png(file = paste0("Results/", sheet_name, "_Funnel.png"), 
        width = 2400, height = 1800, res = 300)
    
    funnel(
      tf,                     # use Trim-and-Fill object
      xlab = "Log Odds Ratio",
      ylab = "Standard Error",
      studlab = TRUE,          # show study labels
      pch = 16,                # dot style
      cex = 0.6,               # smaller dots
      cex.studlab = 0.5,       # smaller font for labels
      back = TRUE
    )
    
    dev.off()
  }
}

##############LOS
data <- read_excel('data/LOS.xlsx')
data$Surgery <- factor(data$Surgery) 
data
# ---- A. Ignore studies with missing SD ----
data_ignore <- na.omit(data)  # remove rows with NA

meta_ignore <- metacont(
  n.e = n_GLP_RA,
  mean.e = mean_GLP_RA,
  sd.e = sd_GLP_RA,
  n.c = n_non_GLP_RA,
  mean.c = mean_non_GLP_RA,
  sd.c = sd_non_GLP_RA,
  data = data_ignore,
  studlab = StudyID,
  sm = "MD",           # Mean Difference
  method.smd = "Cohen",
  random = TRUE,  # Random-effects model
  common = FALSE,  
  method.random.ci = TRUE ,
  subgroup = Surgery,
)

summary(meta_ignore)

# Save forest plot
png(file = paste0("Results/", 'LOS1', "_Forest.png"), 
    width = 2800, height = 1800, res = 300)
forest(meta_ignore,
       main="Meta-analysis ignoring studies with missing SD",
       label.e = 'GLP RA',
       label.c = 'Non GLP RA',
       col.diamond.random ='blue',
       col.common='black',
       col.square ='black',
       col.square.lines='black',
       digits = 1,
       digits.sd = 1,
       layout = "RevMan5")
dev.off()

# ---- B. Impute missing SDs ----
data <- read_excel('data/LOS.xlsx')
data$Surgery <- factor(data$Surgery)

# Compute mean SD for each group (ignoring NA)
mean_sd_GLP_RA <- mean(data$sd_GLP_RA, na.rm = TRUE)
mean_sd_non_GLP_RA <- mean(data$sd_non_GLP_RA, na.rm = TRUE)

# Replace missing SDs with mean SD
data$sd_GLP_RA[is.na(data$sd_GLP_RA)] <- mean_sd_GLP_RA
data$sd_non_GLP_RA[is.na(data$sd_non_GLP_RA)] <- mean_sd_non_GLP_RA

# ---- Meta-analysis with imputed SDs ----
meta_imputed <- metacont(
  n.e = n_GLP_RA,
  mean.e = mean_GLP_RA,
  sd.e = sd_GLP_RA,
  n.c = n_non_GLP_RA,
  mean.c = mean_non_GLP_RA,
  sd.c = sd_non_GLP_RA,
  data = data,
  studlab = StudyID,
  sm = "MD",
  method.smd = "Cohen",
  random = TRUE,
  common = FALSE,
  method.random.ci = TRUE,
  subgroup = Surgery
)

summary(meta_imputed)

# ---- Save forest plot ----
png(file = paste0("Results/", 'LOS_ImputedSD_Forest.png'), 
    width = 2800, height = 1800, res = 300)
forest(meta_imputed,
       main="Meta-analysis with imputed SDs",
       label.e = 'GLP RA',
       label.c = 'Non GLP RA',
       col.diamond.random ='blue',
       col.common='black',
       col.square ='black',
       col.square.lines='black',
       digits = 1,
       digits.sd = 1,
       layout = "RevMan5")
dev.off()
