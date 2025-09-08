library(tidyverse)
library(meta)
library(metafor)
library(readxl)
library(dplyr)

setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda/')

# Get sheet names excluding first two
all_sheets <- excel_sheets('data/GLP_Reformat.xlsx')
sheet_names <- all_sheets[-c(1,2)]

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

# ######################################## Publication bias ######################################## 
# # Create overall model without subgroups
# for(sheet_name in sheet_names){
#   print(sheet_name)
#   
#   data <- read_excel('data/GLP_Reformat.xlsx', sheet = sheet_name)
#   data <- na.omit(data)
#   
#   
#   # Fit model without subgroups (for bias tests)
#   m.overall <- metabin(
#     event.e = n_events_GLP_RA, 
#     n.e = n_total_GLP_RA,
#     event.c = n_events_non_GLP_RA,
#     n.c = n_total_non_GLP_RA,
#     studlab = StudyID,
#     data = data,
#     sm = "OR",
#     method = "MH",
#     MH.exact = FALSE,
#     common = FALSE,
#     random = TRUE,
#     method.tau = "DL",
#     method.ci = "knha",
#     incr = 0.5,
#     allstudies = TRUE
#   )
#   
#   # Publication bias (Egger’s test)
#   if (nrow(data) >= 7) {
#     print(metabias(m.overall, k.min = 7, method.bias = "linreg"))
#     
#     # Funnel plot
#     png(file = paste0("Results/", sheet_name, "_Funnel.png"), 
#         width = 2400, height = 1800, res = 300)
#     funnel(
#       m.overall,
#       xlab = "Log Odds Ratio",
#       ylab = "Standard Error",
#       studlab = TRUE,
#       pch = 16,                     # dot style
#       cex = 0.6,                     # smaller dots
#       cex.studlab = 0.5,  # smaller font for labels     
#       # contour = c(0.95),
#       # col.contour = c("#d6f5f5"),
#       back = TRUE,
#       
#       
#     )
#     # legend("topright", legend = c("p < 0.05"),
#     #        fill = c("#d6f5f5"))
#     
#     dev.off()
#   }}
