# ==============================
#  GLP-1 RA Meta-Analysis Pipeline
#  Forest + Funnel (Subgroup-Colored) + Bias Tests
# ==============================

library(tidyverse)
library(meta)
library(metafor)
library(readxl)
library(dplyr)
library(ggplot2)

# Set working directory
setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda/')

# Create Results folder if not exists
if (!dir.exists("Results")) dir.create("Results")

# ------------------------------------------------------------------
# 1) Forest plots (fixed/random based on I² > 50%)
# ------------------------------------------------------------------
all_sheets  <- excel_sheets('data/GLP-1 Stats w Ranson.xlsx')
sheet_names <- setdiff(all_sheets[-c(1,2)], "LOS")

cat("Starting forest plot generation...\n")

for (sheet_name in sheet_names) {
  cat("\n=== Processing sheet:", sheet_name, "===\n")
  
  data <- read_excel('data/GLP-1 Stats w Ranson.xlsx', sheet = sheet_name) 
  
  data <-data %>%filter(StudyID != "Ranson 2025") %>%  # Exclude Ranson 2025
    na.omit()
  

  
  if (nrow(data) == 0) {
    cat("  No data after na.omit(). Skipping.\n")
    next
  }
  
  ## ---- Fixed-effect model to check heterogeneity ----
  m.fixed <- metabin(
    event.e = n_events_GLP_RA,   n.e = n_total_GLP_RA,
    event.c = n_events_non_GLP_RA, n.c = n_total_non_GLP_RA,
    studlab = StudyID, data = data,
    sm = "OR", method = "MH", MH.exact = FALSE,
    common = TRUE, random = FALSE,
    subgroup = Surgery, title = sheet_name,
    incr = 0.5, allstudies = TRUE
  )
  
  I2   <- m.fixed$I2
  Tau2 <- m.fixed$tau^2
  
  ## ---- Choose model: Random if I² > 50% ----
  if (I2 > 0.50) {
    cat("  -> Using RANDOM-EFFECTS (I² =", round(I2*100, 1), "%, τ² =", round(Tau2, 4), ")\n")
    m.bin <- metabin(
      event.e = n_events_GLP_RA,   n.e = n_total_GLP_RA,
      event.c = n_events_non_GLP_RA, n.c = n_total_non_GLP_RA,
      studlab = StudyID, data = data,
      sm = "OR", method = "MH", MH.exact = FALSE,
      common = FALSE, random = TRUE,
      method.tau = "DL", method.ci = "knha",
      subgroup = Surgery, title = sheet_name,
      incr = 0.5, allstudies = TRUE
    )
  } else {
    cat("  -> Using FIXED-EFFECT (I² =", round(I2*100, 1), "%, τ² =", round(Tau2, 4), ")\n")
    m.bin <- m.fixed
  }
  
  ## ---- Save Forest Plot ----
  png(file = file.path("Results", paste0(sheet_name, "_Forest.png")),
      width = 2800, height = 2200, res = 300)
  forest(m.bin,
         label.e = "GLP-1 RA", label.c = "Non GLP-1 RA",
         fontsize = 10,
         col.square = "black", col.diamond.common = "#F5EB27",
         col.diamond.random = "blue", layout = "RevMan5",
         print.tau2 = TRUE, print.I2 = TRUE)
  dev.off()
  
  cat("  Forest plot saved.\n")
}

# ------------------------------------------------------------------
# 2) Publication Bias: Egger's, Trim-and-Fill + Funnel Plot (Subgroup-Colored)
# ------------------------------------------------------------------
cat("\nStarting publication bias analysis...\n")

for (sheet_name in sheet_names) {
  cat("\n======================================\n")
  cat("Publication Bias Analysis:", sheet_name, "\n")
  cat("======================================\n")
  
  data <- read_excel('data/GLP-1 Stats w Ranson.xlsx', sheet = sheet_name) %>%
    na.omit()
  
  if (nrow(data) == 0) {
    cat("  No data. Skipping.\n")
    next
  }
  
  ## ---- Overall Random-Effects Model (no subgroups) ----
  m.overall <- metabin(
    event.e = n_events_GLP_RA,   n.e = n_total_GLP_RA,
    event.c = n_events_non_GLP_RA, n.c = n_total_non_GLP_RA,
    studlab = StudyID, data = data,
    sm = "OR", method = "MH", MH.exact = FALSE,
    common = FALSE, random = TRUE,
    method.tau = "DL", method.ci = "knha",
    incr = 0.5, allstudies = TRUE
  )
  
  ## ---- Egger’s Test ------------------------------------------------
  cat("\nEgger's Test for Funnel Plot Asymmetry\n")
  cat("--------------------------------------\n")
  if (nrow(data) >= 7) {
    egger_test <- metabias(m.overall, k.min = 7, method.bias = "linreg")
    print(egger_test)
  } else {
    cat("  < 7 studies → Egger’s test skipped.\n")
  }
  
  ## ---- Trim-and-Fill + Funnel Plot ---------------------------------
  if (nrow(data) >= 3) {  # Need at least 3 for funnel
    cat("\nTrim-and-Fill Analysis\n")
    cat("--------------------------------------\n")
    tf <- trimfill(m.overall)
    print(tf)
    
    ## Extract observed studies
    n_obs <- nrow(data)
    obs_df <- data.frame(
      TE      = tf$TE[1:n_obs],
      seTE    = tf$seTE[1:n_obs],
      studlab = tf$studlab[1:n_obs],
      Surgery = factor(data$Surgery),
      type    = "Observed"
    )
    
    ## Extract imputed studies (if any)
    if (tf$k0 > 0) {
      imp_df <- data.frame(
        TE      = tf$TE[(n_obs + 1):length(tf$TE)],
        seTE    = tf$seTE[(n_obs + 1):length(tf$seTE)],
        studlab = paste("Imputed", 1:tf$k0),
        Surgery = NA,
        type    = "Imputed"
      )
      funnel_df <- rbind(obs_df, imp_df)
    } else {
      funnel_df <- obs_df
    }
    
    ## Pooled estimate
    pooled_log <- tf$TE.random
    
    ## ---- Build 95% pseudo-confidence funnel lines ----
    max_se <- max(funnel_df$seTE, na.rm = TRUE) * 1.2
    contour_se <- seq(0.01, max_se, length.out = 200)
    lower_ci <- pooled_log - 1.96 * contour_se
    upper_ci <- pooled_log + 1.96 * contour_se
    
    contour_lines <- rbind(
      data.frame(se = contour_se, bound = lower_ci, side = "lower"),
      data.frame(se = contour_se, bound = upper_ci, side = "upper")
    )
    
    ## ---- ggplot2 Funnel Plot 
    p <- ggplot() +
      
      # 95% pseudo-confidence contour
      geom_line(data = contour_lines, 
                aes(x = bound, y = se, group = side),
                color = "gray60", linetype = "dashed", size = 0.7) +
      
      # Observed studies (colored by Surgery)
      geom_point(data = subset(funnel_df, type == "Observed"),
                 aes(x = TE, y = seTE, color = Surgery),
                 size = 3.2, alpha = 0.9) +
      
      # Study labels
      geom_text(data = subset(funnel_df, type == "Observed"),
                aes(x = TE, y = seTE, label = studlab),
                hjust = 0, nudge_x = 0.05, size = 3, color = "black", fontface = "plain") +
      
      # Imputed studies (red, if any)
      geom_point(data = subset(funnel_df, type == "Imputed"),
                 aes(x = TE, y = seTE), color = "red", size = 3, shape = 5) +
      
      # Pooled effect line
      geom_vline(xintercept = pooled_log, color = "black", size = 0.8) +
      
      # Scales
      scale_y_reverse(breaks = scales::pretty_breaks(n = 6),
                      limits = c(max_se, 0)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      
      # Labels
      labs(
        x = "Log Odds Ratio (GLP-1 RA vs. Non-GLP-1 RA)",
        y = "Standard Error",
        title = paste0("Funnel Plot – ", sheet_name),
        color = "Surgery Type",
        caption = if(tf$k0 > 0) paste(tf$k0, "studies imputed (red diamonds)") else "No studies imputed"
      ) +
      
      # Theme
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 10, color = "gray50"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(override.aes = list(size = 4)))
    
    ## Save
    ggsave(filename = file.path("Results", paste0(sheet_name, "_Funnel_Subgroup.png")),
           plot = p, width = 9, height = 7, dpi = 300)
    
    cat("  Funnel plot saved (with correct funnel shape & subgroup colors)\n")
    
  } else {
    cat("  < 3 studies → Funnel plot skipped.\n")
  }
}

cat("\n=== ALL ANALYSES COMPLETED SUCCESSFULLY ===\n")

  
##############LOS
setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda/')
data <- read_excel('data/GLP-1 Stats w Ranson.xlsx', sheet='LOS (with SD)')
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
data <- read_excel('data/GLP-1 Stats w Ranson.xlsx', sheet='LOS (with SD)')
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
