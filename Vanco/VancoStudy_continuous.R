library(meta)
library(metafor)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

############################### Continuous outcome meta-analysis
setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda')

file_path <- "data/Vanco Stats.xlsx"
data <- read_excel(file_path, sheet = 'Serum Vanco Levels 2hr')

out_dir <- "Results_Vanco"
dir.create(out_dir, showWarnings = FALSE)




# Continuous outcome meta-analysis
m <- metamean(
  n = n_total,
  mean = vanco_mean,
  sd = vanco_sd,
  studlab = StudyID,
  data = data,
  sm = "MRAW",          # pooled raw mean
  method.tau = "REML",  # random-effects estimator
  common = FALSE,
  random = TRUE,
  subgroup = Group
)

# View results
summary(m)

# Forest plot with subgroup analysis
png(file.path(out_dir, paste0("Serum Vanco Levels 2hr", "_Forest.png")), width = 3000, height = 3400, res = 300)
forest(
  m,
  bylab = "Group",
  print.subgroup.name = TRUE,
  print.I2 = TRUE,
  print.tau2 = TRUE,
  layout = "RevMan5",
  col.diamond.random = "blue",
  col.diamond.lines = "white",
  col.square = "black",
  col.square.lines = "white",  
  col.study.labels = "black",
  colgap.left = "1cm",
  fontsize = 12,
  subgroup.hetstat = FALSE
)
dev.off()

# ---- Leave-one-out ----
loo <- metainf(m)

png(file.path(out_dir, paste0("Serum_Vanco_Levels_2hr", "_LeaveOneOut.png")),
    width = 3000, height = 2500, res = 300)

forest(loo,
       layout = "RevMan5",
       # xlim = c(0, 0.04),
       col.square = "black",
       col.diamond = "blue",
       col.square.lines = "white",  
       fontsize = 10)

dev.off()

metabias(m, method.bias = "linreg")
trimfill(m)

# ---- Funnel plot data ----
funnel_df <- data.frame(
  study = m$studlab,
  effect = m$TE,
  se = m$seTE,
  group = data$Group
)

pooled <- m$TE.random
se_seq <- seq(0, max(funnel_df$se, na.rm = TRUE), length.out = 100)

funnel_limits <- data.frame(
  se = se_seq,
  upper = pooled + 1.96 * se_seq,
  lower = pooled - 1.96 * se_seq
)

# ---- Funnel plot ----
png(file.path(out_dir, paste0('Serum_Vanco_Levels_2hr', "_Funnel.png")),
    width = 2000, height = 2000, res = 300)

print(
  ggplot(funnel_df, aes(x = effect, y = se)) +
    
    geom_line(data = funnel_limits, aes(x = upper, y = se),
              linetype = "dashed", linewidth = 0.5) +
    geom_line(data = funnel_limits, aes(x = lower, y = se),
              linetype = "dashed", linewidth = 0.5) +
    
    geom_vline(xintercept = pooled, linetype = "dashed", linewidth = 0.5) +
    
    geom_point(aes(color = group), size = 1, alpha = 0.8) +
    
    scale_y_reverse(sec.axis = dup_axis()) +
    scale_x_continuous(sec.axis = dup_axis()) +
    
    labs(
      x = "Mean (µg/ml)",
      y = "Standard Error",
      title = paste0("Funnel Plot - ", sheet)
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      axis.text.x.top = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank()
    ) +
    
    geom_text_repel(aes(label = study), size = 2, max.overlaps = 10)
)

dev.off()

funnel(m)
