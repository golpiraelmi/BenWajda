library(meta)
library(metafor)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd('/Users/golpira/Python/University of Calgary/UofC-Git/Ben Wajda')
file_path <- "data/Vanco Stats_formatted.xlsx"
sheets <- excel_sheets(file_path)  

# exclude_sheets <- c("VTE Events 90d", "Reoperation Rate 90d","SSI 90d")
# sheets <- setdiff(excel_sheets(file_path), exclude_sheets)


out_dir <- "Results_Vanco"
dir.create(out_dir, showWarnings = FALSE)

for (sheet in sheets) {
  
  message("Processing sheet: ", sheet)
  
  # ---- Read sheet ----
  data <- read_excel(file_path, sheet = sheet)
  
  # Clean sheet name for filenames
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", sheet)
  
  # ---- Meta-analysis ----
  meta_analysis <- metaprop(
    event = data$n_events,
    n = data$n_total,
    studlab = data$StudyID,
    subgroup = data$Group,
    method = "Inverse",
    sm = "PLOGIT", 
    incr = 0.5,
    method.incr = "if0all",
    common = FALSE,
    random = TRUE
  )
  
  # ---- Forest plot ----
  png(file.path(out_dir, paste0(safe_name, "_Forest.png")),
      width = 3000, height = 3400, res = 300)
  
  forest(meta_analysis,
         backtransf = TRUE,
         # xlim = c(0.0, 0.2),
         cex.lab = 0.7,
         layout = "RevMan5",
         mar = c(0, 0, 0, 0),
         col.diamond.random = "blue",
         col.diamond.lines = "white",
         col.square = "black",
         col.square.lines = "white",  
         col.study.labels = "black",
         colgap.left = "1cm",
         fontsize = 12,
         subgroup.hetstat = FALSE)
  
  dev.off()
  
  # ---- Leave-one-out ----
  loo <- metainf(meta_analysis)
  
  png(file.path(out_dir, paste0(safe_name, "_LeaveOneOut.png")),
      width = 3000, height = 2500, res = 300)
  
  forest(loo,
         layout = "RevMan5",
         backtransf = TRUE,
         # xlim = c(0, 0.04),
         col.square = "black",
         col.diamond = "blue",
         col.square.lines = "white",  
         fontsize = 10)
  
  dev.off()
  
  # ---- Bias tests ----
  metabias(meta_analysis, method.bias = "linreg")
  trimfill(meta_analysis)
  
  # ---- Funnel plot data ----
  funnel_df <- data.frame(
    study = meta_analysis$studlab,
    effect = meta_analysis$TE,
    se = meta_analysis$seTE,
    group = data$Group
  )
  
  pooled <- meta_analysis$TE.random
  se_seq <- seq(0, max(funnel_df$se, na.rm = TRUE), length.out = 100)
  
  funnel_limits <- data.frame(
    se = se_seq,
    upper = pooled + 1.96 * se_seq,
    lower = pooled - 1.96 * se_seq
  )
  
  # ---- Funnel plot ----
  png(file.path(out_dir, paste0(safe_name, "_Funnel.png")),
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
        x = "Logit Proportion (PLO scale)",
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
      
      geom_text_repel(aes(label = study), size = 2, max.overlaps = Inf)
  )
  
  dev.off()
}

