# Run here: Rscript make_forecasting_results.R
# Defaults read the frozen scores; optional CLI args select another score file/output directory.
# Tables report LPS differences and MASE ratios relative to EXH.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
score_file <- if (length(args) >= 1) args[1] else "inputs/origin_scores.csv"
output_dir <- if (length(args) >= 2) args[2] else "outputs"

origin <- read.csv(score_file)
origin <- origin[order(origin$model, origin$origin), ]
regimes <- 2:15
ms_models <- c(paste0("sthmsh", regimes), paste0("stmsh", regimes))
models <- c(ms_models, "sv", "exh")
label <- function(m) ifelse(m == "sv", "SV", ifelse(m == "exh", "EXH",
  sub("^sthmsh([0-9]+)$", "HMSH(\\1)", sub("^stmsh([0-9]+)$", "MSH(\\1)", m))))
family_of <- function(m) ifelse(startsWith(m, "sthmsh"), "HMSH",
  ifelse(startsWith(m, "stmsh"), "MSH", toupper(m)))
regimes_of <- function(m) suppressWarnings(as.integer(sub("^[a-z]+", "", m)))

# All models use the same 117 one-month-ahead forecasts.
ref <- origin[origin$model == "exh", ]
stopifnot(nrow(origin) == 3510,
          min(ref$target) == "2016-02-01", max(ref$target) == "2025-10-01")
for (m in models) {
  d <- origin[origin$model == m, ]
  stopifnot(nrow(d) == 117, identical(d$origin, ref$origin), identical(d$target, ref$target))
}

summary_all <- do.call(rbind, lapply(models, function(m) {
  x <- origin[origin$model == m, ]
  data.frame(label = label(m), family = family_of(m), M = regimes_of(m), model = m,
             n_origins = nrow(x), LPS = mean(x$lps), MASE = mean(x$mase))
}))
exh <- summary_all[summary_all$model == "exh", ]
summary_all$delta_LPS_EXH <- summary_all$LPS - exh$LPS
summary_all$MASE_ratio_EXH <- summary_all$MASE / exh$MASE
selected <- sapply(c("HMSH", "MSH"), function(f) {
  x <- summary_all[summary_all$family == f, ]; x$model[which.max(x$LPS)]
})
best_mase <- sapply(c("HMSH", "MSH"), function(f) {
  x <- summary_all[summary_all$family == f, ]; x$model[which.min(x$MASE)]
})
summary_all$selected_by_LPS <- summary_all$model %in% selected
summary_all$best_MASE <- summary_all$model %in% best_mase

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
write.csv(summary_all[c("label", "M", "model", "LPS", "MASE",
                        "delta_LPS_EXH", "MASE_ratio_EXH")],
          file.path(output_dir, "results/model_summary.csv"), row.names = FALSE)

f3 <- function(x) formatC(x, format = "f", digits = 3)
row_of <- function(m, tab = summary_all) tab[tab$model == m, ]
table_cells <- function(m) {
  r <- row_of(m)
  values <- c(f3(r$delta_LPS_EXH), f3(r$MASE_ratio_EXH))
  bold <- c(r$selected_by_LPS, r$best_MASE)
  values[bold] <- paste0("\\textbf{", values[bold], "}")
  values
}

# Same caption, note, layout and bold selections as the manuscript table.
tab_main <- c(
  "\\begin{table}[!t]", "\\centering",
  "\\caption{Forecasting performance results: one-period-ahead density and point forecast measures for SVARs",
  "with various volatility specifications relative to that of SVAR with EXH.}",
  "\\label{tab:forecast_metrics}", "\\small", "\\setstretch{1.0}",
  "\\renewcommand{\\arraystretch}{0.92}", "",
  "\\begin{tabular*}{0.94\\textwidth}{@{\\extracolsep{\\fill}}lrrrr@{}}", "\\toprule",
  "& \\multicolumn{2}{c}{HMSH} & \\multicolumn{2}{c}{MSH}\\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(l){4-5}",
  "$M$ & $\\Delta$LPS & MASE & $\\Delta$LPS & MASE\\\\", "\\midrule",
  sapply(regimes, function(M) paste0(sprintf("%-2d", M), " & ",
    paste(c(table_cells(paste0("sthmsh", M)), table_cells(paste0("stmsh", M))), collapse = " & "), " \\\\")),
  "\\midrule", paste0("SV & ", paste(table_cells("sv"), collapse = " & "), " & & \\\\") ,
  "\\bottomrule", "\\end{tabular*}", "", "\\medskip", "\\begin{minipage}{\\textwidth}", "\\footnotesize",
  "\\textit{Note:} The table reports forecasting performance measures for a 1-month-ahead forecast from a recursive",
  "expanding-window exercise. The difference of LPSs is the average log predictive score across all periods",
  "minus that computed for the SVAR with EXH. Its positive values mean that a model performs better than",
  "the benchmark EXH model. Ratios of MASEs are the point forecast measures relative to those",
  "of the benchmark model. Values less than 1 indicate improved performance relative to the benchmark.",
  "\\end{minipage}", "\\end{table}")
writeLines(tab_main, file.path(output_dir, "tables/tab_forecast_metrics.tex"))

plot_lps <- function() {
  par(mar = c(4.5, 4.5, 1, 6), las = 1)
  colours <- c(HMSH = "#ff69b4", MSH = "#ffd700")
  benchmark_colours <- c(SV = "#78A9D1", EXH = "#000000")
  s <- summary_all[summary_all$family %in% names(colours), ]
  # MASE ratios 0.9--1.0 share the plotting height of LPS values 38--42.
  mase_y <- function(ratio) 38 + (ratio - 0.9) * 40
  plot(s$M, s$LPS, type = "n", ylim = c(38, 42), xlim = c(2, 15), xaxt = "n", bty = "l",
       xlab = "M", ylab = "LPS")
  axis(1, at = regimes)
  right_ticks <- seq(0.9, 1, by = 0.025)
  axis(4, at = mase_y(right_ticks), labels = formatC(right_ticks, format = "f", digits = 3))
  mtext("MASE relative to EXH", side = 4, line = 4.3, las = 0)
  sv <- row_of("sv"); exh <- row_of("exh")
  abline(h = c(sv$LPS, exh$LPS), lty = 1, lwd = 2, col = benchmark_colours)
  abline(h = mase_y(c(sv$MASE_ratio_EXH, exh$MASE_ratio_EXH)), lty = 2,
         lwd = 2, col = benchmark_colours)
  for (f in names(colours)) {
    d <- s[s$family == f, ]; d <- d[order(d$M), ]
    lines(d$M, d$LPS, type = "o", pch = 16, col = colours[[f]], lwd = 2, cex = 0.7)
    lines(d$M, mase_y(d$MASE_ratio_EXH), type = "o", pch = 1,
          lty = 2, col = colours[[f]], lwd = 2, cex = 0.7)
  }
  legend(2.2, 39.65, c("LPS (density)", "MASE (point)"), col = "grey30", lty = c(1, 2),
         lwd = 2, pch = c(16, 1), bty = "n", cex = 0.85, seg.len = 3)
  model_colours <- c(colours, benchmark_colours)
  legend(12, 39.65, names(model_colours), col = model_colours,
         lty = 1, lwd = 2, bty = "n", cex = 0.85)
}
pdf(file.path(output_dir, "figures/lps_by_regimes.pdf"), width = 7, height = 4.5)
plot_lps(); invisible(dev.off())
png(file.path(output_dir, "figures/lps_by_regimes.png"), width = 1680, height = 1080, res = 240)
plot_lps(); invisible(dev.off())
cat("Selected by LPS:", label(selected), "\n")
cat("Best MASE:", label(best_mase), "\n")
cat("Wrote LPS/MASE results, tables and figures to", output_dir, "\n")
