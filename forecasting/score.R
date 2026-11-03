# Joint LPS and MASE: adapted from the author's score_final_results.R.
source("estimations/model.R")
models = c(paste0("hmsh", 2:15), paste0("msh", 2:15), "sv", "exh")
rows = list()
k = 1L
for (model in models) for (iteration in 1:117) {
  origin = 504 + iteration
  model_id = match(model, models)
  load(sprintf("forecasting/results/bsvar%02d_%d.rda", model_id, iteration))
  actual = as.numeric(y_evaluation[1, ])
  scale = colMeans(abs(diff(y_full[1:origin, ])))
  point = apply(fore$forecasts[,1,], 1, median)
  mase = mean(abs(actual - point) / scale)
  draws = dim(fore$forecast_mean)[3]
  log_density = vapply(seq_len(draws), function(s)
    mvtnorm::dmvnorm(actual, fore$forecast_mean[,1,s], fore$forecast_covariance[,,1,s], log = TRUE), 0.0)
  maximum = max(log_density)
  lps = maximum + log(mean(exp(log_density - maximum)))
  rows[[k]] = data.frame(model = sub("^hmsh", "sthmsh", sub("^msh", "stmsh", model)),
    origin = format(zoo::as.yearmon(time(y_full)[origin]), "%Y-%m-01"),
    target = format(zoo::as.yearmon(time(y_full)[origin + 1]), "%Y-%m-01"), lps = lps, mase = mase)
  k = k + 1L
}
dir.create("forecasting/outputs", showWarnings = FALSE)
write.csv(do.call(rbind, rows), "forecasting/outputs/origin_scores-rerun.csv", row.names = FALSE)
