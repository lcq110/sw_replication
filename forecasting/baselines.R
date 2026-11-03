# From the repository root: Rscript forecasting/baselines.R sv (or exh).
# Expanding-window workflow follows the ordinary R&R scripts.
source("estimations/model.R")
model = commandArgs(TRUE)[1]
model_id = if (model == "sv") 29L else 30L
set.seed(20241110)
star = model_specification(y_full, model)
if (model == "exh") {
  # Prediction project's EXH default prior; empirical fits use a different prior.
  star$prior$sigma_s = 1
  star$prior$sigma_nu = 3
}
burn = estimate(star, S = 1e6, thin = 1e5)
star = burn$last_draw
dir.create("forecasting/burnin/starting_values", recursive = TRUE, showWarnings = FALSE)
state_file = paste0("forecasting/burnin/starting_values/e_", model, "_star.rda")
save(star, file = state_file)
dir.create("forecasting/results", showWarnings = FALSE)
for (iteration in 1:117) {
  origin = 504 + iteration
  y_evaluation = y_full[origin + 1, , drop = FALSE]
  load(state_file)
  spec = star
  T = origin - spec$p
  spec$data_matrices$Y = spec$data_matrices$Y[,1:T,drop=FALSE]
  spec$data_matrices$X = spec$data_matrices$X[,1:T,drop=FALSE]
  spec$starting_values$lambda = spec$starting_values$lambda[,1:T,drop=FALSE]
  if (model == "sv") {
    spec$starting_values$h = spec$starting_values$h[,1:T,drop=FALSE]
    spec$starting_values$S = spec$starting_values$S[,1:T,drop=FALSE]
  } else {
    spec$starting_values$xi = spec$starting_values$xi[,1:T,drop=FALSE]
  }
  set.seed(20260917 + (model_id - 1) * 10000 + origin)
  burn = estimate(spec, S = 5e3, thin = 5e3)
  post = estimate(burn, S = 1e4)
  fore = forecast(post, horizon = 1)
  save(fore, y_evaluation,
       file = sprintf("forecasting/results/bsvar%02d_%d.rda", model_id, iteration))
}
