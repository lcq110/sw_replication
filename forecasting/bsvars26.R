
# bsvar normal MSH, 13 regimes, full B
################################################

model   = "26"
S       = 10e3
S_burn  = 5e3

forecast_origin_first = 505 - 1

library(bsvars)

# data work
##################################################
load("../data/bpss21.rda")
y_full = ts(
  cbind(diff(bpss21[,1:2], lag = 12), bpss21[-(1:12),6],
        diff(bpss21[,c(3:5,7)], lag = 12), bpss21[-(1:12),8:10]),
  start = c(1974,1), frequency = 12,
  names = colnames(bpss21)[c(1,2,6,3,4,5,7,8,9,10)]
)

dir.create("results", showWarnings = FALSE)

for (iteration in 1:117) {
  forecast_origin = forecast_origin_first + iteration
  y = y_full[1:forecast_origin,]
  y_evaluation = y_full[forecast_origin + 1,,drop = FALSE]

  # bvar estimation and forecasting
  ##################################################
  load("burnin/starting_values/e_stmsh13_star.rda")
  spec = star

  T = nrow(y) - spec$p
  spec$data_matrices$Y = spec$data_matrices$Y[,1:T,drop = FALSE]
  spec$data_matrices$X = spec$data_matrices$X[,1:T,drop = FALSE]
  spec$starting_values$lambda = spec$starting_values$lambda[,1:T,drop = FALSE]
  spec$starting_values$xi = spec$starting_values$xi[,1:T,drop = FALSE]

  set.seed(20260917 + (as.integer(model) - 1) * 10000 + forecast_origin)
  burn = estimate(spec, S_burn, thin = S_burn)
  post = estimate(burn, S)

  fore = forecast(post, horizon = 1)

  save(
    fore, y_evaluation,
    file = paste0("results/bsvar",model,"_",iteration,".rda")
  )
}
