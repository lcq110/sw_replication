library(bsvars)
load("data/bpss21.rda")
load("data/bpss21_regimes.rda")
y_full = ts(cbind(diff(bpss21[,1:2], lag = 12), bpss21[-(1:12),6],
                  diff(bpss21[,c(3:5,7)], lag = 12), bpss21[-(1:12),8:10]),
            start = c(1974,1), frequency = 12,
            names = colnames(bpss21)[c(1,2,6,3,4,5,7,8,9,10)])
# Specifications follow empirical_rr/e_*.R and the ordinary forecast burn-in scripts.
model_specification = function(y, model) {
  N = ncol(y)
  B = matrix(TRUE, N, N)
  stationary = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  if (startsWith(model, "hmsh")) {
    M = as.integer(sub("hmsh", "", model))
    spec = specify_bsvar_hmsh$new(y, p = 10, B = B, M = M,
                                stationary = stationary, finiteM = TRUE)
  } else if (startsWith(model, "msh")) {
    M = as.integer(sub("msh", "", model))
    spec = specify_bsvar_msh$new(y, p = 10, B = B, M = M,
                               stationary = stationary, finiteM = TRUE)
  } else if (model == "exh") {
    regimes = bpss21_regimes[-(1:12)][seq_len(nrow(y))]
    M = length(unique(regimes))
    spec = specify_bsvar_exh$new(y, p = 10, B = B, stationary = stationary,
                                variance_regimes = regimes)
  } else if (model == "sv") {
    return(specify_bsvar_sv$new(y, p = 10, B = B, stationary = stationary))
  }
  spec$prior$sigma_s = M
  spec$prior$sigma_nu = M + 2
  spec
}
