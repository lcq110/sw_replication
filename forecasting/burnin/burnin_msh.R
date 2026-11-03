library(bsvars)

Mselect   = 2:15
burn_S    = 1e6
burn_thin = 1e5

load("../data/bpss21.rda")

y = ts(
  cbind(
    diff(bpss21[,1:2],lag=12),
    bpss21[-(1:12),6],
    diff(bpss21[,c(3:5,7)],lag=12),
    bpss21[-(1:12),8:10]
  ),
  start = c(1974, 1),
  frequency = 12,
  names = colnames(bpss21)[c(1,2,6,3,4,5,7,8,9,10)]
)

N     = dim(y)[2]
B     = matrix(TRUE, N, N)

dir.create("burnin/starting_values", recursive = TRUE, showWarnings = FALSE)

for (M in Mselect) {
  model = paste0("stmsh", M)

  # Step 1: specify and estimate
  set.seed(123)
  spec  = specify_bsvar_msh$new(
    data = y,
    stationary = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    p = 10,
    M = M,
    B = B,
    finiteM = TRUE
    # distribution = "t",
  )
  spec$prior$sigma_s  = M * 1
  spec$prior$sigma_nu = spec$prior$sigma_s + 2

  burn  = estimate(spec, S = burn_S, thin = burn_thin)
  star  = burn$last_draw

  save(
    star,
    file = paste0("burnin/starting_values/e_", model, "_star.rda")
  )
}
