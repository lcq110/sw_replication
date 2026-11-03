# From the replication root: Rscript estimations/estimate.R hmsh4 [prior_case]
# Prior cases 1/2/3 follow e_{hmsh,msh}_rob.R; paper tab:sddr uses 1 and 3.
source("estimations/model.R")
args = commandArgs(TRUE)
model = args[1]
prior_case = if (length(args) > 1) as.integer(args[2]) else 0L
set.seed(123)
spec = model_specification(y_full, model)
if (prior_case == 0L) {
  burn = estimate(spec, S = 1e6, thin = 1e5)
  if (model != "sv") {
    # One pass of the author's Step 2. Historical repeat count is not recorded.
    set.seed(1234)
    burn = estimate(burn, S = 2e5, thin = 1e4)
  }
} else {
  load(paste0("estimations/results/e_", model, ".rda"))
  spec$starting_values = post$last_draw$starting_values
  M = if (startsWith(model, "hmsh")) as.integer(sub("hmsh", "", model)) else as.integer(sub("msh", "", model))
  if (prior_case == 1L) spec$prior$sigma_s = 1
  if (prior_case == 2L) spec$prior$sigma_s = 8
  if (prior_case == 3L) spec$prior$PR_TR = matrix(1, M, M) + 9 * diag(M)
  spec$prior$sigma_nu = spec$prior$sigma_s + 2
  burn = estimate(spec, S = 2e5, thin = 1e4)
}
set.seed(1234)
post = estimate(burn, S = 5e4, thin = 1)
dir.create("estimations/results", recursive = TRUE, showWarnings = FALSE)
name = if (prior_case == 0L) model else paste0(model, "_", prior_case)
save(post, file = paste0("estimations/results/e_", name, ".rda"))
