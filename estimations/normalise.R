# Supply a reference row explicitly, after checking its economic identity.
# The author's indices for the original posteriors were 8/3/8/10 (HMSH/MSH/EXH/SV).
library(bsvars)
args = commandArgs(TRUE)
model = args[1]
reference_row = as.integer(args[2])
load(paste0("estimations/results/e_", model, ".rda"))
source("estimations/align_reference.R")
N = nrow(post$posterior$B)
Bhat = post$last_draw$starting_values$B[c(reference_row, (1:N)[-reference_row]), ]
Bhat[1, ] = sign(Bhat[1, 3]) * Bhat[1, ]
family = if (startsWith(model, "hmsh")) "hmsh" else if (startsWith(model, "msh")) "msh" else model
alignment = align_reference(post$posterior, Bhat)
normalised = alignment$draws
constructor = get(paste0("specify_posterior_bsvar_", family), asNamespace("bsvars"))
post = constructor$new(post$last_draw, normalised)
post$set_normalised()
save(post, Bhat, reference_row, alignment, file = paste0("estimations/results/e_", model, "_norm.rda"))
sddr = verify_volatility(post)
csds = compute_conditional_sd(post)
summ = summary(post)
save(sddr, csds, summ, file = paste0("estimations/results/e_", model, "_post.rda"))
irfs = compute_impulse_responses(post, horizon = 60)
save(irfs, file = paste0("estimations/results/e_", model, "_irfs.rda"))
