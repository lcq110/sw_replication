# Paper tab:sddr: baseline, s=1, and duration sensitivity (author cases 1/3).
models = c("hmsh4", "msh10", "exh", "sv")
sddrs = matrix(NA_real_, 4, 3,
  dimnames = list(c("HMSH(4)", "MSH(10)", "EXH", "SV"),
                  c("baseline", "s=1", "duration")))
for (i in seq_along(models)) {
  load(paste0("estimations/results/e_", models[i], "_post.rda"))
  sddrs[i,1] = sddr$logSDDR[1]
}
for (i in 1:2) for (j in 1:2) {
  load(paste0("estimations/results/e_", models[i], "_", c(1,3)[j], "_post.rda"))
  sddrs[i,j+1] = sddr$logSDDR[1]
}
print(round(sddrs, 2))
print(xtable::xtable(sddrs, digits = 2), file = "estimations/out_sddr.tex")
