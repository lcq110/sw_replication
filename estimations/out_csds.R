


output = "csds"

which_hmsh = 1
which_msh = 1

models_hmsh = paste0("hmsh",c(4))
models_msh = paste0("msh",c(10))

mps_hmsh = 1
mps_msh = 1
which_exh = 1
which_sv = 1

models  = c(models_hmsh[which_hmsh],models_msh[which_msh],"exh","sv")
model_names = c("HMSH(4)","MSH(10)","EXH","SV")
mps_id  = c(mps_hmsh[which_hmsh],mps_msh[which_msh],which_exh,which_sv)
M       = length(models)

pdf(
  file = paste0("estimations/out_",output,".pdf"),
  height = 8,
  width = 8
)
par(
  mfrow = c(M,1),
  # mfrow = c(1,1),
  mar = c(3,3,1,1),
  cex.lab = 2,
  cex.axis = 2
)
for (model in 1:M) {
  
  cat("Model: ", models[model], "\n")
  
  load(paste0("estimations/results/e_",models[model],"_post.rda"))
  
  if (model == M) {
    ylim = c(0,20)
  } else {
    ylim = c(0,4)
  }
  
  bsvars::plot_ribbon(
    csds[mps_id[model],,],
    ylim = ylim,
    start_at = 1,
    axes = FALSE
  )
  if (model == M){
    axis(
      1,
      at = which(as.numeric(colnames(csds)) %% 5 == 0),
      labels = seq(from = 1975, to = 2025, by = 5)
    )
    axis(
      2,
      at = c(0,1,4,20),
      labels = c(NA,1,4,NA)
    )
    axis(
      2,
      at = 4,
      labels = 4
    )
    text(
      550,
      16,
      model_names[model],
      cex = 2.5
    )
  } else {
    axis(
      2,
      at = c(0,1,4),
      labels = c(NA,1,4)
    )
    text(
      550,
      3,
      model_names[model],
      cex = 2.5
    )
  }
}
dev.off()



