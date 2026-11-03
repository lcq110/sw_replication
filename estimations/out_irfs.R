


output = "irfs"

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

load("data/bpss21.rda")
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
var_order = c(1:5,9,10,8,6,7)
cny       = c("ip","p","R","hhc","bc","GZ","ES","TS","m","pcm")

ylim = rbind(
  c(-0.5,.5),
  c(-0.14,.14),
  c(-0.2,0.5),
  c(-0.58,.5),
  c(-0.8,.8),
  c(-0.35,0.35),
  c(-1.2,0.8),
  c(-0.3,0.1),
  c(-0.035,.035),
  c(-0.03,0.07)
)
ylab = list(
  c(-0.5,0,.5),
  c(-0.1,0,.1),
  c(0,0.25,.5),
  c(-0.5,0,0.5),
  c(-0.7,0,.7),
  c(-.3,0,.3),
  c(-1.2,0,.5),
  c(-0.3,0,0.1),
  c(-0.03,0,.03),
  c(0,0.07)
)

pdf(
  file = paste0("estimations/out_",output,"1.pdf"),
  height = 9,
  width = 9
)
par(
  mfcol = c(5,M),
  # mfrow = c(1,1),
  mar = c(2.8,4.2,1,0.1),
  cex.lab = 1.5,
  cex.axis = 1.5
)
for (model in 1:M) {
  
  cat("Model: ", models[model], "\n")
  
  load(paste0("estimations/results/e_",models[model],"_irfs.rda"))
  
  mps_irfs = (0.25 / mean(irfs[3,mps_id[model],1,])) * irfs[,mps_id[model],,]
  
  for (i in 1:5) {
    
    if (i == 1) {
      main = model_names[model]
    } else {
      main = ""
    }
    
    bsvars::plot_ribbon(
      mps_irfs[var_order[i],,],
      start_at = 0,
      ylim = ylim[var_order[i],],
      axes = FALSE,
      main = main
    )
    abline(h = 0)
    if (i == 5){
      axis(
        1,
        at = c(0,seq(from = 13, to = 61, by = 12)),
        labels = c(NA,1:5)
      )
    }
    if (model == 1) {
      axis(
        2,
        at = ylab[[var_order[i]]],
        labels = ylab[[var_order[i]]]
      )
      title(ylab = cny[i])
    }
  }
}
dev.off()






pdf(
  file = paste0("estimations/out_",output,"2.pdf"),
  height = 9,
  width = 9
)
par(
  mfcol = c(5,M),
  # mfrow = c(1,1),
  mar = c(2.8,4.2,1,0.1),
  cex.lab = 1.5,
  cex.axis = 1.5
)
for (model in 1:M) {
  
  cat("Model: ", models[model], "\n")
  
  load(paste0("estimations/results/e_",models[model],"_irfs.rda"))
  
  mps_irfs = (0.25 / mean(irfs[3,mps_id[model],1,])) * irfs[,mps_id[model],,]
  
  for (i in 6:10) {
    
    if (i == 6) {
      main = model_names[model]
    } else {
      main = ""
    }
    
    bsvars::plot_ribbon(
      mps_irfs[var_order[i],,],
      start_at = 0,
      ylim = ylim[var_order[i],],
      axes = FALSE,
      main = main
    )
    abline(h = 0)
    if (i == 10){
      axis(
        1,
        at = c(0,seq(from = 13, to = 61, by = 12)),
        labels = c(NA,1:5)
      )
    }
    if (model == 1) {
      axis(
        2,
        at = ylab[[var_order[i]]],
        labels = ylab[[var_order[i]]]
      )
      title(ylab = cny[i])
    }
  }
}
dev.off()




