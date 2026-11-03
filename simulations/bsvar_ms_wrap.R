
# read spartan outputs
######################################################
# phrase      = c("exh", "hmsh", "hmsh2", "msh", "msh2", "sv")
phrase      = c("exh", "hmsh2", "msh2", "sv")
samples     = c(260, 780)
iterations  = 1100

sddrs       = array(NA, c(length(phrase), length(samples), iterations))
BB          = array(NA, c(2, length(phrase), length(samples), iterations))
irfs        = array(NA, c(2, length(phrase), length(samples), iterations))
sigma1      = array(NA, c(2, samples[1], length(phrase), iterations))
sigma2      = array(NA, c(2, samples[2], length(phrase), iterations))
i           = 1

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    # if (s == 2 && p ==3) next
    for (iteration in 1:1100) {
      
      id              = paste0(paste0(rep(0, 4 - nchar(as.character(iteration))), collapse = ""), as.character(iteration))
      load(paste0("results/", phrase[p],"_", samples[s], "_", id, ".rda"))
      
      sddrs[p,s,iteration]    = sddr
      BB[,p,s,iteration]      = Bm
      irfs[,p,s,iteration]    = ir1
      if (s == 1) {
        sigma1[,,p,iteration] = sm
      } else {
        sigma2[,,p,iteration] = sm
      }
      
    }
  }
}

save(
  phrase, samples, iterations, sddrs, BB, sigma1, sigma2, irfs,
  file = "results/bsvar_ms_wrap.rda"
)


sddrs       = array(NA, c(length(phrase), length(samples), iterations))
BB          = array(NA, c(2, length(phrase), length(samples), iterations))
irfs        = array(NA, c(2, length(phrase), length(samples), iterations))
sigma1      = array(NA, c(2, samples[1], length(phrase), iterations))
sigma2      = array(NA, c(2, samples[2], length(phrase), iterations))
i           = 1

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    # if (s == 2 && p ==3) next
    for (iteration in 1:1100) {
      
      id              = paste0(paste0(rep(0, 4 - nchar(as.character(iteration))), collapse = ""), as.character(iteration))
      load(paste0("results/", phrase[p],"_", samples[s], "_", id, "_n2.rda"))
      
      sddrs[p,s,iteration]    = sddr
      BB[,p,s,iteration]      = Bm
      irfs[,p,s,iteration]    = ir1
      if (s == 1) {
        sigma1[,,p,iteration] = sm
      } else {
        sigma2[,,p,iteration] = sm
      }
      
    }
  }
}

save(
  phrase, samples, iterations, sddrs, BB, sigma1, sigma2, irfs,
  file = "results/bsvar_ms_wrap_n2.rda"
)






