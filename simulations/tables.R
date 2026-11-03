# process spartan output
######################################################
load("results/bsvar_ms_wrap.rda")
# load("results/bsvar_ms_wrap_n2.rda")


# l- and p-value for homoskedasticity verification
######################################################
row_names = c("homo", "homo1", "hetero1", "hetero")
col_names = c("SV", "GARCH", "MS", "HMS")
l_value   = array(NA, c(4,4,length(phrase),length(samples)))
q_value   = array(NA, c(3,4,length(phrase),length(samples)))

size_adj = function(x) {
  cp = mean(x[,1][order(x[,1])][95:96])
  pv_tmp = apply(x > cp, 2, mean)
  return(pv_tmp)
}

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {

    sddr_sp_tmp = matrix(sddrs[p,s,], nrow = 100)
    lv_tmp = apply(sddr_sp_tmp > 0, 2, mean)
    l_value[,,p,s] = array(
      cbind(
        lv_tmp[1:4],
        c(lv_tmp[1],lv_tmp[5:7]),
        c(lv_tmp[1],lv_tmp[8:10]),
        c(lv_tmp[1],lv_tmp[c(8:9,11)])
      ),
      c(4, 4),
      dimnames = list(row_names, col_names)
    )

    q_value[,,p,s] = array(
      cbind(
        size_adj(sddr_sp_tmp[,2:4]),
        size_adj(sddr_sp_tmp[,5:7]),
        size_adj(sddr_sp_tmp[,8:11])[1:3],
        size_adj(sddr_sp_tmp[,8:11])[c(1,2,4)]
      ),
      c(3, 4),
      dimnames = list(row_names[2:4], col_names)
    )

  }
}
xtable::xtable(
  cbind(
    rbind(
      l_value[,,2,1],
      l_value[,,2,2]
    ),
    rbind(
      l_value[,,3,1],
      l_value[,,3,2]
    ),
    rbind(
      l_value[,,1,1],
      l_value[,,1,2]
    ),
    rbind(
      l_value[,,4,1],
      l_value[,,4,2]
    )
  )
)
phrase[c(2,3,1,4)]

xtable::xtable(
  cbind(
    rbind(
      q_value[,,2,1],
      q_value[,,2,2]
    ),
    rbind(
      q_value[,,3,1],
      q_value[,,3,2]
    ),
    rbind(
      q_value[,,1,1],
      q_value[,,1,2]
    ),
    rbind(
      q_value[,,4,1],
      q_value[,,4,2]
    )
  )
)
phrase[c(2,3,1,4)]


# estimation mean for B (bias)
######################################################
bia_tmp     = apply(BB, 2:4, mean)
bia_tmp     = BB[1,,,]
only_hetero = c(4,7,10,11)
dgps_hetero = c("dgp_sv","dgp_garch","dgp_ms","dgp_hms")
bias_B      = array(NA,
                    c(4,length(phrase),length(samples)),
                    dimnames = list(dgps_hetero,phrase,samples)
)
trim_mean = function(x) {
  codahdi = coda::HPDinterval(coda::as.mcmc(x), prob = 0.9)
  return(mean(x[x > codahdi[1] & x < codahdi[2]]))
}


for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    Bmse_sp_tmp = matrix(bia_tmp[p,s,], nrow = 100)
    bias_B[,p,s] = apply(Bmse_sp_tmp[,only_hetero],2,trim_mean)
  }
}
bias_B


# estimation rmse for B (efficiency)
######################################################
rmse = function(x) {mean(x^2)}
mse_tmp     = apply(irfs, 2:4, rmse)
only_hetero = c(4,7,10,11)
dgps_hetero = c("dgp_sv","dgp_garch","dgp_ms","dgp_hms")
rmse_B      = array(NA,
                    c(4,length(phrase),length(samples)),
                    dimnames = list(dgps_hetero,phrase,samples)
)
trim_mean_pos = function(x) {
  return(mean(x[x < sort(x)[80]], na.rm = TRUE))
}

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    Bmse_sp_tmp = matrix(mse_tmp[p,s,], nrow = 100)
    rmse_B[,p,s] = sqrt(apply(Bmse_sp_tmp[,only_hetero],2,trim_mean_pos))
  }
}
rmse_B = rmse_B[,c(2,3,1,4),]
round((rmse_B[,,1] / rmse_B[,1,1]), 3) #T=260 samples
round((rmse_B[,,2] / rmse_B[,1,2]), 3) #T=780 samples

xtable::xtable(
  rbind(
    round((rmse_B[,,1] / rmse_B[,1,1]), 3), #T=260 samples
    round((rmse_B[,,2] / rmse_B[,1,2]), 3) #T=780 samples
  )
)




# estimation rmse for IRFs
######################################################
rmse = function(x) {mean(x^2)}
mse_tmp     = apply(irfs, 2:4, rmse)
only_hetero = c(4,7,10,11)
dgps_hetero = c("dgp_sv","dgp_garch","dgp_ms","dgp_hms")
rmse_B      = array(NA,
                    c(4,length(phrase),length(samples)),
                    dimnames = list(dgps_hetero,phrase,samples)
)
trim_mean_pos = function(x) {
  return(mean(x[x < sort(x)[80]], na.rm = TRUE))
}

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    Bmse_sp_tmp = matrix(mse_tmp[p,s,], nrow = 100)
    rmse_B[,p,s] = sqrt(apply(Bmse_sp_tmp[,only_hetero],2,trim_mean_pos))
  }
}
rmse_B = rmse_B[,c(2,3,1,4),]
round((rmse_B[,,1] / rmse_B[,1,1]), 3) #T=260 samples
round((rmse_B[,,2] / rmse_B[,1,2]), 3) #T=780 samples

xtable::xtable(
  rbind(
    round((rmse_B[,,1] / rmse_B[,1,1]), 3), #T=260 samples
    round((rmse_B[,,2] / rmse_B[,1,2]), 3) #T=780 samples
  )
)




# estimation rmse for sigma (efficiency)
######################################################
rmse = function(x) {mean(x^2)}
mse_s1_tmp  = apply(sigma1, 3:4, rmse)
mse_s2_tmp  = apply(sigma2, 3:4, rmse)
mse_s_tmp   = aperm(
  array(c(mse_s1_tmp, mse_s2_tmp), c(dim(mse_s1_tmp), 2)),
  c(1,3,2)
)
rmse_s      = array(NA,
                    c(4,length(phrase),length(samples)),
                    dimnames = list(dgps_hetero,phrase,samples)
)

for (s in 1:length(samples)) {
  for (p in 1:length(phrase)) {
    smse_sp_tmp = matrix(mse_s_tmp[p,s,], nrow = 100)
    rmse_s[,p,s] = sqrt(apply(smse_sp_tmp[,only_hetero],2,trim_mean_pos))
  }
}
rmse_s = rmse_s[,c(2,3,1,4),]
round((rmse_s[,,1] / rmse_s[,1,1]), 3) #T=260 samples
round((rmse_s[,,2] / rmse_s[,1,2]), 3) #T=780 samples

xtable::xtable(
  rbind(
    round((rmse_s[,,1] / rmse_s[,1,1]), 3), #T=260 samples
    round((rmse_s[,,2] / rmse_s[,1,2]), 3) #T=780 samples
  )
)

