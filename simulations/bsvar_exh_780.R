
Rcpp::sourceCpp("src_cpp/bsvar_msh.cpp")
Rcpp::sourceCpp("src_cpp/normalise_jaro.cpp")

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
iteration <- as.integer(args[1])
rm(args)

# iteration = 1

library(bsvars)
set.seed(123 + iteration)

S_burn          = 1e4
S               = 1e4
T               = 780

load("dgps/dgp_B.rda")
B_dgp           = B
rm("B")
id              = paste0(paste0(rep(0, 4 - nchar(as.character(iteration))), collapse = ""), as.character(iteration))

yy              = as.matrix(read.csv(paste0("dgps/dgp", T, "_", id, ".csv"), header = FALSE))
N               = ncol(yy)
T               = nrow(yy)

variance_regimes = rep(1, 130)
for (i in 2:6) {
  variance_regimes = c(variance_regimes, rep(i, 130))
}

B               = matrix(TRUE, N, N)
spec            = specify_bsvar_exh$new(
  yy, 
  B = B,
  variance_regimes = variance_regimes
)
prior           = spec$prior$get_prior()

starting_values = spec$starting_values$get_starting_values()
starting_values$B = B_dgp
starting_values$xi = cbind(starting_values$xi[,1], starting_values$xi)
starting_values$hyper[,1] = 500
VB              = spec$identification$VB

burn  = bsvar_exh_cpp ( S_burn, t(yy), prior, VB, starting_values, thin = 1)
post  = bsvar_exh_cpp ( S, t(yy), prior, VB, burn$last_draw, thin = 1)

P     = pick_permutation_sse (post$last_draw$B, B_dgp,
                              allPermutations(N) - 1, 
                              allSigns(N), 
                              communication_matrix(N, N))
B_hat = P %*% post$last_draw$B
post$posterior = normalise_jaro_bsvar_exh (post$posterior, B_hat = B_hat)

sddr            = -as.numeric(verify_volatility_msh_cpp( post$posterior, prior, t(yy) ))

load("dgps/dgp_B.rda")
load(paste0("dgps/dgp_sds", T, "_", id, ".rda"))

Bm1 = apply(post$posterior$B, 1:2, mean)
Bm1 = diag(1/diag(Bm1)) %*% Bm1
B1 = diag(1/diag(B)) %*% B
Bm = (Bm1 - B1)[2:3]

IR1 = solve(B)
IR1 = diag(1 / diag(IR1)) %*% IR1
ir1 = sapply(
  1:dim(post$posterior$B)[3],
  \(s){
    ir = solve(post$posterior$B[,,s])
    ir = ir %*% diag(1 / diag(ir))
    ir[2:3]
  }
)
ir1 = apply(ir1, 1, mean) - IR1[2:3]

Bm1 = diag(1/diag(Bm1)) %*% Bm1
B1 = diag(1/diag(B)) %*% B
Bm = (Bm1 - B1)[2:3]

sm1 = t(apply(post$posterior$sigma, 1:2, mean))
sm1 = (sm1 - matrix(colMeans(sm1), nrow(sm1), ncol(sm1), byrow = TRUE)) / matrix(apply(sm1, 2, sd), nrow(sm1), ncol(sm1), byrow = TRUE)
sds1 = sds
for (n in 1:2) {
  if (all(sds[,n]==1)) {
    sds1[,n] = 0
  } else {
    sds1[,n] = (sds[,n] - mean(sds[,n])) / sd(sds[,n])
  }
}
sm = sm1 - sds1 

save(
  sddr, Bm, sm, ir1,
  file = paste0("results/exh_", T, "_", id, ".rda")
)

# norm2
B_hat = B_dgp
post$posterior = normalise_jaro_bsvar_exh (post$posterior, B_hat = B_hat)

sddr            = -as.numeric(verify_volatility_msh_cpp( post$posterior, prior, t(yy) ))

load("dgps/dgp_B.rda")
load(paste0("dgps/dgp_sds", T, "_", id, ".rda"))

Bm1 = apply(post$posterior$B, 1:2, mean)
Bm1 = diag(1/diag(Bm1)) %*% Bm1
B1 = diag(1/diag(B)) %*% B
Bm = (Bm1 - B1)[2:3]

IR1 = solve(B)
IR1 = diag(1 / diag(IR1)) %*% IR1
ir1 = sapply(
  1:dim(post$posterior$B)[3],
  \(s){
    ir = solve(post$posterior$B[,,s])
    ir = ir %*% diag(1 / diag(ir))
    ir[2:3]
  }
)
ir1 = apply(ir1, 1, mean) - IR1[2:3]

sm1 = t(apply(post$posterior$sigma, 1:2, mean))
sm1 = (sm1 - matrix(colMeans(sm1), nrow(sm1), ncol(sm1), byrow = TRUE)) / matrix(apply(sm1, 2, sd), nrow(sm1), ncol(sm1), byrow = TRUE)
sds1 = sds
for (n in 1:2) {
  if (all(sds[,n]==1)) {
    sds1[,n] = 0
  } else {
    sds1[,n] = (sds[,n] - mean(sds[,n])) / sd(sds[,n])
  }
}
sm = sm1 - sds1 

save(
  sddr, Bm, sm, ir1,
  file = paste0("results/exh_", T, "_", id, "_n2.rda")
)
