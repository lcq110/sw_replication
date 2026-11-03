# First-row minimum-SSE alignment intended by empirical_rr/src/normalise_jaro.cpp.
# The old source compares C++ row(1) while its candidates move the shock to row(0),
# and leaves HMSH state arrays unpermuted. This local adapter uses the selected
# first reference row and applies the same signed permutation to every shock field.
align_reference = function(draws, Bhat) {
  N = dim(draws$B)[1]
  S = dim(draws$B)[3]
  permutations = matrix(NA_integer_, N, S)
  signs = numeric(S)
  for (s in seq_len(S)) {
    B = draws$B[,,s]
    positive = rowSums((B - matrix(Bhat[1,], N, N, byrow = TRUE))^2)
    negative = rowSums((-B - matrix(Bhat[1,], N, N, byrow = TRUE))^2)
    candidate = which.min(c(positive, negative))
    chosen = (candidate - 1L) %% N + 1L
    sign = if (candidate <= N) 1 else -1
    permutation = c(chosen, (1:N)[-chosen])
    permutations[,s] = permutation
    signs[s] = sign
    draws$B[,,s] = B[permutation,]
    draws$B[1,,s] = sign * draws$B[1,,s]
    for (field in intersect(c('sigma','lambda','sigma2','h','S'), names(draws)))
      draws[[field]][,,s] = draws[[field]][permutation,,s]
    for (field in intersect(c('df','rho','omega','sigma2v','sigma2_omega','s_'), names(draws)))
      draws[[field]][,s] = draws[[field]][permutation,s]
    # Column 1 belongs to B/shocks; column 2 belongs to unchanged reduced-form A.
    draws$hyper[1:N,1,s] = draws$hyper[permutation,1,s]
    draws$hyper[N + 1:N,1,s] = draws$hyper[N + permutation,1,s]
    if ('xi_cpp' %in% names(draws)) {
      draws$xi_cpp[[s]] = draws$xi_cpp[[s]][,,permutation,drop=FALSE]
      draws$PR_TR_cpp[[s]] = draws$PR_TR_cpp[[s]][,,permutation,drop=FALSE]
      draws$xi[,,,s] = draws$xi[,,permutation,s,drop=FALSE]
      draws$PR_TR[,,,s] = draws$PR_TR[,,permutation,s,drop=FALSE]
      draws$pi_0[,,s] = draws$pi_0[,permutation,s]
    }
  }
  list(draws = draws, permutations = permutations, signs = signs)
}
