#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Option 1: Jarocinski normalisation of B
// // [[Rcpp::export]]
// double normalise_jaro_distance (
//     mat&          B,              // NxN
//     const mat&    B_hat,          // NxN
//     const mat&    sigma2,         // NxT
//     const int     T,
//     const mat&    communication   // NNxNN communication matrix
// ) {
//   const int N     = B.n_rows; 
//   vec   vecBt     = vectorise(B.t());
//   vec   vecB_hatt = vectorise(B_hat.t());
//   
//   mat   B_hat_inv = inv(B_hat);
//   vec   sigma2_sum = sum(sigma2, 1);
//   
//   mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
//   mat   covariance_inv  = T * Finsher_info;
//     
//   double distance = as_scalar(trans(vecBt - vecB_hatt) * covariance_inv * (vecBt - vecB_hatt));
//   
//   return distance;
// } // END normalise_jaro_distance
//     
// 
//     
// Option 1: Jarocinski normalisation of B
// // [[Rcpp::export]]
// mat covariance_inv_bsvars (
//     const mat&    B_hat,          // NxN
//     const mat&    sigma2,         // NxT
//     const int     T,
//     const mat&    communication   // NNxNN communication matrix
// ) {
//   const int N     = B_hat.n_rows; 
//   vec   vecB_hatt = vectorise(B_hat.t());
//   mat   B_hat_inv = inv(B_hat);
//   vec   sigma2_sum = sum(sigma2, 1);
//   
//   mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
//   mat   covariance_inv  = T * Finsher_info;
//   
//   return covariance_inv;
// } // END covariance_inv_bsvars
// 
// 
// Option 1: Jarocinski normalisation of B
// // [[Rcpp::export]]
// mat pick_permutation_bsvars (
//     arma::mat&              B,
//     const arma::mat&        B_hat,
//     const int               T,              // sample size
//     const arma::mat&        cov_inv,
//     const arma::umat&       perms,
//     const arma::mat&        signs,
//     const arma::mat&        communication   // NNxNN communication matrix
// ) {
// 
//   const int   no_perms    = perms.n_rows;
//   const int   no_signs    = signs.n_rows;
//   const int   N           = B.n_rows;
//   
//   mat   I                 = eye(N,N);
//   vec   vecB_hatt         = vectorise(B_hat.t());
//   vec   vecBt(vecB_hatt.n_elem);
//   
//   mat   distance(no_perms, no_signs);
//   
//   for (int p=0; p<no_perms; p++) {
//     for (int s=0; s<no_signs; s++) {
//       vecBt               = vectorise(trans(I.rows(perms.row(p)) * diagmat(signs.row(s)) * B));
//       distance(p, s)      = as_scalar(trans(vecBt - vecB_hatt) * cov_inv * (vecBt - vecB_hatt));
//     } // END s loop
//   } // END p loop
//   uvec sign_min_inds      = index_min(distance, 1);
//   vec sign_mins           = min(distance, 1);
//   int perm_min_ind        = index_min(sign_mins);
//   int sign_min_ind        = sign_min_inds(perm_min_ind);
//   
//   mat P = I.rows(perms.row(perm_min_ind));
//   mat D = diagmat(trans(signs.row(sign_min_ind)));
//   
//   return P * D;
//   // return distance;
// } // END pick_permutation_bsvars



// Option 2: normalise the first row of B
// [[Rcpp::export]]
double normalise_row1_distance (
    mat&          B,              // NxN
    const mat&    B_hat          // NxN
) {
  // the distance is the square error of the first row of B less B_hat
  double distance = accu(square(B.row(1) - B_hat.row(1)));
  return distance;
} // END normalise_row1_distance




// Option 2: normalise the first row of B
// [[Rcpp::export]]
mat pick_permutation_bsvars (
    arma::mat&              B,
    const arma::mat&        B_hat,
    const arma::umat&       perms,
    const arma::mat&        signs
) {
  
  const int   no_perms    = perms.n_rows;
  const int   no_signs    = signs.n_rows;
  const int   N           = B.n_rows;
  
  mat   I                 = eye(N,N);
  mat   distance(no_perms, no_signs);
  
  for (int p=0; p<no_perms; p++) {
    for (int s=0; s<no_signs; s++) {
      mat Bt              = I.rows(perms.row(p)) * diagmat(signs.row(s)) * B;
      distance(p, s)      = normalise_row1_distance( Bt, B_hat );
    } // END s loop
  } // END p loop
  uvec sign_min_inds      = index_min(distance, 1);
  vec sign_mins           = min(distance, 1);
  int perm_min_ind        = index_min(sign_mins);
  int sign_min_ind        = sign_min_inds(perm_min_ind);
  
  mat P = I.rows(perms.row(perm_min_ind));
  mat D = diagmat(trans(signs.row(sign_min_ind)));
  
  return P * D;
  // return distance;
} // END pick_permutation_bsvars



        
    
    
    
// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_sv (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {

  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_lambda  = as<cube>(posterior["lambda"]);
  mat    posterior_df      = as<mat>(posterior["df"]);
  
  cube   posterior_h       = as<cube>(posterior["h"]);
  mat    posterior_rho     = as<mat>(posterior["rho"]);
  mat    posterior_omega   = as<mat>(posterior["omega"]);
  mat    posterior_sigma2v = as<mat>(posterior["sigma2v"]);
  cube   posterior_S       = as<cube>(posterior["S"]);
  mat    posterior_sigma2_omega = as<mat>(posterior["sigma2_omega"]);
  mat    posterior_s_      = as<mat>(posterior["s_"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;

  // Option 1: Jarocinski normalisation of B
  // Function AllPerms("allPermutations");
  // Function AllSigns("allSigns");
  // Function communication_mat("communication_matrix");
  // mat communication = as<mat>(communication_mat(N, N));
  // mat cov_inv(pow(N,2), pow(N,2));

  // Option 2: normalise the first row of B
  Function AllPerms("allPermutations1");
  Function AllSigns("allSigns1");
  
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // Option 1: Jarocinski normalisation of B
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    
    // Option 2: normalise the first row of B
    P               = pick_permutation_bsvars( posterior_B.slice(s), B_hat, allP, allS );
    
    Pabs            = abs(P);

    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_lambda.slice(s) = Pabs * posterior_lambda.slice(s);
    posterior_df.col(s)       = Pabs * posterior_df.col(s);
    
    posterior_h.slice(s)      = Pabs * posterior_h.slice(s);
    posterior_rho.col(s)      = Pabs * posterior_rho.col(s);
    posterior_omega.col(s)    = Pabs * posterior_omega.col(s);
    posterior_sigma2v.col(s)  = Pabs * posterior_sigma2v.col(s);
    posterior_S.slice(s)      = Pabs * posterior_S.slice(s);
    posterior_sigma2_omega.col(s) = Pabs * posterior_sigma2_omega.col(s);
    posterior_s_.col(s)       = Pabs * posterior_s_.col(s);
  } // END s loop

  return List::create(
    _["B"]            = posterior_B,
    _["A"]            = posterior_A,
    _["hyper"]        = posterior_hyper,
    _["sigma"]        = posterior_sigma,
    _["lambda"]       = posterior_lambda,
    _["df"]           = posterior_df,
    _["h"]            = posterior_h,
    _["rho"]          = posterior_rho,
    _["omega"]        = posterior_omega,
    _["sigma2v"]      = posterior_sigma2v,
    _["S"]            = posterior_S,
    _["sigma2_omega"] = posterior_sigma2_omega,
    _["s_"]           = posterior_s_
  );
} // END normalise_jaro_bsvar_sv




// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_exh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_lambda  = as<cube>(posterior["lambda"]);
  mat    posterior_df      = as<mat>(posterior["df"]);
  cube   posterior_sigma2  = as<cube>(posterior["sigma2"]);
  cube   posterior_xi      = as<cube>(posterior["xi"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  // Option 1: Jarocinski normalisation of B
  // Function AllPerms("allPermutations");
  // Function AllSigns("allSigns");
  // Function communication_mat("communication_matrix");
  // mat communication = as<mat>(communication_mat(N, N));
  // mat cov_inv(pow(N,2), pow(N,2));
  
  // Option 2: normalise the first row of B
  Function AllPerms("allPermutations1");
  Function AllSigns("allSigns1");
  
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // Option 1: Jarocinski normalisation of B
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    
    // Option 2: normalise the first row of B
    P               = pick_permutation_bsvars( posterior_B.slice(s), B_hat, allP, allS );
    
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_lambda.slice(s) = Pabs * posterior_lambda.slice(s);
    posterior_df.col(s)       = Pabs * posterior_df.col(s);
    posterior_sigma2.slice(s) = Pabs * posterior_sigma2.slice(s);
    
  } // END s loop
  
  return List::create(
    _["B"]            = posterior_B,
    _["A"]            = posterior_A,
    _["hyper"]        = posterior_hyper,
    _["sigma"]        = posterior_sigma,
    _["lambda"]       = posterior_lambda,
    _["df"]           = posterior_df,
    _["sigma2"]       = posterior_sigma2,
    _["xi"]           = posterior_xi
  );
} // END normalise_jaro_bsvar_exh





// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_msh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_lambda  = as<cube>(posterior["lambda"]);
  mat    posterior_df      = as<mat>(posterior["df"]);
  cube   posterior_sigma2  = as<cube>(posterior["sigma2"]);
  cube   posterior_xi      = as<cube>(posterior["xi"]);
  cube   posterior_PR_TR   = as<cube>(posterior["PR_TR"]);
  mat    posterior_pi_0    = as<mat>(posterior["pi_0"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  // Option 1: Jarocinski normalisation of B
  // Function AllPerms("allPermutations");
  // Function AllSigns("allSigns");
  // Function communication_mat("communication_matrix");
  // mat communication = as<mat>(communication_mat(N, N));
  // mat cov_inv(pow(N,2), pow(N,2));
  
  // Option 2: normalise the first row of B
  Function AllPerms("allPermutations1");
  Function AllSigns("allSigns1");
  
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // Option 1: Jarocinski normalisation of B
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    
    // Option 2: normalise the first row of B
    P               = pick_permutation_bsvars( posterior_B.slice(s), B_hat, allP, allS );
    
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_lambda.slice(s) = Pabs * posterior_lambda.slice(s);
    posterior_df.col(s)       = Pabs * posterior_df.col(s);
    posterior_sigma2.slice(s) = Pabs * posterior_sigma2.slice(s);
    
  } // END s loop
  
  return List::create(
    _["B"]            = posterior_B,
    _["A"]            = posterior_A,
    _["hyper"]        = posterior_hyper,
    _["sigma"]        = posterior_sigma,
    _["lambda"]       = posterior_lambda,
    _["df"]           = posterior_df,
    _["sigma2"]       = posterior_sigma2,
    _["xi"]       = posterior_xi,
    _["PR_TR"]    = posterior_PR_TR,
    _["pi_0"]         = posterior_pi_0
  );
} // END normalise_jaro_bsvar_msh




// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_hmsh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_lambda  = as<cube>(posterior["lambda"]);
  mat    posterior_df      = as<mat>(posterior["df"]);
  cube   posterior_sigma2  = as<cube>(posterior["sigma2"]);
  field<cube> posterior_PR_TR_cpp = as<field<cube>>(posterior["PR_TR_cpp"]);
  field<cube> posterior_xi_cpp = as<field<cube>>(posterior["xi_cpp"]);
  cube   posterior_pi_0    = as<cube>(posterior["pi_0"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  // Option 1: Jarocinski normalisation of B
  // Function AllPerms("allPermutations");
  // Function AllSigns("allSigns");
  // Function communication_mat("communication_matrix");
  // mat communication = as<mat>(communication_mat(N, N));
  // mat cov_inv(pow(N,2), pow(N,2));
  
  // Option 2: normalise the first row of B
  Function AllPerms("allPermutations1");
  Function AllSigns("allSigns1");
  
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // Option 1: Jarocinski normalisation of B
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    
    // Option 2: normalise the first row of B
    P               = pick_permutation_bsvars( posterior_B.slice(s), B_hat, allP, allS );
    
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_lambda.slice(s) = Pabs * posterior_lambda.slice(s);
    posterior_df.col(s)       = Pabs * posterior_df.col(s);
    posterior_sigma2.slice(s) = Pabs * posterior_sigma2.slice(s);

  } // END s loop
  
  return List::create(
    _["B"]            = posterior_B,
    _["A"]            = posterior_A,
    _["hyper"]        = posterior_hyper,
    _["sigma"]        = posterior_sigma,
    _["lambda"]       = posterior_lambda,
    _["df"]           = posterior_df,
    _["sigma2"]       = posterior_sigma2,
    _["xi_cpp"]       = posterior_xi_cpp,
    _["PR_TR_cpp"]    = posterior_PR_TR_cpp,
    _["pi_0"]         = posterior_pi_0
  );
} // END normalise_jaro_bsvar_hmsh



// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_hmsh20 (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
    // arma::cube              posterior_xi
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_A       = as<cube>(posterior["A"]);
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_sigma2  = as<cube>(posterior["sigma2"]);
  field<cube> posterior_xi_cpp = as<field<cube>>(posterior["xi_cpp"]);
  // cube   posterior_hyper   = as<cube>(posterior["hyper"]);  // 7x2xS (gamma_0, gamma_+, s_0, s_+, s_)
  // cube   posterior_lambda  = as<cube>(posterior["lambda"]);
  // mat    posterior_df      = as<mat>(posterior["df"]);
  // cube  posterior_pi_0    = as<cube>(posterior["pi_0"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  // Option 1: Jarocinski normalisation of B
  // Function AllPerms("allPermutations");
  // Function AllSigns("allSigns");
  // Function communication_mat("communication_matrix");
  // mat communication = as<mat>(communication_mat(N, N));
  // mat cov_inv(pow(N,2), pow(N,2));
  
  // Option 2: normalise the first row of B
  Function AllPerms("allPermutations1");
  Function AllSigns("allSigns1");
  
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  cube xi_tmp;
  
  for (int s=0; s<S; s++) {
    Rcout << " s: " << s << endl;
    // Option 1: Jarocinski normalisation of B
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    
    // Option 2: normalise the first row of B
    P               = pick_permutation_bsvars( posterior_B.slice(s), B_hat, allP, allS );
    
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_sigma2.slice(s) = Pabs * posterior_sigma2.slice(s);
    // posterior_hyper.slice(s).rows(0, N - 1) = Pabs * posterior_hyper.slice(s).rows(0, N - 1);
    // posterior_hyper.slice(s).rows(N, 2 * N - 1) = Pabs * posterior_hyper.slice(s).rows(N, 2 * N - 1);
    // posterior_lambda.slice(s) = Pabs * posterior_lambda.slice(s);
    // posterior_df.col(s)       = Pabs * posterior_df.col(s);
    // posterior_xi.slice(s)     = posterior_xi(s) * Pabs.t();
    
    uvec permute_s  = index_max(Pabs, 1);
    for (int i=0; i<N; i++) {
      xi_tmp.slice(i) = posterior_xi_cpp(s).slice(permute_s(i));
    }
    posterior_xi_cpp(s) = xi_tmp;
  } // END s loop
  
  // return List::create(
  //   _["B"]            = posterior_B,
  //   _["A"]            = posterior_A,
  //   _["hyper"]        = posterior_hyper,
  //   _["sigma"]        = posterior_sigma,
  //   _["lambda"]       = posterior_lambda,
  //   _["df"]           = posterior_df,
  //   _["sigma2"]       = posterior_sigma2,
  //   _["xi_cpp"]       = posterior_xi_cpp,
  //   _["pi_0"]         = posterior_pi_0
  // );
  return List::create(
    _["B"]            = posterior_B,
    _["A"]            = posterior_A,
    _["sigma"]        = posterior_sigma,
    _["sigma2"]       = posterior_sigma2,
    _["xi_cpp"]       = posterior_xi_cpp
  );
} // END normalise_jaro_bsvar_hmsh20





/*** R

# // Option 1: Jarocinski normalisation of B
# allPermutations = function(N) {
#   return(do.call(rbind, combinat::permn(1:N)))
# }
# 
# allSigns = function(N) {
#   diag.signs    = matrix(NA, 2^N, N)
#   for (n in 1:N) {
#     diag.signs[,n]    = kronecker(c(-1,1), rep(1, 2^(n - 1)))
#   }
#   return(diag.signs)
# }
# 
# communication_matrix = function(N, M) {
#   stopifnot("N must be a positive integer." = N %% 1 == 0 & N > 0)
#   stopifnot("M must be a positive integer." = M %% 1 == 0 & M > 0)
#   
#   NM = N * M
#   
#   return(diag(NM)[as.vector(t(matrix(1:NM, N, M))),])
# }



# // Option 2: normalise the first row of B
allPermutations1 = function(N) {
  out = matrix(NA, N, N)
  out[1,] = 1:N
  for (i in 2:N) {
    out[i,] = c(i, (1:N)[-i])
  }
  return(out)
}

allSigns1 = function(N) {
  as.matrix(
    cbind(
      expand.grid(c(-1,1)),
      matrix(1, 2, N - 1)
    )
  )
}

*/
