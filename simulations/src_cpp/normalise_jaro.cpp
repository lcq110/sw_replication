#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double normalise_jaro_distance (
    mat&          B,              // NxN
    const mat&    B_hat,          // NxN
    const mat&    sigma2,         // NxT
    const int     T,
    const mat&    communication   // NNxNN communication matrix
) {
  const int N     = B.n_rows; 
  vec   vecBt     = vectorise(B.t());
  vec   vecB_hatt = vectorise(B_hat.t());
  
  mat   B_hat_inv = inv(B_hat);
  vec   sigma2_sum = sum(sigma2, 1);
  
  mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
  mat   covariance_inv  = T * Finsher_info;
    
  double distance = as_scalar(trans(vecBt - vecB_hatt) * covariance_inv * (vecBt - vecB_hatt));
  
  return distance;
} // END normalise_jaro_distance
    

    
// [[Rcpp::export]]
mat covariance_inv_bsvars (
    const mat&    B_hat,          // NxN
    const mat&    sigma2,         // NxT
    const int     T,
    const mat&    communication   // NNxNN communication matrix
) {
  const int N     = B_hat.n_rows; 
  vec   vecB_hatt = vectorise(B_hat.t());
  mat   B_hat_inv = inv(B_hat);
  vec   sigma2_sum = sum(sigma2, 1);
  
  mat   Finsher_info    = communication * kron(B_hat_inv, B_hat_inv.t()) + kron(eye(N, N), B_hat_inv * diagmat(sigma2_sum / T) * B_hat_inv.t());
  mat   covariance_inv  = T * Finsher_info;
  
  return covariance_inv;
} // END covariance_inv_bsvars


// [[Rcpp::export]]
mat pick_permutation_bsvars (
    arma::mat&              B,
    const arma::mat&        B_hat,
    const int               T,              // sample size
    const arma::mat&        cov_inv,
    const arma::umat&       perms,
    const arma::mat&        signs,
    const arma::mat&        communication   // NNxNN communication matrix
) {

  const int   no_perms    = perms.n_rows;
  const int   no_signs    = signs.n_rows;
  const int   N           = B.n_rows;
  
  mat   I                 = eye(N,N);
  vec   vecB_hatt         = vectorise(B_hat.t());
  vec   vecBt(vecB_hatt.n_elem);
  
  mat   distance(no_perms, no_signs);
  
  for (int p=0; p<no_perms; p++) {
    for (int s=0; s<no_signs; s++) {
      vecBt               = vectorise(trans(I.rows(perms.row(p)) * diagmat(signs.row(s)) * B));
      distance(p, s)      = as_scalar(trans(vecBt - vecB_hatt) * cov_inv * (vecBt - vecB_hatt));
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
mat pick_permutation_sse (
    arma::mat&              B,
    const arma::mat&        B_hat,
    const arma::umat&       perms,
    const arma::mat&        signs,
    const arma::mat&        communication   // NNxNN communication matrix
) {
  
  const int   no_perms    = perms.n_rows;
  const int   no_signs    = signs.n_rows;
  const int   N           = B.n_rows;
  mat   I                 = eye(N,N);
  mat   Bt;
  mat   distance(no_perms, no_signs);
  
  for (int p=0; p<no_perms; p++) {
    for (int s=0; s<no_signs; s++) {
      Bt               = I.rows(perms.row(p)) * diagmat(signs.row(s)) * B;
      distance(p, s)      = accu(square(Bt - B_hat));
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
Rcpp::List normalise_jaro_bsvar_msh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {

  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_sigma2   = as<cube>(posterior["sigma2"]);
  cube   posterior_PR_TR   = as<cube>(posterior["PR_TR"]);
  cube   posterior_xi   = as<cube>(posterior["xi"]);
  mat    posterior_pi_0       = as<mat>(posterior["pi_0"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;

  Function AllPerms("allPermutations");
  Function AllSigns("allSigns");
  Function communication_mat("communication_matrix");

  mat communication = as<mat>(communication_mat(N, N));
  
  mat cov_inv(pow(N,2), pow(N,2));
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    P               = pick_permutation_sse (posterior_B.slice(s), B_hat, allP, allS, communication);
    Pabs            = abs(P);

    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_sigma2.slice(s)  = Pabs * posterior_sigma2.slice(s);
    
  } // END s loop

  
  return List::create(
    _["B"]            = posterior_B,
    _["sigma"]        = posterior_sigma,
    _["sigma2"]        = posterior_sigma2,
    _["PR_TR"]        = posterior_PR_TR,
    _["xi"]        = posterior_xi,
    _["pi_0"]           = posterior_pi_0
  );
} // END normalise_jaro_bsvar_sv






// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_hmsh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_sigma2   = as<cube>(posterior["sigma2"]);
  field<cube> posterior_xi  = as<field<cube>>(posterior["xi_cpp"]);
  field<cube> posterior_PR_TR = as<field<cube>>(posterior["PR_TR_cpp"]);
  cube posterior_pi_0 = as<cube>(posterior["pi_0"]);
  
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  Function AllPerms("allPermutations");
  Function AllSigns("allSigns");
  Function communication_mat("communication_matrix");
  
  mat communication = as<mat>(communication_mat(N, N));
  
  mat cov_inv(pow(N,2), pow(N,2));
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    P               = pick_permutation_sse (posterior_B.slice(s), B_hat, allP, allS, communication);
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_sigma2.slice(s)  = Pabs * posterior_sigma2.slice(s);
    
    vec ind = cumsum(ones(N)) - 1;
    vec ord = Pabs * ind;
    cube post_xi      = posterior_xi(s);
    cube post_target  = post_xi;
    for (int n=0; n<N; n++) {
      post_target.slice(n) = post_xi.slice(ord(n));
    }
    posterior_xi(s)   = post_target;
  } // END s loop
  
  
  return List::create(
    _["B"]            = posterior_B,
    _["sigma"]        = posterior_sigma,
    _["sigma2"]        = posterior_sigma2,
    _["xi_cpp"]        = posterior_xi,
    _["PR_TR_cpp"]  = posterior_PR_TR,
    _["pi_0"]     = posterior_pi_0
  );
} // END normalise_jaro_bsvar_sv







// [[Rcpp::export]]
Rcpp::List normalise_jaro_bsvar_exh (
    Rcpp::List&             posterior,
    const arma::mat&        B_hat
) {
  
  cube   posterior_B       = as<cube>(posterior["B"]);
  cube   posterior_sigma   = as<cube>(posterior["sigma"]);
  cube   posterior_sigma2   = as<cube>(posterior["sigma2"]);
  cube   posterior_xi   = as<cube>(posterior["xi"]);
  
  const int   N     = posterior_B.n_rows;
  const int   T     = posterior_sigma.n_cols;
  const int   S     = posterior_B.n_slices;
  
  Function AllPerms("allPermutations");
  Function AllSigns("allSigns");
  Function communication_mat("communication_matrix");
  
  mat communication = as<mat>(communication_mat(N, N));
  
  mat cov_inv(pow(N,2), pow(N,2));
  umat allP         = as<umat>(AllPerms(N)) - 1;
  mat allS          = as<mat>(AllSigns(N));
  
  mat P(N, N);
  mat Pabs(N, N);
  
  for (int s=0; s<S; s++) {
    // mat sigma2      = square(posterior_sigma.slice(s));
    // cov_inv         = covariance_inv_bsvars(B_hat, sigma2, T, communication);
    // P               = pick_permutation_bsvars(posterior_B.slice(s), B_hat, T, cov_inv, allP, allS, communication);
    P               = pick_permutation_sse (posterior_B.slice(s), B_hat, allP, allS, communication);
    Pabs            = abs(P);
    
    posterior_B.slice(s)      = P * posterior_B.slice(s);
    posterior_sigma.slice(s)  = Pabs * posterior_sigma.slice(s);
    posterior_sigma2.slice(s)  = Pabs * posterior_sigma2.slice(s);
  
  } // END s loop
  
  
  return List::create(
    _["B"]            = posterior_B,
    _["sigma"]        = posterior_sigma,
    _["sigma2"]        = posterior_sigma2,
    _["xi"]        = posterior_xi
  );
} // END normalise_jaro_bsvar_sv












/*** R


allPermutations = function(N) {
  return(rbind(1:N, permute::allPerms(N)))
}

allSigns = function(N) {
  diag.signs    = matrix(NA, 2^N, N)
  for (n in 1:N) {
    diag.signs[,n]    = kronecker(c(-1,1), rep(1, 2^(n - 1)))
  }
  return(diag.signs)
}

communication_matrix = function(N, M) {
  stopifnot("N must be a positive integer." = N %% 1 == 0 & N > 0)
  stopifnot("M must be a positive integer." = M %% 1 == 0 & M > 0)
  
  NM = N * M
  
  return(diag(NM)[as.vector(t(matrix(1:NM, N, M))),])
}

*/
