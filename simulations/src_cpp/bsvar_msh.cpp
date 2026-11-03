
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppTN)]]
#include <RcppTN.h>

using namespace Rcpp;
using namespace arma;




// [[Rcpp::export]]
arma::mat orthogonal_complement_matrix_TW (const arma::mat& x) {
  // # x is a mxn matrix and m>n
  // # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  // # t(x)%*%out = 0 and det(cbind(x,out))!=0
  int n_nrow     = x.n_rows;
  int n_ncol     = x.n_cols;
  mat Q;
  mat R;
  qr(Q, R, x);
  mat ocm = Q.tail_cols(n_nrow-n_ncol);
  return ocm;
} // END orthogonal_complement_matrix_TW



// [[Rcpp::export]]
arma::vec log_mean (
    arma::mat     log_density     // n x s matrix with log density ordinates
) {
  int S               = log_density.n_cols;
  vec c_log_density   = max(log_density, 1);
  vec log_numerator   = c_log_density - log(S) + log( sum( exp(log_density.each_col() - c_log_density), 1) );
  return log_numerator;
} // log_mean 



// [[Rcpp::export]]
int csample_num1 (
    Rcpp::NumericVector x,
    Rcpp::NumericVector prob = NumericVector::create()
) {
  bool replace = false;
  NumericVector ret = Rcpp::RcppArmadillo::sample(x, 1, replace, prob);
  int out           = ret(0);
  return out;
} // END csample_num1







// [[Rcpp::export]]
arma::vec Ergodic_PR_TR (
    const arma::mat&  PR_TR               // MxM
) {
  const int   M = PR_TR.n_rows;
  
  rowvec  i_mat(M, fill::ones);
  mat     A   = join_cols(eye(M,M)-PR_TR.t(),i_mat);
  mat     aj  = eye(M+1,M+1);
  
  vec     ergodic   = solve(A.t() * A, A.t() * aj.col(M));
  return ergodic;
} // END Ergodic_PR_TR




// [[Rcpp::export]]
arma::mat count_regime_transitions (
    const arma::mat& xi,
    const int M
) {
  const int T = xi.n_cols;
  round(xi);
  mat count(M, M);
  for (int t=1; t<T; t++) {
    count( xi(0,t-1), xi(0,t))++;
  }
  return count;
} // END count_regime_transitions




// [[Rcpp::export]]
arma::rowvec rDirichlet1 (
    const arma::rowvec&   alpha     // Kx1
) {
  const int K   = alpha.size();
  rowvec    draw(K);
  for (int k=0; k<K; k++) {
    draw(k)     = randg(distr_param(alpha(k), 1.0));
  }
  return draw/sum(draw);
} // END rDirichlet1




// [[Rcpp::export]]
arma::rowvec rIG2_Dirichlet1 (
    const arma::rowvec&  s,           // 1xM
    const arma::rowvec&  nu           // 1xM
) {
  const int   M     = s.n_cols;
  rowvec      draw  = s;
  for (int m=0; m<M; m++) {
    draw(m)        /= chi2rnd(nu(m));
  }
  return draw/sum(draw);
} // END rIG2_Dirichlet1




// [[Rcpp::export]]
arma::mat filtering_msh (
    const arma::mat&  U,                  // NxT
    const arma::mat&  sigma,              // NxM - in fact, a matrix with regime-specific VARIANCES!!!
    const arma::mat&  PR_TR,              // MxM
    const arma::vec&  pi_0                // Mx1
) {
  const int   T = U.n_cols;
  const int   N = U.n_rows;
  const int   M = PR_TR.n_rows;
  
  mat         log_eta_t(M, T);
  mat         xi_t_t(M, T);
  
  // This loop evaluates mvnormal pdf at U - simplified operations for zero-mean diagonal-covariance case
  for (int m=0; m<M; m++) {
    rowvec log_d    = -0.5 * sum(pow( pow(sigma.col(m), -0.5) % U.each_col(), 2), 0);
    log_d          += -0.5 * N * log(2*M_PI) - 0.5 * log(prod(sigma.col(m)));
    log_eta_t.row(m) = log_d;
  } // END m loop
  
  vec xi_tm1_tm1    = pi_0;
  
  for (int t=0; t<T; t++) {
    vec log_weight  = log_eta_t.col(t) + log(PR_TR.t() * xi_tm1_tm1);
    log_weight     -= max(log_weight);
    vec weight      = exp(log_weight);
    xi_t_t.col(t)   = weight / sum(weight);
    xi_tm1_tm1      = xi_t_t.col(t);
  } // END t loop
  
  return xi_t_t;
} // END filtering_msh




// [[Rcpp::export]]
arma::mat smoothing_msh (
    const arma::mat&  U,                  // NxT
    const arma::mat&  PR_TR,              // MxM
    const arma::mat&  filtered            // MxT
) {
  const int   T = U.n_cols;
  const int   M = PR_TR.n_rows;
  
  mat   smoothed(M, T);
  smoothed.col(T-1)   = filtered.col(T-1);
  
  for (int t=T-2; t>=0; --t) {
    smoothed.col(t)   = (PR_TR * (smoothed.col(t+1)/(PR_TR.t() * filtered.col(t)) )) % filtered.col(t);
    if (any(smoothed.col(t) < 0) || any(smoothed.col(t) > 1)) {
      for (int m=0; m<M; m++) {
        if (smoothed(m,t) > 1) {smoothed(m,t) = 1;}
        if (smoothed(m,t) < 0) {smoothed(m,t) = 0;}
      }
      smoothed.col(t) = smoothed.col(t) / accu(smoothed.col(t));
    }
  } // END t loop
  
  return smoothed;
} // smoothing_msh





// [[Rcpp::export]]
arma::mat sample_Markov_process_msh (
    arma::mat&        aux_xi,             // 1xT
    const arma::mat&  U,                  // NxT
    const arma::mat&  aux_sigma2,         // NxM
    const arma::mat&  aux_PR_TR,         // MxM
    const arma::vec&  aux_pi_0,          // Mx1
    const bool        finiteM = true
) {
  // the function changes the value of aux_xi by reference (filling it with a new draw)
  
  // const int minimum_regime_occurrences = finiteM ? 3 : 0;
  const int minimum_regime_occurrences = 0;
  const int max_iterations = finiteM ? 10 : 1;
  
  const int   T   = U.n_cols;
  const int   M   = aux_PR_TR.n_rows;
  mat filtered    = filtering_msh(U, aux_sigma2, aux_PR_TR, aux_pi_0);
  mat smoothed    = smoothing_msh(U, aux_PR_TR, filtered);
  
  for (int iteration=0; iteration<max_iterations; iteration++) {
    mat aux_xi_tmp(1, T);
    int draw = csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
    aux_xi_tmp(0,T-1) = draw-1;
    
    for (int t=T-2; t>=0; --t) {
      vec xi_Tmj = aux_PR_TR.col(aux_xi_tmp(0,t+1)) % filtered.col(t);
      xi_Tmj /= sum(xi_Tmj);
      draw = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
      aux_xi_tmp(0,t) = draw-1;
    }
    
    rowvec occurrences(M, fill::zeros);
    for (int t=0; t<T; t++) occurrences(aux_xi_tmp(0,t))++;
    if (!finiteM || min(occurrences) >= minimum_regime_occurrences) {
      aux_xi = aux_xi_tmp;
      break;
    }
  }
  
  return aux_xi;
} // END sample_Markov_process_msh











// [[Rcpp::export]]
arma::cube sample_Markov_process_hmsh (
    arma::cube&       aux_xi,             // 1xTxN
    const arma::mat&  U,                  // NxT
    const arma::mat&  aux_sigma2,         // NxM
    const arma::cube& aux_PR_TR,          // MxMxN
    const arma::mat&  aux_pi_0,           // MxN
    const bool        finiteM = true
) {
  
  // const int minimum_regime_occurrences = finiteM ? 3 : 0;
  const int minimum_regime_occurrences = finiteM ? 3 : 0;
  const int max_iterations = finiteM ? 10 : 1;
  
  const int   T   = U.n_cols;
  const int   N   = U.n_rows;
  const int   M   = aux_PR_TR.n_rows;
  // Rcout << "M: " << M << endl;
  for (int n=0; n<N; n++) {
    mat filtered    = filtering_msh(U.row(n), aux_sigma2.row(n), aux_PR_TR.slice(n), aux_pi_0.col(n));
    mat smoothed    = smoothing_msh(U.row(n), aux_PR_TR.slice(n), filtered);
    
    for (int iteration=0; iteration<max_iterations; iteration++) {
      mat aux_xi_tmp(1, T);
      int draw = csample_num1(wrap(seq_len(M)), wrap(smoothed.col(T-1)));
      aux_xi_tmp(0,T-1) = draw-1;
      
      for (int t=T-2; t>=0; --t) {
        vec xi_Tmj = aux_PR_TR.slice(n).col(aux_xi_tmp(0,t+1)) % filtered.col(t);
        xi_Tmj /= sum(xi_Tmj);
        draw = csample_num1(wrap(seq_len(M)), wrap(xi_Tmj));
        aux_xi_tmp(0,t) = draw-1;
      }
      
      rowvec occurrences(M, fill::zeros);
      for (int t=0; t<T; t++) occurrences(aux_xi_tmp(0,t))++;
      if (!finiteM || min(occurrences) >= minimum_regime_occurrences) {
        aux_xi.slice(n) = aux_xi_tmp;
        break;
      }
    }
  } // END n loop
  
  return aux_xi;
} // END sample_Markov_process_hmsh







// [[Rcpp::export]]
Rcpp::List sample_transition_probabilities (
    arma::mat           aux_PR_TR,    // MxM 
    arma::vec           aux_pi_0,     // Mx1
    const arma::mat&    aux_xi,       // 1xT
    const Rcpp::List&   prior,         // a list of priors - original dimensions
    const bool          MSnotMIX = true
) {
  // Rcout << "t" << endl;
  // the function changes the value of aux_PR_TR and aux_pi_0 by reference (filling it with a new draw)
  const int   M           = aux_PR_TR.n_rows;
  const mat   prior_PR_TR = as<mat>(prior["PR_TR"]);
  // Rcout << "t1" << endl;
  if ( MSnotMIX ) {
    // Rcout << "t2" << endl;
    vec prob_xi1          = aux_pi_0 % aux_PR_TR.col(aux_xi(0,0));
    prob_xi1             /= sum(prob_xi1);
    int S0_draw           = csample_num1(wrap(seq_len(M)), wrap(prob_xi1));
    
    // Rcout << "t3" << endl;
    // Rcout << "M: " << M << endl;
    mat transitions       = count_regime_transitions(aux_xi, M);
    // Rcout << "transitions(S0_draw-1, aux_xi(0,0)): " << transitions(S0_draw-1, aux_xi(0,0)) << endl;
    transitions(S0_draw-1, aux_xi(0,0))++;
    // Rcout << "transitions: " << transitions << endl;
    // Rcout << "prior_PR_TR: " << prior_PR_TR << endl;
    mat posterior_alpha   = transitions + prior_PR_TR;
    // Rcout << "t4" << endl;
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = rDirichlet1(posterior_alpha.row(m));
    }
    // Rcout << "t5" << endl;
    rowvec posterior_alpha_0(M, fill::value((1.0)));
    posterior_alpha_0(S0_draw-1)++;
    aux_pi_0              = trans(rDirichlet1(posterior_alpha_0));
  } else {
    // Rcout << "t6" << endl;
    rowvec occurrences(M, fill::zeros);
    for (uword t=0; t<aux_xi.n_cols; t++) occurrences(aux_xi(0,t))++;
    rowvec posterior_alpha= occurrences + prior_PR_TR.row(0);
    aux_pi_0              = trans(rDirichlet1(posterior_alpha));
    // Rcout << "t7" << endl;
    for (int m=0; m<M; m++) {
      aux_PR_TR.row(m)    = aux_pi_0.t();
    }
  }
  // Rcout << "t8" << endl;
  return List::create(
    _["PR_TR"]        = aux_PR_TR,
    _["pi_0"]         = aux_pi_0
  );
} // END sample_transition_probabilities




// [[Rcpp::export]]
arma::mat sample_variances_msh (
    arma::mat&          aux_sigma2, // NxM
    const arma::mat&    U,          // NxT
    const arma::mat&    aux_xi,     // 1xT state variables
    const Rcpp::List&   prior       // a list of priors - original dimensions
) {
  // the function changes the value of aux_sigma2 by reference (filling it with a new draw)
  const int   M     = aux_sigma2.n_cols;
  const int   N     = U.n_rows;
  const int   T     = U.n_cols;
  const double MM   = M;
  
  rowvec posterior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
  mat posterior_s(N, M);
  posterior_s.fill(prior["sigma_s"]);
  mat sq_resid = square(U); // NxT
  for (int m=0; m<M; m++) {
    for (int t=0; t<T; t++) {
      if (aux_xi(0,t)==m) {
        posterior_nu(m)++;
        posterior_s.col(m) += sq_resid.col(t);
      }
    }
  }
  // This is the version with restriction sum(aux_sigma2,0) = M
  for (int n=0; n<N; n++) {
    aux_sigma2.row(n)     = MM*rIG2_Dirichlet1( posterior_s.row(n), posterior_nu);
  }
  
  return aux_sigma2;
} // END sample_variances_msh



// [[Rcpp::interfaces(cpp, r)]]
// [[Rcpp::export]]
arma::mat sample_variances_hmsh (
    arma::mat&          aux_sigma2, // NxM
    const arma::mat&    U,          // NxT 
    const arma::cube&   aux_xi,     // 1xTxN state variables
    const Rcpp::List&   prior       // a list of priors - original dimensions
) {
  // the function changes the value of aux_sigma2 by reference (filling it with a new draw)
  const int   M     = aux_sigma2.n_cols;
  const int   N     = U.n_rows;
  const int   T     = U.n_cols;
  const double MM   = M;
  
  mat sq_resid = square(U); // NxT
  for (int n=0; n<N; n++) {
    rowvec posterior_s(M);
    posterior_s.fill(prior["sigma_s"]);
    rowvec posterior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
    for (int m=0; m<M; m++) {
      for (int t=0; t<T; t++) {
        if (aux_xi(0,t,n)==m) {
          posterior_nu(m)++;
          posterior_s(m) += sq_resid(n,t);
        }
      }
    }
    
    aux_sigma2.row(n)     = MM*rIG2_Dirichlet1( posterior_s, posterior_nu);
  }
  
  return aux_sigma2;
} // END sample_variances_hmsh








// [[Rcpp::export]]
arma::mat sample_B_heterosk1 (
    arma::mat         aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  aux_sigma,      // NxT conditional STANDARD DEVIATIONS
    const arma::mat&  Y,              // NxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B0 and aux_Bplus by reference (filling it with a new draw)
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  const int radial_degrees  = posterior_nu - N + 1;
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat shocks                = Y;
  
  
  for (int n=0; n<N; n++) {
    
    // set scale matrix
    mat shocks_sigma        = shocks.each_row() / aux_sigma.row(n);
    mat posterior_SS_inv    = (pow(aux_hyper(n,0), -1) * prior_SS_inv) + shocks_sigma * shocks_sigma.t();
    mat posterior_S_inv     = VB(n) * posterior_SS_inv * VB(n).t();
    posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
    
    // sample B
    mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u(radial_degrees, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = pow(as_scalar(sum(pow(u,2))), 0.5);
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)           = b0n * VB(n);
  } // END n loop
  
  return aux_B;
} // END sample_B_heterosk1







// [[Rcpp::export]]
Rcpp::List bsvar_msh_cpp (
    int&              S,              // No. of posterior draws
    arma::mat&        Y,              // NxT dependent variables
    Rcpp::List&       prior,          // a list of priors - original dimensions
    arma::field<arma::mat>& VB,       // restrictions on B0
    Rcpp::List&       starting_values,
    const int               thin = 100,     // introduce thinning
    const bool              finiteM = true,
    const bool              MSnotMIX = true,
    const bool              show_progress = true
) {
  
  // Rcout<<"0";
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const vec adptive_alpha_gamma = as<vec>(NumericVector::create(0.44, 0.6));
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  // const int   K     = X.n_rows;
  // Rcout<<"1";
  mat   aux_B       = as<mat>(starting_values["B"]);
  // mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_sigma2  = as<mat>(starting_values["sigma2"]);
  mat   aux_sigma(N, T);
  mat   aux_PR_TR   = as<mat>(starting_values["PR_TR"]);
  vec   aux_pi_0    = as<vec>(starting_values["pi_0"]);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  // mat   aux_lambda  = as<mat>(starting_values["lambda"]);
  // mat   aux_lambda_sqrt(N, T, fill::ones);
  mat   aux_hetero(N, T, fill::ones);
  // vec   aux_df      = as<vec>(starting_values["df"]);
  mat   U;
  
  const int   M     = aux_PR_TR.n_rows;
  
  const int   SS     = floor(S / thin);
  // Rcout<<"2";
  cube  posterior_B(N, N, SS);
  // cube  posterior_A(N, K, SS);
  cube  posterior_sigma2(N, M, SS);
  cube  posterior_PR_TR(M, M, SS);
  mat   posterior_pi_0(M, SS);
  cube  posterior_xi(1, T, SS);
  // cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  posterior_sigma(N, T, SS);
  // cube  posterior_lambda(N, T, SS);
  // mat   posterior_df(N, SS);
  // Rcout<<"3";
  int   ss = 0;
  for (int t=0; t<T; t++) {
    aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi(0,t)) , 0.5 );
  }
  
  aux_hetero = aux_sigma;
  
  // Rcout<<"4";
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    aux_B             = sample_B_heterosk1(aux_B, aux_hyper, aux_hetero, Y, prior, VB);
    
    U                 = aux_B * Y / aux_sigma;
    aux_xi            = sample_Markov_process_msh(aux_xi, U, aux_sigma2, aux_PR_TR, aux_pi_0, finiteM);
    // Rcout<<"8";
    // sample aux_PR_TR
    List aux_PR_tmp   = sample_transition_probabilities(aux_PR_TR, aux_pi_0, aux_xi, prior, MSnotMIX);
    aux_PR_TR         = as<mat>(aux_PR_tmp["PR_TR"]);
    aux_pi_0          = as<vec>(aux_PR_tmp["pi_0"]);
    // Rcout<<"9";
    // sample aux_sigma2
    U                 = aux_B * Y;
    aux_sigma2        = sample_variances_msh(aux_sigma2, U, aux_xi, prior);
    
    for (int t=0; t<T; t++) {
      aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi(0,t)) , 0.5 );
    }
    aux_hetero      = aux_sigma;
    
    // Rcout<<"10";
    if (s % thin == 0) {
      posterior_B.slice(ss)      = aux_B;
      posterior_sigma2.slice(ss) = aux_sigma2;
      posterior_PR_TR.slice(ss)  = aux_PR_TR;
      posterior_pi_0.col(ss)     = aux_pi_0;
      posterior_xi.slice(ss)     = aux_xi;
      posterior_sigma.slice(ss)  = aux_sigma;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["sigma2"]   = aux_sigma2,
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi,
      _["hyper"]    = aux_hyper,
      _["sigma"]    = aux_sigma
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["sigma2"]   = posterior_sigma2,
      _["PR_TR"]    = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi"]       = posterior_xi,
      _["sigma"]    = posterior_sigma
    )
  );
} // END bsvar_msh



// [[Rcpp::export]]
Rcpp::List bsvar_hmsh_cpp (
    const int&              S,              // No. of posterior draws
    const arma::mat&        Y,              // NxT dependent variables
    const Rcpp::List&       prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB,       // restrictions on B0
    const Rcpp::List&       starting_values,
    const int               thin = 100,     // introduce thinning
    const bool              finiteM = true,
    const bool              MSnotMIX = true,
    const bool              show_progress = true
) {
  Rcout << "0" << endl;
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const vec adptive_alpha_gamma = as<vec>(NumericVector::create(0.44, 0.6));
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  // const int   K     = X.n_rows;
  // Rcout << "1" << endl;
  mat   aux_B       = as<mat>(starting_values["B"]);
  // mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_sigma2  = as<mat>(starting_values["sigma2"]);     // NxM
  mat   aux_sigma(N, T);
  cube  aux_PR_TR   = as<cube>(starting_values["PR_TR"]);     // MxMxN
  mat   aux_pi_0    = as<mat>(starting_values["pi_0"]);       // MxN
  cube  aux_xi      = as<cube>(starting_values["xi"]);        // MxTxN
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  // mat   aux_lambda  = as<mat>(starting_values["lambda"]);
  // mat   aux_lambda_sqrt(N, T, fill::ones);
  mat   aux_hetero(N, T, fill::ones);
  // vec   aux_df      = as<vec>(starting_values["df"]);
  mat   U;
  // Rcout << "2" << endl;
  const int   M     = aux_PR_TR.n_rows;
  
  const int   SS     = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  cube  posterior_sigma2(N, M, SS);
  field<cube> posterior_PR_TR(SS);                // (S)(M, M, N) 
  cube  posterior_pi_0(M, N, SS);
  field<cube> posterior_xi(SS);                   // (S)(M, T, N)
  cube  posterior_sigma(N, T, SS);
  // Rcout << "3" << endl;
  int   ss = 0;
  for (int t=0; t<T; t++) {
    for (int n=0; n<N; n++) {
      aux_sigma(n,t)  = sqrt( aux_sigma2(n, aux_xi(0,t,n)) );
    }
  }
  
  aux_hetero = aux_sigma;
  // Rcout << "4" << endl;
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_B
    // Rcout << "5" << endl;
    aux_B             = sample_B_heterosk1(aux_B, aux_hyper, aux_hetero, Y, prior, VB);
    
    // sample aux_xi
    // Rcout << "5" << endl;
    U                 = aux_B * Y / aux_sigma;
    aux_xi            = sample_Markov_process_hmsh(aux_xi, U, aux_sigma2, aux_PR_TR, aux_pi_0, finiteM);
    
    // sample aux_PR_TR
    // Rcout << "7" << endl;
    for (int n=0; n<N; n++) {
      List aux_PR_tmp     = sample_transition_probabilities(aux_PR_TR.slice(n), aux_pi_0.col(n), aux_xi.slice(n), prior, MSnotMIX);
      // Rcout << "7a" << endl;
      aux_PR_TR.slice(n)  = as<mat>(aux_PR_tmp["PR_TR"]);
      aux_pi_0.col(n)     = as<vec>(aux_PR_tmp["pi_0"]);
    }
    
    // sample aux_sigma2
    // Rcout << "8" << endl;
    U                 = aux_B * Y ;
    aux_sigma2        = sample_variances_hmsh(aux_sigma2, U, aux_xi, prior);
    // Rcout << "9" << endl;
    for (int t=0; t<T; t++) {
      for (int n=0; n<N; n++) {
        aux_sigma(n,t)  = sqrt( aux_sigma2(n, aux_xi(0,t,n)) );
      }
    }
    aux_hetero      = aux_sigma;
    // Rcout << "10" << endl;
    if (s % thin == 0) {
      posterior_B.slice(ss)      = aux_B;
      posterior_sigma2.slice(ss) = aux_sigma2;
      posterior_PR_TR(ss)        = aux_PR_TR;
      posterior_pi_0.slice(ss)   = aux_pi_0;
      posterior_xi(ss)           = aux_xi;
      posterior_sigma.slice(ss)  = aux_sigma;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["sigma2"]   = aux_sigma2,
      _["PR_TR"]    = aux_PR_TR,
      _["pi_0"]     = aux_pi_0,
      _["xi"]       = aux_xi,
      _["hyper"]    = aux_hyper,
      _["sigma"]    = aux_sigma
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["sigma2"]   = posterior_sigma2,
      _["PR_TR_cpp"]  = posterior_PR_TR,
      _["pi_0"]     = posterior_pi_0,
      _["xi_cpp"]   = posterior_xi,
      _["sigma"]    = posterior_sigma
    )
  );
} // END bsvar_msh











// [[Rcpp::export]]
Rcpp::List bsvar_exh_cpp (
    int&              S,              // No. of posterior draws
    arma::mat&        Y,              // NxT dependent variables
    Rcpp::List&       prior,          // a list of priors - original dimensions
    arma::field<arma::mat>& VB,       // restrictions on B0
    Rcpp::List&       starting_values,
    const int               thin = 100,     // introduce thinning
    const bool              show_progress = true
) {
  
  // Rcout<<"0";
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  
  const int   T     = Y.n_cols;
  const int   N     = Y.n_rows;
  // const int   K     = X.n_rows;
  // Rcout<<"1";
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_sigma2  = as<mat>(starting_values["sigma2"]);
  mat   aux_sigma(N, T);
  mat   aux_xi      = as<mat>(starting_values["xi"]);
  mat   aux_hetero(N, T, fill::ones);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  mat   U;
  
  const int   M     = aux_xi.max() + 1;
  const int   SS     = floor(S / thin);
  // Rcout<<"2";
  cube  posterior_B(N, N, SS);
  // cube  posterior_A(N, K, SS);
  cube  posterior_sigma2(N, M, SS);
  // cube  posterior_PR_TR(M, M, SS);
  // mat   posterior_pi_0(M, SS);
  cube  posterior_xi(1, T, SS);
  // cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  posterior_sigma(N, T, SS);
  
  // Rcout<<"3";
  int   ss = 0;
  for (int t=0; t<T; t++) {
    aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi(0,t)) , 0.5 );
  }
  
  aux_hetero = aux_sigma;
  
  // Rcout<<"4";
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_B
    // Rcout<<"5";
    aux_B             = sample_B_heterosk1(aux_B, aux_hyper, aux_hetero, Y, prior, VB);
    
    // Rcout<<"6";
    U                 = aux_B * Y;
    mat PR_TR         = count_regime_transitions(aux_xi, M);
    // Rcout<<"7";
    PR_TR             = diagmat(sum(PR_TR, 1)) * PR_TR;
    aux_sigma2        = sample_variances_msh(aux_sigma2, U, aux_xi, prior);
    
    // Rcout<<"8";
    for (int t=0; t<T; t++) {
      aux_sigma.col(t)    = pow( aux_sigma2.col(aux_xi(0,t)) , 0.5 );
    }
    aux_hetero      = aux_sigma;
    
    // Rcout<<"10";
    if (s % thin == 0) {
      posterior_B.slice(ss)      = aux_B;
      posterior_sigma2.slice(ss) = aux_sigma2;
      posterior_xi.slice(ss)     = aux_xi;
      posterior_sigma.slice(ss)  = aux_sigma;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["sigma2"]   = aux_sigma2,
      _["xi"]       = aux_xi,
      _["hyper"]    = aux_hyper,
      _["sigma"]    = aux_sigma
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["sigma2"]   = posterior_sigma2,
      _["xi"]       = posterior_xi,
      _["sigma"]    = posterior_sigma
    )
  );
} // END bsvar_exh











// [[Rcpp::export]]
double dig2dirichlet (
    const arma::rowvec&  x,                // M-vector of positive rv summing up to 1
    const arma::rowvec&  a,                // M-vector of positive parameters
    const arma::rowvec&  b,                // M-vector of positive parameters
    const bool  logarithm = true  
) {
  // density ordinate of ig2-based Dirichlet distribution
  double constant_ig2d    = lgamma(accu(a) / 2) - accu(lgamma(a  / 2)) - accu(log(b));
  double kernel_ig2d      = accu(0.5 * (a + 2) % (log(b) - log(x))) - ( 0.5 * accu(a) * log(accu(b / x)) );
  double out              = constant_ig2d + kernel_ig2d;
  
  if ( !logarithm ) {
    out =  exp(out);
  }
  return out;
} // END dig2dirichlet




// [[Rcpp::export]]
double ddirichlet (
    const arma::rowvec&  x,                // M-vector of positive rv summing up to 1
    const arma::rowvec&  a,                // M-vector of positive parameters
    const bool  logarithm = true  
) {
  // density ordinate of Dirichlet distribution
  double constant_d       = lgamma(accu(a)) - accu(lgamma(a));
  double kernel_d         = accu((a - 1) % log(x));
  double out              = constant_d + kernel_d;
  
  if ( !logarithm ) {
    out =  exp(out);
  }
  return out;
} // END ddirichlet




// [[Rcpp::export]]
arma::vec verify_volatility_msh_cpp (
    const Rcpp::List&       posterior,  // a list of posteriors
    const Rcpp::List&       prior,      // a list of priors - original dimensions
    const arma::mat&        Y
) {
  
  cube  posterior_sigma2    = posterior["sigma2"];
  cube  posterior_B         = posterior["B"];
  cube  posterior_xi        = posterior["xi"];
  
  const int   M             = posterior_sigma2.n_cols;
  double      MM            = posterior_sigma2.n_cols;
  const int   N             = posterior_B.n_rows;
  const int   T             = Y.n_cols;
  const int   S             = posterior_B.n_slices;
  
  rowvec      homoskedasticity_hypothesis(M, fill::value(1 / MM));
  
  // compute denominator
  
  rowvec  prior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
  rowvec  prior_s(M, fill::value(as<double>(prior["sigma_s"])));
  
  double  log_denominator = dig2dirichlet( homoskedasticity_hypothesis, prior_nu, prior_s );
  
  // compute numerator
  mat     log_numerator_s(1, S);
  for (int s = 0; s < S; s++) {
    // for (int n = 0; n < N; n++) {
    int n = 0;
    rowvec posterior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
    
    mat posterior_s(N, M);
    posterior_s.fill(prior["sigma_s"]);
    for (int m=0; m<M; m++) {
      for (int t=0; t<T; t++) {
        if (posterior_xi(0,t,s)==m) {
          posterior_nu(m)++;
          posterior_s.col(m) += square(posterior_B.slice(s) * Y.col(t));
        }
      }
    }
    
    log_numerator_s(n,s)  = dig2dirichlet( homoskedasticity_hypothesis, posterior_nu, posterior_s.row(n) );
    
    // } // END n loop
  } // END s loop
  
  // compute the log of the mean numerator exp(log_numerator)
  vec log_numerator   = log_mean(log_numerator_s);
  vec out             = log_numerator - log_denominator;
  
  // compute the standard error 
  return out;
} // END verify_volatility_msh_cpp








// [[Rcpp::export]]
arma::vec verify_volatility_hmsh_cpp (
    const Rcpp::List&       posterior,  // a list of posteriors
    const Rcpp::List&       prior,      // a list of priors - original dimensions
    const arma::mat&        Y
) {
  
  cube  posterior_sigma2    = as<cube>(posterior["sigma2"]);
  cube  posterior_B         = as<cube>(posterior["B"]);
  field<cube> posterior_xi  = as<field<cube>>(posterior["xi_cpp"]);
  
  const int   M             = posterior_sigma2.n_cols;
  double      MM            = posterior_sigma2.n_cols;
  const int   N             = posterior_B.n_rows;
  const int   T             = Y.n_cols;
  const int   S             = posterior_B.n_slices;
  
  rowvec      homoskedasticity_hypothesis(M, fill::value(1 / MM));
  
  // compute denominator
  
  rowvec  prior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
  rowvec  prior_s(M, fill::value(as<double>(prior["sigma_s"])));
  
  double  log_denominator = dig2dirichlet( homoskedasticity_hypothesis, prior_nu, prior_s );
  
  // compute numerator
  mat     log_numerator_s(1, S);
  for (int s = 0; s < S; s++) {
    // for (int n = 0; n < N; n++) {
    int n = 0;
    mat post_xi           = posterior_xi(s).slice(n);
    rowvec posterior_nu(M, fill::value(as<double>(prior["sigma_nu"])));
    
    mat posterior_s(N, M);
    posterior_s.fill(prior["sigma_s"]);
    for (int m=0; m<M; m++) {
      for (int t=0; t<T; t++) {
        if (post_xi(0,t)==m) {
          posterior_nu(m)++;
          posterior_s.col(m) += square(posterior_B.slice(s) * Y.col(t) );
        }
      }
    }
    
    log_numerator_s(n,s)  = dig2dirichlet( homoskedasticity_hypothesis, posterior_nu, posterior_s.row(n) );
    
    // } // END n loop
  } // END s loop
  
  // compute the log of the mean numerator exp(log_numerator)
  vec log_numerator   = log_mean(log_numerator_s);
  vec out             = log_numerator - log_denominator;
  // compute the standard error 
  return out;
} // END verify_volatility_hmsh_cpp

