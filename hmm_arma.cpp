// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//  [[Rcpp::export]]
double C_CalcNegLlk(const arma::vec wpar, 
                  const arma::vec data, 
                  const int n_states) {
  
  // unpack and transform parameters  
  arma::vec lambda(wpar.rows(0, n_states - 1));
  lambda = exp(lambda); 
  // lambda is computed cumulatively to prevent label-switching 
  lambda = cumsum(lambda);
  arma::mat tpm(n_states, n_states);
  int cur = n_states; 
  for (int i = 0; i < n_states; ++i) {
    tpm(i, i) = 1; 
    for (int j = 0; j < n_states; ++j) {
       if (i != j) {
         tpm(i, j) = exp(wpar(cur)); 
         ++cur; 
       }
    }
    tpm.row(i) /= accu(tpm.row(i)); 
  } 
  // compute stationary distribution
  arma::rowvec delta(n_states); 
  arma::mat I = arma::eye<arma::mat>(n_states, n_states); 
  arma::mat tpminv = I; 
  tpminv -= tpm; 
  tpminv += 1; 
  arma::rowvec ivec = arma::ones<arma::rowvec>(n_states); 
  // if tpm is ill-conditioned then just use uniform initial distribution 
  try {
    tpminv = inv(tpminv);
    delta = ivec * tpminv;
  } catch(...) {
    delta.ones(); 
    delta /= n_states; 
  }
  // compute observation probabilities 
  int n = data.n_rows; 
  arma::mat prob(n, n_states); 
  for (int s = 0; s < n_states; ++s) {
    for (int i = 0; i < n; ++i) {
      prob(i, s) = R::dpois(data(i), lambda(s), false); 
    }
  }
  // compute log-likelihood 
  double llk = 0; 
  arma::rowvec phi(delta); 
  double sumphi = 0; 
  for (int i = 0; i < n; ++i) {
    phi %= prob.row(i); 
    phi *= tpm;
    sumphi = accu(phi); 
    llk += log(sumphi);  
    phi /= sumphi; 
  }
  double nll = -llk; 
  return nll; 
}