#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in data 
  DATA_VECTOR(data);
  DATA_INTEGER(n_states);
  // Read in parameters
  PARAMETER_VECTOR(wpar);
  // Unpack and transform parameters
  vector<Type> lambda(wpar.head(n_states));
  lambda = exp(lambda);
  for (int i = 0; i < n_states - 1; ++i) lambda(i + 1) += lambda(i);  
  matrix<Type> tpm(n_states, n_states);
  int cur = n_states; 
  for (int i = 0; i < n_states; ++i) {
    tpm(i, i) = 1; 
   for (int j = 0; j < n_states; ++j) {
    if (i != j) {
      tpm(i, j) = exp(wpar[cur]); 
      ++cur; 
    } 
   } 
   tpm.row(i) /= tpm.row(i).sum();
  } 
  // Compute stationary distribution 
  matrix<Type> I = matrix<Type>::Identity(n_states, n_states);
  matrix<Type> tpminv = I; 
  tpminv -= tpm; 
  tpminv = (tpminv.array() + 1).matrix(); 
  matrix<Type> ivec(1, n_states); for (int i = 0; i < n_states; ++i) ivec(0, i) = 1; 
  tpminv = tpminv.inverse();
  matrix<Type> delta = ivec * tpminv;
  // compute observation probabilities
  int n = data.rows();  
  matrix<Type> prob(n, n_states); 
  for (int s = 0; s < n_states; ++s) {
    for (int i = 0; i < n; ++i) {
      prob(i, s) = dpois(data(i), lambda[s]); 
    }
  } 
  // compute log-likelihood 
  Type llk = 0;
  matrix<Type> phi(delta);
  Type sumphi = 0;
  for (int i = 0; i < n; ++i) {
    phi = (phi.array() * prob.row(i).array()).matrix();  
    phi = phi * tpm;
    sumphi = phi.sum();
    llk += log(sumphi);
    phi /= sumphi;
  }
  Type nll = -llk; 
  //std::cout << "llk: " << llk << std::endl; 
  return nll;
}
