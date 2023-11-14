#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type inv_logit(Type x){return 1/(1 + exp(- x));}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(z);
  DATA_VECTOR(x);

  PARAMETER(m50);
  PARAMETER(m95);
  PARAMETER(c);
  
  vector<Type> X = log(19) * ((x - m50) / (m95 - m50));
  vector<Type> p = c * inv_logit(X);
    
    
  Type nll=-sum(dbinom(z,Type(1.0),p,true));
  
  SIMULATE {
    z = rbinom(Type(1.0), p);
    REPORT(z);
  }
  
  return nll;
}

