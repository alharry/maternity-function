#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type inv_logit(Type x){return 1/(1 + exp(- x));}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  
  PARAMETER(l50);
  PARAMETER(l95);
  
  vector<Type> X = log(19) * ((x - l50) / (l95 - l50));
  vector<Type> p = inv_logit(X);
  
  
  Type nll=-sum(dbinom(y,Type(1.0),p,true));
  
  SIMULATE {
    y = rbinom(Type(1.0), p);
    REPORT(y);
  }

  
  return nll;
}

