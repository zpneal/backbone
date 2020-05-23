#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double scobit_loglike_cpp(NumericVector x1, NumericVector x2, NumericVector y,
                          NumericVector params) {
  int n = x1.length();
  double b0 = params[0];
  double b1 = params[1];
  double b2 = params[2];
  // double b3 = params[3];
  double a  = params[3];
  double f = 0;

  for(int i=0; i<n;i++){
    // f+= (1 - y[i])*log(pow(1 + exp(b0+b1*x1[i]+b2*x2[i]+b3*x1[i]*x2[i]),-a))+y[i]*log(1-pow(1+exp(b0+b1*x1[i]+b2*x2[i]+b3*x1[i]*x2[i]),-a));
    f+= (1 - y[i])*log(pow(1 + exp(b0+b1*x1[i]+b2*x2[i]),-a))+y[i]*log(1-pow(1+exp(b0+b1*x1[i]+b2*x2[i]),-a));
  }
  if(f!=f){
    f=-10000000000;
  }
  return -f;
}

// [[Rcpp::export]]
NumericVector scobit_loglike_gr_cpp(NumericVector x1, NumericVector x2, NumericVector y,
                                    NumericVector params){
  int n = x1.length();
  double b0 = params[0];
  double b1 = params[1];
  double b2 = params[2];
  // double b3 = params[3];
  // double a  = params[4];
  double a  = params[3];
  // double f = 0;
  NumericVector Fexp(n);
  NumericVector Fexp1a(n);
  for(int i=0; i<n;i++){
    // Fexp[i]   = exp(b0+b1*x1[i]+b2*x2[i]+b3*x1[i]*x2[i]);
    // Fexp1a[i] = pow(1+exp(b0+b1*x1[i]+b2*x2[i]+b3*x1[i]*x2[i]),-a);

    Fexp[i]   = exp(b0+b1*x1[i]+b2*x2[i]);
    Fexp1a[i] = pow(1+exp(b0+b1*x1[i]+b2*x2[i]),-a);
  }

  NumericVector derivFb0(n);
  NumericVector derivFb1(n);
  NumericVector derivFb2(n);
  // NumericVector derivFb3(n);
  NumericVector derivFa(n);
  for(int i=0; i<n;i++){
    derivFb0[i] = -a*Fexp[i]*pow(1+Fexp[i],-(a+1));
    derivFb1[i] = -a*Fexp[i]*pow(1+Fexp[i],-(a+1))*x1[i];
    derivFb2[i] = -a*Fexp[i]*pow(1+Fexp[i],-(a+1))*x2[i];
    // derivFb3[i] = -a*Fexp[i]*pow(1+Fexp[i],-(a+1))*x1[i]*x2[i];
    derivFa[i]  = -log(1+Fexp[i])*Fexp1a[i];
  }
  NumericVector derivLb0(n);
  NumericVector derivLb1(n);
  NumericVector derivLb2(n);
  NumericVector derivLb3(n);
  NumericVector derivLa(n);
  for(int i=0; i<n;i++){
    derivLb0[i] =-1 * (-y[i]/(1-Fexp1a[i]) + (1-y[i])/(Fexp1a[i])) * derivFb0[i];
    derivLb1[i] =-1 * (-y[i]/(1-Fexp1a[i]) + (1-y[i])/(Fexp1a[i])) * derivFb1[i];
    derivLb2[i] =-1 * (-y[i]/(1-Fexp1a[i]) + (1-y[i])/(Fexp1a[i])) * derivFb2[i];
    // derivLb3[i] =-1 * (-y[i]/(1-Fexp1a[i]) + (1-y[i])/(Fexp1a[i])) * derivFb3[i];
    derivLa[i]  =-1 * (-y[i]/(1-Fexp1a[i]) + (1-y[i])/Fexp1a[i]) * derivFa[i];
  }
  NumericVector res(4);
  res[0] = sum(derivLb0);
  res[1] = sum(derivLb1);
  res[2] = sum(derivLb2);
  // res[3] = sum(derivLb3);
  // res[4] = sum(derivLa);
  res[3] = sum(derivLa);

  return res;
}
