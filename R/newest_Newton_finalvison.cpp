#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <math.h>
#include <iostream>

  
  
using namespace std;
using namespace Rcpp;
using namespace arma;

int modM(int a, int n)
{
  return a - (int)(floor(a / n) * n);
}

arma::vec findsign(arma::mat X, arma::vec gfunction, arma::vec XTbeta, double h, arma::vec XTalpha, bool usealpha){
  int n = X.n_cols;
  //int p = X.n_rows;
  arma::vec resultvector = arma::zeros(n);
  arma::vec u = arma::ones(n);
  if (usealpha) {
    for (int i = 0; i < n; i++) {
      double xalpha = XTalpha(i);
      arma::vec kernalvec = arma::normpdf( (xalpha*u -XTalpha) / h) / h;
      resultvector(i) = dot(kernalvec, gfunction) / sum(kernalvec);
    }
  }else{
    for (int i = 0; i < n; i++) {
      double xbeta = XTbeta(i);
      arma::vec kernalvec = arma::normpdf( (xbeta*u -XTbeta) / h) / h;
      resultvector(i) = dot(kernalvec, gfunction) / sum(kernalvec);
    }
  }
  return arma::sign(resultvector);
}


arma::vec findsign_test(arma::mat Xtest, arma::vec gfunction, arma::vec XTbeta, double h, arma::vec XTalpha, bool usealpha, arma::vec alpha, arma::vec beta){
  int n1 = Xtest.n_cols;
  int n = gfunction.n_elem;
  arma::vec resultvector = arma::zeros(n1);
  arma::vec u = arma::ones(n);
  arma::vec XTbeta_or_alpha;
  if (usealpha) {
    XTbeta_or_alpha = Xtest.t() * alpha;
  } else{
    XTbeta_or_alpha = Xtest.t() * beta;
  }
  if (usealpha) {
    for (int i = 0; i < n1; i++) {
      double xalpha = XTbeta_or_alpha(i);
      arma::vec kernalvec = arma::normpdf( (xalpha*u -XTalpha) / h) / h;
      resultvector(i) = dot(kernalvec, gfunction) / sum(kernalvec);
    }
  }else{
    for (int i = 0; i < n1; i++) {
      double xbeta = XTbeta_or_alpha(i);
      arma::vec kernalvec = arma::normpdf( (xbeta*u -XTbeta) / h) / h;
      resultvector(i) = dot(kernalvec, gfunction) / sum(kernalvec);
    }
  }
  return arma::sign(resultvector);
}


// [[Rcpp::export]]
arma::vec weight_defunc1_c(arma::mat X, arma::vec alpha, arma::vec Y, double n) {
  arma::vec result = exp(X.t() * alpha);
  result = X * (arma::pow(Y,2) % (arma::pow(result, 4) - 1) % result / (arma::pow(result, 3)));
  return result;
}
// [[Rcpp::export]]
arma::vec weight_defunc2_c(arma::mat X, arma::vec alpha, arma::vec Y, double n) {
  arma::vec tmp = exp(X.t() * alpha);
  int p = alpha.n_elem;
  int n2 = Y.n_elem;
  arma::vec result = arma::zeros(p);
  for(int j = 0; j < p; j++) {
    for (int i = 0; i < n2; i++) {
      result(j) += tmp(i) * X(j,i) * ((pow(tmp(i), 4)) * pow(Y(i) -1,2) - pow(Y(i), 2)) / pow(tmp(i), 3);
    }
  }
  return result;
}


// [[Rcpp::export]]
arma::mat find_hessian(arma::mat X, arma::vec alpha, arma::vec Y, arma::uvec index1, arma::uvec index2, int p){
  arma::mat hessian = arma::zeros(p,p);
  int n1 = index1.n_elem;
  int n2 = index2.n_elem;
  arma::mat X1 = X.cols(index1);
  arma::mat X2 = X.cols(index2);
  arma::vec tmp = exp(X.t() * alpha);
  for (int  j = 0; j < p; j++){
    for (int k = 0; k < p; k++) {
      for (int l = 0; l < n2; l++) {
        int i = index2(l);
        hessian(j,k) += X2(j,l) * X2(k,l) * 2*pow(tmp(i), 2) * (pow((Y(i) - 1),2) + 2 * pow(Y(i),2)/(pow(tmp(i), 2))); 
      }
      for (int l = 0; l < n1; l++) {
        int i = index1(l);
        hessian(j,k) += X1(j,l) * X1(k,l) * pow(Y(i),2) * (2*pow(tmp(i), 2) + 2*  pow(tmp(i), -2)); 
      }
    }
  }
  return hessian;
}


double weight_function_c(arma::mat X, arma::vec alpha, arma::vec Y, double n, arma::uvec index) {
  arma::vec result = exp(X.t() * alpha);
  arma::vec newvec = arma::zeros(n, 1);
  newvec.rows(index) = result.rows(index);
  result = Y % (1 + arma::pow(result, 2)) / result - newvec;
  return 0.5 * pow(arma::norm(result), 2);
}


arma::vec newton_c(double stoptol, arma::mat X, arma::vec alpha, arma::vec Y, double n, arma::uvec index1, arma::uvec index2) {
  arma::mat X1 = X.cols(index1);
  arma::mat X2 = X.cols(index2);
  int p = alpha.n_elem;
  arma::mat hessian = arma::zeros(p,p);
  arma::vec gradient1 = weight_defunc1_c(X1, alpha, Y.rows(index1), n);
  arma::vec gradient2 = weight_defunc2_c(X2, alpha, Y.rows(index2), n);
  arma::vec gradient = gradient1 + gradient2;
  double error = arma::norm(gradient) / (1 + arma::norm(gradient1) + arma::norm(gradient2));
  int iter = 0;


  while(error > stoptol && iter < 1000){
    hessian = find_hessian(X, alpha, Y, index1, index2, p);
    //vec shessian = arma::svd(hessian);
    //cout << " max = " << max(shessian) << " min = " << min(shessian) << endl;
    arma::vec direction = arma::solve(hessian, gradient, solve_opts::allow_ugly);

    double rho = 1;
    double obj = weight_function_c(X, alpha, Y, n, index2);
    for(int i = 0; i < 50; i++) {
      if (weight_function_c(X, alpha - rho * direction, Y, n, index2) <= obj - 1e-4 * rho * dot(gradient, direction)){
        alpha = alpha - rho * direction;
        //cout << " out in iter = " << i << endl;
        
        break;
      }
      rho = rho * 0.5;
    }
    
    iter = iter + 1;
    

    
    gradient1 = weight_defunc1_c(X1, alpha, Y.rows(index1), n);
    gradient2 = weight_defunc2_c(X2, alpha, Y.rows(index2), n);
    gradient = gradient1 + gradient2;
    error = arma::norm(gradient) / (1 + arma::norm(gradient1) + arma::norm(gradient2));
    int printyes = 0;
    if ((modM(iter, 10) == 0 || iter < 5 || error < stoptol) && printyes == 1) {
      cout << "      iter = " << iter << " obj = " << obj << "  error = " << error << endl;
    }

  }
  
  return alpha;
}


// [[Rcpp::export]]
Rcpp::List mytestnewest(arma::mat X, arma::vec Y, arma::vec beta, double h, double stoptol, int maxiter, bool printyes, arma::vec alpha, bool usealpha, arma::mat Xtest, arma::vec Ytest, bool usetest){
  //int p = X.n_rows;
  int n = X.n_cols;
  arma::vec Ynew = Y;
  double objtest = 0;
  

  arma::uvec Index_great = find(Ynew != 0);
  arma::uvec Index_zero = find(Ynew == 0);
  arma::mat Xnew = X.cols(Index_great);
  arma::mat XnewXnewT = Xnew * Xnew.t();
  arma::mat XnewYnew = Xnew * Ynew(Index_great);
  arma::mat XTalpha = X.t() * alpha;
  
  
  
  // no beta initial 
 // beta = arma::ones(p);
  arma::vec betaold = beta;
  // find beta
  beta = newton_c(stoptol, X, beta, Ynew, n, Index_zero, Index_great);

  /*
  for (int i = 0; i < p; i++){
    cout << beta(i) << endl;
  }
  */
  arma::vec Yhat = arma::pow(exp(X.t() * beta), 2);
  Yhat = Yhat / (1 + Yhat);
  Yhat.rows(Index_zero) *= 0;
  
  
  double objtrain = pow(norm(Yhat -Y, 2),2) / n;
  arma::vec XTbeta = X.t() * beta;
  arma::vec XTbeta_exp = exp(XTbeta);
  double diffbeta = norm(beta - betaold, 2);
  double diffbeta_relative = diffbeta / (1 + norm(betaold));
  int iter = 1;
  double obj = pow(norm(Ynew(Index_great) - XTbeta_exp(Index_great),2), 2) / n;
  if (printyes){
    cout << " iter     betadiff      betadiff_relative     obj       objtrain         Nindex" << endl;
    cout << " " << dec << iter << "     " << scientific << diffbeta << "          " << diffbeta_relative << "    " << obj << "    "  << objtrain << "   " << Index_great.n_elem << endl;
  }
  arma::vec objtestall = arma::zeros(5);

  arma::vec gfunction = - (XTbeta_exp - 2*(1 + pow(XTbeta_exp, 2)) % Y / XTbeta_exp) % XTbeta_exp;
  arma::vec gfunction_old;
  arma::vec signvalue = findsign(X,gfunction,XTbeta,h,XTalpha,usealpha);
  Index_great = find(signvalue >= 0);
  Index_zero = find(signvalue < 0);
  double objteststep1 = 0;
  arma::vec Yhatteststep1 = arma::pow(exp(Xtest.t() * beta),2);
  if (usetest) {
    arma::vec signteststep1 = findsign_test(Xtest, gfunction,XTbeta, h, XTalpha, usealpha, alpha, beta);
    Yhatteststep1 = Yhatteststep1/(1 + Yhatteststep1);
    arma::uvec index_zero_test_step1 = find(signteststep1 < 0);
    Yhatteststep1.rows(index_zero_test_step1) *= 0;
    int n1 = Ytest.n_elem;
    objteststep1 = pow(norm(Yhatteststep1 - Ytest, 2), 2) / n1;
    objtestall(0) =  objteststep1;
  }
  arma::vec betastep1 = beta;
  //double objtrain_old = 1000000000000;
  
  
  while (iter < 5) {
    iter++;
      
    
    betaold = beta;
    // find beta
    beta = newton_c(stoptol, X, beta, Ynew, n, Index_zero, Index_great);
    XTbeta = X.t() * beta;
    XTbeta_exp = exp(XTbeta);
    diffbeta = norm(beta - betaold, 2);
    diffbeta_relative = diffbeta / (1 + norm(betaold));
    obj = pow(norm(Ynew(Index_great) - XTbeta_exp(Index_great),2), 2) / n;
    gfunction_old = gfunction;
    
    gfunction = - (XTbeta_exp - 2*(1 + pow(XTbeta_exp, 2)) % Y / XTbeta_exp) % XTbeta_exp;
    signvalue = findsign(X,gfunction,XTbeta,h,XTalpha,usealpha);
    Index_great = find(signvalue >= 0);
    Index_zero = find(signvalue < 0);
    
    
    Yhat = arma::pow(exp(X.t() * beta), 2);
    Yhat = Yhat / (1 + Yhat);
    Yhat.rows(Index_zero) *= 0;
    objtrain = pow(norm(Yhat -Y, 2),2) / n;
    
    
    
    
    if(printyes && (modM(iter,10)==0 || iter < 5 || iter == maxiter || diffbeta < stoptol)){
      cout << " " << dec << iter << "     " << scientific << diffbeta << "          " << diffbeta_relative << "    " << obj << "    "  << objtrain << "   " << Index_great.n_elem << endl;
    }
    XTbeta = X.t() * beta;
    
    arma::vec Yhattest = arma::pow(exp(Xtest.t() * beta),2);
    if (usetest) {
      arma::vec signtest = findsign_test(Xtest, gfunction,XTbeta, h, XTalpha, usealpha, alpha, beta);
      Yhattest = Yhattest/(1 + Yhattest);
      arma::uvec index_zero_test = find(signtest < 0);
      Yhattest.rows(index_zero_test) *= 0;
      int n1 = Ytest.n_elem;
      objtest = pow(norm(Yhattest - Ytest, 2), 2) / n1;
      objtestall(iter-1) = objtest;
    }
  }
  

  //cout << " objtest = " << objtest<< " objteststep1 = " << objteststep1 << endl;
  
  if (usetest) {
    return Rcpp::List::create(Named("beta") = beta,
                              Named("objtrain") = objtrain,
                              Named("objtest") = objtest,
                              Named("betastep1") = betastep1,
                              Named("objteststep1") = objteststep1,
                              Named("objtestall") = objtestall
    );
  } else{
    
    return Rcpp::List::create(Named("beta") = beta,
                              Named("betastep1") = betastep1,
                              Named("objtrain") = objtrain,
                              Named("objtestall") = objtestall);
  }
  
  
  
  
}


