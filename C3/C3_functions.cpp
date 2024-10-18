//#include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector matrix2vector(arma::mat m, const bool byrow=false){
  if (byrow) {
    return as<NumericVector>(wrap(m.as_row()));
  } else {
    return as<NumericVector>(wrap(m.as_col()));
  }
}


// [[Rcpp::export]]
IntegerVector rank_rcpp(NumericVector x) {
  NumericVector sorted = clone(x).sort();
  return match(x, sorted);
}

// [[Rcpp::export]]
NumericVector emp_chieta_rcpp(NumericVector x, NumericVector y, double u){
  if (x.size() != y.size()){
    stop("x and y must have the same length");
  }
  int N = x.size();
  LogicalVector peaks = as<NumericVector>(rank_rcpp(x))/(double)N > u & as<NumericVector>(rank_rcpp(y))/(double)N > u;
  double chi = sum(peaks)/(double)N/(1.0-u);
  double eta = log(1.0-u)/log(sum(peaks)/(double)N);
  NumericVector extdep = NumericVector::create(Named("chi")=chi, _["eta"]=eta, _["N"]=N, _["e"]=sum(peaks));
  return(extdep);
}

// [[Rcpp::export]]
List chieta_matrix(NumericMatrix x, double u){
  int n = x.ncol();
  NumericMatrix chi(n,n);
  NumericMatrix eta(n,n);
  std::fill( chi.begin(), chi.end(), NumericVector::get_na() ) ;
  std::fill( eta.begin(), eta.end(), NumericVector::get_na() ) ;
  for (int i=0; i<n-1; ++i){ 
    for (int j=i+1; j<n; ++j){
      NumericVector chieta = emp_chieta_rcpp(x(_,i), x(_,j), u);
      chi(i,j) = chieta["chi"];
      eta(i,j) = chieta["eta"];
    }
  }
  List extdep = List::create(Named("chi") = chi , _["eta"] = eta);
  return(extdep);
}

