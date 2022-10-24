
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame harpSpatial_basic_scores(NumericMatrix obfield, NumericMatrix fcfield) {
  int i,j, ni=obfield.nrow(), nj=obfield.ncol() ;
  double bias=0., mse=0., mae=0., tmp ;
  for (j=0 ; j < nj ; j++) {
    for (i=0 ; i < ni ; i++) {
      bias += (tmp = fcfield(i,j) - obfield(i,j)) ;
      mse  += tmp * tmp ;
      mae  += abs(tmp) ;
    }
  }
  bias /= ni*nj ;
  mse  /= ni*nj ;
  mae  /= ni*nj ;
  return Rcpp::DataFrame::create(Named("mse")  = mse,
                                 Named("mae")  = mae,
                                 Named("bias") = bias);

}
