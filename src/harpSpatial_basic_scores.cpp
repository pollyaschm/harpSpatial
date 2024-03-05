
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame harpSpatial_basic_scores(NumericMatrix obfield, NumericMatrix fcfield) {
// NOTE: "abs" would give integer abs, so make sure to write std::abs
  int i,j, ni=obfield.nrow(), nj=obfield.ncol() ;
  double bias=0., mse=0., mae=0., tmp ;
  for (j=0 ; j < nj ; j++) {
    for (i=0 ; i < ni ; i++) {
      tmp = fcfield(i,j) - obfield(i,j);
      bias += tmp ;
      mse  += tmp * tmp ;
      mae  += std::abs(tmp) ;
    }
  }
//  Rcout << "bias: " << bias << " mse: " << mse << " mae: " << mae << std::endl;
//  Rcout << "ni: " << ni << " nj: " << nj << std::endl;
  bias /= ni*nj ;
  mse  /= ni*nj ;
  mae  /= ni*nj ;
  
//  Rcout << "bias: " << bias << " mse: " << mse << " mae: " << mae << std::endl;

  return Rcpp::DataFrame::create(Named("mse")  = mse,
                                 Named("mae")  = mae,
                                 Named("bias") = bias);

}
