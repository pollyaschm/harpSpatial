
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame harpSpatial_basic_scores(NumericMatrix obfield, NumericMatrix fcfield) {
// NOTE: "abs" would give integer abs, so make sure to write std::abs
  int i,j, ni=obfield.nrow(), nj=obfield.ncol() ;
  double bias=0., mse=0., mae=0., rmse=0., tmp ;
  double obmean, fcmean, corrnum=0., corrdenom1=0., corrdenom2=0., Rpearson=0. ;
  obmean = mean(obfield) ;
  fcmean = mean(fcfield) ;
  for (j=0 ; j < nj ; j++) {
    for (i=0 ; i < ni ; i++) {
      tmp = fcfield(i,j) - obfield(i,j);
      bias += tmp ;
      mse  += tmp * tmp ;
      mae  += std::abs(tmp) ;
      corrnum    += (obfield(i,j) - obmean) * (fcfield(i,j) - fcmean) ;
      corrdenom1 += (obfield(i,j) - obmean) * (obfield(i,j) - obmean) ;
      corrdenom2 += (fcfield(i,j) - fcmean) * (fcfield(i,j) - fcmean) ;
    }
  }
//  Rcout << "bias: " << bias << " mse: " << mse << " mae: " << mae << std::endl;
//  Rcout << "ni: " << ni << " nj: " << nj << std::endl;
  bias /= ni*nj ;
  mse  /= ni*nj ;
  mae  /= ni*nj ;

  rmse  = sqrt(mse) ;

  corrnum /= ni*nj ;
  corrdenom1 /= ni*nj ;
  corrdenom2 /= ni*nj ;
  Rpearson = corrnum / (sqrt(corrdenom1) * sqrt(corrdenom2)) ;
  
//  Rcout << "bias: " << bias << " mse: " << mse << " mae: " << mae << std::endl;
//  Rcout << "rmse: " << rmse << std::endl;
//  Rcout << "corrnum: " << corrnum << std::endl;
//  Rcout << "corrdenom1: " << corrdenom1 << "corrdenom2: " << corrdenom2 << std::endl;
//  Rcout << "Rpearson: " << Rpearson << std::endl;

  return Rcpp::DataFrame::create(Named("mse")  = mse,
                                 Named("mae")  = mae,
                                 Named("bias") = bias,
                                 Named("rmse") = rmse,
				 Named("Rpearson") = Rpearson);

}
