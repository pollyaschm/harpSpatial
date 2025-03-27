#include <Rcpp.h>
#include <harpCore.h>
#include "spatial_shared.h"
using namespace Rcpp;
using namespace harpCore;


// [[Rcpp::export]]
NumericMatrix cumsum2d(NumericMatrix indat) {
  // result(I,J) = sum_{i<=I, j<=J} indat(i,j)
  // Not used in FSS, but for "Agreement Scale"
  // For FSS, we do in-line thresholding (cumsum2d_bin)
  int i,j, ni=indat.nrow(), nj=indat.ncol() ;
  NumericMatrix result(ni,nj);

  for (j=0; j< nj ; j++) {
    result(0, j) = indat(0, j);
    for (i=1; i< ni; i++) {
      result(i, j) = indat(i, j) + result(i-1, j);
    }
  }
  for (i=0; i< ni; i++) {
    for (j=1; j< nj ; j++) {
      result(i, j) += result(i, j-1);
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix cumsum2d_bin(NumericMatrix indat, float threshold) {
  // result(I,J) = sum_{i<=I, j<=J} indat(i,j)
  // BUT: here we modify indat to 0/1 (>= threshold)
  int i,j, ni=indat.nrow(), nj=indat.ncol() ;
  NumericMatrix result(ni,nj);

  for (j=0; j< nj ; j++) {
    result(0, j) = (indat(0, j) >= threshold) ;
    for (i=1; i< ni; i++) {
      result(i, j) = (indat(i, j) >= threshold) + result(i-1, j);
    }
  }
  for (i=0; i< ni; i++) {
    for (j=1; j< nj ; j++) {
      result(i, j) += result(i, j-1);
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix window_mean_from_cumsum(NumericMatrix indat, int rad) {
  // windowed average
  // input matrix is output from cumsum2d[_bin]
  // rad is an integer (>=0), window size is 2*rad+1
  // zero-padding
  // TODO : other boundary options:
  //    - periodic or mirror
  //    - reduce rad close to border
  int i, j, ni=indat.nrow(), nj=indat.ncol() ;
  int imax, jmax;
  NumericMatrix result(ni, nj);

  for (i=0 ; i < ni; i++) {
    imax = std::min(i+rad, ni-1) ;
    for(j=0 ; j < nj; j++) {
      jmax = std::min(j+rad, nj-1) ;
      result(i, j) = indat(imax, jmax) ;
      if (i > rad) {
        result(i,j) -= indat(i-rad-1, jmax);
        if (j > rad)
          result(i,j) += indat(i-rad-1, j-rad-1) - indat(imax, j-rad-1);
      }
      else if (j > rad) result(i,j) -= indat(imax, j-rad-1);
      result(i, j) /= (2*rad+1)*(2*rad+1);
    }
  }
  return result;
}


// [[Rcpp::export]]
Rcpp::IntegerVector vector_to_bin(NumericVector indat, float threshold) {

  int ni = indat.length();
  IntegerVector  result(ni);

  for (int i = 0; i < ni; i++) {
    result(i) = indat(i) >= threshold?1:0;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector window_sum_from_cumsum_for_ij(NumericMatrix indat, int rad, NumericMatrix indices) {

  //Rcout << "dims2: " << indat.nrow() << " " << indat.ncol() << "\n";
  // windowed average
  // input matrix is output from cumsum2d[_bin]
  // rad is an integer (>=0), window size is 2*rad+1
  // zero-padding
  // TODO : other boundary options:
  //    - periodic or mirror
  //    - reduce rad close to border

  int i, j, ni = indat.nrow(), nj = indat.ncol();
  int no = indices.ncol();

  int imax, jmax;
  NumericVector result(no);
  for (int k = 0; k < no; k++) {
    i = (int) indices(0, k);
    j = (int) indices(1, k);

    imax = std::min(i + rad, ni - 1);
    jmax = std::min(j + rad, nj - 1);
    result(k) = indat(imax, jmax);
    if (i > rad) {
      result(k) -= indat(i - rad - 1, jmax);
      if (j > rad)
        result(k) += indat(i - rad - 1, j - rad - 1) - indat(imax, j - rad - 1);
    } else if (j > rad)
      result(k) -= indat(imax, j - rad - 1);
  }

  return result;
}
