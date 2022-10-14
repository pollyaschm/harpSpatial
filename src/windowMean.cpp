// for fast FSS & other "fuzzy" score calculation
// reference: N. Faggian et al. "Fast calculation of the franctions skill score"
//            MAUSAM, 66 (2015) 457-466
// 1. get the 'summed_area_table' for a given threshold (cumsum2d)
// 2. use this for fast calculation for Fractions tables

#include <Rcpp.h>
using namespace Rcpp;

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
  int i, j, N, ni=indat.nrow(), nj=indat.ncol() ;
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

// a fast windowing call that can be used with spatialVx
// No longer used
// [[Rcpp::export]]
NumericMatrix windowMean(NumericMatrix indat, NumericVector radius) {
  int rad=(int) radius[0] ;
  if (rad==0) return indat;
  return window_mean_from_cumsum(cumsum2d(indat), rad);
}

// [[Rcpp::export]]
double fss_from_fractions(NumericMatrix m1, NumericMatrix m2) {
  int ni=m1.ncol(), nj=m1.nrow();
  int i, j;
  double fss1=0., fss2=0.;

  for (i=0 ; i < ni ; i++) {  
    for (j=0 ; j < nj ; j++) {  
      fss1 += (m1(i,j)-m2(i,j))*(m1(i,j)-m2(i,j)) ;
      fss2 += m1(i,j)*m1(i,j) + m2(i,j)*m2(i,j) ;
    }
  }
  if (fss2 < 1.0E-3) return 0. ;
  return (1. - fss1/fss2) ;
}

// [[Rcpp::export]]
DataFrame score_fss(NumericMatrix obfield, NumericMatrix fcfield,
                    NumericVector thresholds, NumericVector window_sizes) {
  int i, j, k;
  int n_thresholds=thresholds.length(), n_sizes=window_sizes.length();
  int ni=fcfield.ncol(), nj=fcfield.nrow();
  NumericVector res_fss(n_thresholds * n_sizes);
  NumericVector res_thresh(n_thresholds * n_sizes);
  NumericVector res_size(n_thresholds * n_sizes);
  NumericMatrix frac_fc(ni,nj), frac_ob(ni,nj);
  NumericMatrix cum_fc(ni,nj), cum_ob(ni,nj);

  // TODO:
  // if (ob.nrow() != ni || ob.ncol != nj) ERROR
  //
  for (i=0 ; i < n_thresholds ; i++) {
    // calculate cumsum matrices
    cum_fc = cumsum2d_bin(fcfield, thresholds[i]);
    cum_ob = cumsum2d_bin(obfield, thresholds[j]);
    for (j=0 ; j < n_sizes ; j++) {
      k = i*n_sizes + j;
      res_thresh(k) = thresholds(i);
      res_size(k) = window_sizes(j);
      // fraction matrices
      frac_fc = window_mean_from_cumsum(cum_fc, (int) window_sizes[j]);
      frac_ob = window_mean_from_cumsum(cum_ob, (int) window_sizes[j]);
 
      res_fss(k) = fss_from_fractions(frac_fc, frac_ob) ;
        //mean( (frac_fc(_,_)-frac_ob(_,_))^2) /
        //mean(frac_fc^2 + frac_ob^2);
      // other "fuzzy" scores: ETS, ...
    }
  }

  return Rcpp::DataFrame::create(Named("threshold")=res_thresh,
                                 Named("scale")=res_size,
                                 Named("value")=res_fss);
}
  

