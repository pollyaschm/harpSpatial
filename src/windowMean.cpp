// for fast FSS & other "neighborhood" ("fuzzy") score calculation
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
  // obsolete: use neighborhood_scores
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
DataFrame harpSpatial_neighborhood_scores(NumericMatrix obfield, NumericMatrix fcfield,
                    NumericVector thresholds, NumericVector scales) {
  int i, j, k, th, sc;
  double a, b, c, dd ;
  double fss1, fss2 ;
  int n_thresholds=thresholds.length(), n_scales=scales.length();
  int ni=fcfield.nrow(), nj=fcfield.ncol();

  NumericMatrix frac_fc(ni,nj), frac_ob(ni,nj);
  NumericMatrix cum_fc(ni,nj), cum_ob(ni,nj);

  // numeric vectors for the result
  NumericVector res_thresh(n_thresholds * n_scales);
  NumericVector res_size(n_thresholds * n_scales);
  NumericVector res_fss(n_thresholds * n_scales);
  NumericVector res_a(n_thresholds * n_scales);
  NumericVector res_b(n_thresholds * n_scales);
  NumericVector res_c(n_thresholds * n_scales);
  NumericVector res_d(n_thresholds * n_scales);
  // NOTE: we calculate several scores together
  //       it would be redundant to calculate the cumulated matrices twice...
  // TODO:
  // if (ob.nrow() != ni || ob.ncol != nj) ERROR
  //
  for (th=0 ; th < n_thresholds ; th++) {
    // calculate cumsum matrices for given threshold
    cum_fc = cumsum2d_bin(fcfield, thresholds[th]);
    cum_ob = cumsum2d_bin(obfield, thresholds[th]);
    for (sc=0 ; sc < n_scales ; sc++) {
      k = th*n_scales + sc;
      res_thresh(k) = thresholds(th);
      res_size(k) = scales(sc);
      // fraction matrices
      frac_fc = window_mean_from_cumsum(cum_fc, (int) scales[sc]);
      frac_ob = window_mean_from_cumsum(cum_ob, (int) scales[sc]);
      fss1 = fss2 = a = b = c = 0.;
      for (j=0 ; j < nj ; j++) {
        for (i=0 ; i < ni ; i++) {
          // FSS
          fss1 += (frac_fc(i,j)-frac_ob(i,j))*(frac_fc(i,j)-frac_ob(i,j)) ;
          fss2 += frac_fc(i,j)*frac_fc(i,j) + frac_ob(i,j)*frac_ob(i,j) ;

          // Neighborhood Adapted Contingency Table
          // ref Stein & Stoop 2019
          // method:
          // f_fc = a + b of unadapted contingency tab
          // f_ob = a + c of unadapted contingency tab
          // if the difference (b-c) is neg., min(b,c)=b
          if ((dd = frac_fc(i,j)-frac_ob(i,j)) < 0) {
            a += frac_fc(i,j) ;
            c -= dd ;
          } else {
            a += frac_ob(i,j) ;
            b += dd ;
          }
          // TODO: Other scores
          
        } //i
      } //j
      res_fss[k] = (fss2 < 1.0E-3) ? 0. : 1. - fss1/fss2 ;
      res_a[k]   = a / (ni*nj) ;
      res_b[k]   = b / (ni*nj) ;
      res_c[k]   = c / (ni*nj) ;
      res_d[k]   = 1. - res_a[k] - res_b[k] - res_c[k] ;
    } //sc
  } //th

  return Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                 Named("scale")     = res_size,
                                 Named("fss")       = res_fss,
                                 Named("a")         = res_a,
                                 Named("b")         = res_b,
                                 Named("c")         = res_c,
                                 Named("d")         = res_d);
}
  

