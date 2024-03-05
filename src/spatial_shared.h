#ifndef SPATIAL_HEADER_H
#define SPATIAL_HEADER_H

#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix cumsum2d(NumericMatrix indat);
NumericMatrix cumsum2d_bin(NumericMatrix indat, float threshold);
NumericMatrix window_mean_from_cumsum(NumericMatrix indat, int rad);
IntegerVector vector_to_bin(NumericVector indat, float threshold);
NumericVector window_sum_from_cumsum_for_ij(NumericMatrix indat, int rad, NumericMatrix indices);

#endif // SPATIAL_HEADER_H
