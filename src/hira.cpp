#include <Rcpp.h>
#include "spatial_shared.h"
using namespace Rcpp;





float TS(float a, float b, float c) {
  float denum = a + b + c;
  
  float ts = denum > 0 ? a / denum : 0;
  
}

// [[Rcpp::export]]
Rcpp::List hira_hira(NumericVector obsvect, NumericMatrix indices, NumericMatrix fcfield,NumericVector thresholds,
NumericVector scales, NumericVector strategies) {
	   
	   Rcout << "dims: " << strategies.length();
	   return NULL;
}

// [[Rcpp::export]]
Rcpp::List get_hira_scores(NumericVector obsvect, NumericMatrix indices, NumericMatrix fcfield,NumericVector thresholds,
NumericVector scales, NumericVector strategies) {
  // indices is Nx2 matrix the location of the grid point in the model passed grid.
  // indices(0,:) is the first index and indices(1,:) is the second index of fcfield which represents the model grid.


  int n_thresholds = thresholds.length();
  int n_scales = scales.length();
  int ni = fcfield.nrow(), nj = fcfield.ncol();
  int no = obsvect.length();
  
  int nstrat = strategies.length(); // Number of Strategies
  
  NumericVector sum_fc(no);
  NumericVector sum_ob(no);
  LogicalVector bin_ob(no);
  NumericMatrix cum_fc(ni, nj);
  NumericMatrix cum_ob(ni, nj);
  NumericMatrix obsongrid(ni, nj); 
    
  //Rcout << "dims: " << cum_fc.nrow() << " " << cum_fc.ncol() << "\n";


  
  // numeric vectors for the result
  NumericVector res_thresh(n_thresholds * n_scales);
  NumericVector res_size(n_thresholds * n_scales);
  
  // Multi Event
  NumericVector res_me_a(n_thresholds * n_scales);
  NumericVector res_me_b(n_thresholds * n_scales);
  NumericVector res_me_c(n_thresholds * n_scales);
  NumericVector res_me_d(n_thresholds * n_scales);
  // Pragramtic
  NumericVector res_pra_bss(n_thresholds * n_scales);
  NumericVector res_pra_bs(n_thresholds * n_scales);
  //Conditional Square Root RPS
  
  // Theat Detection by Daniel 
  NumericVector res_td_a(n_thresholds * n_scales);
  NumericVector res_td_b(n_thresholds * n_scales);
  NumericVector res_td_c(n_thresholds * n_scales);
  NumericVector res_td_d(n_thresholds * n_scales);
  
  
  // NumericMatrix res_cdf(n_thresholds , n_scales);
  NumericVector res_csrr_pre_prs(n_thresholds * n_scales);
  NumericVector res_csrr_pre_px(n_thresholds * n_scales);
  
  bool is_multi_event = false; // 0
  bool is_pragmatic = false; // 1 
  bool is_pph = false; // 2 Practically Perfext Hindcast
  bool is_csrr = false; // 3 Conditional square root for RPS
  
  for (int is = 0; is < nstrat; is++) {
    int stra = (int) strategies(is);
    
    switch (stra) {
    case 0:
      is_multi_event = true;
      break;
    case 1:
      // Pragmatic Stratigy
      is_pragmatic = true;
      break;
    case 2:
      is_pph = true;
      break;
    case 3:
      is_csrr = true;
      break;
      //default:
      // code block
    }
  }
  
  
  if (is_pph) {
    for (int j=0 ; j < nj ; j++) {
      for (int i=0 ; i < ni ; i++) {
        obsongrid(i,j) = 0;
      }
      }
    for (int j = 0; j < no; j++) { 
      obsongrid(indices(j,0),indices(j,1)) = obsvect[j];
    }
    
  }
  
  int k = 0;
  for (int th = 0; th < n_thresholds; th++) {
    
    bin_ob = vector_to_bin(obsvect, thresholds[th]);
    cum_fc = cumsum2d_bin(fcfield, thresholds[th]);
    
	if (is_pph) {
      cum_ob = cumsum2d_bin(obsongrid, thresholds[th]);
    }
     
    for (int sc = 0; sc < n_scales; sc++) {
      ;
      res_thresh(k) = thresholds(th);
      res_size(k) = scales(sc);
      
      int rad = (int) scales[sc];
      

  
      sum_fc = window_sum_from_cumsum_for_ij(cum_fc, rad, indices);
	  
	  if (is_pph) {
	     sum_ob = window_sum_from_cumsum_for_ij(cum_ob, rad, indices);
      }
	  
       
      if (is_multi_event) {
        res_me_a[k] = 0;
        res_me_b[k] = 0;
        res_me_c[k] = 0;
        for (int j = 0; j < no; j++) {
          //TODO: It could be esier to use only the bitwise operations. the empysise of logical values to is only just to make sure
          bool is_fc = sum_fc[j] > 0;
          res_me_a[k] += (bin_ob[j]) && (is_fc);
          res_me_b[k] += (~bin_ob[j]) && (is_fc);
          res_me_c[k] += (bin_ob[j]) && (~is_fc);
         
        }
         res_me_d[k] = no - res_me_a[k] - res_me_b[k] - res_me_c[k];
      }
     
      if (is_pragmatic) {
        float nume = 0;
        float px_ave = 0;
        for (int j = 0; j < no; j++) {
          float diff = sum_fc(j) / ((2 * rad + 1) * (2 * rad + 1)) - (bin_ob[j]);
          nume += diff * diff;
          px_ave += bin_ob[j];
        }
        px_ave = px_ave / no;
        
        float denume = 0;
        for (int j = 0; j < no; j++) {
          float diff = px_ave - (bin_ob[j]);
          denume += diff * diff;
        }
        
        res_pra_bss[k] = denume > 0.001 ? 1 - nume / denume : -9999;
        res_pra_bs[k] = nume / no;
        
      }
      
      if (is_pph) {
        // Practically Perfect hindcast
        // Here I implement a slightly modified method.
        // Look for all grid points that fall in the same neighborhood of some scale.
        
    
        res_td_a[k] = 0;  
        res_td_b[k] = 0;
        res_td_c[k] = 0;
        res_td_d[k] = 0;
        for(int iob=0; iob<no;iob++){
          int i = indices(iob,0);
          int j = indices(iob,1);
          int co = sum_ob(i,j);
          int cf = sum_fc(i,j);
          res_td_a[k] += co > 0 && cf >= co;  
          res_td_b[k] += co == 0 && cf !=0;
          res_td_c[k] += co > 0 && cf < co;
        }

        res_td_d[k] = no - res_td_a[k] - res_td_b[k] - res_td_c[k];
        
      }
      
      if (is_csrr) {
        // Conditional Square root of RPS
        float rps = 0;
        float csrr_Ix = 0;
        float csrr_Py = 0;
        res_csrr_pre_px[k] = 0;
        for (int j = 0; j < no; j++) {
          csrr_Ix = (bin_ob[j]);
          csrr_Py = sum_fc[j] / ((2 * rad + 1) * (2 * rad + 1));
          float diff = csrr_Ix - csrr_Py;
          rps += diff * diff;
          res_csrr_pre_px[k] += csrr_Ix;
        }
        
        // this quantity should also be normalized over the fraction of observations the excceded the lowest threshold.
        res_csrr_pre_prs[k] = rps / no; // This quantity should be summed over all thresholds  and divided  over ( number of thresholds -1) to be presented later
        res_csrr_pre_px[k] = res_csrr_pre_px[k] / no;
      }
      k++; 
    } // sc
  } // th
  
   Rcpp::List resultList;
   

  if (is_multi_event) {
	Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                                 Named("scale") = res_size);
    df["hit"] = res_me_a;
    df["fa"] = res_me_b;
    df["miss"] = res_me_c;
    df["cr"] = res_me_d;
	
	resultList["me"] = df;
  }
  if (is_pragmatic) {
	  
	Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                                 Named("scale") = res_size);
    df["bss"] = res_pra_bss;
    df["bs"] = res_pra_bs;
	resultList["pragm"] = df;
  }
  
  if (is_csrr) {
	  
	Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                                 Named("scale") = res_size);
												 
    df["prs"] = res_csrr_pre_prs;
    df["px"] = res_csrr_pre_px;
	
	resultList["csrr"] = df;
	
  }
  
  if (is_pph) {
	Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("threshold") = res_thresh,
                                                 Named("scale") = res_size);
												 
    df["hit"] = res_td_a;
    df["fa"] = res_td_b;
    df["miss"] = res_td_c;
    df["cr"] = res_td_d;
	
	resultList["pph"] = df;
	
  }
  return resultList; 
}
  