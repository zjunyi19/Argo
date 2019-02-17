//
//  main.cpp
//  STATS
//
//  Created by Junyi Zhang on 2/13/19.
//  Copyright © 2019 Junyi Zhang. All rights reserved.
//

#include<Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericMatrix calcB(NumericVector vec_x, NumericVector knots, int M) {
  int K = knots.size();
  
  // initialize tau
  NumericVector tau(K+2*M-2);
  int tau_size = tau.size();
  for (int i = M-1; i < M-1+K; ++i) {
    tau[i] = knots[i+1-M];
  }
  for (int i = M-1+K; i < tau_size; ++i) {
    tau[i] = 1;
  }
  
  NumericMatrix B(vec_x.size(), K + M - 2);
  
  // loop through all x
  for (int k = 0; k < vec_x.size(); ++k) {
    
    double x = vec_x[k];
    NumericVector dp(tau_size-1);
    
    // initialize B[i][1](x) <-- the matrix in the book
    // case1: m = 1
    for (int i = 0; i < tau_size-1; ++i) {
      if (tau[i] <= x && x <= tau[i + 1]) {
        dp[i] = 1;
      }
    }
    
    // case2: m = 2,3,....,M-1
    for (int m = 2; m <= M - 1; ++m) {
      for (int i = 0; i < tau_size-m; ++i) {
        double first_para = 0.0, second_para = 0.0;
        if (tau[i+m-1] != tau[i]) {
          first_para = (x - tau[i]) / (tau[i+m-1] - tau[i]);
        }
        if (tau[i+m] != tau[i+1]) {
          second_para = (tau[i+m] - x) / (tau[i+m] - tau[i+1]);
        }
        dp[i] = first_para * dp[i] + second_para * dp[i+1];
      }
    }
    
    // case3, m = M
    for (int i = 0; i < tau_size-M; ++i) {
      double first_para = 0, second_para = 0;
      if (tau[i+M-1] != tau[i]) {
        first_para = (x - tau[i]) / (tau[i+M-1] - tau[i]);
      }
      if (tau[i+M] != tau[i+1]) {
        second_para = (tau[i+M] - x) / (tau[i+M] - tau[i+1]);
      }
      B(k,i) = first_para * dp[i] + second_para * dp[i+1];
    }
    
  }
  return B;
}
