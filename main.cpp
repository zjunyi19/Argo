//
//  main.cpp
//  STATS
//
//  Created by Junyi Zhang on 2/13/19.
//  Copyright © 2019 Junyi Zhang. All rights reserved.
//

#include <RcppArmadillo.h>
#include <algorithm>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix calc_Omega(NumericVector knots, int M) {

  int K = knots.size();
  vec arma_knots = as<arma::vec>(knots);

  // initialize allknots
  vec allknots(K+2*M-2);
  allknots.zeros();
  int knots_size = allknots.size();
  for (int i = M - 1; i < M-1+K; ++i) {
    allknots[i] = arma_knots[i-(M-1)];
  }
  for (int i = K+M-1; i < knots_size; ++i) {
    allknots[i] = 1;
  }

  // initialize xtilde
  vec xtilde(3*(K+2*M-3));
  for (int i = 0; i < xtilde.size(); i += 3) {
    xtilde[i] = allknots[i/3];
    xtilde[i+1] = (allknots[i/3] + allknots[i/3+1]) / 2;
    xtilde[i+2] = allknots[i/3+1];
  }

  // second derivative of B matrix
  mat B2(3*(K+2*M-3), K+M-2);

  // loop through all the X
  for (int k = 0; k < xtilde.size(); ++k) {

    double x = xtilde[k];
    vec dp1(knots_size-1);
    dp1.zeros();

    // calculating the first derivative
    // case1: m = 1
    for (int i = 0; i < knots_size-1; ++i) {
      if (allknots[i] <= x && x <= allknots[i + 1]) {
        dp1[i] = 1;
      }
    }

    // case2: m = 2,3,....,M-2
    int m = 2;
    while (m <= M-2) {
      for (int i = 0; i < knots_size-m; ++i) {
        double first_para = 0, second_para = 0;
        if (allknots[i+m-1] != allknots[i]) {
          first_para = (x - allknots[i]) / (allknots[i+m-1] - allknots[i]);
        }
        if (allknots[i+m] != allknots[i+1]) {
          second_para = (allknots[i+m] - x) / (allknots[i+m] - allknots[i+1]);
        }
        dp1[i] = first_para * dp1[i] + second_para * dp1[i+1];
      }
      ++m;
    }

    // case3, m = M-1
    for (int i = 0; i < knots_size-m; ++i) {
      double first_para = 0, second_para = 0;
      if (allknots[i+m-1] != allknots[i]) {
        first_para = 1 / (allknots[i+m-1] - allknots[i]);
      }
      if (allknots[i+m] != allknots[i+1]) {
        second_para = -1 / (allknots[i+m] - allknots[i+1]);
      }
      dp1[i] = (m-1) * (first_para * dp1[i] + second_para * dp1[i+1]);
    }

    // calculating the second derivative
    // case4: m = M
    for (int i = 0; i < knots_size-M; ++i) {
      double first_para = 0, second_para = 0;
      if (allknots[i+M-1] != allknots[i]) {
        first_para = 1 / (allknots[i+M-1] - allknots[i]);
      }
      if (allknots[i+M] != allknots[i+1]) {
        second_para = 1 / (allknots[i+M] - allknots[i+1]);
      }
      B2(k,i) = (M-1) * (first_para * dp1[i] - second_para * dp1[i+1]);
    }
  }

  vec w_tilde(3*(K+2*M-3));
  for (int i = 0; i < w_tilde.size(); i += 3) {
    double delta_k = allknots[i/3 + 1] - allknots[i/3];
    w_tilde[i] = w_tilde[i+2] = 1/6.0 * delta_k;
    w_tilde[i+1] = 4.0/6.0 * delta_k;
  }
  mat omega = B2.t() * diagmat(w_tilde) * B2;
  return wrap(omega);
}

