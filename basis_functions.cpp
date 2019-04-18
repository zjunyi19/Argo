// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;
using namespace arma;

NumericVector bsplvb(NumericVector t, int jhigh, int index, double xq, int left, NumericVector biatx) {
  int jmax = 20;
  NumericVector deltar(jhigh); // deltar[0] is dummy variable
  NumericVector deltal(jhigh); // deltal[0] is dummy variable

  int j;
  if (index == 1) {
    j = 1;
    biatx[1] = 1;
    if (j >= jhigh) {return biatx;}
  }

  while (j < jhigh) {
    int jp1 = j+1;
    deltar[j] = t[left+j] - xq;
    deltal[j] = xq - t[left+1-j];
    double saved = 0;
    for (int i = 1; i <= j; ++i) {
      double term = 0;
      if (deltar[i] + deltal[jp1-i] != 0) {
        term = biatx[i]/(deltar[i] + deltal[jp1-i]);
      }
      biatx[i] = saved + deltar[i] * term;
      saved = deltal[jp1-i] * term;
    }
    biatx[jp1] = saved;
    j = jp1;
  }
  return biatx;
}


NumericVector bsplvd(NumericVector t, int k, double xp, int left, int nderiv) {
  // mhigh is usually equal to nderiv
  int mhigh = max(min(nderiv,k), 1);
  int kp1 = k+1;
  NumericVector temp(5);
  return(bsplvb(t, kp1-mhigh, 1, xp, left, temp));
}

// [[Rcpp::export]]
sp_mat basis_mat(NumericVector x, NumericVector knots, int order) {
  NumericVector ileft(x.size() + 1); // ileft[0] is dummy variable
  int K = knots.size();
  NumericVector t(K + 6 + 1); //t[0] is dummy variable
  t[1] = t[2] = t[3] = knots[0];
  for (int i = 4; i <= (K + 3); ++i) {
    t[i] = knots[i-3-1];
  }
  t[K+4] = t[K+5] = t[K+6] = knots[K-1];

  // same as function interv in source.R file
  int i = 1;
  for (int k = 1; k <= x.size(); ++k) {
    if (x[k-1] < t[i]) {
      ileft[k] = i - 1;
    } else {
      if (i == t.size()) {
        ileft[k] = i - 4 - 1;
      } else {
        ++i;
        --k;
      }
    }
  }
  //////////////////////////////////////////////

  vec vnikx(x.size() * 4); // vnikx[0] is dummy variable
  for(int i = 1; i <= x.size(); ++i) {
    NumericVector temp(5); // temp[0] is dummy variable
    temp = bsplvd(t, order, x[i-1], ileft[i], 1);
    vnikx[(i-1)*4] = temp[1];
    vnikx[(i-1)*4+1] = temp[2];
    vnikx[(i-1)*4+2] = temp[3];
    vnikx[(i-1)*4+3] = temp[4];
  }
  umat locations(2, vnikx.size());
  for (int i = 0; i < x.size(); ++i) {
    locations(0, i*4) = locations(0, i*4+1) = locations(0, i*4+2) = locations(0, i*4+3) = i;
    locations(1, i*4) = ileft[i+1]-4;
    locations(1, i*4+1) = ileft[i+1]-3;
    locations(1, i*4+2) = ileft[i+1]-2;
    locations(1, i*4+3) = ileft[i+1]-1;
  }
  sp_mat G1(locations, vnikx);

  return G1;
}
