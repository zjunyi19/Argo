// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace std;
using namespace arma;

void bsplvb(NumericVector t, int jhigh, int index, float x, int left,
            mat& biatx, int& j, vec& deltar, vec& deltal) {
  if (index == 1) {
    j = 1;
    biatx.zeros();
    deltar.zeros();
    deltal.zeros();
    biatx(1,1) = 1.0;
    if (j >= jhigh) {return;}
  }

  while(j < jhigh) {
    int jp1 = j+1;
    deltar[j] = t[left+j] - x;
    deltal[j] = x - t[left+1-j];
    float saved = 0;
    for (int i = 1; i <= j; ++i) {
      float term = 0;
      if (deltar[i] + deltal[jp1-i] != 0) {
        term = biatx(i,1)/(deltar[i] + deltal[jp1-i]);
      }
      biatx(i,1) = saved + deltar[i] * term;
      saved = deltal[jp1-i] * term;
    }
    biatx(jp1,1) = saved;
    j = jp1;
  }
  return;
}

void bsplvd(NumericVector t, int k, float x, int left, mat& a, mat& dbiatx,
            int nderiv, vec& deltar, vec& deltal) {
  int mhigh = max(min(nderiv,k), 1);
  int kp1 = k+1;
  int j_set = 1;
  bsplvb(t, kp1-mhigh, 1, x, left, dbiatx, j_set, deltar, deltal);

  if (mhigh == 1) {return;}
  int ideriv = mhigh;
  for (int m = 2; m <= mhigh; ++m) {
    int jp1mid = 1;
    for (int j = ideriv; j <= k; ++j) {
      dbiatx(j, ideriv) = dbiatx(jp1mid, 1);
      jp1mid++;
    }
    ideriv--;
    bsplvb(t, kp1-ideriv, 2, x, left, dbiatx, j_set, deltar, deltal);

  }

  int jlow = 1;
  for (int i = 1; i <= k; ++i) {
    for (int j = jlow; j <= k; ++j) {
      a(j,i) = 0.0;
    }
    jlow = i;
    a(i,i) = 1.0;
  }

  for (int m = 2; m <= mhigh; ++m) {
    int kp1mm = kp1 - m;
    double fkp1mm = (double)kp1mm;
    int il = left;
    int i = k;

    for (int ldummy = 1; ldummy <= kp1mm; ++ldummy) {
      double factor = 0;
      if (t[il+kp1mm] - t[il] != 0) {
        factor = fkp1mm / (t[il+kp1mm] - t[il]);
      }
      for (int j = 1; j <= i; ++j) {
        a(i,j) = (a(i,j) - a(i-1,j)) * factor;
      }
      il--;
      i--;
    }

    for (int i = 1; i <= k; ++i) {
      double sum = 0.0;
      int jlow = max(i,m);
      for (int j = jlow; j <= k; ++j) {
        sum += a(j,i) * dbiatx(j,m);
      }
      dbiatx(i,m) = sum;
    }
  }
  return;
}

// [[Rcpp::export]]
sp_mat second_deriv_mat(NumericVector knots) {
  mat vnikx(5,4); // vnikx column 0 and row 0 are dummy variables
  vec deltar(21);
  vec deltal(21);
  vec yw1(5);
  vec yw2(5);
  mat a(5,5);

  NumericVector t(knots.size() + 6 + 1); //t[0] is dummy variable
  t[1] = t[2] = t[3] = knots[0];
  for (int i = 4; i <= (knots.size() + 3); ++i) {
    t[i] = knots[i-3-1];
  }
  t[knots.size()+4] = t[knots.size()+5] = t[knots.size()+6] = knots[knots.size()-1];
  int nb = knots.size() + 2;
  NumericVector sg0(nb+1); // sg0[0] is dummy variable
  NumericVector sg1(nb+1); // sg1[0] is dummy variable
  NumericVector sg2(nb+1); // sg2[0] is dummy variable
  NumericVector sg3(nb+1); // sg3[0] is dummy variable

  int ileft = 4;
  for (int i = 1; i <= nb; ++i) {
    yw1.zeros();
    yw2.zeros();
    if(i <= 4) {
      ileft = 4;
    } else {
      ileft = i;
    }
    //ileft = find_interval(t[i],t,knots.size(),ileft);
    // left end second derivatives
    bsplvd(t, 4, t[i], ileft, a, vnikx, 3, deltar, deltal);

    for (int ii=1; ii <= 4; ++ii) {
      yw1[ii] = vnikx(ii,3);
    }
    // right end second derivative
    bsplvd(t, 4, t[i+1],ileft, a, vnikx, 3, deltar, deltal);

    for (int ii=1; ii <= 4; ++ii) {
      yw2[ii] = vnikx(ii,3) - yw1[ii];
    }
    double wpt = t[i+1]-t[i];
    if(ileft >= 4) {
      for(int ii = 1; ii <= 4; ++ii) {
        int jj = ii;
        sg0[ileft-4+ii] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        jj = ii + 1;
        if (jj <= 4) {
          sg1[ileft+ii-4] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
        jj = ii + 2;
        if (jj <= 4) {
          sg2[ileft+ii-4] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
        jj = ii + 3;
        if (jj <= 4) {
          sg3[ileft+ii-4] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
      }
    } else if (ileft == 3) {
      for (int ii = 1; ii <= 3; ++ii) {
        int jj = ii;
        sg0[ileft-3+ii] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        jj = ii + 1;
        if (jj <= 3) {
          sg1[ileft+ii-3] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
        jj = ii + 2;
        if (jj <= 3) {
          sg2[ileft+ii-3] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
      }
    } else if (ileft == 2) {
      for (int ii = 1; ii <= 2; ++ii) {
        int jj = ii;
        sg0[ileft-2+ii] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        jj = ii + 1;
        if (jj <= 2) {
          sg1[ileft+ii-2] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
        }
      }
    } else if (ileft == 1) {
      int ii = 1;
      int jj = ii;
      sg0[ileft-1+ii] += wpt*(yw1[ii]*yw1[jj]+ (yw2[ii]*yw1[jj] + yw2[jj]*yw1[ii])*0.5 + yw2(ii)*yw2(jj) / 3.0);
    }
  }

  mat result(nb,nb);
  result.zeros();
  for (int i = 1; i <= nb; ++i) {
    result(i-1,i-1) = sg0[i];
    if(i < nb) {
      result(i-1,i) = sg1[i];
      result(i,i-1) = sg1[i];
    }
    if(i+1 < nb) {
      result(i-1,i+1) = sg2[i];
      result(i+1,i-1) = sg2[i];
    }
    if(i+2 < nb) {
      result(i-1,i+2) = sg3[i];
      result(i+2,i-1) = sg3[i];
    }
  }
  sp_mat temp(result);
  return temp;

}