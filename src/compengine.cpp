//
// C++ implementations of slow-to-evaluate features in the "CompEngine" set
// to speed up computational efficiency
//
// Author: Trent Henderson, 12 June 2026
//

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ----------------
// Internal helpers
// ----------------

// Sum of squared residuals from a simple linear regression of yv on xv over the inclusive index range [a, b] (0-based).
static double lin_sse(const std::vector<double>& xv,
                      const std::vector<double>& yv,
                      int a, int b) {
  int n = b - a + 1;
  double sx = 0.0, sy = 0.0;
  for (int i = a; i <= b; ++i) { sx += xv[i]; sy += yv[i]; }
  double xbar = sx / n, ybar = sy / n;
  double sxx = 0.0, sxy = 0.0;
  for (int i = a; i <= b; ++i) {
    double dx = xv[i] - xbar;
    sxx += dx * dx;
    sxy += dx * (yv[i] - ybar);
  }
  double slope = (sxx > 0.0) ? sxy / sxx : 0.0;
  double intercept = ybar - slope * xbar;
  double sse = 0.0;
  for (int i = a; i <= b; ++i) {
    double r = yv[i] - (intercept + slope * xv[i]);
    sse += r * r;
  }
  return sse;
}

// First zero crossing of the ACF
static int firstzero_ac_raw(const double* y, int N) {
  if (N < 2) return 0;
  double mu = 0.0;
  for (int i = 0; i < N; ++i) mu += y[i];
  mu /= N;
  double denom = 0.0;
  for (int i = 0; i < N; ++i) { double d = y[i] - mu; denom += d * d; }
  if (denom <= 0.0) return 0;
  for (int tau = 1; tau <= N - 1; ++tau) {
    double num = 0.0;
    for (int i = 0; i < N - tau; ++i)
      num += (y[i] - mu) * (y[i + tau] - mu);
    if (num < 0.0) return tau;
  }
  return 0;
}

// First zero crossing of the autocorrelation function
// @param y the input time series
// [[Rcpp::export]]
int firstzero_ac_cpp(NumericVector y) {
  return firstzero_ac_raw(REAL(y), y.size());
}

// ----------------------
// sampenc / sampen_first
// ----------------------

// Sample entropy counting
// @param y the input time series
// @param M embedding dimension
// @param r threshold
// [[Rcpp::export]]
double sampenc_cpp(NumericVector y, int M = 6, double r = 0.3) {
  int N = y.size();
  std::vector<double> lastrun(N, 0.0), run(N, 0.0);
  std::vector<double> A(M, 0.0), B(M, 0.0);

  for (int i = 0; i < N - 1; ++i) {
    double y1 = y[i];
    int nj = N - 1 - i;
    for (int jj = 0; jj < nj; ++jj) {
      int j = i + jj + 1;
      double d = y[j] - y1;
      if (!ISNAN(d) && std::fabs(d) < r) {
        run[jj] = lastrun[jj] + 1.0;
        int M1 = (int)std::min((double)M, run[jj]);
        for (int m = 0; m < M1; ++m) A[m] += 1.0;
        if (j < N - 1)
          for (int m = 0; m < M1; ++m) B[m] += 1.0;
      } else {
        run[jj] = 0.0;
      }
    }
    for (int jj = 0; jj < nj; ++jj) lastrun[jj] = run[jj];
  }

  double p = A[1] / B[0];
  return -std::log(p);
}

// C++ replacement for sampen_first
// @param y the input time series
// [[Rcpp::export]]
double sampen_first_cpp(NumericVector y) {
  return sampenc_cpp(y, 6, 0.3);
}

// ----------------
// walker_propcross
// ----------------

// Hypothetical walker crossing fraction
// @param y the input time series
// @param p the particle gap
// [[Rcpp::export]]
double walker_propcross_cpp(NumericVector y, double p = 0.1) {
  int N = y.size();
  if (N < 2) return NA_REAL;
  std::vector<double> w(N, 0.0);
  for (int i = 1; i < N; ++i)
    w[i] = w[i - 1] + p * (y[i - 1] - w[i - 1]);

  double cnt = 0.0;
  for (int i = 0; i < N - 1; ++i) {
    double prod = (w[i] - y[i]) * (w[i + 1] - y[i + 1]);
    if (!ISNAN(prod) && prod < 0.0) cnt += 1.0;
  }
  return cnt / (double)(N - 1);
}

// ------------------
// localsimple_taures
// ------------------

// Local simple forecasting residual tau
// @param y the input time series
// @param forecastMeth the forecasting method, defaults to \code{mean}.
// @param trainLength the training length. Defaults to \code{-1}
// [[Rcpp::export]]
int localsimple_taures_cpp(NumericVector y,
                           std::string forecastMeth = "mean",
                           int trainLength = -1) {
  int N = y.size();
  int lp;
  if (trainLength > 0) {
    lp = trainLength;
  } else if (forecastMeth == "mean") {
    lp = 1;
  } else if (forecastMeth == "lfit") {
    lp = firstzero_ac_raw(REAL(y), N);
  } else {
    stop("Unknown forecast method");
  }

  if (lp >= N)
    stop("Time series too short for forecasting in `localsimple_taures`");
  if (lp < 1)
    stop("Training window length is %d; series has no zero crossing of the ACF", lp);

  int nev = N - lp;
  std::vector<double> res(nev);

  if (forecastMeth == "mean") {
    double s = 0.0;
    for (int k = 0; k < lp; ++k) s += y[k];
    for (int i = 0; i < nev; ++i) {
      res[i] = s / lp - y[lp + i];
      s += y[lp + i] - y[i];
    }
  } else {
    double xbar = (lp + 1) / 2.0;
    double Sxx  = lp * ((double)lp * lp - 1.0) / 12.0;
    for (int i = 0; i < nev; ++i) {
      double sy = 0.0, sxy = 0.0;
      for (int k = 0; k < lp; ++k) {
        sy  += y[i + k];
        sxy += (k + 1) * y[i + k];
      }
      double ybar  = sy / lp;
      double slope = (Sxx > 0.0) ? (sxy - xbar * sy) / Sxx : 0.0;
      double intercept = ybar - slope * xbar;
      res[i] = (intercept + slope * (lp + 1)) - y[lp + i];
    }
  }

  return firstzero_ac_raw(res.data(), nev);
}

// -----------------
// fluctanal_prop_r1
// -----------------

// Fluctuation analysis
// @param x the input time series
// [[Rcpp::export]]
double fluctanal_prop_r1_cpp(NumericVector x) {
  const int tauStep = 50;
  int N = x.size();

  // cumulative sum with NA -> 0
  std::vector<double> y(N);
  double s = 0.0;
  for (int i = 0; i < N; ++i) {
    double v = ISNAN(x[i]) ? 0.0 : x[i];
    s += v;
    y[i] = s;
  }

  // taur = unique(round(exp(seq(log(5), log(floor(N/2)), length.out = 50))))

  std::vector<int> taur;
  double lo = std::log(5.0), hi = std::log(std::floor(N / 2.0));
  for (int i = 0; i < tauStep; ++i) {
    double v = lo + (hi - lo) * i / (double)(tauStep - 1);
    int t = (int)std::nearbyint(std::exp(v));
    if (taur.empty() || t != taur.back()) taur.push_back(t);
  }
  int ntau = (int)taur.size();
  if (ntau < 8)
    stop("This time series is too short to analyse using this fluctuation analysis");

  std::vector<double> Fl(ntau);
  for (int i = 0; i < ntau; ++i) {
    int tau  = taur[i];
    int nseg = N / tau;
    double xbar = (tau + 1) / 2.0;
    double Sxx  = tau * ((double)tau * tau - 1.0) / 12.0;

    double sumdq = 0.0;
    for (int seg = 0; seg < nseg; ++seg) {
      int off = seg * tau;
      double sy = 0.0, sxy = 0.0;
      for (int k = 0; k < tau; ++k) {
        sy  += y[off + k];
        sxy += (k + 1) * y[off + k];
      }
      double ybar  = sy / tau;
      double slope = (Sxx > 0.0) ? (sxy - xbar * sy) / Sxx : 0.0;
      double intercept = ybar - slope * xbar;

      double rmin =  R_PosInf, rmax = R_NegInf;
      for (int k = 0; k < tau; ++k) {
        double r = y[off + k] - (intercept + slope * (k + 1));
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
      }
      double dr = rmax - rmin;
      sumdq += dr * dr;
    }
    Fl[i] = std::sqrt(sumdq / nseg);
  }

  // Two-regime breakpoint search on the log-log fluctuation function
  
  std::vector<double> logtt(ntau), logFF(ntau);
  for (int i = 0; i < ntau; ++i) {
    logtt[i] = std::log((double)taur[i]);
    logFF[i] = std::log(Fl[i]);
  }

  const int minPoints = 6;
  double best = R_PosInf;
  int breakPt = NA_INTEGER;
  for (int i = minPoints; i <= ntau - minPoints; ++i) {
    double e = std::sqrt(lin_sse(logtt, logFF, 0, i - 1)) +
               std::sqrt(lin_sse(logtt, logFF, i - 1, ntau - 1));
    if (e < best) { best = e; breakPt = i; }
  }

  return (double)breakPt / (double)ntau;
}

// -----------
// firstmin_ac
// -----------

// First local minimum of the ACF
// @param y the input time series
// [[Rcpp::export]]
double firstmin_ac_cpp(NumericVector y) {
  int N = y.size();
  if (N < 2) return NA_REAL;
  const double* yp = REAL(y);

  double mu = 0.0;
  for (int i = 0; i < N; ++i) mu += yp[i];
  mu /= N;
  double denom = 0.0;
  for (int i = 0; i < N; ++i) { double d = yp[i] - mu; denom += d * d; }
  if (denom <= 0.0) { Rcpp::warning("No minimum was found."); return NA_REAL; }

  std::vector<double> ac;
  ac.reserve(64);
  for (int i = 1; i <= N - 1; ++i) {
    double num = 0.0;
    for (int k = 0; k < N - i; ++k)
      num += (yp[k] - mu) * (yp[k + i] - mu);
    ac.push_back(num / denom);

    if (i == 2 && ac[1] > ac[0]) return 1.0;
    if (i > 2 && ac[i - 3] > ac[i - 2] && ac[i - 2] < ac[i - 1])
      return (double)(i - 1);
  }
  return (double)(N - 1);
}

// --------------------------
// spreadrandomlocal_meantaul
// --------------------------

// Bootstrap stationarity measure
// @param y the input time series
// @param l the length of local time-series segments to analyse as a positive integer. Can also be a specified character string: "ac2": twice the first zero-crossing of the autocorrelation function
// @param starts start indices
// [[Rcpp::export]]
double spreadrandomlocal_meantaul_cpp(NumericVector y, int l, IntegerVector starts) {
  int N = y.size();
  if (l < 1)
    stop("Segment length `l` invalid for this series");
  const double* yp = REAL(y);
  int numSegs = starts.size();
  
  double total = 0.0;
  for (int j = 0; j < numSegs; ++j) {
    int ist = starts[j];
    if (ist < 1 || ist + l - 1 > N)
      stop("Start index %d out of range for segment length %d", ist, l);
    total += firstzero_ac_raw(yp + (ist - 1), l);
  }
  return total / (double)numSegs;
}
