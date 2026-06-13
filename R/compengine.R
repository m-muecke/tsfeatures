#' CompEngine feature set
#'
#' Calculate the features that have been used in CompEngine database, using method introduced in package
#' \code{hctsa}.
#'
#' The features involved can be grouped as \code{autocorrelation},
#' \code{prediction}, \code{stationarity}, \code{distribution}, and \code{scaling}.
#'
#' @param x the input time series
#' @return a vector with CompEngine features
#' @seealso \code{\link{autocorr_features}}
#' @seealso \code{\link{pred_features}}
#' @seealso \code{\link{station_features}}
#' @seealso \code{\link{dist_features}}
#' @seealso \code{\link{scal_features}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
compengine <- function(x) {
  c(
    autocorr_features(x),
    pred_features(x),
    station_features(x),
    dist_features(x),
    scal_features(x)
  )
}

#' The autocorrelation feature set from software package \code{hctsa}
#'
#' Calculate the features that grouped as autocorrelation set,
#' which have been used in CompEngine database, using method introduced in package \code{hctsa}.
#'
#' Features in this set are \code{embed2_incircle_1},
#' \code{embed2_incircle_2},
#' \code{ac_9},
#' \code{firstmin_ac},
#' \code{trev_num},
#' \code{motiftwo_entro3},
#' and \code{walker_propcross}.
#'
#' @param x the input time series
#' @return a vector with autocorrelation features
#' @seealso \code{\link{embed2_incircle}}
#' @seealso \code{\link{ac_9}}
#' @seealso \code{\link{firstmin_ac}}
#' @seealso \code{\link{trev_num}}
#' @seealso \code{\link{motiftwo_entro3}}
#' @seealso \code{\link{walker_propcross}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
autocorr_features <- function(x) {
  acfv <- stats::acf(x, length(x) - 1, plot = FALSE, na.action = na.pass)
  output <- c(
    embed2_incircle_1 = embed2_incircle(x, 1, acfv = acfv),
    embed2_incircle_2 = embed2_incircle(x, 2, acfv = acfv),
    ac_9 = ac_9(x, acfv),
    firstmin_ac = firstmin_ac(x, acfv),
    trev_num = trev_num(x),
    motiftwo_entro3 = motiftwo_entro3(x),
    walker_propcross = walker_propcross(x)
  )
  return(output)
}

#' The prediction feature set from software package \code{hctsa}
#'
#' Calculate the features that grouped as prediction set,
#' which have been used in CompEngine database, using method introduced in package \code{hctsa}.
#'
#' Features in this set are \code{localsimple_mean1},
#' \code{localsimple_lfitac},
#' and \code{sampen_first}.
#'
#' @param x the input time series
#' @return a vector with prediction features
#' @seealso \code{\link{localsimple_taures}}
#' @seealso \code{\link{sampen_first}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
pred_features <- function(x) {
  output <- c(
    localsimple_mean1 = localsimple_taures(x, "mean"),
    localsimple_lfitac = localsimple_taures(x, "lfit"),
    sampen_first = sampen_first(x)
  )
  return(output)
}

#' The stationarity feature set from software package \code{hctsa}
#'
#' Calculate the features that grouped as stationarity set,
#' which have been used in CompEngine database, using method introduced in package \code{hctsa}.
#'
#' Features in this set are \code{std1st_der},
#' \code{spreadrandomlocal_meantaul_50},
#' and \code{spreadrandomlocal_meantaul_ac2}.
#'
#' @param x the input time series
#' @return a vector with stationarity features
#' @seealso \code{\link{std1st_der}}
#' @seealso \code{\link{spreadrandomlocal_meantaul}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
station_features <- function(x) {
  output <- c(
    std1st_der = std1st_der(x),
    spreadrandomlocal_meantaul_50 = spreadrandomlocal_meantaul(x, 50),
    spreadrandomlocal_meantaul_ac2 = spreadrandomlocal_meantaul(x, "ac2")
  )
  return(output)
}

#' The distribution feature set from software package \code{hctsa}
#'
#' Calculate the features that grouped as distribution set,
#' which have been used in CompEngine database, using method introduced in package \code{hctsa}.
#'
#' Features in this set are \code{histogram_mode_10}
#' and \code{outlierinclude_mdrmd}.
#'
#' @param x the input time series
#' @return a vector with distribution features
#' @seealso \code{\link{histogram_mode}}
#' @seealso \code{\link{outlierinclude_mdrmd}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
dist_features <- function(x) {
  output <- c(
    histogram_mode_10 = histogram_mode(x),
    outlierinclude_mdrmd = outlierinclude_mdrmd(x)
  )
  return(output)
}

#' The scaling feature set from software package \code{hctsa}
#'
#' Calculate the features that grouped as scaling set,
#' which have been used in CompEngine database, using method introduced in package \code{hctsa}.
#'
#' Feature in this set is \code{fluctanal_prop_r1}.
#'
#' @param x the input time series
#' @return a vector with scaling features
#' @seealso \code{\link{fluctanal_prop_r1}}
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
scal_features <- function(x) {
  output <- c(fluctanal_prop_r1 = fluctanal_prop_r1(x))
  return(output)
}

# autocorr ----------------------------------------------------------------

# CO_Embed2_Basic_tau_incircle_1
# CO_Embed2_Basic_tau_incircle_1
#' Points inside a given circular boundary in a 2-d embedding space from software package \code{hctsa}
#'
#' The time lag is set to the first zero crossing of the autocorrelation function.
#'
#' @param y the input time series
#' @param boundary the given circular boundary, setting to 1 or 2 in CompEngine. Default to 1.
#' @param acfv vector of autocorrelation, if exist, used to avoid repeated computation.
#' @return the proportion of points inside a given circular boundary
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
embed2_incircle <- function(
    y,
    boundary = NULL,
    acfv = stats::acf(y, length(y) - 1, plot = FALSE, na.action = na.pass)
) {
  if (is.null(boundary)) {
    warning(
      "`embed2_incircle()` using `boundary = 1`. Set value with `boundary`."
    )
    boundary <- 1
  }
  tau <- firstzero_ac(y, acfv)
  xt <- y[1:(length(y) - tau)] # part of the time series
  xtp <- y[(1 + tau):length(y)] # time-lagged time series
  N <- length(y) - tau # Length of each time series subsegment
  
  # CIRCLES (points inside a given circular boundary)
  return(sum(xtp^2 + xt^2 < boundary, na.rm = TRUE) / N)
}

# CO_firstzero_ac
#' The first zero crossing of the autocorrelation function from software package \code{hctsa}
#'
#' Search up to a maximum of the length of the time series
#'
#' Accelerated implementation: for a complete (no-NA) series with no
#' user-supplied \code{acfv}, dispatches to \code{firstzero_ac_cpp()}, which
#' computes autocorrelations lag-by-lag and stops at the first crossing.
#' Otherwise the original pure-R logic is used.
#'
#' @param y the input time series
#' @param acfv vector of autocorrelation, if exist, used to avoid repeated computation.
#' @return The first zero crossing of the autocorrelation function
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
firstzero_ac <- function(
    y,
    acfv = stats::acf(y, N - 1, plot = FALSE, na.action = na.pass)
) {
  N <- length(y)
  # Fast path: lazy default `acfv` is never evaluated here, so the full
  # N-1 lag ACF is skipped entirely.
  if (missing(acfv) && !anyNA(y)) {
    return(firstzero_ac_cpp(as.numeric(y)))
  }
  tau <- which(acfv$acf[-1] < 0)
  if (length(tau) == 0L) {
    # Nothing to see here
    return(0)
  } else if (all(is.na(tau))) {
    # All missing
    return(0)
  } else if (!any(tau)) {
    # No negatives, so set output to sample size
    return(N)
  } else {
    # Return lag of first negative
    return(tau[1])
  }
}

# ac_9
#' Autocorrelation at lag 9. Included for completion and consistency.
#'
#' @param y the input time series
#' @param acfv vector of autocorrelation, if exist, used to avoid repeated computation.
#' @return autocorrelation at lag 9
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
ac_9 <- function(
    y,
    acfv = stats::acf(y, 9, plot = FALSE, na.action = na.pass)
) {
  acfv$acf[10]
}

# CO_firstmin_ac
#' Time of first minimum in the autocorrelation function from software package \code{hctsa}
#'
#' Accelerated implementation: for a complete (no-NA) series with no user-supplied \code{acfv}, dispatches to \code{firstmin_ac_cpp()}
#'
#' @param x the input time series
#' @param acfv vector of autocorrelation, if exist, used to avoid repeated computation.
#' @return The lag of the first minimum
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' firstmin_ac(WWWusage)
#' @export
firstmin_ac <- function(
    x,
    acfv = stats::acf(x, lag.max = N - 1, plot = FALSE, na.action = na.pass)
) {
  # hctsa uses autocorr in MatLab to calculate autocorrelation
  N <- length(x)
  # Fast path: skips computing the full N-1 lag ACF.
  if (missing(acfv) && !anyNA(x)) {
    return(firstmin_ac_cpp(as.numeric(x)))
  }
  # getting acf for all lags
  # possible delay when sample size is too big
  autoCorr <- numeric(N - 1)
  autoCorr[1:(N - 1)] <- acfv$acf[-1]
  for (i in 1:length(autoCorr)) {
    if (is.na(autoCorr[i])) {
      warning("No minimum was found.")
      return(NA)
    }
    if (i == 2 && autoCorr[2] > autoCorr[1]) {
      return(1)
    } else if (
      i > 2 &&
      autoCorr[i - 2] > autoCorr[i - 1] &&
      autoCorr[i - 1] < autoCorr[i]
    ) {
      return(i - 1)
    }
  }
  return(N - 1)
}

# CO_trev_1_num
#' Normalized nonlinear autocorrelation, the numerator of the trev function of a time series from software package \code{hctsa}
#'
#' Calculates the numerator of the trev function, a normalized nonlinear autocorrelation,
#' The time lag is set to 1.
#'
#'
#' @param y the input time series
#' @return the numerator of the trev function of a time series
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' trev_num(WWWusage)
#' @export
trev_num <- function(y) {
  yn <- y[1:(length(y) - 1)]
  yn1 <- y[2:length(y)]
  mean((yn1 - yn)^3, na.rm = TRUE)
}

# SB_MotifTwo_mean_hhh
#' Local motifs in a binary symbolization of the time series from software package \code{hctsa}
#'
#'
#' Coarse-graining is performed. Time-series values above its mean are given 1,
#' and those below the mean are 0.
#'
#' @param y the input time series
#' @return Entropy of words in the binary alphabet of length 3.
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @examples
#' motiftwo_entro3(WWWusage)
#' @export
#'
motiftwo_entro3 <- function(y) {
  yBin <- binarize_mean(y)
  N <- length(yBin)
  if (N < 5) {
    warning("Time series too short")
  }
  
  r1 <- yBin == 1
  r0 <- yBin == 0
  
  r1 <- r1[1:(length(r1) - 1)]
  r0 <- r0[1:(length(r0) - 1)]
  
  r00 <- r0 & yBin[2:N] == 0
  r01 <- r0 & yBin[2:N] == 1
  r10 <- r1 & yBin[2:N] == 0
  r11 <- r1 & yBin[2:N] == 1
  
  r00 <- r00[1:(length(r00) - 1)]
  r01 <- r01[1:(length(r01) - 1)]
  r10 <- r10[1:(length(r10) - 1)]
  r11 <- r11[1:(length(r11) - 1)]
  
  r000 <- r00 & yBin[3:N] == 0
  r001 <- r00 & yBin[3:N] == 1
  r010 <- r01 & yBin[3:N] == 0
  r011 <- r01 & yBin[3:N] == 1
  r100 <- r10 & yBin[3:N] == 0
  r101 <- r10 & yBin[3:N] == 1
  r110 <- r11 & yBin[3:N] == 0
  r111 <- r11 & yBin[3:N] == 1
  
  out.ddd <- mean(r000)
  out.ddu <- mean(r001)
  out.dud <- mean(r010)
  out.duu <- mean(r011)
  out.udd <- mean(r100)
  out.udu <- mean(r101)
  out.uud <- mean(r110)
  out.uuu <- mean(r111)
  ppp <- c(
    out.ddd,
    out.ddu,
    out.dud,
    out.duu,
    out.udd,
    out.udu,
    out.uud,
    out.uuu
  )
  out.hhh <- f_entropy(ppp)
  return(out.hhh)
}

# BF_BF_binarize_mean
#' Converts an input vector into a binarized version from software package \code{hctsa}
#'
#' @param y the input time series
#' @return Time-series values above its mean are given 1, and those below the mean are 0.
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export

binarize_mean <- function(y) {
  y <- y - mean(y)
  Y <- numeric(length(y))
  Y[y > 0] <- 1
  return(Y)
}

f_entropy <- function(x) {
  # entropy of a set of counts, log(0)=0
  -sum(x[x > 0] * log(x[x > 0]))
}

# PH_Walker_prop_01_sw_propcross
#' Simulates a hypothetical walker moving through the time domain from software package \code{hctsa}
#'
#' The hypothetical particle (or 'walker') moves in response to values of the
#' time series at each point.
#' The walker narrows the gap between its value and that
#' of the time series by 10%.
#'
#' Accelerated implementation: the sequential walker recurrence is computed in C++ via \code{walker_propcross_cpp()}
#'
#' @param y the input time series
#' @return fraction of time series length that walker crosses time series
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
#'
#'
walker_propcross <- function(y) {
  walker_propcross_cpp(as.numeric(y))
}

# pred --------------------------------------------------------------------

# FC_localsimple_mean1_taures
# FC_localsimple_lfit_taures
#' The first zero crossing of the autocorrelation function of the residuals from Simple local time-series forecasting from software package \code{hctsa}
#'
#' Simple predictors using the past trainLength values of the time series to
#' predict its next value.
#'
#' Accelerated implementation: for a complete (no-NA) series, dispatches to \code{localsimple_taures_cpp()}
#'
#' @param y the input time series
#' @param forecastMeth the forecasting method, default to \code{mean}.
#' \code{mean}: local mean prediction using the past trainLength time-series values.
#' \code{lfit}: local linear prediction using the past trainLength time-series values.
#' @param trainLength the number of time-series values to use to forecast the next value.
#' Default to 1 when using method \code{mean} and 3 when using method \code{lfit}.
#' @return The first zero crossing of the autocorrelation function of the residuals
#' @export
localsimple_taures <- function(
    y,
    forecastMeth = c("mean", "lfit"),
    trainLength = NULL
) {
  forecastMeth <- match.arg(forecastMeth)
  # Fast path
  if (!anyNA(y)) {
    return(localsimple_taures_cpp(
      as.numeric(y),
      forecastMeth,
      if (is.null(trainLength)) -1L else as.integer(trainLength)
    ))
  }
  if (is.null(trainLength)) {
    lp <- switch(forecastMeth, mean = 1, lfit = firstzero_ac(y))
  }
  
  N <- length(y)
  evalr <- (lp + 1):N
  
  if (lp >= length(y)) {
    stop("Time series too short for forecasting in `localsimple_taures`")
  }
  
  res <- numeric(length(evalr))
  if (forecastMeth == "mean") {
    for (i in 1:length(evalr)) {
      res[i] <- mean(y[(evalr[i] - lp):(evalr[i] - 1)]) - y[evalr[i]]
    }
  }
  if (forecastMeth == "lfit") {
    for (i in 1:length(evalr)) {
      # Fit linear
      a <- 1:lp
      b <- y[(evalr[i] - lp):(evalr[i] - 1)]
      lm.ab <- lm(b ~ a, data = data.frame(a, b))
      res[i] <- predict(lm.ab, newdata = data.frame(a = lp + 1)) - y[evalr[i]]
      # p = polyfit((1:lp)',y(evalr(i)-lp:evalr(i)-1),1)
      #       res(i) = polyval(p,lp+1) - y(evalr(i)); % prediction - value
    }
  }
  out.taures <- firstzero_ac(res)
  return(out.taures)
}

# EN_SampEn_5_03_sampen1
#' Second Sample Entropy of a time series from software package \code{hctsa}
#'
#' Modified from the Ben Fulcher's \code{EN_SampEn} which uses code from PhysioNet.
#' The publicly-available PhysioNet Matlab code, sampenc (renamed here to
#' RN_sampenc) is available from:
#' http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
#'
#' Embedding dimension is set to 5.
#' The threshold is set to 0.3.
#'
#' Accelerated implementation: computed in C++ via \code{sampenc_cpp()}.
#'
#' @param y the input time series
#' @references cf. "Physiological time-series analysis using approximate entropy and sample
#' entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
#' Physiol., 278(6) H2039 (2000)
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
sampen_first <- function(y) {
  M <- 5
  r <- 0.3
  sampEn <- sampenc(y, M + 1, r)
  return(sampEn)
}

# PN_sampenc
#' Second Sample Entropy from software package \code{hctsa}
#'
#' Modified from the Ben Fulcher version of original code sampenc.m from
#' http://physionet.org/physiotools/sampen/
#' http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
#' Code by DK Lake (dlake@virginia.edu), JR Moorman and Cao Hanqing.
#'
#' Accelerated implementation: the loops are in C++ via \code{sampenc_cpp()}
#'
#' @param y the input time series
#' @param M embedding dimension
#' @param r threshold
#'
#' @references cf. "Physiological time-series analysis using approximate entropy and sample
#' entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
#' Physiol., 278(6) H2039 (2000)
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
sampenc <- function(y, M = 6, r = 0.3) {
  sampenc_cpp(as.numeric(y), as.integer(M), r)
}

# stationarity ------------------------------------------------------------

# SY_StdNthDer_1
#' Standard deviation of the first derivative of the time series from software package \code{hctsa}
#'
#' Modified from \code{SY_StdNthDer} in \code{hctsa}. Based on an idea by Vladimir Vassilevsky.
#'
#' @param y the input time series. Missing values will be removed.
#' @return Standard deviation of the first derivative of the time series.
#' @references cf. http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
std1st_der <- function(y) {
  if (length(y) < 2) {
    stop("Time series is too short to compute differences")
  }
  yd <- diff(y)
  return(sd(yd, na.rm = TRUE))
}

# SY_SpreadRandomLocal_50_100_meantaul
# SY_SpreadRandomLocal_ac2_100_meantaul
#'  Bootstrap-based stationarity measure from software package \code{hctsa}
#'
#' 100 time-series segments of length \code{l} are selected at random from the time series and
#' the mean of the first zero-crossings of the autocorrelation function in each segment is calculated.
#'
#' Accelerated implementation: for a complete (no-NA) series, dispatches to \code{spreadrandomlocal_meantaul_cpp()}
#'
#' @param y the input time series
#' @param l the length of local time-series segments to analyse as a positive integer. Can also be a specified character string: "ac2": twice the first zero-crossing of the autocorrelation function
#' @return mean of the first zero-crossings of the autocorrelation function
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
spreadrandomlocal_meantaul <- function(y, l = 50) {
  if (is.character(l) && "ac2" %in% l) {
    l <- 2 * firstzero_ac(y)
  }
  if (!is.numeric(l)) {
    stop("Unknown specifier `l`")
  }
  numSegs <- 100
  N <- length(y)
  if (l > 0.9 * N) {
    warning(
      "This time series is too short. Specify proper segment length in `l`"
    )
    return(NA_real_)
  }
  
  if (!anyNA(y)) {
    ists <- sample(N - 1 - l, numSegs, replace = TRUE)
    return(spreadrandomlocal_meantaul_cpp(
      as.numeric(y), as.integer(l), as.integer(ists)
    ))
  }
  
  qs <- numeric(numSegs)
  
  for (j in 1:numSegs) {
    ist <- sample(N - 1 - l, 1)
    ifh <- ist + l - 1
    rs <- ist:ifh
    ysub <- y[rs]
    taul <- firstzero_ac(ysub)
    qs[j] <- taul
  }
  return(mean(qs, na.rm = TRUE))
}

# distribution ------------------------------------------------------------

# DN_histogram_mode_10
#' Mode of a data vector from software package \code{hctsa}
#'
#' Measures the mode of the data vector using histograms with a given number of bins as suggestion.
#' The value calculated is different from \code{hctsa} and \code{CompEngine} as the histogram edges are calculated differently.
#'
#' @param y the input data vector
#' @param numBins the number of bins to use in the histogram.
#' @return the mode
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
#' @importFrom graphics hist
#' @importFrom stats predict

histogram_mode <- function(y, numBins = 10) {
  # Compute the histogram from the data:
  if (is.numeric(numBins)) {
    histdata <- hist(y, plot = FALSE, breaks = numBins)
    binCenters <- histdata$mids
  } else {
    stop("Unknown format for numBins")
  }
  # Compute bin centers from bin edges:
  # binCenters <- mean([binEdges(1:end-1) binEdges(2:end)])
  # Mean position of maximums (if multiple):
  out <- mean(binCenters[which.max(histdata$counts)])
  return(out)
}

# DN_OutlierInclude_abs_001_mdrmd
#' How median depend on distributional outliers from software package \code{hctsa}
#'
#' Measures median as more and
#' more outliers are included in the calculation according to a specified rule,
#' of outliers being furthest from the mean.
#'
#' The threshold for including time-series data points in the analysis increases
#' from zero to the maximum deviation, in increments of 0.01*sigma (by default),
#' where sigma is the standard deviation of the time series.
#'
#' At each threshold,  proportion of time series points
#' included and median are calculated, and outputs from the
#' algorithm measure how these statistical quantities change as more extreme
#' points are included in the calculation.
#'
#' Outliers are defined as furthest from the mean.
#'
#' @param y the input time series (ideally z-scored)
#' @param zscored Should y be z-scored before computing the statistic. Default: TRUE
#' @return median  of the median of range indices
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
#' @importFrom stats ts tsp sd

outlierinclude_mdrmd <- function(y, zscored = TRUE) {
  if (length(unique(y)) == 1L) {
    stop("The time series is a constant!")
  }
  if (zscored) {
    tmp <- ts(c(scale(y)))
    tsp(tmp) <- tsp(y)
    y <- tmp
    isd <- 1
  } else {
    isd <- sd(y, na.rm = TRUE) # Modified to fit the 0.01*sigma increment in description
  }
  N <- length(y)
  inc <- 0.01 * isd
  # inc <- 0.01
  thr <- seq(from = 0, to = max(abs(y), na.rm = TRUE), by = inc)
  tot <- N
  if (length(thr) == 0) {
    stop("peculiar time series")
  }
  
  msDt <- numeric(length(thr))
  msDtp <- numeric(length(thr))
  for (i in 1:length(thr)) {
    th <- thr[i] # the threshold
    # Construct a time series consisting of inter-event intervals for parts
    # of the time serie exceeding the threshold, th
    r <- which(abs(y) >= th)
    
    Dt_exc <- diff(r) # Delta t (interval) time series exceeding threshold
    msDt[i] <- median(r) / (N / 2) - 1
    msDtp[i] <- length(Dt_exc) / tot * 100
    # this is just really measuring the distribution:
    # the proportion of possible values
    # that are actually used in
    # calculation
  }
  
  # Trim off where the statistic power is lacking: less than 2% of data
  # included
  trimthr <- 2 # percent
  mj <- which(msDtp > trimthr)[length(which(msDtp > trimthr))]
  if (length(mj) != 0) {
    msDt <- msDt[1:mj]
    msDtp <- msDtp[1:mj]
    thr <- thr[1:mj]
  } else {
    stop("the statistic power is lacking: less than 2% of data included")
  }
  
  out.mdrmd <- median(msDt)
  return(out.mdrmd)
}

# scaling ----------------------------------------------------------------

# SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1
#' Implements fluctuation analysis from software package \code{hctsa}
#'
#' Fits a polynomial of order 1 and then returns the
#' range. The order of fluctuations is 2, corresponding to root mean
#' square fluctuations.
#'
#' Accelerated implementation: computed in C++ via \code{fluctanal_prop_r1_cpp()}
#'
#' @param x the input time series (or any vector)
#' @references B.D. Fulcher and N.S. Jones. hctsa: A computational framework for automated time-series phenotyping using massive feature extraction. Cell Systems 5, 527 (2017).
#' @references B.D. Fulcher, M.A. Little, N.S. Jones Highly comparative time-series analysis: the empirical structure of time series and their methods. J. Roy. Soc. Interface 10, 83 (2013).
#' @author Yangzhuoran Yang
#' @export
fluctanal_prop_r1 <- function(x) {
  fluctanal_prop_r1_cpp(as.numeric(x))
}
