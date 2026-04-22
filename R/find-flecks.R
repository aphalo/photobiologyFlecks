#' Find flecks
#'
#' Detects and characterizes "flecks" in a time series of irradiances.
#'
#' @details A sunfleck is characterized by an increase followed by a decrease
#' in irradiance. As a first step zero crossings of the derivative are located
#' using the same code as in \code{\link{find_zeros}()}. In a second step,
#' flecks are searched and once found, are checked for asymmetry
#' between baselines. If found, tries to extend baseline a bit. If still
#' asymmetric, behaviour is as defined by \code{asmMethod}. Then fleck is
#' compared againts criteria given by \code{minTime}, \code{minAmp} and
#' \code{minPdiff}. If passed, the baselines are trimmed. Conditions are
#' checked once more, and trimming is reversed if conditions are not passed
#' any more. In the last step overlapping flecks are removed and the time
#' interval between successive flecks is computed. Finally, the fleck properties
#' are returned in a data frame.
#'
#' @inheritParams find_zeros
#' @param minTime numeric Flecks found with duration below \code{minTime} will
#'   be discarded. Keep at 0 to keep all flecks.
#' @param minAmp numeric Flecks found with amplitude below \code{minAmp} will be
#'   discarded. Keep at 0 to keep all flecks.
#' @param minPdiff numeric Flecks found with a percent difference between peak
#'   and baseline that is below \code{minPdiff} will be discarded. Keep at 0 to
#'   keep all flecks.
#' @param asymmetry numeric Threshold value to qualify fleck as asymeric.
#'   Asymetry happens when the baseline have large difference in their value.
#' @param trimCV Control trimming threshold. Fleck baselines are trimmed
#'   iteratively based on the coefficient of variation between two points  at
#'   each baseline side.
#' @param asmMethod character One of \code{"mean"}, \code{"max"}, \code{"rm"}.
#'   Decides what to do with asymetric flecks \code{"mean"} averages the two
#'   baselines, \code{"max"} keeps the largest baseline, \code{"rm"} discard the
#'   asymetric fleck.
#' @param bounds numeric vector of length 2 For relative amplitude calculations,
#'   normalize between 0 and 1 by default.
#' @param timeSplit integer Increase time-series data frequency by linear
#'   interpolation. A value of 10 usually guarantee accuracy.
#' @param shadeflecks logical If true, run the function in \emph{shadefleck}
#'   mode instead of in the default \emph{sunfleck} mode, i.e., the function
#'   will find troughs instead of peaks in the data.
#' @param time.digits,var.digits integer Argument passed to parameter
#'   \code{digits} in internal calls to \code{\link{round}()} for \code{time}
#'   and \code{var} values.
#' @param verbose logical If \code{TRUE}, provides more information while
#'   running.
#'
#' @references
#'
#' Durand M, Matule B, Burgess AJ, Robson TM. 2021. Sunfleck properties from
#'   time series of fluctuating light. \emph{Agricultural and Forest Meteorology}
#'   \strong{308-309}, 108554. \doi{10.1016/j.agrformet.2021.108554}
#'
#' @export
#'
find_flecks <- function(time,
                        var,
                        zero.lim = 0.0005,
                        minTime = 0,
                        minAmp = 0,
                        minPdiff = 0,
                        asymmetry = 1/4,
                        trimCV = 0.05,
                        asmMethod = c("mean", "max", "rm"),
                        bounds = c(0, 1),
                        timeSplit = 10,
                        shadeflecks = FALSE,
                        time.digits = 3,
                        var.digits = 2,
                        verbose = TRUE) {

  # find zeros here so that interpolation is done only once
  # linear interpolation of whole series is use to increase time resolution
  # (note: interpolation is strictly needed only at the zeros!)
  timeStep <- time[2] - time[1]
  x_out <- seq(min(time), max(time), timeStep / timeSplit)
  int.df <- stats::approx(x = time, y = var, xout = x_out)
  sgf <- diff(int.df[["y"]]) / diff(int.df[["x"]])
  n1n2 <- sgf[1:(length(sgf) - 1)] * sgf[2:length(sgf)]
  zeros <- which(n1n2 < 0 & abs(n1n2) > zero.lim) + 1L
  names(zeros) <- ifelse(sgf[zeros] < 0, "up", "low")

  # Annotate with baseline higher or lower than peak
  if (shadeflecks) {
    zz <- zeros[names(zeros) == "up"]
  } else {
    zz <- zeros[names(zeros) == "low"]
  }

  # input
  timeInt <- int.df$x
  varInt <- int.df$y

  # Initialize output
  l0 <- matrix(nrow = length(zeros), ncol = 17)
  l0 <- as.data.frame(l0)
  l0 <- as.list(l0)

  # index
  no <- 0

  for (i in 1:length(zz)) {
    flag <- 0
    newi <- which(zz[i] == zeros)

    # Skip the last data points
    if ((newi + 2) > length(zeros)) {
      next()
    }

    b1 <- zeros[newi]
    b2 <- zeros[newi + 2]
    p <- zeros[newi + 1]

    Xb1 <- round(timeInt[b1], time.digits)
    Xb2 <- round(timeInt[b2], time.digits)
    Xp <- round(timeInt[p], time.digits)

    Yb1 <- round(varInt[b1], var.digits)
    Yb2 <- round(varInt[b2], var.digits)
    Yp <- round(varInt[p], var.digits)

    ## Check if left and right baseline are similar, if not look for points before or after
    # Left side too small
    if (abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry) {
      flag <- 1
      pi <- 0
      while (abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry & pi < 3) {
        pi <- pi + 1
        b1 <- zeros[newi - pi]
        Xb1 <- round(timeInt[b1], time.digits)
        Yb1 <- round(varInt[b1], var.digits)

        # For first point
        if (length(b1) == 0) {
          b1 <- zeros[newi]
          Xb1 <- round(timeInt[b1], time.digits)
          Yb1 <- round(varInt[b1], var.digits)
          break
        }
      }
    }
    # Right side too small
    if (abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry) {
      flag <- 1
      pi <- 0
      while(abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry & pi < 3)
      {
        pi = pi + 1
        b2 <- zeros[newi + 2 + pi]
        Xb2 <- round(timeInt[b2], time.digits)
        Yb2 <- round(varInt[b2], var.digits)

        # For first point
        if (is.na(b2)) {
          b2 <- zeros[newi + 2]
          Xb2 <- round(timeInt[b2], time.digits)
          Yb2 <- round(varInt[b2], var.digits)
          break}
      }
    }

    # If left and right baseline are still not similar, use return to original values and flag
    if(abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry |
       abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry) {
      ASYMMETRIC <- TRUE
      if (asmMethod != "rm") {
        b1 <- zeros[newi]
        b2 <- zeros[newi + 2]
        Xb1 <- round(timeInt[b1], time.digits)
        Xb2 <- round(timeInt[b2], time.digits)
        Yb1 <- round(varInt[b1], time.digits)
        Yb2 <- round(varInt[b2], time.digits)
      } else {next()}
    } else {
      ASYMMETRIC <- FALSE
    }

    # Skip points where baseline is higher than peak
    if (shadeflecks == FALSE) {
      if (Yb1 > Yp | Yb2 > Yp) {
        next()
      }
    } else {
      if (Yb1 < Yp | Yb2 < Yp) {
        next()
      }
    }

    # Initial calculation of fleck data
    if (ASYMMETRIC == TRUE) {
      if (asmMethod == "max" & shadeflecks == TRUE) {
        Ybaseline <-  max(c(Yb1, Yb2))
      }
      if (asmMethod == "max" & shadeflecks == FALSE) {
        Ybaseline <-  min(c(Yb1, Yb2))
      }
      if (asmMethod == "mean") {
        Ybaseline <- mean(c(Yb1, Yb2))
      }
    } else {
      Ybaseline <- mean(c(Yb1, Yb2))
    }

    # Get closest baseline to run condition check on that
    if (shadeflecks == TRUE) {
      YCloseBsl <- min(c(Yb1, Yb2))
    } else {
      YCloseBsl <- max(c(Yb1, Yb2))
    }

    amp <- round(abs(Yp - YCloseBsl), 2)
    pdiff <- round(amp / YCloseBsl, 4)
    duration <- round(abs(Xb2 - Xb1), 3)
    irrad <- sum(varInt[seq(b1,b2,timeSplit)]) * timeStep
    interpBaseline <- stats::approx(x = c(timeInt[b1], timeInt[b2]),
                                    y = c(varInt[b1], varInt[b2]),
                                    xout = timeInt[seq(b1, b2, timeSplit)])
    irradBaseline <- sum(interpBaseline$y) * timeStep
    irradSunfleck <- irrad - irradBaseline

    # If fleck pass criteria, then trim and recalculate
    if (duration > minTime & amp > minAmp & pdiff > minPdiff) {
      # Left trimming
      it <- 0
      trimLeft <- 0
      old_trimLeft = 1
      while (old_trimLeft != trimLeft) {
        it = it + 1

        if(b1 + it * timeSplit < length(varInt)){
          Yk1 = round(varInt[b1 + it * timeSplit], 2)
          Xk1 = round(timeInt[b1 + it * timeSplit], time.digits)

          old_trimLeft <- trimLeft
          if (stats::sd(c(Yb1, Yk1)) / mean(c(Yb1, Yk1)) < trimCV) {
            trimLeft = trimLeft + 1
          }

          # If trimming goes too far, abandon trimming
          if (Xk1 == Xp) {
            trimLeft <- 0
            old_trimLeft <- 0
          }
        } else {
          trimLeft <- 0
          old_trimLeft <- 0
        }
      }

      # Right trimming
      it <- 0
      trimRight <- 0
      old_trimRight <- 1
      while (old_trimRight != trimRight) {
        it = it + 1
        if (b2 - it * timeSplit > 0) {
          Yk2 = round(varInt[b2 - it * timeSplit], 2)
          Xk2 = round(timeInt[b2 - it * timeSplit], 3)

          old_trimRight <- trimRight
          if (stats::sd(c(Yb2,Yk2)) / mean(c(Yb2,Yk2)) < trimCV) {
            trimRight = trimRight + 1
          }

          # If trimming goes too far, abandon trimming
          if (Xk2 == Xp) {
            trimRight <- 0
            old_trimRight <- 0
          }
        } else {
          trimRight <- 0
          old_trimRight <- 0
        }
      }

      b1_new <- b1 + trimLeft * timeSplit
      b2_new <- b2 - trimRight * timeSplit

      # Another safeguard for trimming
      if (b2_new > b1_new) {
        b1 <- b1_new
        b2 <- b2_new
      }

      Xb1 <- round(timeInt[b1], time.digits)
      Xb2 <- round(timeInt[b2], time.digits)
      Yb1 <- round(varInt[b1], 2)
      Yb2 <- round(varInt[b2], 2)

      # Final calculation of fleck data
      if(ASYMMETRIC == TRUE) {
        if (asmMethod == "max" & shadeflecks == TRUE) {
          Ybaseline <-  max(c(Yb1, Yb2))
        }
        if (asmMethod == "max" & shadeflecks == FALSE) {
          Ybaseline <-  min(c(Yb1, Yb2))
        }
        if (asmMethod == "mean") {
          Ybaseline <- mean(c(Yb1, Yb2))
        }
      } else {
        Ybaseline <- mean(c(Yb1, Yb2))
      }
      amp <- round(abs(Yp - Ybaseline), 2)
      relAmp <- round(abs(((Yp - bounds[1]) / (bounds[2] - bounds[1])) -
                            ((Ybaseline - bounds[1]) / (bounds[2] - bounds[1]))),
                      4)
      pdiff <- round(amp / Ybaseline, 4)
      duration <- round(abs(Xb2 - Xb1), 4)

      irrad <- sum(varInt[seq(b1,b2,timeSplit)]) * timeStep
      interpBaseline <- stats::approx(x = c(timeInt[b1], timeInt[b2]),
                                      y = c(varInt[b1], varInt[b2]),
                                      xout = timeInt[seq(b1, b2, timeSplit)])
      irradBaseline <- sum(interpBaseline$y) * timeStep
      irradSunfleck <- irrad - irradBaseline

      relIrrad <- sum((varInt[seq(b1, b2, timeSplit)] -
                         bounds[1]) / (bounds[2] - bounds[1]))
      relBaseline  <- sum((interpBaseline$y - bounds[1]) /
                            (bounds[2] - bounds[1]))
      relIrradSunfleck <- relIrrad - relBaseline

      # Log data (May have to reverse trimming if conditions is not satisfied
      # anymore, will see later)
      no <- no + 1
      row <- c(no, Xp, Yp, Yb1, Yb2, Ybaseline, amp, relAmp, pdiff, duration,
               irradSunfleck, relIrradSunfleck, irrad, Xb1, Xb2,
               ifelse(ASYMMETRIC == T, "asymmetric", "symmetric"),
               zeros[newi + 1])

      for (j in 1:length(row)) {
        l0[[j]][no] <- row[j]
      }

      if (flag == 1 & verbose == TRUE) {
        cat("asymmetry found in no", no, "\n")
      }
    }
  }

  df0 <- as.data.frame(l0)
  df0 <- df0[!(is.na(df0[ , 1])),]

  for(i in 1:15) {
    df0[ , i] <- as.numeric(df0[ , i])
  }

  message("Found ", nrow(df0), ifelse(shadeflecks, " shadeflecks", " sunflecks"))
  if (nrow(df0) == 0) {
    message("You may want to try with different 'min...' values.")
    return(data.frame())
  }
  colnames(df0) <-  c("no", "peakTime", "peak",
                      "baseline1", "baseline2", "baseline",
                      "amplitude", "relAmp",
                      "percDiff", "duration",
                      "irradGain", "relIrrad", "irradTot",
                      "baselineTime1", "baselineTime2",
                      "symmetry", "z")

  # Remove overlapping sunflecks
  for(iROW in 1:nrow(df0))
  {
    if((iROW+1) > nrow(df0)){next()}
    endTime <- df0[iROW,"baselineTime2"]
    startTime <- df0[iROW+1,"baselineTime1"]

    if(startTime < endTime)
    {
      if(df0[iROW,"peak"] < df0[iROW+1,"peak"]){iRM <- df0[iROW,"no"]} else {iRM <- df0[iROW+1,"no"]}
      df0 <- df0[!df0$no == iRM,]
      if(verbose == TRUE){cat("Removed fleck no", iRM, "because of overlap", "\n")}
    }
  }

  # Calculate time between two sunflecks
  if(nrow(df0) > 1)
  {
    df0$timeInterval <- NA
    for (iROW in 2:nrow(df0)) {
      df0$timeInterval[iROW] <- (df0$peakTime[iROW] - df0$peakTime[iROW-1])
    }
  }
  return(df0)
}
