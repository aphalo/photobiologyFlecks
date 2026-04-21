#' Find flecks
#'
#' Detects flecks from time series and zero vector returned by
#' \code{\link{find_zeros}(()}.
#'
#' @details A sunfleck is caracterized by an increased followed by a decrease.
#' Function search zero vector for such events. Once found, checks for asymetry
#' between baselines. If found, tries to extend baseline a bit. If still
#' asymeric, behaviour is as defined by \code{asmMethod}. Then fleck is
#' checks for criteria as defined by \code{minTime}, \code{minAmp} and
#' \code{minPdiff}. If passed, the baselines are trimmed. Conditions are
#' checked once more, and trimming is reversed if conditions are not passed
#' anymore. Finally, log the fleck in a table returned at the end. At the end,
#' remove flecks overlapping, and calculate interval in time between two flecks.
#'
#' @inheritParams find_zeros
#'
#' @param zeroes  Vector of zeroes returned by \code{find_zeroes()} on the same
#'   \code{time} and \code{var} arguments.
#' @param minTime numeric Flecks found with duration below `minTime` will be
#'   discarded. Keep at 0 to keep all flecks.
#' @param minAmp numeric Flecks found with amplitude below `minAmp` will be
#'   discarded. Keep at 0 to keep all flecks.
#' @param minPdiff numeric Flecks found with a percent difference between peak
#'   and baseline that is below `minPdiff` will be discarded. Keep at 0 to keep
#'   all flecks.
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
#'   mode instead of in the default \emph{sunfleck} mode,
#'   i.e., the function will find troughs instead of peaks in the data.
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
                        zeroes = find_zeros(time, var),
                        minTime = 0,
                        minAmp = 0,
                        minPdiff = 0,
                        asymmetry = 1/4,
                        trimCV = 0.05,
                        asmMethod = c("mean", "max", "rm"),
                        bounds = c(0, 1),
                        timeSplit = 10,
                        shadeflecks = FALSE,
                        verbose = TRUE) {
  no = 0

  # Define if baseline value is higher or lower than peak
  if(shadeflecks == FALSE) {upOrLow <- "low"} else {upOrLow <- "up"}
  zz <- zeroes[names(zeroes) == upOrLow]

  # Zeros are on the interpolated curve so
  # we have to redefine it here.
  timeStep <- time[2] - time[1]
  x_out <- seq(0, max(time), timeStep / timeSplit)
  int <- stats::approx(x = time, y = var, xout = x_out)

  timeInt <- int$x
  varInt <- int$y

  # Initialize output
  l0 <- matrix(nrow = length(zeroes), ncol = 17)
  l0 <- as.data.frame(l0)
  l0 <- as.list(l0)

  for (i in 1:length(zz))
  {
    flag = 0
    newi <- which(zz[i] == zeroes)

    # Skip the last data points
    if((newi+2) > length(zeroes)){next()}

    b1 <- zeroes[newi]
    b2 <- zeroes[newi + 2]
    p <- zeroes[newi + 1]
    Xb1 <- round(timeInt[b1], 3)
    Xb2 <- round(timeInt[b2], 3)
    Xp <- round(timeInt[p], 3)
    Yb1 <- round(varInt[b1], 2)
    Yb2 <- round(varInt[b2], 2)
    Yp <- round(varInt[p], 2)

    ## Check if left and right baseline are similar, if not look for points before or after
    # Left side too small
    if(abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry)
    {
      flag = 1
      pi = 0
      while(abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry & pi < 3)
      {
        pi = pi + 1
        b1 <- zeroes[newi - pi]
        Xb1 <- round(timeInt[b1], 3)
        Yb1 <- round(varInt[b1], 2)

        # For first point
        if(length(b1) == 0){b1 <- zeroes[newi] ; Xb1 <- round(timeInt[b1], 3) ; Yb1 <- round(varInt[b1], 2) ; break}
      }
    }
    # Right side too small
    if(abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry)
    {
      flag = 1
      pi = 0
      while(abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry & pi < 3)
      {
        pi = pi + 1
        b2 <- zeroes[newi + 2 + pi]
        Xb2 <- round(timeInt[b2], 3)
        Yb2 <- round(varInt[b2], 2)

        # For first point
        if(is.na(b2)){b2 <- zeroes[newi + 2] ; Xb2 <- round(timeInt[b2], 3) ; Yb2 <- round(varInt[b2], 2) ; break}
      }
    }

    # If left and right baseline are still not similar, use return to original values and flag
    if(abs(Yp - Yb1) < abs(Yp - Yb2) * asymmetry | abs(Yp - Yb2) < abs(Yp - Yb1) * asymmetry)
    {
      ASYMMETRIC <- TRUE
      if(asmMethod != "rm")
      {
        b1 <- zeroes[newi]
        b2 <- zeroes[newi + 2]
        Xb1 <- round(timeInt[b1], 3)
        Xb2 <- round(timeInt[b2], 3)
        Yb1 <- round(varInt[b1], 2)
        Yb2 <- round(varInt[b2], 2)
      } else {next()}
    } else {
      ASYMMETRIC <- FALSE
    }

    # Skip points where baseline is higher than peak
    if(shadeflecks == FALSE){if(Yb1 > Yp | Yb2 > Yp){next()}} else {if(Yb1 < Yp | Yb2 < Yp){next()}}

    # Initial calculation of fleck data
    if(ASYMMETRIC == TRUE)
    {
      if(asmMethod == "max" & shadeflecks == TRUE){Ybaseline <-  max(c(Yb1, Yb2))}
      if(asmMethod == "max" & shadeflecks == FALSE){Ybaseline <-  min(c(Yb1, Yb2))}
      if(asmMethod == "mean"){Ybaseline <- mean(c(Yb1, Yb2))}
    } else {
      Ybaseline <- mean(c(Yb1, Yb2))
    }

    # Get closest baseline to run condition check on that
    if(shadeflecks == TRUE){YCloseBsl <- min(c(Yb1, Yb2))} else {YCloseBsl <- max(c(Yb1, Yb2))}

    amp <- round(abs(Yp - YCloseBsl), 2)
    pdiff <- round(amp / YCloseBsl, 4)
    duration <- round(abs(Xb2 - Xb1), 3)
    irrad <- sum(varInt[seq(b1,b2,timeSplit)]) * timeStep
    interpBaseline <- stats::approx(x = c(timeInt[b1], timeInt[b2]), y = c(varInt[b1], varInt[b2]), xout = timeInt[seq(b1,b2,timeSplit)])
    irradBaseline <- sum(interpBaseline$y) * timeStep
    irradSunfleck <- irrad - irradBaseline

    # If fleck pass criteria, then trim and recalculate
    if(duration > minTime & amp > minAmp & pdiff > minPdiff)
    {
      # Left trimming
      it = 0 ; trimLeft = 0 ; old_trimLeft = 1
      while(old_trimLeft != trimLeft)
      {
        it = it + 1

        if(b1 + it * timeSplit < length(varInt)){
          Yk1 = round(varInt[b1 + it * timeSplit],2)
          Xk1 = round(timeInt[b1 + it * timeSplit],3)

          old_trimLeft <- trimLeft
          if(stats::sd(c(Yb1,Yk1)) / mean(c(Yb1,Yk1)) < trimCV){trimLeft = trimLeft + 1}

          # If trimming goes too far, abandon trimming
          if(Xk1 == Xp){trimLeft = 0 ; old_trimLeft = 0}
        } else {
          trimLeft = 0 ; old_trimLeft = 0
        }
      }

      # Right trimming
      it = 0 ; trimRight = 0 ; old_trimRight = 1
      while(old_trimRight != trimRight)
      {
        it = it + 1
        if(b2 - it * timeSplit > 0){
          Yk2 = round(varInt[b2 - it * timeSplit],2)
          Xk2 = round(timeInt[b2 - it * timeSplit],3)

          old_trimRight <- trimRight
          if(stats::sd(c(Yb2,Yk2)) / mean(c(Yb2,Yk2)) < trimCV){trimRight = trimRight + 1}

          # If trimming goes too far, abandon trimming
          if(Xk2 == Xp){trimRight = 0 ; old_trimRight = 0}
        } else {
          trimRight = 0 ; old_trimRight = 0
        }
      }

      b1_new <- b1 + trimLeft * timeSplit
      b2_new <- b2 - trimRight * timeSplit

      # Another safeguard for trimming
      if(b2_new > b1_new){b1 <- b1_new ; b2 <- b2_new}

      Xb1 <- round(timeInt[b1], 3)
      Xb2 <- round(timeInt[b2], 3)
      Yb1 <- round(varInt[b1], 2)
      Yb2 <- round(varInt[b2], 2)

      # Final calculation of fleck data
      if(ASYMMETRIC == TRUE)
      {
        if(asmMethod == "max" & shadeflecks == TRUE){Ybaseline <-  max(c(Yb1, Yb2))}
        if(asmMethod == "max" & shadeflecks == FALSE){Ybaseline <-  min(c(Yb1, Yb2))}
        if(asmMethod == "mean"){Ybaseline <- mean(c(Yb1, Yb2))}

      } else {
        Ybaseline <- mean(c(Yb1, Yb2))
      }
      amp <- round(abs(Yp - Ybaseline), 2)
      relAmp <- round(abs(((Yp - bounds[1]) / (bounds[2] - bounds[1])) - ((Ybaseline - bounds[1]) / (bounds[2] - bounds[1]))),4)
      pdiff <- round(amp / Ybaseline, 4)
      duration <- round(abs(Xb2 - Xb1), 3)

      irrad <- sum(varInt[seq(b1,b2,timeSplit)]) * timeStep
      interpBaseline <- stats::approx(x = c(timeInt[b1], timeInt[b2]), y = c(varInt[b1], varInt[b2]), xout = timeInt[seq(b1,b2,timeSplit)])
      irradBaseline <- sum(interpBaseline$y) * timeStep
      irradSunfleck <- irrad - irradBaseline

      relIrrad <- sum((varInt[seq(b1,b2,timeSplit)] - bounds[1]) / (bounds[2] - bounds[1]))
      relBaseline  <- sum((interpBaseline$y - bounds[1]) / (bounds[2] - bounds[1]))
      relIrradSunfleck <- relIrrad - relBaseline

      # Log data (May have to reverse trimming if conditions is not satisfied anymore, will see later)
      no = no + 1
      row <- c(no, Xp, Yp, Yb1, Yb2, Ybaseline, amp, relAmp, pdiff, duration, irradSunfleck, relIrradSunfleck, irrad, Xb1, Xb2, ifelse(ASYMMETRIC==T, "asymmetric", "symmetric"), zeroes[newi+1])

      for(j in 1:length(row))
      {
        l0[[j]][no] <- row[j]
      }

      if(flag == 1 & verbose == TRUE){
        cat("asymmetry found in no", no, "\n")
      }
    }
  }
  df0 <- as.data.frame(l0)
  df0 <- df0[!(is.na(df0[,1])),]
  for(i in 1:15)
  {
    df0[,i] <- as.numeric(df0[,i])
  }

  if(nrow(df0) == 0){stop("No sun/shade-fleck found. Try with different min parameters.")}
  colnames(df0) <-  c("no", "peakTime", "peak", "baseline1", "baseline2", "baseline", "amplitude", "relAmp", "percDiff", "duration", "irradGain", "relIrrad", "irradTot", "baselineTime1", "baselineTime2", "symmetry", "z")

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
    for(iROW in 2:nrow(df0))
    {
      df0$timeInterval[iROW] <- (df0$peakTime[iROW] - df0$peakTime[iROW-1])
    }
  }
  return(df0)
}
