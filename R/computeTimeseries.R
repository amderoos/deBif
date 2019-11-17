computeTimeseries <- function(session, model, state, parms, clist, pointid, nopts) {

  curvescomputed <- as.numeric(clist[['TotalCurves']])

  if (pointid > 0) {
    ind1 <- round(pointid/1000000)
    ind2 <- round((pointid-ind1*1000000)/1000)
    ind3 <- round(pointid-ind1*1000000-ind2*1000)
    cln1 <- (c('Orbits', 'BifurcationCurves', 'BifurcationBounds'))[ind1]
    ii <- ifelse((ind1 == 3), 2, 1)   # 2 parameter bifurcation points have 2 columns before the state, otherwise 1 only

    initstate <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (ii + (1:length(state)))])
    initparms <- as.numeric(clist[[cln1]][[ind2]]$parameters)
    inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
    if (inittype != "TS")
      initparms[as.numeric(clist[[cln1]][[ind2]]$bifpars)] <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (1:ii)])
    names(initstate) <- names(state)
    names(initparms) <- names(parms)
  } else {
    initstate <- state
    initparms <- parms
    inittype <- "US"
  }

  newcurvenr <- length((clist[['Orbits']]))+1

  times <- seq(0, nopts$tmax, by=abs(nopts$tstep))
  if (nopts$tstep < 0.0) times <- nopts$tmax - times
  nsol <- as.data.frame(do.call('ode', c(list(times=times, func=model, y=initstate, parms=initparms), method=nopts$odemethod)))

  names(nsol) <- c("Time", names(state))

  startPnt <- c("Type" = "TS",
                "Description" = paste(paste0('T=', times[1]),
                                      unlist(lapply(1:length(state),
                                                    function(i) {paste0(names(state[i]), "=", round(nsol[1, (1+i)], 3))})),
                                      collapse=', '))
  endPnt <- c("Type" = "TS",
              "Description" = paste(paste0('T=', times[length(times)]),
                                    unlist(lapply(1:length(state),
                                                  function(i) {paste0(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 3))})),
                                    collapse=', '))

  updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
  curvescomputed <- curvescomputed + 1

  lbl <- paste0("TS", sprintf("%02d", curvescomputed),": ", startPnt["Description"])
  startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])
  endPnt["Description"] <- paste0(sprintf("%04d: ", nrow(nsol)), endPnt["Description"])

  newcurve <- list(label = lbl, type = "TS", initstate = initstate, parameters = initparms, points = nsol,
                   special.points = rbind(c(nsol[1,]), c(nsol[nrow(nsol),])), special.tags = rbind(startPnt, endPnt))

  clist$Orbits[[newcurvenr]] <- newcurve
  clist$TotalCurves <- curvescomputed

  return(clist)
}
