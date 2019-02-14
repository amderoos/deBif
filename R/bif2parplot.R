bif2parplot <- function(curvelist = NULL, output, session, popts) {
  # Plot the bifurcation curves

  par(mar = as.numeric(popts$plotmar))
  logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
  plot(1, 1, type='n', xlab=popts$xlab, ylab=popts$ylab, xaxs = "i", yaxs = "i",
       xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy,
       cex.main=popts$cex.main, cex.lab=popts$cex.lab, cex.axis=popts$cex.axis,
       font.main=popts$font.main, font.sub=popts$font.sub)

  if (!is.null(curvelist) && (length(curvelist) > 0)) {
    lapply((1:length(curvelist)), function(i) {
      cnames <- colnames(curvelist[[i]]$points)
      if (cnames[1] != popts$xlab) {
        msg <- paste0("Curve plotting skipped: parameter '", popts$xlab, "' not one of the curve variables\n")
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        return(NA)
      }
      if (cnames[2] != popts$ylab) {
        msg <- paste0("Curve plotting skipped: parameter '", popts$ylab, "' not one of the curve variables\n")
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        return(NA)
      }
      colindx <- match(curvelist[[i]]$type, c("BP", "HP", "LP"))
      lines(curvelist[[i]]$points[,1], curvelist[[i]]$points[,2], col=popts$colors[colindx], lwd=popts$lwd)

      if (!is.null(curvelist[[i]]$special.points)) {
        lbls <- c(curvelist[[i]]$special.tags[,1])
        bps <- (lbls %in% c("BT", "CP"))
        if (any(bps)) {
          x <- curvelist[[i]]$special.points[bps,1]
          y <- curvelist[[i]]$special.points[bps,2]
          points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
          text(x, y, labels=lbls[bps], pos = popts$biflblpos)
        }
      }
    })
    legend("topright", legend=c("BP", "HP", "LP"), col=popts$colors[c(1, 2, 3)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
  }
}
