bif2parplot <- function(curvelist = NULL, popts) {
  # Plot the bifurcation curves

  par(mar = as.numeric(popts$plotmar))
  logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
  plot(1, 1, type='n', xlab=popts$xlab, ylab=popts$ylab, xaxs = "i", yaxs = "i",
       xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy,
       cex.main=popts$cex.main, cex.lab=popts$cex.lab, cex.axis=popts$cex.axis,
       font.main=popts$font.main, font.sub=popts$font.sub)

  if (!is.null(curvelist) && (length(curvelist) > 0)) {
    lapply((1:length(curvelist)), function(i) {
      colindx <- match(curvelist[[i]]$type, c("BP", "HP", "LP"))
      lines(curvelist[[i]]$points[,1], curvelist[[i]]$points[,2], col=popts$colors[colindx], lwd=popts$lwd)
    })
    legend("topright", legend=c("BP", "HP", "LP"), col=popts$colors[c(1, 2, 3)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
  }
}
