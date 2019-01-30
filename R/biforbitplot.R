biforbitplot <- function(curvelist = NULL, popts) {
  # Plot the computed time series
  if ((popts$ycol > 1) && (popts$y2col > 1)) {
    par(mar = (as.numeric(popts$plotmar) + c(0, 0, 0, 1)))
    logxy <- ifelse(popts$logx == 1, ifelse(popts$logy2 == 1, "xy", "x"), ifelse(popts$logy2 == 1, "y", ""))
    plot(1, 1, type='n', xlab="", ylab="", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
         xlim=c(popts$xmin,popts$xmax), ylim=c(popts$y2min,popts$y2max), log=logxy)
    axis(4, cex.axis=popts$cex.axis)
    mtext(popts$y2lab, side = 4, line = 3, cex=popts$cex.lab)
    par(new = TRUE)
  } else {
    par(mar = as.numeric(popts$plotmar))
  }

  logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
  plot(1, 1, type='n', xlab=popts$xlab, ylab=popts$ylab, xaxs = "i", yaxs = "i",
       xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy,
       cex.main=popts$cex.main, cex.lab=popts$cex.lab, cex.axis=popts$cex.axis,
       font.main=popts$font.main, font.sub=popts$font.sub)

  if (!is.null(curvelist) && (length(curvelist) > 0)) {
    if (popts$ycol == 1) {
      lapply((1:length(curvelist)), function(i) {
        lapply(2:ncol(curvelist[[i]]$curve), function(j) {
          lines(curvelist[[i]]$curve[,1], curvelist[[i]]$curve[,j], col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
        })
      })
      legend("topright", legend=names(curvelist[[1]]$curve)[2:ncol(curvelist[[1]]$curve)], col=popts$colors[1:(ncol(curvelist[[1]]$curve)-1)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
    } else {
      lapply((1:length(curvelist)), function(i) {
        lines(curvelist[[i]]$curve[,1], curvelist[[i]]$curve[,popts$ycol], col=popts$colors[1], lwd=popts$lwd)
      })
      if (popts$y2col > 1) {
        lapply((1:length(curvelist)), function(i) {
          lines(curvelist[[i]]$curve[,1],
                converty2y(curvelist[[i]]$curve[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2),
                col=popts$colors[2], lwd=popts$lwd)
        })
        legend("topright", legend=names(curvelist[[1]]$curve)[c(popts$ycol, popts$y2col)], col=popts$colors[c(1, 2)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
      }
    }
  }
}
