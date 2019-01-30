bif1parplot <- function(curvelist = NULL, popts) {
  # Plot the bifurcation curves
  par(mar = popts$plotmar)
  logxy <- ifelse(popts$logx == 1, ifelse(popts$logy == 1, "xy", "x"), ifelse(popts$logy == 1, "y", ""))
  plot(1, 1, type='n', xlab=popts$xlab, ylab=popts$ylab,
       xlim=c(popts$xmin,popts$xmax), ylim=c(popts$ymin,popts$ymax), log=logxy,
       cex.main=popts$cex.main, cex.lab=popts$cex.lab, cex.axis=popts$cex.axis,
       font.main=popts$font.main, font.sub=popts$font.sub)

  # if (!is.null(curvelist) && (length(curvelist) > 0)) {
  #   lapply((1:length(curvelist)), function(i) {
  #     lapply(2:ncol(curvelist[[i]]$curve), function(j) {
  #       lines(curvelist[[i]]$curve[,1], curvelist[[i]]$curve[,j], col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
  #     })
  #   })
  #   legend("topright", legend=names(curvelist[[1]]$curve)[2:ncol(curvelist[[1]]$curve)], col=popts$colors[1:(ncol(curvelist[[1]]$curve)-1)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
  # }
}
