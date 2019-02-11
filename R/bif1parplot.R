bif1parplot <- function(curvelist = NULL, popts) {
  # Plot the bifurcation curves
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
        lapply(2:ncol(curvelist[[i]]$points), function(j) {
          evmax <- Re(curvelist[[i]]$eigvals[,1])
          sp <- (evmax < 0)
          x <- curvelist[[i]]$points[,1]
          y <- curvelist[[i]]$points[,j]
          x[!sp] <- NA
          y[!sp] <- NA
          lines(x, y, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
          x <- curvelist[[i]]$points[,1]
          y <- curvelist[[i]]$points[,j]
          x[sp] <- NA
          y[sp] <- NA
          lines(x, y, lty=popts$unstablelty, col=popts$colors[min(j-1, length(popts$colors))], lwd=popts$lwd)
          if (!is.null(curvelist[[i]]$special.points)) {
            lbls <- c(curvelist[[i]]$special.tags[,1])
            bps <- (lbls %in% c("BP", "HP", "LP"))
            if (any(bps)) {
              x <- curvelist[[i]]$special.points[bps,1]
              y <- curvelist[[i]]$special.points[bps,j]
              points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
              text(x, y, labels=lbls[bps], pos = popts$biflblpos)
            }
          }
        })
      })
      legend("topright", legend=colnames(curvelist[[1]]$points)[2:ncol(curvelist[[1]]$points)], col=popts$colors[1:(ncol(curvelist[[1]]$points)-1)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
    } else {
      lapply((1:length(curvelist)), function(i) {
        evmax <- Re(curvelist[[i]]$eigvals[,1])
        sp <- (evmax < 0)
        x <- curvelist[[i]]$points[,1]
        y <- curvelist[[i]]$points[,popts$ycol]
        x[!sp] <- NA
        y[!sp] <- NA
        lines(x, y, col=popts$colors[1], lwd=popts$lwd)
        x <- curvelist[[i]]$points[,1]
        y <- curvelist[[i]]$points[,popts$ycol]
        x[sp] <- NA
        y[sp] <- NA
        lines(x, y, lty=popts$unstablelty, col=popts$colors[1], lwd=popts$lwd)
        if (!is.null(curvelist[[i]]$special.points)) {
          lbls <- c(curvelist[[i]]$special.tags[,1])
          bps <- (lbls %in% c("BP", "HP", "LP"))
          if (any(bps)) {
            x <- curvelist[[i]]$special.points[bps,1]
            y <- curvelist[[i]]$special.points[bps,popts$ycol]
            points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
            text(x, y, labels=lbls[bps], pos = popts$biflblpos)
          }
        }
      })
      if (popts$y2col > 1) {
        lapply((1:length(curvelist)), function(i) {
          evmax <- Re(curvelist[[i]]$eigvals[,1])
          sp <- (evmax < 0)
          x <- curvelist[[i]]$points[,1]
          y <- converty2y(curvelist[[i]]$points[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
          x[!sp] <- NA
          y[!sp] <- NA
          lines(x, y, col=popts$colors[2], lwd=popts$lwd)
          x <- curvelist[[i]]$points[,1]
          y <- converty2y(curvelist[[i]]$points[,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
          x[sp] <- NA
          y[sp] <- NA
          lines(x, y, lty=popts$unstablelty, col=popts$colors[2], lwd=popts$lwd)
          if (!is.null(curvelist[[i]]$special.points)) {
            lbls <- c(curvelist[[i]]$special.tags[,1])
            bps <- (lbls %in% c("BP", "HP", "LP"))
            if (any(bps)) {
              x <- curvelist[[i]]$special.points[bps,1]
              y <- converty2y(curvelist[[i]]$special.points[bps,popts$y2col], popts$ymin, popts$ymax, popts$logy, popts$y2min, popts$y2max, popts$logy2)
              points(x, y, pch=popts$bifsym, cex=popts$cex.sym, lwd=2)
              text(x, y, labels=lbls[bps], pos = popts$biflblpos)
            }
          }
        })
        legend("topright", legend=colnames(curvelist[[1]]$points)[c(popts$ycol, popts$y2col)], col=popts$colors[c(1, 2)], lty=1, lwd=popts$lwd, cex=popts$sizeLegend)
      }
    }
  }
}
