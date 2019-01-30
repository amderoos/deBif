converty2y <- function(y, ymin1, ymax1, logy1, ymin2, ymax2, logy2) {
  # Convert values on the second y axis to corresponding values on the first y-axis
  if (logy2 == 0) {
    ytmp <- pmax((y - ymin2)/(ymax2 - ymin2), -0.5)
  } else {
    ytmp <- pmax(log(pmax(y, 1.0E-15)/ymin2)/(log(ymax2/ymin2)), -0.5)
  }

  if (logy1 == 0) {
    yval <- ymin1 + ytmp*(ymax1 - ymin1)
  } else {
    yval <- exp(log(ymin1) + ytmp*(log(ymax1/ymin1)))
  }
  return(yval)
}
