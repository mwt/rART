clm.fits <- function(cmfi, mt, contrasts, singular.ok, ...) {
  # get x,y
  cxi <- model.matrix(mt, cmfi, contrasts)
  cyi <- model.response(cmfi, "numeric")
  # get offset
  coffseti <- model.offset(cmfi)
  if(!is.null(offset)) {
    coffseti <- as.vector(coffseti)
  }
  czi <- lm.fit(cxi, cyi, offset = coffseti, singular.ok=singular.ok, ...)
  czi$coefficients
}

clm.wfits <- function(cmfi, mt, contrasts, singular.ok, ...) {
  # get x,y,w
  cxi <- model.matrix(mt, cmfi, contrasts)
  cyi <- model.response(cmfi, "numeric")
  cwi <- as.vector(model.weights(cmfi))
  # get offset
  coffseti <- model.offset(cmfi)
  if(!is.null(offset)) {
    coffseti <- as.vector(coffseti)
  }
  czi <- lm.wfit(cxi, cyi, cwi, offset = coffseti, singular.ok=singular.ok, ...)
  czi$coefficients
}
