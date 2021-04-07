clm.fits <- function(cmfi, mt, contrasts, singular.ok, ...) {
  cxi <- model.matrix(mt, cmfi, contrasts)
  cyi <- model.response(cmfi, "numeric")
  coffseti <- model.offset(cmfi)
  if(!is.null(offset)) {
    coffseti <- as.vector(coffseti)
  }
  czi <- lm.fit(cxi, cyi, offset = coffseti, singular.ok=singular.ok, ...)
  czi$coefficients
}
