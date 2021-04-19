artiv <- function (formula, data, cluster, subset, weights,
                   na.action, model = TRUE, singular.ok = TRUE,
                   contrasts = NULL, offset, ...)
{
  # from lm()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights",
               "na.action", "offset", "cluster"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)


  # from ivreg()
  ## we use Formula package to extract IV
  formula <- Formula::as.Formula(formula)
  mf$formula <- formula
  mf <- eval(mf, parent.frame())

  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf, contrasts)
  if(length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- NULL
  } else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
  }

  # back to regular lm() stuff
  ## get weights
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")

}
