#' ART Linear Model
#'
#' @param formula an object of class [formula] (or one that
#'   can be coerced to that class): a symbolic description of the
#'   model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible
#'   by [as.data.frame] to a data frame) containing the variables in the model.
#'   If not found in `data`, the variables are taken from `environment(formula)`
#'   , typically the environment from which [lm()] is called.
#' @param cluster the name for the cluster variable contained in `data` or
#'   vector of clusters
#' @param select an optional string or list of strings for the parameter(s) of
#'   interest. That is, the parameters you want to compute tests and CI for.
#'   Parameters as passed to `subset`.
#' @param subset an optional vector specifying a subset of observations
#'   to be used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'   process.  Should be `NULL` or a numeric vector.
#'   If non-NULL, weighted least squares is used with weights
#'   `weights` (that is, minimizing `sum(w*e^2)`); otherwise
#'   ordinary least squares is used.  See also \sQuote{Details},
#' @param na.action a function which indicates what should happen
#'   when the data contain `NA`s.  The default is set by
#'   the `na.action` setting of [options], and is
#'   [na.fail] if that is unset.  The \sQuote{factory-fresh}
#'   default is [na.omit].  Another possible value is
#'   `NULL`, no action.  Value [na.exclude] can be useful.
#' @param method the method to be used; for fitting, currently only
#'   `method = "qr"` is supported; `method = "model.frame"` returns the model
#'   frame (the same as with `model = TRUE`, see below).
#' @param singular.ok logical. If `FALSE` (the default in S but
#'   not in \R) a singular fit is an error.
#' @param contrasts an optional list. See the `contrasts.arg`
#'   [model.matrix.default].
#' @param offset this can be used to specify an *a priori* known
#'   component to be included in the linear predictor during fitting.
#'   This should be `NULL` or a numeric vector or matrix of extents
#'   matching those of the response.  One or more [offset] terms can
#'   be included in the formula instead or as well, and if more than one are
#'   specified their sum is used.  See [model.offset].
#' @param ... additional arguments to be passed to the low level
#'   regression fitting functions (see below).
#'
#' @return an object of [class] `c("artlm", "lm")`.
#' @export
artlm <- function (formula, data, cluster, select = NULL, subset, weights,
                na.action, method = "qr", model = TRUE, x = FALSE,
                y = FALSE, qr = TRUE, singular.ok = TRUE,
                contrasts = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights",
               "na.action", "offset", "cluster"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
            domain = NA)
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- model.offset(mf)
  if(is.matrix(y))
    stop("'y' is a matrix but a vector is required")
  ny <- length(y)
  if(!is.null(offset)) {
    offset <- as.vector(offset)
    if(NROW(offset) != ny)
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    NROW(offset), ny), domain = NA)
  }

  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = numeric(),
              residuals = y,
              fitted.values = 0 * y, weights = w, rank = 0L,
              df.residual = if(!is.null(w)) sum(w != 0) else ny)
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    ## make a list of model frames by cluster
    cmfs <- split(mf, model.extract(mf, "cluster"))
    if(is.null(w)) {
      z <- lm.fit(x, y, offset = offset, singular.ok=singular.ok, ...)
      clbetas <- vapply(cmfs, clm.fits, z$coefficients, mt = mt,
                        contrasts = contrasts, singular.ok=singular.ok, ...)
    } else {
      z <- lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
      clbetas <- vapply(cmfs, clm.wfits, z$coefficients, mt = mt,
                        contrasts = contrasts, singular.ok=singular.ok, ...)
    }
  }
  class(z) <- c("artlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr) z$qr <- NULL
  ## record the number of clusters
  z$ncluster <- length(cmfs)
  ## use all betas if none selected
  if(is.null(select)) {
    z$clbetas <- clbetas
  }
  else {
    # make unselected betas NA
    z$clbetas <- clbetas * NA
    z$clbetas[select,] <- clbetas[select,]
    z$select <- select
  }
  z
}

#' Dummy vcov for ART
#'
#' ART does not generate standard errors or a variance covariance matrix. So,
#' this function will always return `NA`.
#'
#' @param ... variables that are passed to [vcov.lm()]
#'
#' @return A matrix of `NA` values in the same shape as a variance
#'   covariance matrix.
#' @export
vcov.artlm <- function(...) {
  warning("ART does not provide variances")
  stats:::vcov.lm(...)*NA
}

#' Summarizing Linear Models with ART
#'
#' @param object an object of class `"artlm"`, usually, a result of a
#'   call to [artlm()].
#' @param nrg an `integer` that represents the number of random
#'   permutations to calculate when there are more than 10 groups. Defaults
#'   to 10,000.
#' @param ... Other parameters that can be passed to the base
#'   [summary.lm()].
#'
#' @return A summary object as in [summary.lm()].
#' @export
summary.artlm <- function(object, nrg = 10000, ...) {
  raw_lmsum <- summary.lm(object, ...)
  clbetas <- object$clbetas
  q <- object$ncluster
  raw_lmsum$Gs <- random.G(q = q, B = nrg)
  test <- apply(clbetas, 1, CRS.test, G = raw_lmsum$Gs)
  raw_lmsum$coefficients <- na.omit(
    cbind("Estimate" = object$coefficients, t(test)))
  raw_lmsum$aliased <- NULL # if coef are aliased, then it displays wrong
  raw_lmsum
}

#' Confidence Intervals for Linear ART Parameters
#'
#' @param an object of class `"artlm"`, usually, a result of a call to
#'   [artlm()].
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param nrg an `integer` that represents the number of random permutations to
#'   calculate when there are more than 10 groups. Defaults to 10,000.
#' @param level the confidence level required. Defaults to 95%
#'
#' @return  A matrix (or vector) with columns giving lower and upper confidence
#'   limits for each parameter. These will be labeled as (1-level)/2 and
#'   1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#' @export
confint.artlm <- function(object, parm, nrg = 10000, level = 0.95) {
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames
  else if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  q <- object$ncluster
  clbetas <- object$clbetas[parm,]
  G <- random.G(q = q, B = nrg)
  # if there is only one variable, then don't use apply
  if (length(parm) == 1L) {
    ci[] <- t(CRS.CI(clbetas, G = G, alpha = (1 - level)))
  } else {
    ci[] <- t(apply(clbetas, 1, CRS.CI, G = G, alpha = (1 - level)))
  }
  # remove NAs
  ci[complete.cases(ci),]
}
