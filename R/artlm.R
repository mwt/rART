#' ART Linear Model
#'
#' @param formula an object of class \code{"\link{formula}"} (or one that
#'   can be coerced to that class): a symbolic description of the
#'   model to be fitted.
#' @param data an optional data frame, list or environment (or object
#'   coercible by \code{\link{as.data.frame}} to a data frame) containing
#'   the variables in the model.  If not found in \code{data}, the
#'   variables are taken from \code{environment(formula)},
#'   typically the environment from which \code{lm} is called.
#' @param subset an optional vector specifying a subset of observations
#'   to be used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'   process.  Should be \code{NULL} or a numeric vector.
#'   If non-NULL, weighted least squares is used with weights
#'   \code{weights} (that is, minimizing \code{sum(w*e^2)}); otherwise
#'   ordinary least squares is used.  See also \sQuote{Details},
#' @param na.action a function which indicates what should happen
#'   when the data contain \code{NA}s.  The default is set by
#'   the \code{na.action} setting of \code{\link{options}}, and is
#'   \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
#'   default is \code{\link{na.omit}}.  Another possible value is
#'   \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.
#' @param method the method to be used; for fitting, currently only
#'   \code{method = "qr"} is supported; \code{method = "model.frame"} returns
#'   the model frame (the same as with \code{model = TRUE}, see below).
#' @param model,x,y,qr logicals. If \code{TRUE} the corresponding
#'   components of the fit (the model frame, the model matrix, the
#'   response, the QR decomposition) are returned.
#' @param singular.ok logical. If \code{FALSE} (the default in S but
#'   not in \R) a singular fit is an error.
#' @param contrasts an optional list. See the \code{contrasts.arg}
#'   \code{\link{model.matrix.default}}.
#' @param offset this can be used to specify an \emph{a priori} known
#'   component to be included in the linear predictor during fitting.
#'   This should be \code{NULL} or a numeric vector or matrix of extents
#'   matching those of the response.  One or more \code{\link{offset}} terms can be
#'   included in the formula instead or as well, and if more than one are
#'   specified their sum is used.  See \code{\link{model.offset}}.
#' @param ... additional arguments to be passed to the low level
#'   regression fitting functions (see below).
#'
#' @return \code{lm} returns an object of \code{\link{class}} \code{"lm"} or for
#'   multiple responses of class \code{c("mlm", "lm")}.
#' @export
artlm <- function (formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE,
                qr = TRUE, singular.ok = TRUE, contrasts = NULL,
                offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
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
  mlm <- is.matrix(y)
  ny <- if(mlm) nrow(y) else length(y)
  if(!is.null(offset)) {
    if(!mlm) offset <- as.vector(offset)
    if(NROW(offset) != ny)
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    NROW(offset), ny), domain = NA)
  }

  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if(mlm) matrix(NA_real_, 0, ncol(y))
              else numeric(),
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
    z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                               singular.ok=singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
  }
  class(z) <- c(if(mlm) "mlm", "lm")
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
  z
}
