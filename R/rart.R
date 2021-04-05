# Author: Ivan Canay | Northwestern University | November 2020
# Copyright (c) 2020 Northwestern University. All rights reserved.
#-------------------------------------------------------------------
# X-plain:  This file has several function in R that are used in the
#           paper ``A Guide To Practice for Approximate Randomization
#           Tests in Regressions with a Small Number of Clusters''
#-------------------------------------------------------------------
# List of functions included in this file:
#  Functions for the CRS test and confidence interval
#   (1) random.G  : Get groups of transformations
#   (2) CRS.test  : Computes randomization test for t-test
#   (3) CRS.CI    : Computes confidence intervals by test inversion using bisection
#  Functions for computing Clustered Standard Errors (used in applications)
#   (3) cl        : Computes clustered standard errors (BCH and Stata)
#   (4) cce.stata : Computes clustered se (BCH and Stata). Faster than cl
#   (5) cce.brl   : Computes clustered se with Bias Reduced Linearization
#-------------------------------------------------------------------

#------------------------------------------------------------------------
#--- FUNCTIONS FOR THE RANDOMIZATION TEST -------------------------------
#------------------------------------------------------------------------

#-------------------------------------------------------------------
pre.random.G <- function(q, B = 10000) {
  #-------------------------------------------------------------------
  # X-plain: Computes the group G of transformations for the CRS test
  #-------------------------------------------------------------------
  # INPUTS: - q: the dimension of the data / # of clusters (q point estimators)
  #         - B: the number of random draws if q>10.
  #------------------------------------------------------------------
  # RETURNS: - G: a B \times q matrix.
  #------------------------------------------------------------------
  if (q < 11) {
    e <- c(1, -1);
    B <- 2 ^ q;
    perm <- t(expand.grid(rep(list(e), q)));
  }
  # Draw B transformations from G otherwise
  else {
    perm <- matrix(2 * ((runif(q * B) > 1 / 2) - 1 / 2), q, B);
  }
  return(perm)
}
#-------------------------------------------------------------------

#-------------------------------------------------------------------
pre.CRS.test <- function(c.beta, G, lambda = 0, alpha = 0.05, nj = 1) {
  #-------------------------------------------------------------------
  # X-plain: Computes randomization critical value
  #-------------------------------------------------------------------
  # INPUTS:  - c.beta: vector of parameters entering the null hypothesis
  #                    c'beta - with one estimator per cluster (q x 1).
  #          - G: the group of all trasformations (use random.G to get)
  #          - lambda : scalar for the null hypothesis c'beta = lambda (default 0)
  #          - alpha: significance level (dafault 0.05)
  #          - nj : qx1 vector of sample sizes in each cluster (for alternative weighting)
  #------------------------------------------------------------------
  # RETURNS: - rule: binary decision of the test
  #          - Nrule: binary decision of the non-randomized test
  #          - cv: the critival value of the randomization test.
  #          - pv and pv2: p-values according to folmulae 15.5 and 15.7
  #------------------------------------------------------------------
  q = length(c.beta);
  # # Number of clusters/estimators
  M = dim(G)[2];
  # # of elements in G
  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  Sn = sqrt(q) * sqrt(nj) * (c.beta - lambda);

  ObsT = abs(mean(Sn)) / sd(Sn);
  # observed Test stat
  NewX = G * as.vector(Sn);
  # transformed data
  # Compute New Test Stat over transformed data
  NewT = abs(apply(NewX, 2, mean) / apply(NewX, 2, sd));
  NewT = sort(NewT);
  # sort the vector of Test Stat
  k = M - floor(M * alpha);
  # find index for quantile
  Mplus = sum(NewT > NewT[k]);
  # M+ from LR book
  M0 = sum(NewT == NewT[k]);
  # M0 from LR book
  prob = (M * alpha - Mplus) / M0;
  # term a(x) from LR book

  # randomized Test
  rule = (ObsT > NewT[k]) + (ObsT == NewT[k]) * (runif(1) <= prob);
  # Non non-randomized Test
  Nrule = (ObsT > NewT[k]);
  # Compute the p-value (formula (15.5) from LR book
  p.value = sum(NewT >= ObsT) / M;
  # Compute the p-value (formula (15.7) from LR book
  p.value2 = (sum(NewT >= ObsT) + 1) / (M + 1);

  # List of returns
  list(rule = rule, Nrule = Nrule, cv = NewT[k], pv = p.value, pv2 = p.value2)
}
#-------------------------------------------------------------------

#-------------------------------------------------------------------
CRS.CI <- function(c.beta, G, alpha = 0.05, nj = 1) {
  #-------------------------------------------------------------------
  # X-plain: Computes confidence intervals by test inversion using bisection
  #-------------------------------------------------------------------
  # INPUTS:  - c.beta: vector of parameters entering the null hypothesis
  #                    c'beta - with one estimator per cluster (q x 1).
  #          - G: the group of all trasformations (use random.G to get)
  #          - alpha: significance level (dafault 0.05)
  #          - nj : qx1 vector of sample sizes in each cluster (for alternative weighting)
  #-------------------------------------------------------------------
  # RETURNS: - Cn : 1-alpha confidence interval for c'beta
  #-------------------------------------------------------------------
  max_ite = 100
  tolerance = mean(c.beta) / 1000;
  #---------------------------------------------------------------
  q = length(c.beta);
  if (length(nj) != 1 & length(nj) != q) { nj = 1 }

  # Random variable Sn as in Algorithm 2.1 plus initial values
  center = mean(c.beta);
  distance = 2 * sd(c.beta);
  L = center - distance;
  U = center + distance;
  #---------------------------------------------------------------
  # Find lower an upper bounds that are ``rejected''
  lo.test <- CRS.test(c.beta, G, L, alpha);
  ite = 1;
  while (lo.test$Nrule == 0 & ite < 10) {
    L = L - distance;
    lo.test <- CRS.test(c.beta, G, L, alpha, nj);
    ite = ite + 9;
    if (ite == 10) { stop("Could not find proper lower bound") }
  }

  hi.test <- CRS.test(c.beta, G, U, alpha);
  ite = 1;
  while (hi.test$Nrule == 0 & ite < 10) {
    U = U - distance;
    hi.test <- CRS.test(c.beta, G, U, alpha, nj);
    ite = ite + 1;
    if (ite == 10) { stop("Could not find proper upper bound") }
  }

  ite = 1;
  Rej.L = L;
  Acc.L = center;

  while (ite < max_ite) {
    # Check that the numbers are not too close
    if (abs(Acc.L - Rej.L) / 2 < tolerance) { new.L = (Rej.L + Acc.L) / 2; break; }
    # compute new mid point
    new.L = (Rej.L + Acc.L) / 2;
    lo.test <- CRS.test(c.beta, G, new.L, alpha, nj);
    if (lo.test$Nrule == 0) { Acc.L = new.L; } else { Rej.L = new.L; }
    ite = ite + 1;
  }
  print(paste("iterations to find lower bound: ", ite))

  ite = 1;
  Rej.U = U;
  Acc.U = center;

  while (ite < max_ite) {
    # Check that the numbers are not too close
    if (abs(Acc.U - Rej.U) / 2 < tolerance) { new.U = (Rej.U + Acc.U) / 2; break; }
    # compute new mid point
    new.U = (Rej.U + Acc.U) / 2;
    hi.test <- CRS.test(c.beta, G, new.U, alpha, nj);
    if (hi.test$Nrule == 0) { Acc.U = new.U; } else { Rej.U = new.U; }
    ite = ite + 1;
  }
  print(paste("iterations to find upper bound: ", ite))
  Cn = c(new.L, new.U)

  return(Cn)
}

#-------------------------------------------------------------------
#--------- FUNCTIONS FOR Clustered Standard Errors  ----------------
#-------------------------------------------------------------------
pre.cl <- function(dat, reg, cluster, dof = "bch") {
  #-------------------------------------------------------------------
  # X-plain: Computes clustered standard errors for a regression coeff.
  #          Returns estimate of asymptotic var-covariance matrix / N.
  # See http://people.su.se/~ma/clustering.pdf for references.
  #-------------------------------------------------------------------
  # INPUTS:  - dat: data in (Y,X) format.
  #          - reg: the group of all trasformations (from random.G)
  #          - cluster: clustering variable (integers)
  #          - dof: degrees of freedom correction. "bch" (default). Also, "stata".
  #------------------------------------------------------------------
  # WARNING: This function is slower than cce.stata. Use cce.stata.
  #------------------------------------------------------------------
  # RETURNS: - cl: Estimate of asymptotic var-cov matrix divided by N
  #------------------------------------------------------------------
  M <- length(unique(cluster));
  N <- length(cluster);
  L <- reg$rank;
  # DFC as in BCH11 pg. 142. or as on p.16 in NBER WP 18478 by Imbens and Kolesar
  dfc <- (M / (M - 1)) * (dof == "bch") + (N / (N - L)) * (dof == "HC1") + ((N - 1) / (N - L)) * (M / (M - 1)) * (dof == "stata");
  uj <- apply(estfun(reg), 2, function(x) tapply(x, cluster, sum));
  CCE <- dfc * sandwich(reg, meat = crossprod(uj) / N);
  return(CCE);
  # Comment on this function. This function returns the same as the
  # following two steps.
  # 1) run plm using only "firms" as clusters and pooling.
  # > pool <- plm(y ~ x, data=data, model='pooling', index=c('firmid'))
  # then use the plm class "pool" with fuction cvovHC as follows.
  # dfc*vcovHC(pool, type = "HC0", cluster = "group", adjust="T")
}
#-------------------------------------------------------------------

#------------------------------------------------------------------
pre.cce.stata <- function(X, reg, cluster, dof = "stata") {
  # Computes the robust cluster variance estimator implemented in Stata
  # ... as on p.16 in NBER WP 18478 by Imbens and Kolesar
  #------------------------------------------------------------------
  # INPUTS:  - X:       N-by-L matrix of regressors
  #          - reg:     outcome of a model fit (to obtain residuals and for sandwich)
  #          - cluster: N-by-1 vector denoting cluster membership of observations
  #                     clusters can be denoted by anything, not just consecutive numbers
  #          - dof:     degrees of freedom correction: "stata" default. Also "bch".
  #------------------------------------------------------------------
  # WARNING: if constant is included in the regression 'reg', it must be supplied
  #          ... as part of X
  #------------------------------------------------------------------
  # RETURNS: - Estimate of asymptotic var-cov matrix divided by N
  #------------------------------------------------------------------
  N <- length(cluster) # Number of observations
  S <- length(unique(cluster)) # Number of clusters
  L <- dim(X)[2] # Number of regressors
  ## In R, a vector is dimensionless. As such, X[1,] returns an error.
  ## Use the following to introduce consistency with matrix X references.
  if (is.null(L)) {
    L <- 1
    X <- as.matrix(X)
  }
  res <- as.vector(residuals(reg)) # Model residuals

  ## Construct meat matrix
  mymeat <- matrix(rep(0, L ^ 2), L, L)
  for (s in unique(cluster)) {
    tag <- cluster == s
    if (sum(tag) == 1) {
      Xs <- t(X[tag,])
    } else {
      Xs <- X[tag,]
    }
    mymeat <- mymeat + crossprod(crossprod(res[tag], Xs))
  }
  mymeat <- mymeat / N
  # DFC as in BCH11 pg. 142. or as on p.16 in NBER WP 18478 by Imbens and Kolesar
  dfc <- (S / (S - 1)) * (dof == "bch") + (N / (N - L)) * (dof == "HC1") + ((N - 1) / (N - L)) * (S / (S - 1)) * (dof == "stata");

  Vstata <- dfc * sandwich(reg, meat. = mymeat)
  return(Vstata)
}
#-------------------------------------------------------------------

#-------------------------------------------------------------------
pre.cce.brl <- function(X, reg, cluster, method = "eigen") {
  # Computes the robust cluster variance estimator of Bell and McCafrey (2002)
  # ... as on pp.16-17 in NBER WP 18478 by Imbens and Kolesar (\hat{V}_{lz2})
  #----------------------------------------------------------------------------------
  # INPUTS:  - X:       N-by-L matrix of regressors
  #          - reg:     outcome of a model fit (to obtain residuals and for sandwich)
  #          - cluster: N-by-1 vector denoting cluster membership of observations
  #                        clusters can be denoted by anything, not just consecutive numbers
  #          - method: specifies how the square root matrix is obtained
  #            . "eigen" obtains symmetric square root via eigendecomposition. Slow.
  #            . anything else uses Cholesky decomposition which might be faster, but
  #              ... is not symmetric
  #----------------------------------------------------------------------------------
  # RETURNS: list with the following values:
  #          - Vlz2:   L-by-L estimator of the variance-covariance matrix (div. by N)
  #          - Kbm:    L-by-1 BM degrees-of-freedom correction for all coefficients
  #          - Kik:    L-by-1 IK degrees-of-freedom correction for all coefficients
  #----------------------------------------------------------------------------------
  # WARNING: if constant is included in the regression 'reg', it must be supplied
  #          ... as part of X
  # WARNING: due to the need to get eigenvectors to construct symmetric square root,
  #          ... this function gets really slow with large number of observations
  #          ... when method=="eigen". "chol" is faster but not that fast either.
  # WARNING: results can vary greatly between "eigen" and "chol"
  #----------------------------------------------------------------------------------
  N <- length(cluster) # Number of observations
  S <- length(unique(cluster)) # Number of clusters
  L <- dim(X)[2] # Number of regressors
  ## In R, a vector is dimensionless. As such, X[1,] returns an error.
  ## Use the following to introduce consistency with matrix X references.
  if (is.null(L)) {
    L <- 1
    X <- as.matrix(X)
  }
  res <- as.vector(residuals(reg)) # Model residuals

  XXinv <- solve(crossprod(X))
  M <- diag(N) - tcrossprod(X, tcrossprod(X, XXinv))

  ## Construct meat matrix and degrees-of-freedom corrections
  mymeat <- matrix(rep(0, L ^ 2), L, L)
  G <- array(dim = c(L, N, S))
  Omega <- matrix(rep(0, N ^ 2), N, N)
  crossressum <- 0 # Stores the sum of same-cluster residual cross-products
  numcrossres <- 0 # Stores the number of those cross-products (to get mean)
  i <- 0
  for (s in unique(cluster)) {
    i <- i + 1
    tag <- cluster == s
    Ns <- sum(tag)
    if (sum(tag) == 1) {
      Xs <- t(X[tag,])
    } else {
      Xs <- X[tag,]
    }
    Mss <- diag(Ns) - tcrossprod(Xs, tcrossprod(Xs, XXinv))
    restag <- res[tag]

    ## Compute relevant cross-products to produce \Omega as on p.17 (for Kik)
    crossterms <- tcrossprod(restag)
    diag(crossterms) <- NA
    crossressum <- crossressum + sum(crossterms, na.rm = TRUE)
    numcrossres <- numcrossres + Ns ^ 2 - Ns
    Omega <- Omega + tcrossprod(tag)

    if (method == "eigen") {
      ## Eigendecompose Mss and construct symmetric square root
      r <- eigen(Mss, symmetric = TRUE)
      v <- r$vectors
      if (sum(tag) == 1) {
        rootMss <- tcrossprod(v, tcrossprod(v, sqrt(r$values)))
      } else {
        rootMss <- tcrossprod(v, tcrossprod(v, sqrt(diag(r$values))))
      }
      negrootMss <- solve(rootMss)
    } else {
      ## Do Cholesky decomposition to construct square root
      negrootMss <- chol(solve(Mss), pivot = FALSE)
    }

    projres <- crossprod(t(negrootMss), restag)
    mymeat <- mymeat + crossprod(crossprod(projres, Xs))

    ## Construct G matrices for all coefficients to get dof corrections later

    if (sum(tag) == 1) {
      Mtag <- t(M[tag,])
    } else {
      Mtag <- M[tag,]
    }
    protoG <- tcrossprod(crossprod(Mtag, negrootMss), tcrossprod(XXinv, Xs))
    for (k in seq(1, L)) {
      G[k,, i] <- protoG[, k]
    }
    # Close coefficient loop
  }
  # Close cluster loop

  mymeat <- mymeat / N
  Vlz2 <- sandwich(reg, meat. = mymeat)

  ## Construct estimator of E[ee'|X] as suggested on p.17
  sigmadiag2 <- mean(res ^ 2)
  sigmanu2 <- crossressum / numcrossres
  Omega <- sigmanu2 * Omega
  diag(Omega) <- sigmadiag2

  ## Degrees-of-freedom corrections
  Kbm <- vector(length = L)
  Kik <- vector(length = L)
  for (k in seq(1, L)) {
    GG <- crossprod(G[k,,])
    GomG <- crossprod(G[k,,], Omega %*% G[k,,])
    r <- eigen(GG, symmetric = TRUE, only.values = TRUE)
    Kbm[k] <- (sum(r$values) ^ 2) / sum(r$values ^ 2)
    if (sum(tag) != 1) {
      r <- eigen(GomG, symmetric = TRUE, only.values = TRUE)
      Kik[k] <- (sum(r$values) ^ 2) / sum(r$values ^ 2)
    }
  }

  return(list(Vlz2 = Vlz2, Kbm = Kbm, Kik = Kik))
}
#-------------------------------------------------------------------

#-------------------------------------------------------------------
pre.InvSqrt <- function(M) {
  #-------------------------------------------------------------------
  # Computes the inverse of the symmetric square root of M
  #----------------------------------------------------------------------------------
  # INPUTS: - M: An NxN symmetric pos. def. matrix.
  #----------------------------------------------------------------------------------
  ## Eigendecompose Mss and construct symmetric square root
  r <- eigen(M, symmetric = TRUE)
  v <- r$vectors
  rootM <- tcrossprod(v, tcrossprod(v, sqrt(diag(r$values))))
  InvrootM <- solve(rootM)
  #negrootMss <- chol(solve(Mss),pivot=FALSE)
  return(InvrootM)
}
#-------------------------------------------------------------------
