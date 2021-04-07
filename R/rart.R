#' Groups of Transformations for ART
#'
#' This helper function computes the group G of transformations for the CRS test.
#'
#' @param q the dimension of the data / # of clusters (q point estimators)
#' @param B the number of random draws if q>10
#'
#' @return A B x q matrix of transformations (G)
random.G <- function(q, B = 10000) {
  if (q < 11) {
    e <- c(1, -1);
    B <- 2 ^ q;
    perm <- t(expand.grid(rep(list(e), q)));
  }
  # Draw B transformations from G otherwise
  else {
    perm <- matrix(sample(e, q * B, replace = T), q, B);
  }
  return(perm)
}

#-------------------------------------------------------------------

#' Hypothesis Test for ART
#'
#' This helper function computes the randomization critical value for t-tests.
#'
#' @param c.beta vector of parameters entering the null hypothesis c'beta - with one estimator per cluster (q x 1)
#' @param G the group of all trasformations (use random.G to get)
#' @param lambda scalar for the null hypothesis c'beta = lambda (default 0)
#' @param alpha significance level (dafault 0.05)
#' @param nj q x 1 vector of sample sizes in each cluster (for alternative weighting)
#'
#' @return A list contating the following elements:
#' \describe{
#'  \item{`rule`}{binary decision of the test}
#'  \item{`Nrule`}{binary decision of the non-randomized test}
#'  \item{`cv`}{the critical value of the randomization test}
#'  \item{`pv`}{p-value according to folmula 15.5}
#'  \item{`pv2`}{p-value according to folmula 15.7}
#' }
CRS.test <- function(c.beta, G, lambda = 0, alpha = 0.05, nj = 1) {
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

#' Confidence Intervals for ART
#'
#' This helper function computes confidence intervals by test inversion using bisection.
#'
#' @param c.beta vector of parameters entering the null hypothesis c'beta
#'               - with one estimator per cluster (q x 1)
#' @param G the group of all trasformations (use random.G to get)
#' @param alpha significance level (dafault 0.05)
#' @param nj q x 1 vector of sample sizes in each cluster (for alternative weighting)
#'
#' @return The 1-alpha confidence interval for c'beta
CRS.CI <- function(c.beta, G, alpha = 0.05, nj = 1) {
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
