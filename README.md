# rART

The `rART` package provides functions for Approximate Randomization Tests with a Small Number of Clusters introduced in "Randomization tests under an approximate symmetry assumption" (Canay et al. 2017) and further described in "A User's Guide to Approximate Randomization Tests with a Small Number of Clusters" (Cai et al. 2021).

The package is centered around a single function `artlm` which runs the main regression. Several companion functions are provided in order to use the method. In order to demonstrate this, we generate some data.


```r
cn <- 100 # cluster size
nc <- 10   # number of clusters

x1 <- rnorm(cn * nc, sd = 2)
x2 <- rnorm(cn * nc, mean = x1)
x3 <- runif(cn * nc)

eps <- matrix(runif(cn * nc), cn, nc)

# make some correlation within groups
eps <- as.vector(eps + eps[sample(cn),] + eps[sample(cn),])
group = factor(rep(1:nc, each = cn))

# true model
y <- 1 + 10*x1 + x2/20 + eps

df <- data.frame(y, x1, x2, x3, group)

head(df)
#>           y         x1          x2        x3 group
#> 1 -8.990301 -1.1209513 -2.11675002 0.3044642     1
#> 2 -1.943071 -0.4603550 -1.50031002 0.8328188     1
#> 3 34.619987  3.1174166  3.09943639 0.5936475     1
#> 4  3.566540  0.1410168  0.00884165 0.8071966     1
#> 5  4.400315  0.2585755 -2.29076730 0.2940508     1
#> 6 37.156422  3.4301300  4.47070343 0.1410852     1
```

# Regression

The linear ART regression is `artlm`. It is exactly the same as `lm`. However, it requires a cluster variable or vector to be specified. 


```r
(artlm1 <- artlm(y ~ x1 + x2, cluster=group, data=df))
#> 
#> Call:
#> artlm(formula = y ~ x1 + x2, data = df, cluster = group)
#> 
#> Coefficients:
#> (Intercept)           x1           x2  
#>     2.51695     10.02792      0.03328
```

It supports all the features that `lm` supports. For example, you can add fixed effects. 


```r
(artlm2 <- artlm(y ~ x1 + x2 + group - 1, cluster=group, data=df))
#> 
#> Call:
#> artlm(formula = y ~ x1 + x2 + group - 1, data = df, cluster = group)
#> 
#> Coefficients:
#>       x1        x2    group1    group2    group3    group4    group5    group6    group7    group8    group9   group10  
#> 10.02793   0.03229   2.58429   2.53483   2.43459   2.43422   2.51628   2.59488   2.36401   2.55289   2.53291   2.62133
```

You can also specify that you are only interested in a subset of the variables. For example, suppose I am only interested in `x1`. 


```r
(artlm3 <- artlm(y ~ x1 + x2 + group - 1, cluster=group, select = "x1", data=df))
#> 
#> Call:
#> artlm(formula = y ~ x1 + x2 + group - 1, data = df, cluster = group, 
#>     select = "x1")
#> 
#> Coefficients:
#>       x1        x2    group1    group2    group3    group4    group5    group6    group7    group8    group9   group10  
#> 10.02793   0.03229   2.58429   2.53483   2.43459   2.43422   2.51628   2.59488   2.36401   2.55289   2.53291   2.62133
```

Including the `select` option allows you to choose the parameters you are interested in. The result will not display in the regression object itself. However, chosen variables will not de displayed in `summary` etc.

# Summarizing models

You can run the simplest form of ART by running `summary` on any regression object. For example, we can apply it to our weighted fixed effects regression.


```r
summary(artlm2)
#> 
#> Call:
#> artlm(formula = y ~ x1 + x2 + group - 1, data = df, cluster = group)
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.49058 -0.38165  0.01422  0.37043  1.30953 
#> 
#> Coefficients:
#>    Estimate Crit. value t value Pr(>|t|)   
#> x1 10.02793     0.71236 196.106  0.00195 **
#> x2  0.03229     0.71394   0.658  0.07031 . 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.504 on 988 degrees of freedom
#> Multiple R-squared:  0.9994,	Adjusted R-squared:  0.9994 
#> F-statistic: 1.332e+05 on 12 and 988 DF,  p-value: < 2.2e-16
```

This gives us the OLS estimates, the t-statistic, the p-value, and the critical value of the t-test for 95\% confidence.

If we summarize the regression where we selected to view only `x1`, then we will only see this one variable in the summary. Note that no selection is required to ignore the group fixed effects because ART cannot be applied to these coefficients.


```r
summary(artlm3)
#> 
#> Call:
#> artlm(formula = y ~ x1 + x2 + group - 1, data = df, cluster = group, 
#>     select = "x1")
#> 
#> Residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.49058 -0.38165  0.01422  0.37043  1.30953 
#> 
#> Coefficients:
#>    Estimate Crit. value t value Pr(>|t|)   
#> x1  10.0279      0.7124   196.1  0.00195 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.504 on 988 degrees of freedom
#> Multiple R-squared:  0.9994,	Adjusted R-squared:  0.9994 
#> F-statistic: 1.332e+05 on 12 and 988 DF,  p-value: < 2.2e-16
```

# Confidence intervals

Confidence intervals are computed using the `confint` command as usual.


```r
confint(artlm2)
#>           2.5 %     97.5 %
#> x1  9.992375525 10.0690510
#> x2 -0.003316027  0.0752652
```

You can adjust the level of significance using the `level` parameter.


```r
confint(artlm2, level = 0.98)
#>            1 %        99 %
#> x1  9.97959628 10.06905101
#> x2 -0.01448144  0.08521353
```

The command will only display selected variables. So, if we use `artlm3` where we selected only `x1`, then we will only see the confidence interval for this one term. More importantly, the function will not compute intervals for unselected parameters.


```r
confint(artlm3)
#>     2.5 %    97.5 % 
#>  9.992376 10.069051
```

It is also possible to select variables in `confint` instead of in the regression.


```r
confint(artlm2, parm = "x2")
#>        2.5 %       97.5 % 
#> -0.003316027  0.075265196
```

# Linear tests

You can conduct arbitrary linear tests using `ARTHypothesis`. For example, suppose I wanted to test to see if $\beta_1 = \beta_2$. Then, I can run.


```r
ARTHypothesis(artlm2, "x1 = x2")
#> Linear ART hypothesis test
#> 
#> Hypothesis:
#> x1 - x2 = 0
#> 
#>   Crit. value t value Pr(>|t|)   
#> 1     0.71342  96.966 0.001953 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

You can also supply the constraint vector and constant value manually like so.


```r
ARTHypothesis(artlm2, c(1,-1), 0)
#> Linear ART hypothesis test
#> 
#> Hypothesis:
#> x1 - x2 = 0
#> 
#>   Crit. value t value Pr(>|t|)   
#> 1     0.71342  96.966 0.001953 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

You cannot construct a linear test using a parameter that was not selected in the original regression.
