#' @importFrom stats as.formula binomial lm predict sd
NULL
#' rankIC
#'
#' Fit the semiparametric accelerated failure time model with rank regression.
#'
#'
#' @param U left-censoring time, having 0 if left-censored.
#' @param V right-censoring time, having \code{Inf} if right-censored.
#' @param Delta censoring indicator, 1: observed; 0: interval-censored.
#' @param X baseline covariate.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param alpha scalar to adjust informative cluster sizes, 0<alpha<=1. Cluster size will be ignored by using \code{alpha=0}.
#' @param type type of rank estimation, \code{gehan}: gehan estimation; \code{logrank}: log-rank estimation.
#' @param maxit maximum number of iteration for the log-rank estimator, default is 20.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-5.
#' @param nboot number of bootstrapped sample generation for variance estimation of regression estimator (Zeng and Uin, 2008). For stable estimation, \code{nboot=200} is usually recommended.
#'
#' @return \code{aft_rank} returns a data frame containing at least the following components:
#' \itemize{
#'   \item \code{est}: regression estimator.
#'   \item \code{se}: standard error estimates for \code{est}.
#'   \item \code{pvalue}: p-value, which will be appeared when the value of \code{nboot} is larger than one.
#' }
#'
#' @details
#' see Choi et al., (2022+) for detailed method explanation.
#'
#' @references
#' Zeng, D. and Lin, D. (2008).  
#' Efficient resampling methods for nonsmooth estimating functions. Biostatistics 9(2): 355–363.
#' 
#' Choi, T., Choi, S. and Bandyopadhyay, D. (2022+). 
#' Rank estimation for the accelerated failure time model with partially interval-censored data. 
#' Submitted to Biometrics.
#'
#'
#' @examples
#' \dontrun{
#' # Simulations
#' set.seed(111)
#' n = 200
#' x1 = rnorm(n)
#' x2 = ifelse(rnorm(n)>0, 1, 0)
#' X = cbind(x1,x2)
#' T = 2 + x1 + x2 + rnorm(n)
#' U = (1 - 0.25*x1)*runif(n, -6, 3.85)
#' V = U + (1-0.5*x2)*runif(n, 6, 19.6)
#' U = exp(dplyr::case_when(TRUE ~ T, T>V ~ V, T<U ~ -Inf))
#' V = exp(dplyr::case_when(TRUE ~ T, T>V ~ Inf, T<U ~ U))
#' Delta = ifelse(U==V, 1, 0)
#' aft_rank(U,V,X,Delta)
#' 
#' # Data example
#' library(PICBayes)
#' data(mCRC)
#' dt0 = as.data.frame(mCRC)
#' d = with(dt0,
#'          data.frame(U = ifelse(is.na(L), 0, L),
#'                     V = ifelse(is.na(R), Inf, R),
#'                     Delta = 1-IC,
#'                     x1 = TRT_C,
#'                     x2 = KRAS_C,
#'                     id = SITE))
#' U = d$U; V = d$V; X = cbind(d$x1, d$x2); Delta = d$Delta; id = d$id
#' aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
#'          alpha = 1, type = "gehan", nboot = 10)
#' aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
#'          alpha = 1, type = "logrank", nboot = 10)
#' }
#' @export
#'
#'
#'




aft_rank=function(U, V, X, Delta, 
                  id = NULL, alpha = 1, type = c("gehan","logrank"), 
                  maxit = 20, tol = 1e-5, nboot = 0){
  gkern=function(x, h=0.01) 1/(1 + exp(-x/h))

  wgh_fit = function(U, V, X, Delta, 
                     id = NULL, alpha = 1,
                     weight = rep(1, length(U)))
  {
    if (is.null(id)) id = rep(1,length(U))
    options(warn = -1)
    n = length(U)
    y = U
    m = table(id)
    U = log(y); V = log(V)
    n = nrow(X); p = ncol(X)
    d1 = V<Inf; d2 = U>-Inf
    id_i = rep(1:n, each = n)
    id_j = rep(1:n, times = n)
    idd = which(d1[id_i]*d2[id_j]==1)
    m = rep(m, m); mi = m[id_i]; mj = m[id_j]
    zi = weight
    xi = X[id_i,]; xj = X[id_j,]
    yi = ifelse(Delta == 1, log(y), V)[id_i]
    yj = ifelse(Delta == 1, log(y), U)[id_j]
    yy = ((yi - yj)*zi/(mi*mj)^alpha)[idd]
    xx = ((xi - xj)*zi/(mi*mj)^alpha)[idd,]
    yy.new = c(yy, 1e10)
    xx.new = rbind(xx, -colSums(xx))
    quantreg::rq.fit(x = xx.new, y = yy.new)$coef
  }

  wfun = function(U, V, X, Delta, id, alpha, init) 
  {
    n = length(U)
    y = pmax(U, 1e-8);
    U = log(y); V = log(V)
    m = table(id)
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    xi = X[id_i,]; xj = X[id_j,]
    m = rep(m, m); mi = m[id_i]; mj = m[id_j]
    yi = ifelse(Delta == 1, log(y), V)[id_i]
    yj = ifelse(Delta == 1, log(y), U)[id_j]
    eii = c(yi - xi%*%init)
    ejj = c(yj - xj%*%init)
    1/pmax(sapply(split(gkern(ejj - eii)/
                          (mj^alpha), id_i), mean)[id_i], 1e-3)
  }

  aft_rq_se = function(B, est = NULL, type = "gehan", alpha = 1) 
  {
    efun = function(U, V, Delta, X, id, beta, type, alpha) 
    {
      y = pmax(U, 1e-8);
      n = length(U)
      m = table(id)
      id_i = rep(1:n, each = n)
      id_j = rep(1:n, times = n)
      xi = X[id_i,]; xj = X[id_j,]
      m = rep(m, m); mi = m[id_i]; mj = m[id_j]
      d1i = (V<Inf)[id_i]; d2j = (U>-Inf)[id_j]
      idd = which(d1i*d2j == 1)
      yi = ifelse(Delta == 1, log(y), V)[id_i]
      yj = ifelse(Delta == 1, log(y), U)[id_j]
      ei = c(yi - xi%*%beta)
      ej = c(yj - xj%*%beta)
      eta = I(ei <= ej)
      if (type == "gehan") wi = 1
      if (type == "logrank") {
        wi = 1/pmax(sapply(split(gkern(
          ifelse(yj > log(1e-8), ej - ei, 0))/(mj^alpha), 
          id_i),mean)[id_i], 1e-8)
      }
      res = colMeans(d1i*d2j*wi/((mi*mj)^alpha)*(xi - xj)*eta)/n
      res
    }
    p = length(est); n = length(U)
    library(MASS)
    Shat = t(replicate(B,{
      Bid = sample(n, n, replace = TRUE);
      n^(-1/2)*efun(U[Bid], V[Bid], Delta[Bid], 
                    X[Bid,], id[Bid], est, type, alpha)}))
    An = matrix(0, p, p)
    zmat = matrix(stats::rnorm(p*B), B, p)
    Umat = matrix(0, B, p)
    for(b in 1:B){
      Umat[b,] = efun(U, V, Delta, X, id,
                   est+n^(-1/2)*zmat[b,], 
                   type = type, alpha)*n^(-1/2)
    }
    Vi=stats::var(Shat)
    for (i in 1:p) {
      An[i,] = lm(Umat[,i] ~ matrix(zmat, ncol = p))$coef[-1]
    }
    covmat=(solve(An)%*%Vi)%*%solve(An)/n
    covmat=matrix(as.numeric(covmat),p)
    sqrt(diag(covmat))
  }
  X = as.matrix(X)
  init = beta_new = drop(wgh_fit(U, V, X, Delta, id, alpha))

  if (type == "logrank") {
    beta_old = init
    err = 10; iter = 0
    while(iter < maxit & err > tol) {
      wi = wfun(U, V, X, Delta, id, alpha, init = beta_old)
      beta_new = drop(wgh_fit(U, V, X, Delta, 
                              id, alpha = alpha, weight = wi))
      err = sum((beta_new - beta_old)^2)
      iter  = iter + 1
      beta_old = beta_new
    }
  }
  res = cbind(est = beta_new)
  if (nboot > 1) {
    se = aft_rq_se(est = beta_new, type = type, B = nboot, alpha = alpha)
    res = cbind(est = beta_new, 
                se = se, 
                pvalue = 1 - pnorm(abs(beta_new/se)))
  }
  return(round(res, 3))
}




# roxygen2::roxygenize()