#' @importFrom stats as.formula binomial lm predict sd
NULL
#' rankIC
#'
#' Fit the semiparametric accelerated failure time model with rank regression.
#'
#'
#' @param L left-censoring time, having 0 if left-censored.
#' @param R right-censoring time, having \code{Inf} if right-censored.
#' @param delta censoring indicator, 1: observed; 0: interval-censored.
#' @param X baseline covariate.
#' @param id cluster id. If the data does not have clustered structure, set \code{id=NULL}.
#' @param alpha baseline covariate for fitting outcome regression model.
#' @param type type of rank estimation, \code{gehan}: gehan estimation; \code{logrank}: log-rank estimation.
#' @param maxit maximum number of iteration for the log-rank estimator, default is 20.
#' @param tol tolerance of iteration for the log-rank estimator, default is 1e-5.
#' @param nboot number of bootstrapped sample generation for variance estimation of regression estimator (Zeng and Lin, 2008). For stable estimation, \code{nboot=200} is usually recommended.
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
#' Zeng, D. and Lin, D. (2008). Efficient resampling methods for nonsmooth estimat- ing functions, Biostatistics 9(2): 355–363.
#'
#' Choi, T., Choi, S. and Bandyopadhyay, D. (2022+). Rank estimation for the accelerated failure time model with partially interval-censored data. Submitted, \emph{Biometrics}.
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
#' L = (1 - 0.25*x1)*runif(n, -6, 3.85)
#' R = L + (1-0.5*x2)*runif(n, 6, 19.6)
#' L = exp(dplyr::case_when(TRUE ~ T, T>R ~ R, T<L ~ -Inf))
#' R = exp(dplyr::case_when(TRUE ~ T, T>R ~ Inf, T<L ~ L))
#' delta = ifelse(L==R, 1, 0)
#' aft_rank(L,R,X,delta)
#' 
#' # Data example
#' library(PICBayes)
#' data(mCRC)
#' dt0 = as.data.frame(mCRC)
#' d = with(dt0,
#'          data.frame(L = ifelse(is.na(L), 0, L),
#'                     R = ifelse(is.na(R), Inf, R),
#'                     delta = 1-IC,
#'                     x1 = TRT_C,
#'                     x2 = KRAS_C,
#'                     id = SITE))
#' L = d$L; R = d$R; X = cbind(d$x1, d$x2); delta = d$delta; id = d$id
#' aft_rank(L = L, R = R, X = X, delta = delta, id = id, 
#'          alpha = 1, type = "gehan", nboot = 10)
#' aft_rank(L = L, R = R, X = X, delta = delta, id = id, 
#'          alpha = 1, type = "logrank", nboot = 10)
#' }
#' @export
#'
#'
#'




aft_rank=function(L, R, X, delta, 
                  id = NULL, alpha = 1, type = c("gehan","logrank"), 
                  maxit = 20, tol = 1e-5, nboot = 0){
  gkern=function(x, h=0.01) 1/(1 + exp(-x/h))

  wgh_fit = function(L, R, X, delta, 
                     id = NULL, alpha = 1,
                     weight = rep(1, length(L)))
  {
    if (is.null(id)) id = rep(1,length(L))
    options(warn = -1)
    n = length(L)
    y = L
    m = table(id)
    L = log(y); R = log(R)
    n = nrow(X); p = ncol(X)
    d1 = R<Inf; d2 = L>-Inf
    id_i = rep(1:n, each = n)
    id_j = rep(1:n, times = n)
    idd = which(d1[id_i]*d2[id_j]==1)
    m = rep(m, m); mi = m[id_i]; mj = m[id_j]
    zi = weight
    xi = X[id_i,]; xj = X[id_j,]
    yi = ifelse(delta == 1, log(y), R)[id_i]
    yj = ifelse(delta == 1, log(y), L)[id_j]
    yy = ((yi - yj)*zi/(mi*mj)^alpha)[idd]
    xx = ((xi - xj)*zi/(mi*mj)^alpha)[idd,]
    yy.new = c(yy, 1e10)
    xx.new = rbind(xx, -colSums(xx))
    quantreg::rq.fit(x = xx.new, y = yy.new)$coef
  }

  wfun = function(L,R,X,delta,id,alpha,init) 
  {
    n = length(L)
    y = pmax(L, 1e-8);
    L = log(y); R = log(R)
    m = table(id)
    id_i = rep(1:n, each=n)
    id_j = rep(1:n, times=n)
    xi = X[id_i,]; xj = X[id_j,]
    m = rep(m, m); mi = m[id_i]; mj = m[id_j]
    yi = ifelse(delta == 1, log(y), R)[id_i]
    yj = ifelse(delta == 1, log(y), L)[id_j]
    eii = c(yi - xi%*%init)
    ejj = c(yj - xj%*%init)
    1/pmax(sapply(split(gkern(ejj - eii)/
                          (mj^alpha), id_i), mean)[id_i], 1e-3)
  }

  aft_rq_se = function(B, est = NULL, type = "gehan", alpha = 1) 
  {
    efun = function(L, R, delta, X, id, beta, type, alpha) 
    {
      y = pmax(L, 1e-8);
      n = length(L)
      m = table(id)
      id_i = rep(1:n, each = n)
      id_j = rep(1:n, times = n)
      xi = X[id_i,]; xj = X[id_j,]
      m = rep(m, m); mi = m[id_i]; mj = m[id_j]
      d1i = (R<Inf)[id_i]; d2j = (L>-Inf)[id_j]
      idd = which(d1i*d2j == 1)
      yi = ifelse(delta == 1, log(y), R)[id_i]
      yj = ifelse(delta == 1, log(y), L)[id_j]
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
    p = length(est); n = length(L)
    library(MASS)
    Shat = t(replicate(B,{
      Bid = sample(n, n, replace = TRUE);
      n^(-1/2)*efun(L[Bid], R[Bid], delta[Bid], 
                    X[Bid,], id[Bid], est, type, alpha)}))
    An = matrix(0, p, p)
    zmat = matrix(stats::rnorm(p*B), B, p)
    U = matrix(0, B, p)
    for(b in 1:B){
      U[b,] = efun(L, R, delta, X, id,
                   est+n^(-1/2)*zmat[b,], 
                   type = type, alpha)*n^(-1/2)
    }
    Vi=stats::var(Shat)
    for (i in 1:p) {
      An[i,] = lm(U[,i] ~ matrix(zmat, ncol = p))$coef[-1]
    }
    covmat=(solve(An)%*%Vi)%*%solve(An)/n
    covmat=matrix(as.numeric(covmat),p)
    sqrt(diag(covmat))
  }
  X = as.matrix(X)
  init = beta_new = drop(wgh_fit(L, R, X, delta, id, alpha))

  if (type == "logrank") {
    beta_old = init
    err = 10; iter = 0
    while(iter < maxit & err > tol) {
      wi = wfun(L, R, X, delta, id, alpha, init = beta_old)
      beta_new = drop(wgh_fit(L, R, X, delta, 
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
