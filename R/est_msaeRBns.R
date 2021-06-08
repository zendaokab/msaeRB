#' @title EBLUPs Ratio Benchmarking for Non Sampled Area based on a Multivariate Fay Herriot (Model 1)
#'
#' @description This function gives EBLUPs ratio benchmarking for non sampled area based on multivariate Fay-Herriot (Model 1)
#'
#' @param formula an object of class list of formula describe the fitted models
#' @param vardir matrix containing sampling variances of direct estimators. The order is: \code{var1, cov12, ..., cov1r, var2, cov23, ..., cov2r, ..., cov(r-1)(r), var(r)}
#' @param weight matrix containing proportion of units in small areas. The order is: \code{w1, w2, ..., w(r)}
#' @param cluster matrix containing cluster of auxiliary variables. The order is: \code{c1, c2, ..., c(r)}
#' @param samevar logical. If \code{TRUE}, the varians is same. Default is \code{FALSE}
#' @param MAXITER maximum number of iterations for Fisher-scoring. Default is 100
#' @param PRECISION coverage tolerance limit for the Fisher Scoring algorithm. Default value is \code{1e-4}
#' @param data dataframe containing the variables named in formula, vardir, and weight
#'
#' @return This function returns a list with following objects:
#' \item{eblup}{a list containing a value of estimators}
#' \itemize{
#'   \item est.eblup : a dataframe containing EBLUP estimators
#'   \item est.eblupRB : a dataframe containing ratio benchmark estimators
#' }
#'
#' \item{fit}{a list contining following objects:}
#' \itemize{
#'   \item method : fitting method, named "REML"
#'   \item convergence : logical value of convergence of Fisher Scoring
#'   \item iterations : number of iterations of Fisher Scoring algorithm
#'   \item estcoef : a data frame containing estimated model coefficients (\code{beta, std. error, t value, p-value})
#'   \item refvar : estimated random effect variance
#' }
#' \item{random.effect}{a data frame containing values of random effect estimators}
#' \item{agregation}{a data frame containing agregation of direct, EBLUP, and ratio benchmark estimation}
#' @export est_msaeRBns
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm
#' @importFrom MASS mvrnorm
#'
#' @examples
#' ## load dataset
#' data(datamsaeRBns)
#'
#' # Compute EBLUP and Ratio Benchmark using auxiliary variables X1 and X2 for each dependent variable
#'
#' ## Using parameter 'data'
#' Fo = list(f1 = Y1 ~ X1 + X2,
#'           f2 = Y2 ~ X1 + X2,
#'           f3 = Y3 ~ X1 + X2)
#' vardir = c("v1", "v12", "v13", "v2", "v23", "v3")
#' weight = c("w1", "w2", "w3")
#' cluster = c("c1", "c2", "c3")
#'
#' est_msae = est_msaeRBns(Fo, vardir, weight, cluster, data = datamsaeRBns)
#'
#' ## Without parameter 'data'
#' Fo = list(f1 = datamsaeRBns$Y1 ~ datamsaeRBns$X1 + datamsaeRBns$X2,
#'           f2 = datamsaeRBns$Y2 ~ datamsaeRBns$X1 + datamsaeRBns$X2,
#'           f3 = datamsaeRBns$Y3 ~ datamsaeRBns$X1 + datamsaeRBns$X2)
#' vardir = datamsaeRBns[, c("v1", "v12", "v13", "v2", "v23", "v3")]
#' weight = datamsaeRBns[, c("w1", "w2", "w3")]
#' cluster = datamsaeRBns[, c("c1", "c2", "c3")]
#'
#' est_msae = est_msaeRBns(Fo, vardir, weight, cluster)
#'
#' ## Return
#' est_msae$eblup$est.eblupRB # to see the Ratio Benchmark estimators
#'
est_msaeRBns = function (formula, vardir, weight, cluster, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04, data)
{
  r = length(formula)
  if (r <= 1)
    stop("You should use est_saeRB() for univariate")
  R_function = function(vardir, n, r) {
    if (r == 1) {
      R = diag(vardir)
    }
    else {
      R = matrix(rep(0, times = n * r * n * r), nrow = n *
                   r, ncol = n * r)
      k = 1
      for (i in 1:r) {
        for (j in 1:r) {
          if (i <= j) {
            mat0 = matrix(rep(0, times = r * r), nrow = r,
                          ncol = r)
            mat0[i, j] = 1
            matVr = diag(vardir[, k], length(vardir[,
                                                    k]))
            R_hasil = kronecker(mat0, matVr)
            R = R + R_hasil
            k = k + 1
          }
        }
      }
      R = forceSymmetric(R)
      R = R
    }
    return(as.matrix(R))
  }
  if (!missing(data)) {
    formuladata = lapply(formula, function(x) model.frame(x, na.action = na.omit, data))
    y = unlist(lapply(formula, function(x) model.frame(x, na.action = na.omit, data)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x, data)))
    W = as.matrix(data[, weight])
    n = length(y)/r
    if (!all(vardir %in% names(data)))
      stop("Object vardir is not appropriate with data.")
    if (length(vardir) != sum(1:r))
      stop("Length of vardir is not appropriate with data. The length must be ", sum(1:r))
    if (any(is.na(data[, weight])))
      stop("Object weight contains NA values.")
    if (!all(weight %in% names(data)))
      stop("Object weight is not appropriate with data.")
    if (length(weight) != r)
      stop("Length of weight is not appropriate with data. The length must be ", r)
    if (!all(cluster %in% names(data)))
      stop("Object cluster is not appropriate with data.")
    if (length(cluster) != r)
      stop("Length of cluster is not appropriate with data. The length must be ", r)
    vardir = data[, vardir]
    cluster = as.matrix(data[, cluster])
  } else {
    formuladata = lapply(formula, function(x) model.frame(x, na.action = na.omit))
    y = unlist(lapply(formula, function(x) model.frame(x, na.action = na.omit)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x)))
    W = as.matrix(weight)
    n = length(y)/r
    if ((dim(vardir)[1] != n) || (dim(vardir)[2] != sum(1:r)))
      stop("Object vardir is not appropriate with data. It must be ", n ," x ", sum(1:r) ," matrix.")
    if (any(is.na(weight)))
      stop("Object weight contains NA values.")
    if ((dim(weight)[1] != n) || (dim(weight)[2] != r))
      stop("Object weight is not appropriate with data. It must be ", n ," x ", r ," matrix.")
    if (any(is.na(cluster)))
      stop("Object cluster contains NA values.")
    if ((dim(cluster)[1] != n) || (dim(cluster)[2] != r))
      stop("Object cluster is not appropriate with data. It must be ", n ," x ", r ," matrix.")
    cluster = as.matrix(cluster)
  }
  y.matrix = matrix(as.vector(y), n, r)
  y.matrix[y.matrix == 0] = NA
  indexns = unique(which(is.na(y.matrix), arr.ind = TRUE)[, 1])
  indexs = c(1:n)[-indexns]
  n.s = length(indexs)
  n.ns = n - n.s
  y.s = y.matrix[-indexns, ]
  y.s.vec = as.vector(y.s)
  vardir.s = vardir[-indexns, ]
  R = R_function(vardir.s, n.s, r)
  W.s = W[-indexns, ]
  W.s = prop.table(W.s, 2)
  Xindexns = c()
  for (i in 1:r) {
    Xindexns = c(Xindexns, indexns + rep(n, times = n.ns)*(i - 1))
  }
  X.s = X[-Xindexns, ]
  y_names = sapply(formula, "[[", 2)
  Ir = diag(r)
  In = diag(n.s)
  dV = list()
  dV1 = list()
  for (i in 1:r) {
    dV[[i]] = matrix(0, nrow = r, ncol = r)
    dV[[i]][i, i] = 1
    dV1[[i]] = kronecker(dV[[i]], In)
  }
  convergence = TRUE
  if (samevar) {
    Vu = median(diag(R))
    k = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > PRECISION) & (k < MAXITER)) {
      k = k + 1
      Vu1 = Vu
      Gr = kronecker(Vu1, Ir)
      Gn = kronecker(Gr, In)
      V = as.matrix(Gn + R)
      Vinv = solve(V)
      XstVinv = t(Vinv %*% X.s)
      Q = solve(XstVinv %*% X.s)
      P = Vinv - t(XstVinv) %*% Q %*% XstVinv
      Py = P %*% y.s.vec
      s = (-0.5) %*% sum(diag(P)) + 0.5 %*% (t(Py) %*% Py)
      iF =  0.5 %*% sum(diag(P %*% P))
      Vu = Vu1 + solve(iF) %*% s
      diff = abs((Vu - Vu1)/Vu1)
    }
    Vu = as.vector((rep(max(Vu, 0), r)))
    names(Vu) = y_names
    if (k >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    Gn = kronecker(diag(Vu), In)
    V = as.matrix(Gn + R)
    Vinv = solve(V)
    XstVinv = t(Vinv %*% X.s)
    Q = solve(XstVinv %*% X.s)
    P = Vinv - t(XstVinv) %*% Q %*% XstVinv
    Py = P %*% y.s.vec
    beta = Q %*% XstVinv %*% y.s.vec
    res = y.s.vec - X.s %*% beta
    random.effect = data.frame(matrix(Gn %*% Vinv %*% res, n.s, r))
    names(random.effect) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error", "t value", "p-value")
    rownames(coef) = colnames(X)
  }
  else {
    Vu = apply(matrix(diag(R), nrow = n.s, ncol = r), 2, median)
    k = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > rep(PRECISION, r)) & (k < MAXITER)) {
      k = k + 1
      Vu1 = Vu
      if (r == 1) {
        Gr = Vu1
      } else {
        Gr = diag(as.vector(Vu1))
      }
      Gn = kronecker(Gr, In)
      V = as.matrix(Gn + R)
      Vinv = solve(V)
      XstVinv = t(Vinv %*% X.s)
      Q = solve(XstVinv %*% X.s)
      P = Vinv - t(XstVinv) %*% Q %*% XstVinv
      Py = P %*% y.s.vec
      s = sapply(dV1, function(x) (-0.5) * sum(diag(P %*% x)) + 0.5 * (t(Py) %*% x %*% Py))
      iF = matrix(unlist(lapply(dV1, function(x) lapply(dV1, function(y) 0.5 * sum(diag(P %*% x %*% P %*% y))))), r)
      Vu = Vu1 + solve(iF) %*% s
      diff = abs((Vu - Vu1)/Vu1)
    }
    Vu = as.vector(sapply(Vu, max, 0))
    if (k >= MAXITER && diff >= PRECISION){
      convergence = FALSE
    }
    if (r == 1) {
      Gr = Vu
    } else {
      Gr = diag(as.vector(Vu))
    }
    Gn = kronecker(Gr, In)
    V = as.matrix(Gn + R)
    Vinv = solve(V)
    XtVinv = t(Vinv %*% X.s)
    Q = as.matrix(solve(XtVinv %*% X.s))
    P = Vinv - t(XtVinv) %*% Q %*% XtVinv
    Py = P %*% y.s.vec
    beta = Q %*% XtVinv %*% y.s.vec
    res = y.s.vec - X.s %*% beta
    random.effect = data.frame(matrix(Gn %*% Vinv %*% res, n.s, r))
    names(random.effect) = y_names
    se.b = sqrt(diag(Q))
    t.value = beta/se.b
    p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
    coef = as.matrix(cbind(beta, se.b, t.value, p.value))
    colnames(coef) = c("beta", "std. error", "t value", "p-value")
    rownames(coef) = colnames(X)
  }
  random.effect = data.frame(index = indexs, random.effect)
  random.effect.ns = data.frame(index = indexns, matrix(NA, n.ns, r))
  names(random.effect.ns) = names(random.effect)
  random.effect = rbind(random.effect, random.effect.ns)
  random.effect = random.effect[order(random.effect$index), ]
  rownames(random.effect) = random.effect$index
  random.effect = random.effect[ ,-1]
  random.effect.full = matrix(NA, nrow = n, ncol = r)
  colnames(random.effect.full) = y_names
  for (i in 1:r) {
    df = data.frame(cluster[, i], random.effect[, i])
    for (j in 1:n) {
      if (df[j, 2] %in% NA) {
        df[j, 2] = mean(df[df[, 1] == df[j,1], 2], na.rm = TRUE)
      }
    }
    random.effect.full[ ,i] = df[, 2]
  }
  Xbeta = X %*% beta
  eblup = matrix(Xbeta, nrow = n, ncol = r) + random.effect.full
  eblup = as.data.frame(eblup)
  names(eblup) = y_names
  random.effect.full = as.data.frame(random.effect.full)
  eblup.mat = as.matrix(eblup)
  eblup.ratio = matrix(0, nrow = n, ncol = r)
  alfa = c()
  for (i in 1:r) {
    eblup.ratio[, i] = eblup.mat[, i] * (sum(W.s[, i] * y.s[, i])/sum(W[, i] * eblup.mat[, i]))
  }
  eblup.ratio = as.data.frame(eblup.ratio)
  names(eblup.ratio) = y_names
  agregation.direct = diag(t(W.s) %*% y.s)
  agregation.eblup = diag(t(W) %*% eblup.mat)
  agregation.eblup.ratio = diag(t(W) %*% as.matrix(eblup.ratio))
  agregation = as.matrix(rbind(agregation.direct, agregation.eblup, agregation.eblup.ratio))
  colnames(agregation) = y_names
  result = list(eblup = list(est.eblup = NA, est.eblupRB = NA), fit = list(method = NA, convergence = NA, iteration = NA, estcoef = NA, refvar = NA), random.effect = NA, agregation = NA)
  result$eblup$est.eblup = eblup
  result$eblup$est.eblupRB = eblup.ratio
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iteration = k
  result$fit$estcoef = coef
  result$fit$refvar = t(Vu)
  result$random.effect = random.effect.full
  result$agregation = agregation
  return(result)
}
