#' @title Parametric Bootstrap Mean Squared Error Estimators of Ratio Benchmarking for Univariate Small Area Estimation
#'
#' @description Calculates the parametric bootstrap mean squared error estimates of ratio benchmarking for univariate small area estimation
#'
#' @param formula an object of class list of formula describe the fitted model
#' @param vardir vector containing sampling variances of direct estimators
#' @param weight vector containing proportion of units in small areas
#' @param samevar logical. If \code{TRUE}, the varians is same. Default is \code{FALSE}
#' @param B number of bootstrap. Default is 1000
#' @param MAXITER maximum number of iterations for Fisher-scoring. Default is 100
#' @param PRECISION coverage tolerance limit for the Fisher Scoring algorithm. Default value is \code{1e-4}
#' @param data dataframe containing the variables named in formula, vardir, and weight
#'
#' @return
#' \item{mse.eblup}{estimated mean squared errors of the EBLUPs for the small domains based on Prasad Rao}
#' \item{pbmse.eblupRB}{parametric bootstrap mean squared error estimates of the ratio benchmark}
#' \item{running.time}{time for running function}
#'
#' @export mse_saeRB
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \dontrun{
#' # You can still run this example. Because of the long time example, so i added the "dontrun" code
#' ## load dataset
#' data(datamsaeRB)
#'
#' # Compute MSE EBLUP and Ratio Benchmark
#'
#' ## Using parameter 'data'
#' mse_sae = mse_saeRB(Y1 ~ X1 + X2, v1, w1, data = datamsaeRB)
#'
#' ## Without parameter 'data'
#' mse_sae = mse_saeRB(datamsaeRB$Y1 ~ datamsaeRB$X1 + datamsaeRB$X2, datamsaeRB$v1, datamsaeRB$w1)
#'
#' ## Return
#' mse_sae$pbmse.eblupRB # to see the MSE Ratio Benchmark estimators
#' }
mse_saeRB = function(formula, vardir, weight, samevar = FALSE, B = 1000, MAXITER = 100, PRECISION = 1E-04, data) {
  start_time <- Sys.time()
  if (!is.list(formula))
    formula = list(formula)
  r = length(formula)
  if (r > 1)
    stop("You should use mse_msaeRB() for multivariate")
  R_function = function (vardir, n, r) {
    if (r == 1) {
      R = diag(vardir)
    } else {
      R = matrix(rep(0, times = n*r*n*r), nrow = n*r, ncol = n*r)
      k = 1
      for (i in 1:r) {
        for (j in 1:r) {
          if (i <= j) {
            mat0 = matrix(rep(0, times = r*r), nrow = r, ncol = r)
            mat0[i, j] = 1
            matVr = diag(vardir[, k], length(vardir[, k]))
            R_hasil = kronecker(mat0, matVr)
            R = R + R_hasil
            k = k + 1
          }
        }
      }
      R = forceSymmetric(R)
    }
    return(as.matrix(R))
  }
  eblup_inside = function(r, n, samevar, y, X, R, MAXITER = 100, PRECISION = 1E-04) {
    y_names = sapply(formula, "[[", 2)
    Ir = diag(r)
    In = diag(n)
    dV = list()
    dV1 = list()
    for (i in 1:r){
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
        XtVinv = t(Vinv %*% X)
        Q = solve(XtVinv %*% X)
        P = Vinv - t(XtVinv) %*% Q %*% XtVinv
        Py = P %*% y
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
      XtVinv = t(Vinv %*% X)
      Q = solve(XtVinv %*% X)
      P = Vinv - t(XtVinv) %*% Q %*% XtVinv
      Py = P %*% y
      beta = Q %*% XtVinv %*% y
      res = y - X %*% beta
      eblup = data.frame(matrix(X %*% beta + Gn %*% Vinv %*% res, n, r))
      names(eblup) = y_names
      se.b = sqrt(diag(Q))
      t.value = beta/se.b
      p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
      coef = as.matrix(cbind(beta, se.b, t.value, p.value))
      colnames(coef) = c("beta", "std. error", "t value", "p-value")
      rownames(coef) = colnames(X)
      g1 = diag(Gn %*% Vinv %*% R)
      g2 = diag(R %*% Vinv %*% X %*% Q %*% t(X) %*% t(R %*% Vinv))
      dg = Vinv - Gn %*% Vinv %*% Vinv
      gg3 = (dg %*% V %*% t(dg))/iF
      g3 = diag(gg3)
      mse = g1 + g2 + 2 * g3
      mse.df = data.frame(matrix(data = mse, nrow = n, ncol = r))
      names(mse.df) = y_names
    } else {
      Vu = apply(matrix(diag(R), nrow = n, ncol = r), 2, median)
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
        XtVinv = t(Vinv %*% X)
        Q = solve(XtVinv %*% X)
        P = Vinv - t(XtVinv) %*% Q %*% XtVinv
        Py = P %*% y
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
        Gr = Vu1
      } else {
        Gr = diag(as.vector(Vu1))
      }
      Gn = kronecker(Gr, In)
      V = as.matrix(Gn + R)
      Vinv = solve(V)
      XtVinv = t(Vinv %*% X)
      Q = solve(XtVinv %*% X)
      P = Vinv - t(XtVinv) %*% Q %*% XtVinv
      Py = P %*% y
      beta = Q %*% XtVinv %*% y
      res = y - X %*% beta
      eblup = data.frame(matrix(X %*% beta + Gn %*% Vinv %*% res, n, r))
      names(eblup) = y_names
      se.b = sqrt(diag(Q))
      t.value = beta/se.b
      p.value = 2 * pnorm(abs(as.numeric(t.value)), lower.tail = FALSE)
      coef = as.matrix(cbind(beta, se.b, t.value, p.value))
      colnames(coef) = c("beta", "std. error", "t value", "p-value")
      rownames(coef) = colnames(X)
      FI = solve(iF)
      g1 = diag(Gn %*% Vinv %*% R)
      g2 = diag(R %*% Vinv %*% X %*% Q %*% t(X) %*% t(R %*% Vinv))
      dg = lapply(dV1, function(x) x %*% Vinv - Gn %*% Vinv %*% x %*% Vinv)
      gg3 = list()
      for (i in 1:r) {
        for (j in 1:r) {
          gg3[[(i - 1) * r + j]] = FI[i, j] * (dg[[i]] %*% V %*% t(dg[[j]]))
        }
      }
      g3 = diag(Reduce("+", gg3))
      mse = g1 + g2 + 2 * g3
      mse.df = data.frame(matrix(data = mse, nrow = n, ncol = r))
      names(mse.df) = y_names
    }
    result = list(eblup = NA, fit = list(estcoef = NA, refvar = NA), mse = NA, mse_component = list(g1 = NA, g2 = NA, g3 = NA))
    result$eblup = eblup
    result$fit$estcoef = coef
    result$fit$refvar = t(Vu)
    result$mse = mse.df
    result$mse_component$g1 = g1
    result$mse_component$g2 = g2
    result$mse_component$g3 = g3
    return(result)
  }
  namevar = deparse(substitute(vardir))
  nameweight = deparse(substitute(weight))
  if (!missing(data)) {
    formuladata = lapply(formula, function(x) model.frame(x, na.action = na.omit, data))
    y = unlist(lapply(formula, function(x) model.frame(x, na.action = na.omit, data)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x, data)))
    W = as.matrix(data[, nameweight])
    n = length(y)/r
    if (any(is.na(data[, namevar])))
      stop("Object vardir contains NA values.")
    if (any(is.na(data[, nameweight])))
      stop("Object weight contains NA values.")
    R = R_function(data[, namevar], n, r)
    vardir = data[, namevar]
    samevar = samevar
  } else {
    formuladata = lapply(formula, function(x) model.frame(x, na.action = na.omit))
    y = unlist(lapply(formula, function(x) model.frame(x, na.action = na.omit)[[1]]))
    X = Reduce(adiag, lapply(formula, function(x) model.matrix(x)))
    W = as.matrix(weight)
    n = length(y)/r
    if (any(is.na(vardir)))
      stop("Object vardir contains NA values")
    if (any(is.na(weight)))
      stop("Object weight contains NA values.")
    R = R_function(vardir, n, r)
    samevar = samevar
  }
  y_names = sapply(formula, "[[", 2)
  eblup_first = eblup_inside(r = r, n = n, samevar = samevar, y = y, X = X, R = R)
  beta = eblup_first$fit$estcoef[,1]
  A = eblup_first$fit$refvar
  A_mat = kronecker(A, diag(n))
  Vinv = solve(kronecker(A, diag(n)) + R)
  XtVinv = t(Vinv %*% X)
  Q = solve(XtVinv %*% X)
  mse_prasad = eblup_first$mse
  g1d = eblup_first$mse_component$g1
  g2d = eblup_first$mse_component$g2
  sumg1.pb = rep(0, n*r)
  sumg2.pb = rep(0, n*r)
  sumg3.pb = rep(0, n*r)
  cat("\nBootstrap procedure with B =", B, "iterations starts.\n")
  boot <- 1
  while (boot <= B) {
    u.boot = rnorm(n, 0, sqrt(A))
    theta.boot = X %*% beta + u.boot
    e.boot = rnorm(n, 0, sqrt(as.vector(vardir)))
    direct.boot = theta.boot + e.boot
    direct.boot.mat = matrix(direct.boot, nrow = n, ncol = r)
    resultEBLUP = eblup_inside(r = r, n = n, samevar = samevar, y = direct.boot, X = X, R = R)
    sigma2.simula = resultEBLUP$fit$refvar
    beta.simula = resultEBLUP$fit$estcoef[,1]
    Gn.simula = kronecker(as.vector(sigma2.simula), diag(n))
    Vinv.simula = as.matrix(solve(Gn.simula + R))
    Xbeta.simula = X %*% beta.simula
    XtVi.simula = t(Vinv.simula %*% X)
    Q.simula = solve(XtVi.simula %*% X)
    thetaEBLUP.boot1 = Xbeta.simula + Gn.simula %*% Vinv.simula %*% (direct.boot - Xbeta.simula)
    thetaEBLUP.boot1.mat = as.matrix(thetaEBLUP.boot1)
    thetaRATIO.boot1 = thetaEBLUP.boot1 %*% (colSums(W * direct.boot.mat)/colSums(W * thetaEBLUP.boot1.mat))
    g1boot = diag(Gn.simula %*% Vinv.simula %*% R)
    g2boot = diag(R %*% Vinv.simula %*% X %*% Q.simula %*% t(X) %*% t(R %*% Vinv.simula))
    Bstim.eblup = solve(XtVinv %*% X) %*% XtVinv %*% direct.boot
    Xbeta.eblup = X %*% Bstim.eblup
    thetaEBLUP.boot2 = Xbeta.eblup + A_mat %*% Vinv %*% (direct.boot - Xbeta.eblup)
    thetaEBLUP.boot2.mat = as.matrix(thetaEBLUP.boot2)
    thetaRATIO.boot2 = thetaEBLUP.boot2 %*% (colSums(W * direct.boot.mat)/colSums(W * thetaEBLUP.boot2.mat))
    g3boot = (thetaRATIO.boot1 - thetaRATIO.boot2)^2
    sumg1.pb = sumg1.pb + g1boot
    sumg2.pb = sumg2.pb + g2boot
    sumg3.pb = sumg3.pb + g3boot
    boot = boot + 1
  }
  cat("\nBootstrap procedure with B =", B, "iterations ends.\n")
  g1.pb = sumg1.pb/B
  g2.pb = sumg2.pb/B
  g3.pb = sumg3.pb/B
  msebootratio = 2*(g1d + g2d) - g1.pb - g2.pb + g3.pb
  msebootratio.df = data.frame(matrix(data = msebootratio, nrow = n, ncol = r))
  names(msebootratio.df) = y_names
  end_time <- Sys.time()
  running_time = end_time - start_time
  result1 = list(mse.eblup = NA, pbmse.eblupRB = NA, running.time = NA)
  result1$mse.eblup = mse_prasad
  result1$pbmse.eblupRB = msebootratio.df
  result1$running.time = running_time
  return(result1)
}
