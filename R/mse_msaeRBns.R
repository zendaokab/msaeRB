#' @title Parametric Bootstrap Mean Squared Error Estimators of Ratio Benchmarking for Multivariate Non Sampled Area in Small Area Estimation
#'
#' @description Calculates the parametric bootstrap mean squared error estimates of ratio benchmarking for multivariate non sampled area in small area estimation
#'
#' @param formula an object of class list of formula describe the fitted models
#' @param vardir matrix containing sampling variances of direct estimators. The order is: \code{var1, cov12, ..., cov1r, var2, cov23, ..., cov2r, ..., cov(r-1)(r), var(r)}
#' @param weight matrix containing proportion of units in small areas. The order is: \code{w1, w2, ..., w(r)}
#' @param cluster matrix containing cluster of auxiliary variables. The order is: \code{c1, c2, ..., c(r)}
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
#' @export mse_msaeRBns
#'
#' @import abind
#' @importFrom magic adiag
#' @importFrom Matrix forceSymmetric
#' @importFrom stats model.frame na.omit model.matrix median pnorm rnorm
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \donttest{
#' ## load dataset
#' data(datamsaeRBns)
#'
#' # Compute MSE EBLUP and Ratio Benchmark
#' # This is the long running example
#' ## Using parameter 'data'
#' Fo = list(f1 = Y1 ~ X1 + X2,
#'           f2 = Y2 ~ X1 + X2,
#'           f3 = Y3 ~ X1 + X2)
#' vardir = c("v1", "v12", "v13", "v2", "v23", "v3")
#' weight = c("w1", "w2", "w3")
#' cluster = c("c1", "c2", "c3")
#'
#' mse_msae = est_msaeRBns(Fo, vardir, weight, cluster, data = datamsaeRBns)
#'
#' ## Without parameter 'data'
#' Fo = list(f1 = datamsaeRBns$Y1 ~ datamsaeRBns$X1 + datamsaeRBns$X2,
#'           f2 = datamsaeRBns$Y2 ~ datamsaeRBns$X1 + datamsaeRBns$X2,
#'           f3 = datamsaeRBns$Y3 ~ datamsaeRBns$X1 + datamsaeRBns$X2)
#' vardir = datamsaeRBns[, c("v1", "v12", "v13", "v2", "v23", "v3")]
#' weight = datamsaeRBns[, c("w1", "w2", "w3")]
#' cluster = datamsaeRBns[, c("c1", "c2", "c3")]
#'
#' mse_msae = mse_msaeRBns(Fo, vardir, weight, cluster)
#'
#' ## Return
#' mse_msae$pbmse.eblupRB # to see the MSE of Ratio Benchmark
#' }
mse_msaeRBns = function (formula, vardir, weight, cluster, samevar = FALSE, B = 1000, MAXITER = 100, PRECISION = 1e-04, data) {
  start_time = Sys.time()
  r = length(formula)
  if (r <= 1)
    stop("You should use mse_saeRB() for univariate")
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
    }
    return(as.matrix(R))
  }
  temporary_vardir = function(vardir, n, r) {
    dt = list()
    for (h in 1:n) {
      var_mat = matrix(0, nrow = r, ncol = r)
      k = 1
      for (i in 1:r) {
        for (j in 1:r) {
          if (i <= j) {
            var_mat[i, j] = vardir[h, k]
            k = k + 1
          }
        }
      }
      var_mat = forceSymmetric(var_mat)
      dt[[h]] = as.matrix(var_mat)
    }
    return(dt)
  }
  eblup_inside = function(r, n, samevar, y, X, cluster, vardir, MAXITER = 100, PRECISION = 1e-04){
    y.matrix = matrix(as.vector(y), n, r)
    y.matrix[y.matrix==0] = NA
    indexns = unique(which(is.na(y.matrix), arr.ind=TRUE)[,1])
    indexs = c(1:n)[-indexns]
    n.s = length(indexs)
    n.ns = n - n.s
    y.s = y.matrix[-indexns, ]
    y.s.vec = as.vector(y.s)
    corr = cor(y.s)
    hasil_cor = c()
    for (i in 1:r){
      for (j in 1:r){
        if (i < j && i != j){
          hasil_cor = c(hasil_cor, corr[i, j])
        }
      }
    }
    vardir[indexns, ] = NA
    index.vardir = sum(1:r)
    for (i in 2:r){
      index.vardir[i] = index.vardir[1] - sum(2:i)
    }
    index.vardir = sort(index.vardir)
    index.covardir = c(1:sum(1:r))[-index.vardir]
    varians.direct = vardir[, index.vardir]
    covarians.direct = vardir[, -index.vardir]
    for (i in 1:r) {
      df.vardir = data.frame(cluster[, i], varians.direct[, i])
      for (j in 1:n) {
        if (df.vardir[j, 2] %in% NA) {
          df.vardir[j, 2] = mean(df.vardir[df.vardir[, 1] == df.vardir[j,1], 2], na.rm = TRUE)
        }
      }
      varians.direct[ ,i] = df.vardir[, 2]
    }
    index.cor = 0
    if (r > 2){
      for (p in 1:r) {
        for (q in 1:r) {
          if (p < q) {
            index.cor = index.cor + 1
            for (col in 1:n) {
              if (covarians.direct[col, index.cor] %in% NA) {
                covarians.direct[col, index.cor] = hasil_cor[index.cor] * sqrt(varians.direct[col, p] * varians.direct[col, q])
              }
            }
          }
        }
      }
    } else {
      index.cor = index.cor + 1
      for (col in 1:n) {
        if (covarians.direct[col] %in% NA) {
          covarians.direct[col] = hasil_cor[index.cor] * sqrt(varians.direct[col, 1] * varians.direct[col, 2])
        }
      }
    }
    varians.direct = rbind(index.vardir, varians.direct)
    if (r > 2) {
      covarians.direct = rbind(index.covardir, covarians.direct)
    } else {
      covarians.direct = c(index.covardir, covarians.direct)
    }
    vardir.full = cbind(varians.direct, covarians.direct)
    vardir.full = vardir.full[, order(vardir.full[1, ])]
    vardir.full = vardir.full[-1, ]
    R.ns = R_function(vardir.full[indexns, ], n.ns, r)
    vardir.s = vardir[-indexns, ]
    R = R_function(vardir.s, n.s, r)
    W.s = W[-indexns, ]
    W.s = prop.table(W.s, 2)
    Xindexns = c()
    for (i in 1:r) {
      Xindexns = c(Xindexns, indexns + rep(n, times = n.ns)*(i - 1))
    }
    X.s = X[-Xindexns, ]
    X.ns = X[Xindexns, ]
    y_names = sapply(formula, "[[", 2)
    Ir = diag(r)
    In = diag(n.s)
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
      Gn.ns = kronecker(diag(Vu), diag(n.ns))
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
      Gn.ns = kronecker(diag(Vu), diag(n.ns))
      Vinv.ns = solve(Gn.ns + R.ns)
      Q.ns = ginv(t(Vinv.ns %*% X.ns) %*% X.ns)
      g1.s = diag(Gn %*% Vinv %*% R)
      g2.s = diag(R %*% Vinv %*% X.s %*% Q %*% t(X.s) %*% t(R %*% Vinv))
      g1.ns = diag(Gn.ns %*% Vinv.ns %*% R.ns)
      g2.ns = diag(R.ns %*% Vinv.ns %*% X.ns %*% Q.ns %*% t(X.ns) %*% t(R.ns %*% Vinv.ns))
      g12.s = matrix(g1.s + g2.s, n.s, r)
      names(g12.s) = y_names
      g12.s = data.frame(index = indexs, g12.s)
      g12.ns = matrix(g1.ns + g2.ns, n.ns, r)
      names(g12.ns) = y_names
      g12.ns = data.frame(index = indexns, g12.ns)
      g12 = rbind(g12.s, g12.ns)
      g12 = g12[order(g12$index), ]
      rownames(g12) = g12$index
      g12 = g12[, -1]
      dg = Vinv - Gn %*% Vinv %*% Vinv
      gg3 = (dg %*% V %*% t(dg))/iF
      g3 = diag(gg3)
      g3.matrix = matrix(g3, n.s, r)
      g3.s = data.frame(index = indexs, g3.matrix)
      g3.ns = data.frame(index = indexns, matrix(NA, n.ns, r))
      names(g3.ns) = names(g3.s)
      g3.all = rbind(g3.s, g3.ns)
      g3.all = g3.all[order(g3.all$index), ]
      rownames(g3.all) = g3.all$index
      g3.all = g3.all[ ,-1]
      g3.full = matrix(NA, nrow = n, ncol = r)
      colnames(g3.full) = y_names
      for (i in 1:r) {
        df = data.frame(cluster[, i], g3.all[, i])
        for (j in 1:n) {
          if (df[j, 2] %in% NA) {
            df[j, 2] = mean(df[df[, 1] == df[j,1], 2], na.rm = TRUE)
          }
        }
        g3.full[ ,i] = df[, 2]
      }
      names(g3.full) = y_names
      mse.df = as.data.frame(g12 + 2 * g3.full)
    } else {
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
      Q = solve(XtVinv %*% X.s)
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
      FI = solve(iF)
      Gn.ns = kronecker(Gr, diag(n.ns))
      Vinv.ns = solve(Gn.ns + R.ns)
      XtVinv.ns = t(Vinv.ns %*% X.ns)
      Q.ns = ginv(XtVinv.ns %*% X.ns)
      g1.s = diag(Gn %*% Vinv %*% R)
      g2.s = diag(R %*% Vinv %*% X.s %*% Q %*% t(X.s) %*% t(R %*% Vinv))
      g1.ns = diag(Gn.ns %*% Vinv.ns %*% R.ns)
      g2.ns = diag(R.ns %*% Vinv.ns %*% X.ns %*% Q.ns %*% t(X.ns) %*% t(R.ns %*% Vinv.ns))
      g12.s = matrix(g1.s + g2.s, n.s, r)
      names(g12.s) = y_names
      g12.s = data.frame(index = indexs, g12.s)
      g12.ns = matrix(g1.ns + g2.ns, n.ns, r)
      names(g12.ns) = y_names
      g12.ns = data.frame(index = indexns, g12.ns)
      g12 = rbind(g12.s, g12.ns)
      g12 = g12[order(g12$index), ]
      rownames(g12) = g12$index
      g12 = g12[, -1]
      dg = lapply(dV1, function(x) x %*% Vinv - Gn %*% Vinv %*% x %*% Vinv)
      gg3 = list()
      for (i in 1:r) {
        for (j in 1:r) {
          gg3[[(i - 1) * r + j]] = FI[i, j] * (dg[[i]] %*% V %*% t(dg[[j]]))
        }
      }
      g3 = diag(Reduce("+", gg3))
      g3.matrix = matrix(g3, n.s, r)
      g3.s = data.frame(index = indexs, g3.matrix)
      g3.ns = data.frame(index = indexns, matrix(NA, n.ns, r))
      names(g3.ns) = names(g3.s)
      g3.all = rbind(g3.s, g3.ns)
      g3.all = g3.all[order(g3.all$index), ]
      rownames(g3.all) = g3.all$index
      g3.all = g3.all[ ,-1]
      g3.full = matrix(NA, nrow = n, ncol = r)
      colnames(g3.full) = y_names
      for (i in 1:r) {
        df = data.frame(cluster[, i], g3.all[, i])
        for (j in 1:n) {
          if (df[j, 2] %in% NA) {
            df[j, 2] = mean(df[df[, 1] == df[j,1], 2], na.rm = TRUE)
          }
        }
        g3.full[ ,i] = df[, 2]
      }
      names(g3.full) = y_names
      mse.df = as.data.frame(g12 + 2 * g3.full)
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
    result = list(eblup = list(est.eblup = NA, est.eblupRB = NA), fit = list(method = NA, convergence = NA, iteration = NA, estcoef = NA, refvar = NA), random.effect = NA, agregation = NA, mse = NA, g12 = NA, R.s = NA, R.ns = NA)
    result$eblup$est.eblup = eblup
    result$eblup$est.eblupRB = eblup.ratio
    result$fit$method = "REML"
    result$fit$convergence = convergence
    result$fit$iteration = k
    result$fit$estcoef = coef
    result$fit$refvar = t(Vu)
    result$random.effect = random.effect.full
    result$mse = mse.df
    result$g12 = g12
    result$R.s = R
    result$R.ns = R.ns
    return(result)
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
    samevar = samevar
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
    samevar = samevar
  }
  y.matrix = matrix(as.vector(y), n, r)
  y.matrix[y.matrix==0] = NA
  indexns = unique(which(is.na(y.matrix), arr.ind=TRUE)[,1])
  indexs = c(1:n)[-indexns]
  n.s = length(indexs)
  n.ns = n - n.s
  y.s = y.matrix[-indexns, ]
  y.s.vec = as.vector(y.s)
  W.s = W[-indexns, ]
  W.s = prop.table(W.s, 2)
  Xindexns = c()
  for (i in 1:r) {
    Xindexns = c(Xindexns, indexns + rep(n, times = n.ns)*(i - 1))
  }
  X.s = X[-Xindexns, ]
  X.ns = X[Xindexns, ]
  y_names = sapply(formula, "[[", 2)
  temp_vardir = temporary_vardir(vardir[-indexns, ], n = n.s, r = r)
  eblup_first = eblup_inside(r = r, n = n, samevar = samevar, y = y, X = X, cluster = cluster, vardir = vardir)
  R.s = eblup_first$R.s
  R.ns = eblup_first$R.ns
  beta = eblup_first$fit$estcoef[, 1]
  A = eblup_first$fit$refvar
  A_mat = kronecker(diag(as.vector(A)), diag(n.s))
  Vinv = solve(A_mat + R.s)
  XtVinv = t(Vinv %*% X.s)
  Q = solve(XtVinv %*% X.s)
  mse_prasad = eblup_first$mse
  g12 = eblup_first$g12
  sumg12.pb = rep(0, n * r)
  sumg3.pb = rep(0, n * r)
  boot = 1
  while (boot <= B) {
    u.boot = mvrnorm(n = n.s, mu = rep(0, r), Sigma = (diag(as.vector(A), nrow = r, ncol = r)))
    theta.boot = X.s %*% beta + as.vector(u.boot)
    e.boot = matrix(0, nrow = n.s, ncol = r)
    for (i in 1:n.s) {
      e.boot[i, ] = mvrnorm(1, mu = rep(0, r), Sigma = (temp_vardir[[i]]))
    }
    direct.boot = theta.boot + as.vector(e.boot)
    direct.boot.mat = matrix(direct.boot, nrow = n.s, ncol = r)
    direct.boot.mat.s = data.frame(index = indexs, direct.boot.mat)
    direct.boot.mat.ns = data.frame(index = indexns, matrix(0, n.ns, r))
    direct.boot.mat.full = rbind(direct.boot.mat.s, direct.boot.mat.ns)
    direct.boot.mat.full = direct.boot.mat.full[order(direct.boot.mat.full$index), ]
    rownames(direct.boot.mat.full) = direct.boot.mat.full$index
    direct.boot.mat.full = direct.boot.mat.full[, -1]
    resultEBLUP = eblup_inside(r = r, n = n, samevar = samevar, y = unlist(direct.boot.mat.full), X = X, cluster = cluster, vardir = vardir)
    sigma2.simula = resultEBLUP$fit$refvar
    beta.simula = resultEBLUP$fit$estcoef[, 1]
    Gn.simula = kronecker(diag(as.vector(sigma2.simula)), diag(n.s))
    Vinv.simula = as.matrix(solve(Gn.simula + R.s))
    Xbeta.simula = X %*% beta.simula
    Xsbeta.simula = X.s %*% beta.simula
    XtVi.simula = t(Vinv.simula %*% X.s)
    Q.simula = solve(XtVi.simula %*% X.s)
    #randomEffect
    random.effect.simula = Gn.simula %*% Vinv.simula %*% (direct.boot - Xsbeta.simula)
    random.effect.simula = matrix(random.effect.simula, nrow = n.s, ncol = r)
    random.effect.simula = data.frame(index = indexs, random.effect.simula)
    random.effect.simula.ns = data.frame(index = indexns, matrix(NA, n.ns, r))
    names(random.effect.simula.ns) = names(random.effect.simula)
    random.effect.simula = rbind(random.effect.simula, random.effect.simula.ns)
    random.effect.simula = random.effect.simula[order(random.effect.simula$index), ]
    rownames(random.effect.simula) = random.effect.simula$index
    random.effect.simula = random.effect.simula[ ,-1]
    random.effect.simula.full = matrix(NA, nrow = n, ncol = r)
    colnames(random.effect.simula.full) = y_names
    for (i in 1:r) {
      df = data.frame(cluster[, i], random.effect.simula[, i])
      for (j in 1:n) {
        if (df[j, 2] %in% NA) {
          df[j, 2] = mean(df[df[, 1] == df[j,1], 2], na.rm = TRUE)
        }
      }
      random.effect.simula.full[ ,i] = df[, 2]
    }
    thetaEBLUP.boot1 = Xbeta.simula + as.vector(random.effect.simula.full)
    thetaEBLUP.boot1.mat = matrix(thetaEBLUP.boot1, nrow = n, ncol = r)
    thetaRATIO.boot1.mat = matrix(0, nrow = n, ncol = r)
    for (i in 1:r) {
      thetaRATIO.boot1.mat[, i] = thetaEBLUP.boot1.mat[, i] * (sum(W.s[, i] * direct.boot.mat[, i]) / sum(W[, i] * thetaEBLUP.boot1.mat[, i]))
    }
    g1boot.s = diag(Gn.simula %*% Vinv.simula %*% R.s)
    g2boot.s = diag(R.s %*% Vinv.simula %*% X.s %*% Q.simula %*% t(X.s) %*% t(R.s %*% Vinv.simula))
    Gn.simula.ns = kronecker(diag(as.vector(sigma2.simula)), diag(n.ns))
    Vinv.simula.ns = as.matrix(solve(Gn.simula.ns + R.ns))
    XtVinv.simula.ns = t(Vinv.simula.ns %*% X.ns)
    Q.ns = ginv(XtVinv.simula.ns %*% X.ns)
    g1boot.ns = diag(Gn.simula.ns %*% Vinv.simula.ns %*% R.ns)
    g2boot.ns = diag(R.ns %*% Vinv.simula.ns %*% X.ns %*% Q.ns %*% t(X.ns) %*% t(R.ns %*% Vinv.simula.ns))
    g12boot.s = matrix(g1boot.s + g2boot.s, n.s, r)
    names(g12boot.s) = y_names
    g12boot.s = data.frame(index = indexs, g12boot.s)
    g12boot.ns = matrix(g1boot.ns + g2boot.ns, n.ns, r)
    names(g12boot.ns) = y_names
    g12boot.ns = data.frame(index = indexns, g12boot.ns)
    g12boot = rbind(g12boot.s, g12boot.ns)
    g12boot = g12boot[order(g12boot$index), ]
    rownames(g12boot) = g12boot$index
    g12boot = g12boot[, -1]
    Bstim.eblup = solve(XtVinv %*% X.s) %*% XtVinv %*% direct.boot
    Xbeta.eblup = X %*% Bstim.eblup
    Xsbeta.eblup = X.s %*% Bstim.eblup
    random.effect2.simula = A_mat %*% Vinv %*% (direct.boot - Xsbeta.eblup)
    random.effect2.simula = matrix(random.effect2.simula, nrow = n.s, ncol = r)
    random.effect2.simula = data.frame(index = indexs, random.effect2.simula)
    random.effect2.simula.ns = data.frame(index = indexns, matrix(NA, n.ns, r))
    names(random.effect2.simula.ns) = names(random.effect2.simula)
    random.effect2.simula = rbind(random.effect2.simula, random.effect2.simula.ns)
    random.effect2.simula = random.effect2.simula[order(random.effect2.simula$index), ]
    rownames(random.effect2.simula) = random.effect2.simula$index
    random.effect2.simula = random.effect2.simula[ ,-1]
    random.effect2.simula.full = matrix(NA, nrow = n, ncol = r)
    colnames(random.effect2.simula.full) = y_names
    for (i in 1:r) {
      df = data.frame(cluster[, i], random.effect2.simula[, i])
      for (j in 1:n) {
        if (df[j, 2] %in% NA) {
          df[j, 2] = mean(df[df[, 1] == df[j,1], 2], na.rm = TRUE)
        }
      }
      random.effect2.simula.full[ ,i] = df[, 2]
    }
    thetaEBLUP.boot2 = Xbeta.eblup + as.vector(random.effect2.simula.full)
    thetaEBLUP.boot2.mat = matrix(thetaEBLUP.boot2, nrow = n, ncol = r)
    thetaRATIO.boot2.mat = matrix(0, nrow = n, ncol = r)
    for (i in 1:r) {
      thetaRATIO.boot2.mat[, i] = thetaEBLUP.boot2.mat[, i] * (sum(W.s[, i] * direct.boot.mat[, i]) / sum(W[, i] * thetaEBLUP.boot2.mat[, i]))
    }
    g3boot = (thetaRATIO.boot1.mat - thetaRATIO.boot2.mat)^2
    sumg12.pb = sumg12.pb + as.vector(g12boot)
    sumg3.pb = sumg3.pb + as.vector(g3boot)
    boot = boot + 1
  }
  g12.pb = sumg12.pb/B
  g3.pb = sumg3.pb/B
  msebootratio.df = 2 * (g12) - g12.pb + g3.pb
  names(mse_prasad) = y_names
  names(msebootratio.df) = y_names
  end_time = Sys.time()
  running_time = end_time - start_time
  result1 = list(mse.eblup = NA, pbmse.eblupRB = NA, running.time = NA)
  result1$mse.eblup = mse_prasad
  result1$pbmse.eblupRB = msebootratio.df
  result1$running.time = running_time
  return(result1)
}
