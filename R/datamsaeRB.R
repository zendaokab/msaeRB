#' @title Sample Data for Multivariate Small Area Estimation with Ratio Benchmarking
#'
#' @description Dataset to simulate ratio benchmarking of Multivariate Fay-Herriot model
#'
#' This data is generated based on multivariate Fay-Herriot model by these following steps:
#' \enumerate{
#'   \item Generate explanatory variables \code{X1} and \code{X2}. \code{X1 ~ N(10, 1)} and \code{X2 ~ U(9.5, 10.5)}.
#'   \cr Sampling error \code{e} is generated with the following \eqn{\sigma_{e11}}{\sigmae11} = 0.01,  \eqn{\sigma_{e22}}{\sigmae22} = 0.02, \eqn{\sigma_{e33}}{\sigmae33} = 0.03, and \eqn{\rho_{e}}{\rhoe} = 1/2.
#'   \cr For random effect \code{u}, we set \eqn{\sigma_{u11}}{\sigmau11}= 0.02, \eqn{\sigma_{u22}}{\sigmau22}= 0.03, and \eqn{\sigma_{u33}}{\sigmau33}= 0.04.
#'   \cr For the weight, we generate \code{w1, w2, w3} by set {w1, w2, w3 ~ U(10, 20)}
#'   \cr Set beta, \eqn{\beta01} = 10, \eqn{\beta02} = 8, \eqn{\beta03} = 6, \eqn{\beta11} = -0.3, \eqn{\beta12} = 0.2, \eqn{\beta13} = 0.4, \eqn{\beta21} = 0.5, \eqn{\beta22} = -0.1, and \eqn{\beta23} = -0.2.
#'   \cr Calculate direct estimation \code{Y1 Y2 Y3} where \eqn{Y_{i}}{Yi} = \eqn{X * \beta + u_{i} + e_{i}}{X\beta+ui+ei}
#'   \item Then combine the direct estimations \code{Y1 Y2 Y3}, explanatory variables \code{X1 X2}, weight \code{w1 w2 w3}, and sampling varians covarians \code{v1 v12 v13 v2 v23 v3} in a dataframe then named as datamsaeRB
#' }
#'
#' @format A data frame with 30 rows and 14 variables:
#' \describe{
#'   \item{Y1}{Direct Estimation of Y1}
#'   \item{Y2}{Direct Estimation of Y2}
#'   \item{Y3}{Direct Estimation of Y3}
#'   \item{X1}{Auxiliary variable of X1}
#'   \item{X2}{Auxiliary variable of X2}
#'   \item{w1}{Known proportion of units in small areas of Y1}
#'   \item{w2}{Known proportion of units in small areas of Y2}
#'   \item{w3}{Known proportion of units in small areas of Y3}
#'   \item{v1}{Sampling Variance of Y1}
#'   \item{v12}{Sampling Covariance of Y1 and Y2}
#'   \item{v13}{Sampling Covariance of Y1 and Y3}
#'   \item{v2}{Sampling Variance of Y2}
#'   \item{v23}{Sampling Covariance of Y2 and Y3}
#'   \item{v3}{Sampling Variance of Y3}
#' }
#'
"datamsaeRB"
