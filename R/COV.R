#' @name COV
#' @title \bold{Cov}ariance matrix for survival data
#'
#' @include tn.R
#' 
#' @rdname COV
#' @export
#'
COV <- function(x, ...) UseMethod("COV")
#'
#'
#' @param x
#'  For the default method, a \code{numeric} vector of
#'  \emph{number of events}, \eqn{e_t}{e[t]}.
#'  These are assumed to be ordered by discrete times.
#'  \cr
#'  A method is available for objects of \code{class} \code{tne}.
#' @param ... Additional arguments (not implemented).
#'  \cr
#'  The following arguments apply only to the default method.
#' @param n \bold{N}umber at risk (total).
#' @param ncg \bold{N}umber at risk, per \bold{c}ovariate \bold{g}roup.
#'  \cr
#' If there are \eqn{2} groups, this can be given as a \code{vector} with
#' the number at risk for group \eqn{1}.
#' \cr
#' If there are \eqn{\geq 2}{>=2} groups, it is
#' a \code{matrix} with one column for each group.
#'
#' @details Gives variance-covariance matrix for comparing survival
#' data for two or more groups.
#' \cr
#' Inputs are vectors corresponding to observations at a set of discrete
#' time points for right censored data, except for \eqn{n1},
#' the no. at risk by predictor.
#' \cr
#' This should be specified as a vector for one group,
#' otherwise as a matrix with each column corresponding to a group.
#'
#' @return An \code{array}.
#' \cr
#' The first two dimensions = the number of covariate groups \eqn{K},
#' \eqn{k = 1, 2, \ldots K}.
#' This is the square matrix below.
#' \cr
#' The third dimension is the number of observations
#' (discrete time points).
#' \cr \cr
#' To calculate this, we use \code{x} (= \eqn{e_t}{e[t]} below) and
#' \eqn{n_1}{n1}, the number at risk in covariate group \eqn{1}.
#' \cr
#' Where there are \eqn{2} groups, the resulting sparse square matrix
#' (i.e. the non-diagonal elements are \eqn{0})
#' at time \eqn{t} has diagonal elements:
#'  \deqn{cov_t = - \frac{n_{0t} n_{1t} e_t (n_t - e_t)}{n_t^2(n_t-1)}}{
#'        cov[t] = - n0[t] * n1[t] * e[t] * (n[t] - e[t]) /
#'                  (n[t]^2 * (n[t] - 1))}
#' For \eqn{\geq 2}{>=2} groups, the resulting square matrix
#' has diagonal elements given by:
#'  \deqn{cov_{kkt} = \frac{n_{kt}(n_t - n_{kt}) e_t(n_t - e_t)}{
#'                          n_t^2(n_t - 1)}}{
#'    cov[k, k, t] = n[k, t] * (n[t] - n[k, t]) * e[t] * (n[t] - e[t]) /
#'                   (n[t]^2 * (n[t] - 1))}
#' The off diagonal elements are:
#' \deqn{cov_{klt} = \frac{-n_{kt} n_{lt} e_t (n_t-e_t) }{
#'                         n_t^2(n_t-1)}}{
#'       cov[k, l, t] = - n[k, t] * n[l, t] * e[t] * (n[t] - e[t]) /
#'                      n[t]^2 * (n[t] - 1)}
#'
#' @seealso Called by \code{\link{comp}}
#'
#' @keywords survival
#'
#' @rdname COV
#' @aliases COV.tn
#' @method COV tn
#' @export
#'
#' @examples
#' ## Two covariate groups
#' ## K&M 2nd ed. Example 7.2, pg 210.
#' data("kidney", package="KMsurv")
#' k1 <- with(kidney,
#'            tne(Surv(time=time, event=delta) ~ type, shape="long"))
#' ##with(k1, COV(x=e, n=n, n1=n_1))
#' COV(k1)
#' ## Three covariate groups
#' ## K&M 2nd ed. Example 7.4, pp 212-214.
#' data("bmt", package="KMsurv")
#' b1 <- tne(Surv(time=t2, event=d3) ~ group, data=bmt, shape="long")
#' COV(b1)
COV.tne <- function(x, ...){
    ## no. of groups
    g1 <- attr(x, "ncg")
    if (g1 <= 1) stop("Only valid if more than one covariate group")
    ## if 2 groups only
    if (g1==2){
        res2 <- x[, (ncg / n) * (1 - (ncg / n)) * ((n - e) / (n - 1)) * e, by=list(t, cg)]
        res2 <- res2[, sum(V1), by=t]
        res1 <- res2[, V1]
        names(res1) <- res2[, t]
    }
    if (g1 > 2){
        t1 <- x[, sort(unique(t))]
        y1 <- x[, unique(cg)]
        ## get no. at risk for each unique time and covariate group
        n1 <- lapply(t1, FUN=function(t1) (setkey(x[t >= t1, max(ncg), by=cg][y1], cg)))
        n1 <- t(sapply(n1, function(x) x[, V1]))
        n1[is.na(n1)] <- 0
        ## no. events, no. at risk at each time
        x1 <- x[, list("e"=sum(e), "n"=max(n)), by=t]
        data.table::setkey(x1, t)
        ## 'base variance'; used in all calcuations below
        bv1 <- x1[, e * (n - e) / (n^2 * (n - 1))]
        ## diagonal elements
        r1 <- bv1 * t(apply(n1, MARGIN=1, FUN=
                            function(i) (i * (sum(i) - i))))
        ## off-diagonal elements
        r2 <- bv1 * - t(apply(n1, MARGIN=1, FUN=
                              function(i) apply(utils::combn(i, g1 - 1), MARGIN=2, FUN=prod)))
        lt1 <- length(t1)
        res1 <- lapply(seq.int(lt1), FUN=
                       function(i){
                           res1 <- diag(r1[i, ])
                           res1[upper.tri(res1)] <- res1[lower.tri(res1)] <- r2[i, ]
                           return(res1)
                       })
        res1 <- as.array(unlist(res1))
        dim(res1) <- c(g1, g1, lt1)
        dimnames(res1) <- list(x[, unique(cg)], x[, unique(cg)], t1)
    }
    class(res1) <- c("COV", class(res1))
    data.table::setattr(x, "COV", res1)
    return(attr(x, "COV"))
}    
###---------------- 
COV.stratTne <- function(x, ...){
    return(lapply(x, FUN=COV))
}
###---------------- 
COV.numeric <- function(x, ..., n, n1){
  stopifnot(all(sapply(list(x, n, n1), is.numeric)))
  ## ensure all same length
  stopifnot(
    diff(range(sapply(list(x, n), length)))
    < .Machine$double.eps)
  ## no. of groups
  g1 <- ncol(n1)
  if (is.null(g1)) g1 <- 1L
  ## if 2 groups only
  if (g1==1){
    cov1 <- (n1 / n) * (1 - (n1 / n)) * ((n - x) / (n - 1)) * x
    return(cov1)
  }
  ## hold results
  a1 <- array(data=0,  dim=c(g1, g1, length(n)))
  ## diagonal elements
  for (i in seq_len(g1)){
    a1[i, i, ] <-
      (n1[, i] * (n - n1[, i]) * x * (n - x)) / (n^2 * (n - 1))
  }
  ## off-diagonal elements
  for (j in seq_len(g1)){
    for (k in 1:g1){
      if (j==k) next
      a1[j, k, ] <-
        - (n1[, j] * n1[, k] * x * (n - x)) / (n^2 * (n-1))
    }
  }
  dimnames(a1) <- list(1:g1, 1:g1, seq.int(length(x)))
  return(a1)
}
###
## ## diagonal elements
##     d1 <- n1 * (x[, n] - n1) * x[, e * (n - e) / (n^2 * (n - 1L))]
##     ## off-diagonals
##     od1 <- - apply(n1, 1, prod) * x[, e * (n - e) / (n^2 * (n-1L))]
## apply(n1)
## f1 <- function(x){
##     nc1 <- ncol(x)
##     c1 <- t(combn(nc1, (nc1 - 1)))
##     for (i in nrow(c1)){
##         return(x[, c1[1]] * x[, c1[2]])
##     }
## }
## ##     for (i in seq.int(nrow(x))){
##       diag(res1[, , i]) <- d1[seq.int(from=(i + 1 - g1),
##                                       to=g1*i)]
##       upper.tri(res1[, , i]) <- od1[seq.int(from=(i + 1 - g1),
##                                             to=g1*i)]
##       lower.tri(res1[, , i]) <- od1[seq.int(from=(i + 1 - g1),
##                                             to=g1*i)]
##     }   
##  [1] 5.327934e-05 5.406574e-05 5.486968e-05 5.569169e-05 5.653231e-05
##   [6] 5.739210e-05 1.156468e-04 1.192461e-04 6.200012e-05 6.298816e-05
##  [11] 6.400000e-05 6.503642e-05 1.311129e-04 6.830135e-05 6.944444e-05
##  [16] 1.400361e-04 7.305136e-05 7.431629e-05 7.561437e-05 7.694675e-05
##  [21] 7.831467e-05 1.580024e-04 8.264463e-05 8.416800e-05 8.573388e-05
##  [26] 8.734387e-05 8.899964e-05 9.070295e-05 1.831160e-04 9.611688e-05
##  [31] 9.802960e-05 1.000000e-04 1.020304e-04 1.041233e-04 1.062812e-04
##  [36] 1.085069e-04 1.108033e-04 1.131734e-04 1.156203e-04 1.181474e-04
##  [41] 0.000000e+00 1.234568e-04 1.262467e-04 1.291322e-04 1.321178e-04
##  [46] 1.352082e-04 1.384083e-04 1.417234e-04 1.451589e-04 1.487210e-04
##  [51] 1.524158e-04 1.562500e-04 1.602307e-04 1.643655e-04 3.328865e-04
##  [56] 1.777778e-04 1.826150e-04 1.876525e-04 1.929012e-04 1.983733e-04
##  [61] 2.040816e-04 2.100399e-04 2.162630e-04 2.227668e-04 2.295684e-04
##  [66] 2.366864e-04 0.000000e+00 2.519526e-04 2.601457e-04 2.687450e-04
##  [71] 2.777778e-04 2.872738e-04 2.972652e-04 3.077870e-04 3.188776e-04
##  [76] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
##  [81] 0.000000e+00 0.000000e+00 0.000000e+00 4.526935e-04 4.725898e-04
##  [86] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
##  [91] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
##  [96] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## [101] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## [106] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## [111] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## [116] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
## [121] 0.000000e+00 1.234568e-02 0.000000e+00 0.000000e+00 0.000000e+00
## [126] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00          NaN
## >
## (n1[i, ] * (sum(n1[i, ]) - n1[i, ]) * v1[i])
## r1 <- sapply(X=seq.int(130), FUN=function(i) n1[i, ] * (sum(n1[i, ]) - n1[i, ]) * v1[i])
## sapply(X=seq.int(130), FUN=function(i) n1[i, ] * (sum(n1[i, ]) - n1[i, ]))

## r1 * v1
## apply(n1, 1, FUN=function(i) sum(i))

      
## r1 <- sapply(X=seq.int(130), FUN=function(i) apply(combn(n1[i, ], g1 - 1), MARGIN=2, FUN=prod) * -v1[i])

## mapply(FUN=function(i, j, k) - (n1[i, j] * n1[i, k] * v1[i]),
##        i=1:130, j=1:3, k=1:3)

## v1 <- x[, e * (n - e) / (n^2 * (n - 1)), by=t][, V1]


## mapply(function(i, j)
##     n1[i, j] * (sum(n1[, j] - n1[i, j])) *  x[, e * (n - e) / (n^2 * (n - 1))],
## i=1:3, j=1:130)


## dim=c(g1, g1, length(t1)))
## dimnames(res1) <- list(1:g1, 1:g1, t1)
## n1 <- mapply(function(t1, g1) {
##     y <- rep(NaN, g1)
##     r1 <- x[t >= t1, max(ncg), by=cg]
##     y[r1[, cg]] <- r1[, V1]
##     y
## }, t1, g1)
## mapply(function(i, j)
##     n1[i, j] * (sum(n1[, j] - n1[i, j])) *  x[t==t1[i], e * (n - e) / (n^2 * (n - 1))],
## i=1:3, j=1:130)

## for (i in seq_len(g1)){
## res1[i, i, ] <-
##               (n1[, i] * (x[, n] - n1[, i]) * x[, e * (n - e) / (n^2 * (n - 1))])
##       }



## for(i in 1:length(t1)){

## for(j in 1:nrow(n1)){
##     diag(res1[j, j, i]) <- n1[, V1] * (sum(n1[, V1]) - n1[, V1]) * x[t==t1[i], e * (n - e) / (n^2 * (n - 1))]}}


## n1 <- rep(NaN, g1)
## r1 <-
    
## n1[r1[,cg]] <- 
## length(n1) <- 6
## lapply(n1, is.null)
## t2 <- dimnames(res1)[[3]][i]
## n1 <- x[t >= t2, max(ncg), by=cg][, V1]
## l1 <- sapply(t1, function(t1)
             
## y1 <- unlist(lapply(l1, function(n1)
## length(y1) <- g1 * length(t1)
## y1[is.na(y1)] <- NaN
## res1[1,,] <- y1
##              diag(res1)
##     }
##  (n1[, i] * (x[, n] - n1[, i]) * x[t==t2, e * (n - e) / (n^2 * (n - 1))])

## x[(t1 %in% t), ]
## sum(e * (n - e) / (n^2 * (n - 1)))]

## n1 <- x[, ncg, by=cg]
## res1 <- x[, sum(e * (n - e) / (n^2 * (n - 1))), by=t]
## x1 <- res1[1,][, V1]
## x[t<res1[, t], max(ncg), by=cg]

## x1* -(38*54)

##  38  54  45 = 137 


## na1 <- paste0(seq.int(g1))
## x[t==res2[, t], ncg, by=cg]

## res1[, (na1) := mapply(function(x) list(x / y^2),
##                               res1[, .SD, .SDcols=grep("Sv", names(res1))],
##                               res1[, .SD, .SDcols=grep("S[^v]", names(res1))])]
## split(x, x[, cg])


## res1[, sum(V1), by=t]

## x1 <- x[, sorted(unique(t))]
## rc1 <- data.table::CJ(seq.int(3), seq.int(3))

## data.table::setkey(n1, t)
## res1[, V1]
## x[CJ(ncg), (x[, t])]

## ## x1 <- res1[, sum(V1), by=list(t, cg)]
##   ## more than 2 groups?
##   if(g1 > 2){
## ### hold results
##       res1 <- array(data=0,  dim=c(g1, g1, nrow(x)))
##       n1 <- as.matrix(x[, .SD, .SDcols=grep("n_", names(x))])
##       ## diagonal elements
##       for (i in seq_len(g1)){
##           res1[i, i, ] <-
##               (n1[, i] * (x[, n] - n1[, i]) * x[, e * (n - e) / (n^2 * (n - 1))])
##       }
##     ## off-diagonal elements
##       for (j in seq_len(g1)){
##           for (k in 1:g1){
##               if (j==k) next
##               res1[j, k, ] <-
##                   - (n1[, j] * n1[, k] * x[, e * (n - e) / (n^2 * (n - 1))])
##           }
##       }
##       dimnames(res1) <- list(1:g1, 1:g1, x[, t])
##       class(res1) <- c("COV", class(res1))
##       data.table::setattr(x, "COV", res1)
##   }
##   return(attr(x, "COV"))
## }
