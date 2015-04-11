##' @name comp
##' @title Compare survival curves
##' @rdname comp
##' @keywords htest
##'
##' @export comp
##'
comp <- function(x, ...){
    UseMethod("comp")
    }
##' 
##' @include tne.R
##' @include covMatSurv.R
##' @include gastric.R
##' 
##' @param x A \code{survfit} or \code{coxph} object
##' @param ... Additional arguments
##' @param FHp \eqn{p} for Fleming-Harrington test
##' @param FHq \eqn{q} for Fleming-Harrington test
##' @param lim limit used for Renyi tests when generating supremum of
##' absolute value of Brownian motion
##' @param scores scores for tests for trend
##' @return A \code{list} with two elements, \code{tne} and \code{tests}.
##' \cr 
##' The first is a \code{data.table} with one row for each time at
##' which an event occurred.
##' \cr
##' Columns show time, no. at risk and no. events
##' (by stratum and overall).
##' \cr \cr
##' The second contains the tests, as a \code{list}.
##' \cr \cr
##' The first element is the log-rank family of tests.
##' \cr
##' The following additional tests depend on the no. of strata:
##' \cr
##' For a \code{survfit} or a \code{coxph} object with 2 strata,
##' these are the Supremum (Renyi) family of tests.
##' \cr
##' For a \code{survfit} or \code{coxph} object with at least 3 strata,
##' there are tests for trend.
##'
##' @details
##'  The \bold{log-rank} tests are given by the general expression:
##'  \deqn{ Q = \sum{ W_i (e_i - \hat{e}_i)}^T \sum{ W_i \hat{V_i} W_i^{-1}} \sum{ W_i (e_i - \hat{e}_i)} }{
##'         Q = [sum W(e-E)]^T [sum WVW]^-1 [sum W(e-E)] }
##' Where \eqn{W} is the weight, given below,  \eqn{e} is the no. of events,
##' \eqn{\hat{e}}{E} is the no. of expected events for that time and
##' \eqn{\hat{V}}{V} is the variance-covariance matrix given by \code{\link{covMatSurv}}.
##'
##' The sum is taken to the largest observed survival time
##' (i.e. censored observations are excluded).
##' \cr
##' If there are \eqn{K} groups, then \eqn{K-1} are selected (arbitrary).
##' Likewise the corresponding variance-covariance matrix is reduced to the
##' appropriate \eqn{K-1 \times K-1}{K-1 * K-1} dimensions.
##' \eqn{Q} is distributed as chi-square with \eqn{K-1} degrees of freedom.
##' \cr \cr
##' For 2 strata this simplifies to:
##'  \deqn{ Q = \frac{ \sum{ W_i [e1_i - n1_i (\frac{e_i}{n_i})]} }{ \sqrt{ \sum{ W_i^2 \frac{n1_i}{n_i} ( 1- \frac{n1_i}{n_i} ) ( \frac{n_i - e_i}{n_i-1}) e_i }}}}{
##' Q = SUM W [e1 - n1.(e/n)] / (SUM W^2.e1/e.(1-(n1/n)).(n-e/n-1).e)^0.5 }
##' Here \eqn{e} and \eqn{n} refer to
##' the no. events and no. at risk overall and
##' \eqn{e1} and \eqn{n1} refer to the no. events and no. at risk in group \eqn{1}.
##' \cr \cr
##' The weights are given as follows: \describe{
##'   \item{Log-rank}{weight = 1}
##'   \item{Gehan-Breslow generalized Wilcoxon}{weight = \eqn{n},
##' the no. at risk}
##'   \item{Tarone-Ware}{weight = \eqn{\sqrt{n}}{n^0.5}}
##'   \item{Peto-Peto}{weight = \eqn{\bar{S}(t)}{S(t)},
##' a  modified estimator of survival function given by
##' \deqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i+1}}}{
##'  S(t) = cumprod( 1- e/(n+1))}
##' }
##'   \item{modified Peto-Peto (by Andersen)}{
##'   weight = \eqn{\bar{S}(t) \frac{n}{n+1} }{S(t) n/n+1}}
##'   \item{Fleming-Harrington}{weight at \eqn{t_0 = 1} and thereafter is:
##'  \deqn{ \hat{S}(t_{i-1})^p [1-\hat{S}(t_{i-1})^q]}
##'  Here \eqn{\hat{S}} is the Kaplan-Meier (product-limit) estimator. Note that both \eqn{p} and \eqn{q} need to be \eqn{ \geq 0}{>=0}
##'  }
##' }
##' The \bold{Supremum (Renyi)} family of tests are designed to detect differences in survival curves which cross.
##' \cr
##' That is, an early difference in survival in favor of one group
##' is balanced by a later reversal.
##' \cr
##' The same weights as above are used.
##' \cr
##' They are calculated by finding
##' \deqn{ Z(t_i) = \sum_{t_k \leq t_i} W(t_k)[e1_k - n1_k\frac{e_k}{n_k}], \quad i=1,2,...,k}{
##'        Z(t[i]) = SUM W(t[k]) [ e1[k] - n1[k]e[k]/n[k] ]}
##' (which is similar to the numerator used to find \eqn{Q}
##' in the log-rank test for 2 groups above).
##' \cr
##' and it's variance:
##' \deqn{ \sigma^2(\tau) = \sum_{t_k \leq \tau} W(t_k)^2 \frac{n1_k n2_k (n_k-e_k) e_k}{n_k^2 (n_k-1)} }{
##'       simga^2(tau) = SUM(k=1,2,...,tau) W(t[k]) [ n1[k].n2[k].(n[k]-e[k]).e[k] / n[k]^2.(n[k]-1) ] }
##' where \eqn{\tau}{tau} is the largest \eqn{t}
##' where both groups have at least one subject at risk.
##' \cr \cr
##' Then calculate:
##' \deqn{ Q = \frac{ \sup{|Z(t)|}}{\sigma(\tau)}, \quad t<\tau }{
##' Q = sup( |Z(t)| ) / sigma(tau), t<tau}
##' When the null hypothesis is true,
##' the distribution of \eqn{Q} is approximately
##' \deqn{Q \sim \sup{|B(x)|, \quad 0 \leq x \leq 1}}{
##' Q ~ sup( |B(x)|, 0 <= x <= 1)}
##' And for a standard Brownian motion (Wiener) process:
##' \deqn{ Pr[\sup|B(t)|>x] = 1-\frac{4}{\pi} \sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} \exp{ \frac{-\pi^2(2k+1)^2}{8x^2}}}{
##'        Pr[sup|B(t)|>x] = 1 - 4/pi SUM (-1)^k/2k+1 exp(-pi^2 (2k+1)^2/x^2)}
##' \bold{Tests for trend} are designed to detect ordered differences in survival curves.
##' \cr
##' That is, for at least one group:
##' \deqn{S_1(t) \geq S_2(t) \geq ... \geq S_K(t) \quad t \leq \tau}{
##'       S1(t) >= S2(t) >= ... >= SK(t) for t <= tau}
##' where \eqn{\tau}{tau} is the largest \eqn{t} where all groups have at least one subject at risk.
##' The null hypothesis is that
##' \deqn{S_1(t) = S_2(t) = ... = S_K(t) \quad t \leq \tau}{
##'       S1(t) = S2(t) = ... = SK(t) for t <= tau}
##' Scores used to construct the test are typically \eqn{s = 1,2,...,K},
##' but may be given as a vector representing a numeric characteristic of the group.
##' \cr
##' They are calculated by finding
##' \deqn{ Z_j(t_i) = \sum_{t_i \leq \tau} W(t_i)[e_{ji} - n_{ji} \frac{e_i}{n_i}], \quad j=1,2,...,K}{
##'       Z[t(i)] = SUM W[t(i)] [e[j](i) - n[j](i).e(i)/n(i) ]}
##' The test statistic is
##' \deqn{Z = \frac{ \sum_{j=1}^K s_jZ_j(\tau)}{\sqrt{\sum_{j=1}^K \sum_{g=1}^K s_js_g \sigma_{jg}}} }{
##'       Z = SUM(j=1...K) s[j]Z[j] / SUM(j=1..K) SUM(g=1..K) s[j]s[g]sigma[jg]}
##' where \eqn{\sigma}{sigma} is the the appropriate element in the
##' variance-covariance matrix (as in \code{\link{covMatSurv}}).
##' \cr
##' If ordering is present, the statistic \eqn{Z} will be greater than the upper \eqn{\alpha}{alpha}th
##' percentile of a standard normal distribution.
##' 
##' @note Regarding the Fleming-Harrington weights: \itemize{
##' \item \eqn{p = q = 0} gives the log-rank test, i.e. \eqn{W=1}
##' \item \eqn{p=1, q=0} gives a version of the Mann-Whitney-Wilcoxon test
##' (tests if populations distributions are identical)
##' \item \eqn{p=0, q>0} gives more weight to differences later on
##' \item \eqn{p>0, q=0} gives more weight to differences early on
##' }
##' The example using \code{alloauto} data illustrates this.
##' Here the log-rank statistic
##' has a p-value of  around 0.5
##' as the late advantage of allogenic transplants
##' is offset by the high early mortality. However using
##' Fleming-Harrington weights of \eqn{p=0, q=0.5},
##' emphasising differences later in time, gives a p-value of 0.04.
##'
##' @examples
##' ### 2 curves
##' data(kidney,package="KMsurv")
##' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' comp(s1)
##' ### 3 curves
##' data(bmt, package="KMsurv")
##' comp(survfit(Surv(time=t2, event=d3) ~ group, data=bmt))
##' ### see effect of F-H test
##' data(alloauto, package="KMsurv")
##' s3 <- survfit(Surv(time, delta) ~ type, data=alloauto)
##' comp(s3, FHp=0, FHq=1)
##' ### see trend tests
##' data(larynx, package="KMsurv")
##' s4 <- survfit(Surv(time, delta) ~ stage, data=larynx)
##' comp(s4)
##' ### Renyi tests
##' data("gastric", package="survMisc")
##' s5 <- survfit(Surv(time, event) ~ group, data=gastric)
##' comp(s5)
##' 
##' @references Gehan A.
##' A Generalized Wilcoxon Test for Comparing Arbitrarily
##' Singly-Censored Samples.
##' Biometrika 1965 Jun. 52(1/2):203--23.
##' \href{http://www.jstor.org/stable/2333825}{JSTOR}
##' @references Tarone RE, Ware J 1977
##' On Distribution-Free Tests for Equality of Survival Distributions.
##' \emph{Biometrika};\bold{64}(1):156--60.
##' \href{http://www.jstor.org/stable/2335790}{JSTOR}
##' @references Peto R, Peto J 1972
##' Asymptotically Efficient Rank Invariant Test Procedures.
##' \emph{J Royal Statistical Society} \bold{135}(2):186--207.
##' \href{http://www.jstor.org/stable/2344317}{JSTOR}
##' @references Fleming TR, Harrington DP, O'Sullivan M 1987
##' Supremum Versions of the Log-Rank and Generalized Wilcoxon Statistics.
##' \emph{J  American Statistical Association} \bold{82}(397):312--20.
##' \href{http://www.jstor.org/stable/2289169}{JSTOR}
##' @references Billingsly P 1999
##' \emph{Convergence of Probability Measures.}
##' New York: John Wiley & Sons.
##' \href{http://dx.doi.org/10.1002/9780470316962}{Wiley (paywall)}
##' @references Examples are from
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Examples 7.2, 7.4, 7.5, 7.6, 7.9, pp 210-225.
##'
##' @rdname comp
##' @aliases comp.survfit
##' @method comp survfit
##' @export
comp.survfit <- function(x, ..., FHp=1, FHq=1, lim=1e4, scores=NULL) {
    if (!class(x)=="survfit") stop("Only applies to object of class 'survfit'")
    m1 <- "Values for p and q for Fleming-Harrington tests must be >=0"
    if (FHp<0 | FHq<0) stop(m1)
###
    m1 <- tne(x, what="all", eventsOnly=TRUE)
###
### 2 groups only:
    if (length(x$strata)==2){
        res1 <- comp2Surv(n=m1$n, e=m1$e,
                          n1=m1[, grep("n_", colnames(m1))[1], with=FALSE],
                          e1=m1[, grep("e_", colnames(m1))[1], with=FALSE],
                          FHp=FHp, FHq=FHq, lim=lim)
    } else {
        res1 <- compNSurv(t=m1$t, n=m1$n, e=m1$e,
                          n1=as.matrix(m1[, grep("n_",colnames(m1)),
                          with=FALSE]),
                          e1=as.matrix(m1[, grep("e_",colnames(m1)),
                          with=FALSE]),
                          FHp=FHp, FHq=FHq, scores=scores)
    }
###
    res <- list(
        tne = m1,
        tests = res1)
    return(res)
}
##' @rdname comp
##' @aliases comp.coxph
##' @method comp coxph
##' @export
##' 
##' @examples
##' c1 <- coxph(Surv(time=time, event=delta) ~ type, data=kidney )
##' comp(c1)
##' 
comp.coxph <- function(x, ..., FHp=1, FHq=1, scores=NULL, lim=1e4){
     if (!class(x)=="coxph") stop("Only applies to object of class 'coxph'")
     if (FHp<0|FHq<0) stop("Values for p and q for Fleming-Harrington tests must be >=0")
###
    m1 <- tne.coxph(x, what="all", eventsOnly=TRUE)
###
### 2 groups only:
    if (ncol(m1)==7){
        res1 <- comp2Surv(n=m1$n, e=m1$e,
                          n1=m1[, grep("n_", colnames(m1))[1], with=FALSE],
                          e1=m1[, grep("e_", colnames(m1))[1], with=FALSE],
                          FHp=FHp, FHq=FHq, round1=5)
    } else {
        res1 <- compNSurv(t=m1$t, n=m1$n, e=m1$e,
                          n1=as.matrix(m1[, grep("n_",colnames(m1)),
                          with=FALSE]),
                          e1=as.matrix(m1[, grep("e_",colnames(m1)),
                          with=FALSE]),
                          FHp=FHp, FHq=FHq, scores=scores)
    }
     res <- list(
         tne = m1,
         tests = res1)
     return(res)
}
###
###----------------------------------------
###
### compare 2 survival curves
###
comp2Surv <- function(n, e, n1, e1,
                      FHp=1,FHq=1,
                      lim=1e4,round1=5){
    if (FHp<0|FHq<0) stop
    ("Values for p and q for Fleming-Harrington tests must be >=0")
    stopifnot( all.equal(length(n),length(e),length(n1),length(e1)) )
    if(!isTRUE(all(vapply(c(n,e,n1,e1),
                          FUN=is.numeric,
                          FUN.VALUE=TRUE)==TRUE))) stop
    ("All vectors must be numeric")
###
### make observed - expected for group 1
    eME1 <- e1-(n1*e/n)
### variance (2 groups)
    var1 <- as.matrix( (n1/n)*(1-(n1/n))*((n-e)/(n-1))*e )
### remove NaNs -
### first or last no.s may not be possible to calculate
    if (any(is.nan(var1))){
### index
        in1 <- which(is.nan(var1))
        var1 <- var1[-in1]
        eME1 <- eME1[-in1]
        e <- e[-in1]
        n <- n[-in1]
        e1 <- e1[-in1]
        n1 <- n1[-in1]
    }
### display chisq, degrees of freedom and rounded result
    dis1 <- function(chi1, df1=1, rounded=round1){
        c(chi1,df1, round(1-stats::pchisq(chi1,df1), digits=rounded))
    }
###
### WEIGHTS
### log-rank, weight = 1
    chi1 <- sum(eME1)^2/sum(var1)
    LR1 <- dis1(chi1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    chi1 <- sum(n*eME1)^2/sum(n^2*var1)
    GB1 <- dis1(chi1)
### Tarone-Ware, weight = sqrt(n)
    chi1 <- sum(sqrt(n)*eME1)^2/ sum(sqrt(n)^2*var1)
    TW1 <- dis1(chi1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (1- (e/(n+1))) )
    chi1 <- sum(S1*eME1)^2/sum(S1^2*var1)
    PP1 <- dis1(chi1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    S2 <- (S1*n)/(n+1)
    chi1 <- sum(S2*eME1)^2/sum(S2^2*var1)
    mPP1 <- dis1(chi1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    S3 <- cumprod( (n-e)/n )
    S3 <- c(1,S3[1:(length(S3)-1)])
    FHw <- S3^FHp*((1-S3)^FHq)
    chi1 <- sum(FHw*eME1)^2/sum(FHw^2*var1)
    FH1 <- dis1(chi1)
### results
    FHn <- paste("Flem~-Harr~ with p=",FHp,", q=",FHq,sep="")
    res1 <- rbind(LR1,GB1,TW1,PP1,mPP1,FH1)
    dimnames(res1) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHn), c(
            "ChiSq", "df", "p"))
###
### Renyi statistics (analagous to 2-sample Kolmogorov-Smirnov test)
### aka supremum tests
###
### Probability of supremum of absolute value of
### standard Brownian motion process B(t)
    probSupBr <- function(y, limit=lim){
        k <- seq(from=0, to=limit, by=1)
        res <- 1-( (4/pi)* sum(  ((-1)^k) /(2*k+1) * exp( (-(pi^2)*(2*k+1)^2) / (8*y^2) )  ) )
        return(res)
    }
    Z1 <- cumsum(e1-(n1*(e/n)))
### display chisq, degrees of freedom and rounded result
    dis2 <- function(Q1,rounded=round1){
        c(Q1,round(probSupBr(Q1),digits=rounded))
      }
    Q1 <- max(abs(Z1))/sqrt(sum(var1))
### log-rank weights
    RenLR1 <- dis2(Q1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    Z1 <- cumsum( n * (e1-(n1 * (e / n))) )
    Q1 <- max(abs(Z1)) / sqrt(sum(n^2*var1))
    RenGB1 <- dis2(Q1)
### Tarone-Ware, weight = sqrt(n)
    Z1 <- cumsum( sqrt(n)* (e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(n*var1))
    RenTW1 <- dis2(Q1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (n-e)/n )
    Z1 <- cumsum( S1 *(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(S1^2*var1))
    RenPP1 <- dis2(Q1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    Zfun <- function(i)
    Z1 <- cumsum( S2*(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(S2^2*var1))
    RenmPP1 <- dis2(Q1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    Zfun <- function(i)
    Z1 <- cumsum( FHw *(e1-(n1*(e/n))) )
    Q1 <- max(abs(Z1))/sqrt(sum(FHw^2*var1))
    RenFH1 <- dis2(Q1)
### result
    FHnR <- paste("Renyi ",FHn,sep="")
    res2 <- rbind(RenLR1,RenGB1,RenTW1,RenPP1,RenmPP1,RenFH1)
    dimnames(res2) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHnR), c(
            "Q","p"))
### final
    return( list(lrTests=res1,
                 supTests=res2))
}
###
###----------------------------------------
###
### compare N survival curves
###
compNSurv <- function (t, n, e, n1, e1,
                       FHp=FHp, FHq=FHq,
                       scores=NULL, round1=5){
### events observed minus expected
    eME <- e1-(n1*e/n)
### degrees of freedom
    df1 <- ifelse(is.null(dim(eME)),1,ncol(eME))-1
### make covariance matrix (n groups)
    v1 <- covMatSurv(t, n, e, n1)
### remove NaNs -
### last no.s may not be possible to calculate
    if (any(is.nan(v1))){
### index
        in1 <- which(is.nan(v1))
        dim1 <- dim(v1)
        v1 <- v1[-in1]
        dim(v1) <- c(dim1[1], dim1[2], (dim1[3]-1))
        eME <- eME[-nrow(eME), ]
        e <- e[-length(e)]
        n <- n[-length(n)]
        e1 <- e1[-nrow(e1), ]
        n1 <- n1[-nrow(n1), ]
    }
    cov1 <- rowSums(v1, dims=2)
### display chisq, degrees of freedom and rounded result
    dis <- function(chi1, df=df1, rounded=round1){
        c(chi1, df1, round(1-stats::pchisq(chi1, df), digits=rounded))
    }
### WEIGHTS
### log-rank, weight = 1
    eME1 <- colSums(eME)
    chi1 <- eME1[1:df1] %*% solve(cov1[1:df1,1:df1]) %*% matrix((eME1[1:df1]),nrow=df1,ncol=1)
    LR1 <- dis(chi1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    v2 <- sweep(v1, MARGIN=3, STATS=n^2, FUN="*")
    covG <- rowSums(v2, dims=2)
    eMEG <- colSums(sweep(eME,1,n,"*"))
    chi1 <- eMEG[1:df1] %*% solve(covG[1:df1,1:(df1)]) %*% matrix((eMEG[1:df1]),df1,1)
    GB1 <- dis(chi1)
### Tarone-Ware, weight = sqrt(n)
    v2 <- sweep(v1,MARGIN=3,STATS=n,FUN="*")
    covTw <- rowSums(v2,dims=2)
    eMETw <- colSums(sweep(eME,1,sqrt(n),"*"))
    chi1 <- eMETw[1:df1] %*% solve(covTw[1:df1,1:df1]) %*% matrix((eMETw[1:df1]),df1,1)
    TW1 <- dis(chi1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    S1 <- cumprod( (1- (e/(n+1))) )
    v2 <- sweep(v1,MARGIN=3,STATS=S1^2,FUN="*")
    covPP <- rowSums(v2,dims=2)
    eMEPP <- colSums(sweep(eME,1,S1,"*"))
    chi1 <- eMEPP[1:df1] %*% solve(covPP[1:df1,1:df1]) %*% matrix((eMEPP[1:df1]),df1,1)
    PP1 <- dis(chi1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    S2 <- (S1*n)/(n+1)
    v2 <- sweep(v1,MARGIN=3,STATS=S2^2,FUN="*")
    covmPP <- rowSums(v2,dims=2)
    eMEmPP <- colSums(sweep(eME,1,S2,"*"))
    chi1 <- eMEmPP[1:df1] %*% solve(covmPP[1:df1,1:df1]) %*% matrix((eMEmPP[1:df1]),df1,1)
    mPP1 <- dis(chi1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    S1 <- cumprod( (n-e)/n )
    S1 <- c(1,S1[1:(length(S1)-1)])
    FH1 <- S1^FHp*((1-S1)^FHq)
    v2 <- sweep(v1,MARGIN=3,STATS=FH1^2,FUN="*")
    covFH <- rowSums(v2,dims=2)
    eMEFH <- colSums(sweep(eME,1,FH1,"*"))
    chi1 <- eMEFH[1:df1] %*% solve(covFH[1:df1,1:df1]) %*% matrix((eMEFH[1:df1]),df1,1)
    FH1 <- dis(chi1)
### results
    FHn <- paste("Flem~-Harr~ with p=",FHp,", q=",FHq,sep="")
    res1 <- rbind(LR1,GB1,TW1,PP1,mPP1,FH1)
    dimnames(res1) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHn), c(
            "ChiSq", "df", "p"))
###
### Trend tests
    if (is.null(scores)) scores <- 1:(df1+1)
### no. predictors (from degrees of freedom, above)
    lp1 <- df1+1
### display chisq, degrees of freedom and rounded result
    dis2 <- function(Z1, rounded=round1){
        c(Z1, round(1-stats::pnorm(Z1), digits=rounded))
      }
### log-rank weights
    sfun <- function(i, j) scores[i] * scores[j] * cov1[i, j]
    p2 <- as.vector(sapply(seq(lp1), function (x) rep(x, lp1)))
    s1 <- sum(mapply(sfun, rep(seq(lp1), lp1), p2))
### the above is shorthand for:
### get denominator for formula
    .getDenom <- function(lp, cov){
        s1 <- NULL
        for (i in seq(scores)){
            for (j in seq(scores)){
                s1 <- c(s1, cov[i, j] * scores[i] * scores[j])
            }
        }
        return(s1)
    }
###
### eME1 = obs - exp, from log-rank test above
    Z1 <- sum(scores*eME1)/sqrt(s1)
    LRT1 <- dis2(Z1)
### Gehan-Breslow generalized Wilcoxon, weight = n
    sfun <- function(i,j) scores[i]*scores[j]*covG[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEG)/sqrt(s1)
    GBT1 <- dis2(Z1)
### Tarone-Ware, weight = sqrt(n)
    sfun <- function(i,j) scores[i]*scores[j]*covTw[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMETw)/sqrt(s1)
    TWT1 <- dis2(Z1)
### Peto-Peto, weight = S(t) = modified estimator of survival function
    sfun <- function(i,j) scores[i]*scores[j]*covPP[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEPP)/sqrt(s1)
    PPT1 <- dis2(Z1)
### modified Peto-Peto (by Andersen), weight = S(t)n/n+1
    sfun <- function(i,j) scores[i]*scores[j]*covmPP[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEmPP)/sqrt(s1)
    mPPT1 <- dis2(Z1)
### Fleming-Harrington; weight 1st element is 1 as depends on [i-1]
    sfun <- function(i,j) scores[i]*scores[j]*covFH[i,j]
    s1 <- sum(mapply(sfun, rep(seq(lp1),lp1), p2))
    Z1 <- sum(scores*eMEFH)/sqrt(s1)
    FHT1 <- dis2(Z1)
### result
    FHnT <- paste("Trend F-H with p=",FHp,", q=",FHq,sep="")
    res2 <- rbind(LRT1,GBT1,TWT1,PPT1,mPPT1,FHT1)
    dimnames(res2) <- list(c(
        "Log-rank",
        "Gehan-Breslow (mod~ Wilcoxon)",
        "Tarone-Ware",
        "Peto-Peto",
        "Mod~ Peto-Peto (Andersen)",
        FHnT), c(
            "Z", "p"))
### final
    return( list(lrTests=res1,
                 trendTests=res2,
                 scores=scores))
}
