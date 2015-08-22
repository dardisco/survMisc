#' @name comp
#' @title Compare survival curves
#' @rdname comp
#' @keywords htest
#'
#' @export comp
#'
comp <- function(x, ...) UseMethod("comp")
#'
#' @include tne.R
#' @include print.R
#' @include COV.R
#' @include predict.R
#'
#' @param x A \code{tne} object
#' @param p \eqn{p} for Fleming-Harrington test
#' @param q \eqn{q} for Fleming-Harrington test
#' @param scores scores for tests for trend
#' 
#' @return The \code{tne} object is given
#'  additional \code{attributes}.
#'  \cr
#' The following are always added:
#' \item{lrt}{The \bold{l}og-\bold{r}ank family of \bold{t}ests}
#' \item{wlr}{The \bold{w}eights used in calculating the \bold{l}og-\bold{r}ank tests}
#' An additional item depends on the number of covariate groups.
#' \cr
#' If this is \eqn{=2}:
#' \item{sup}{The \bold{sup}remum or Renyi family of tests}
#' and if this is \eqn{>2}:
#' \item{tft}{Tests for trend}
#' 
#' @details
#' The \bold{log-rank} tests are formed from the following elements,
#' with values for each time where there is at least one event:
#' \itemize{
#'  \item \eqn{W_i}{W[i]}, the weights, given below.
#'  \item \eqn{e_i}{e[i]}, the number of events (per time).
#'  \item \eqn{\hat{e_i}}{P[i]}, the number of \emph{predicted} events,
#'                               given by \code{\link{predict}}.
#'  \item \eqn{COV_i}{COV[, , i]}, the covariance matrix for time \eqn{i},
#'                                 given by \code{\link{COV}}.
#' }
#' It is calculated as:
#'  \deqn{Q_i = \sum{W_i (e_i - \hat{e}_i)}^T
#'              \sum{W_i \hat{COV_i} W_i^{-1}}
#'              \sum{W_i (e_i - \hat{e}_i)}}{
#'       Q[i] = sum(W[i] * (e[i] - P[i]))^T *
#'              sum(W[i] * COV[, , i] * W[i])^-1 *
#'              sum(W[i] * (e[i] - P[i]))}
#' 
#' If there are \eqn{K} groups, then \eqn{K-1} are selected (arbitrary).
#'  \cr
#' Likewise the corresponding variance-covariance matrix is reduced to the
#' appropriate \eqn{K-1 \times K-1}{K-1 * K-1} dimensions.
#'  \cr
#' \eqn{Q} is distributed as chi-square with \eqn{K-1} degrees of freedom.
#' \cr \cr
#' For \eqn{2} covariate groups, we can use:
#' \itemize{
#'  \item \eqn{e_i}{e[i]} the number of events (per time).
#'  \item \eqn{n_i}{e[i]} the number at risk overall.
#'  \item \eqn{e1_i}{e1[i]} the number of events in group \eqn{1}.
#'  \item \eqn{n1_i}{n1[i]} the number at risk in group \eqn{1}.
#' } 
#' Then:
#'  \deqn{Q = \frac{\sum{W_i [e1_i - n1_i (\frac{e_i}{n_i})]} }{
#'                 \sqrt{\sum{W_i^2 \frac{n1_i}{n_i}
#'                            (1 - \frac{n1_i}{n_i})
#'                            (\frac{n_i - e_i}{n_i - 1}) e_i }}}}{
#'        Q = sum(W[i] * (e1[i] - n1[i] * e[i] / n[i])) /
#'           sqrt(sum(W[i]^2 * e1[i] / e[i] * (1 - n1[i] / n[i]) * (n[i] - e[i] / (n[i] - 1)) *e[i]))}
#' The weights are given as follows:
#' \tabular{cll}{
#' \eqn{1} \tab log-rank \tab \cr
#' \eqn{n_i}{n[i]} \tab Gehan-Breslow generalized Wilcoxon \tab \cr
#' \eqn{\sqrt{n_i}}{sqrt(n[i])} \tab Tarone-Ware \tab \cr
#' \eqn{S1_i}{S1[i]} \tab Peto-Peto's modified survival estimate \tab
#'                    \eqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i + 1}}}{
#'                             S(t) = cumprod( 1- e/(n+1))} \cr
#' \eqn{S2_i}{S2[i]} \tab modified Peto-Peto (by Andersen) \tab
#'                    \eqn{\bar{S}(t)=\prod{1 - \frac{e_i}{n_i + 1}}}{
#'                              S(t) = cumprod(1 - e/(n + 1))} \cr
#' \eqn{FH_i}{FH[i]} \tab Fleming-Harrington \tab
#'                   weight at \eqn{t_0 = 1} and thereafter is:
#'                   \eqn{\hat{S}(t_{i-1})^p [1-\hat{S}(t_{i-1})^q]}{
#'                        S(t[i - 1])^p * (1 - S(t)[i - 1]^q)}
#' }
#' Here, \eqn{\hat{S}(t)}{S(t)} is the Kaplan-Meier (product-limit) estimator.
#' \cr
#' Note that both \eqn{p} and \eqn{q} need to be \eqn{\geq 0}{>=0}
#' \cr \cr
#' The \bold{supremum (Renyi)} family of tests are designed
#'  to detect differences in survival curves which cross.
#' \cr
#' That is, an early difference in survival in favor of one group
#'  is balanced by a later reversal.
#' \cr
#' The same weights as above are used.
#' \cr
#' They are calculated by finding
#' \deqn{Z(t_i) = \sum_{t_k \leq t_i} W(t_k)[e1_k - n1_k\frac{e_k}{n_k}], \quad i=1,2,...,k}{
#'       Z(t[i]) = SUM W(t[k]) [ e1[k] - n1[k]e[k]/n[k] ]}
#' (which is similar to the numerator used to find \eqn{Q}
#' in the log-rank test for 2 groups above).
#' \cr
#' and it's variance:
#' \deqn{\sigma^2(\tau) = \sum_{t_k \leq \tau} W(t_k)^2 \frac{n1_k n2_k (n_k-e_k) e_k}{n_k^2 (n_k-1)} }{
#'       simga^2(tau) = SUM(k=1,2,...,tau) W(t[k]) [ n1[k].n2[k].(n[k]-e[k]).e[k] / n[k]^2.(n[k]-1) ] }
#' where \eqn{\tau}{tau} is the largest \eqn{t}
#' where both groups have at least one subject at risk.
#' \cr \cr
#' Then calculate:
#' \deqn{ Q = \frac{ \sup{|Z(t)|}}{\sigma(\tau)}, \quad t<\tau }{
#' Q = sup( |Z(t)| ) / sigma(tau), t<tau}
#' When the null hypothesis is true,
#' the distribution of \eqn{Q} is approximately
#'  \deqn{Q \sim \sup{|B(x)|, \quad 0 \leq x \leq 1}}{
#'        Q ~ sup( |B(x)|, 0 <= x <= 1)}
#' And for a standard Brownian motion (Wiener) process:
#'  \deqn{Pr[\sup|B(t)|>x] = 1 - \frac{4}{\pi}
#'                           \sum_{k=0}^{\infty}
#'                           \frac{(- 1)^k}{2k + 1} \exp{\frac{-\pi^2(2k + 1)^2}{8x^2}}}{
#'        Pr[sup|B(t)|>x] = 1 - 4/pi sum((-1)^k / (2 * k + 1) * exp(-pi^2 (2k + 1)^2 / x^2))}
#' \bold{Tests for trend} are designed to detect ordered differences in survival curves.
#' \cr
#' That is, for at least one group:
#' \deqn{S_1(t) \geq S_2(t) \geq ... \geq S_K(t) \quad t \leq \tau}{
#'       S1(t) >= S2(t) >= ... >= SK(t) for t <= tau}
#' where \eqn{\tau}{tau} is the largest \eqn{t} where all groups have at least one subject at risk.
#' The null hypothesis is that
#' \deqn{S_1(t) = S_2(t) = ... = S_K(t) \quad t \leq \tau}{
#'       S1(t) = S2(t) = ... = SK(t) for t <= tau}
#' Scores used to construct the test are typically \eqn{s = 1,2,...,K},
#' but may be given as a vector representing a numeric characteristic of the group.
#' \cr
#' They are calculated by finding
#' \deqn{ Z_j(t_i) = \sum_{t_i \leq \tau} W(t_i)[e_{ji} - n_{ji} \frac{e_i}{n_i}], \quad j=1,2,...,K}{
#'       Z[t(i)] = SUM W[t(i)] [e[j](i) - n[j](i).e(i)/n(i) ]}
#' The test statistic is
#' \deqn{Z = \frac{ \sum_{j=1}^K s_jZ_j(\tau)}{\sqrt{\sum_{j=1}^K \sum_{g=1}^K s_js_g \sigma_{jg}}} }{
#'       Z = SUM(j=1...K) s[j]Z[j] / SUM(j=1..K) SUM(g=1..K) s[j]s[g]sigma[jg]}
#' where \eqn{\sigma}{sigma} is the the appropriate element in the
#' variance-covariance matrix (as in \code{\link{covMatSurv}}).
#' \cr
#' If ordering is present, the statistic \eqn{Z} will be greater than the upper \eqn{\alpha}{alpha}th
#' percentile of a standard normal distribution.
#'
#' @note Regarding the Fleming-Harrington weights: \itemize{
#' \item \eqn{p = q = 0} gives the log-rank test, i.e. \eqn{W=1}
#' \item \eqn{p=1, q=0} gives a version of the Mann-Whitney-Wilcoxon test
#' (tests if populations distributions are identical)
#' \item \eqn{p=0, q>0} gives more weight to differences later on
#' \item \eqn{p>0, q=0} gives more weight to differences early on
#' }
#' The example using \code{alloauto} data illustrates this.
#' Here the log-rank statistic
#' has a p-value of  around 0.5
#' as the late advantage of allogenic transplants
#' is offset by the high early mortality. However using
#' Fleming-Harrington weights of \eqn{p=0, q=0.5},
#' emphasising differences later in time, gives a p-value of 0.04.
#'
#' @references Gehan A.
#' A Generalized Wilcoxon Test for Comparing Arbitrarily
#' Singly-Censored Samples.
#' Biometrika 1965 Jun. 52(1/2):203--23.
#' \href{http://www.jstor.org/stable/2333825}{JSTOR}
#' @references Tarone RE, Ware J 1977
#' On Distribution-Free Tests for Equality of Survival Distributions.
#' \emph{Biometrika};\bold{64}(1):156--60.
#' \href{http://www.jstor.org/stable/2335790}{JSTOR}
#' @references Peto R, Peto J 1972
#' Asymptotically Efficient Rank Invariant Test Procedures.
#' \emph{J Royal Statistical Society} \bold{135}(2):186--207.
#' \href{http://www.jstor.org/stable/2344317}{JSTOR}
#' @references Fleming TR, Harrington DP, O'Sullivan M 1987
#' Supremum Versions of the Log-Rank and Generalized Wilcoxon Statistics.
#' \emph{J  American Statistical Association} \bold{82}(397):312--20.
#' \href{http://www.jstor.org/stable/2289169}{JSTOR}
#' @references Billingsly P 1999
#' \emph{Convergence of Probability Measures.}
#' New York: John Wiley & Sons.
#' \href{http://dx.doi.org/10.1002/9780470316962}{Wiley (paywall)}
#' 
#' @examples
#' ## Two covariate groups
#' ## K&M 2nd ed. Example 7.2, Table 7.2, pp 209--210.
#' data("kidney", package="KMsurv")
#' tne1 <- tne(Surv(time=time, event=delta) ~ type, data=kidney, shape="long")
#' comp(tne1)
#' ## supremum (Renyi) test; two-sided; two covariate groups
#' ## K&M 2nd ed. Example 7.9, pp 223--226.
#' data("gastric", package="survMisc")
#' tne1 <- tne(Surv(time, event) ~ group, data=gastric)
#' comp(tne1)
#' ## Three covariate groups
#' ## K&M 2nd ed. Example 7.4, pp 212-214.
#' data("bmt", package="KMsurv")
#' tne1 <- tne(Surv(time=t2, event=d3) ~ group, data=bmt)
#' comp(tne1)
#' ## Tests for trend
#' ## K&M 2nd ed. Example 7.6, pp 217-218.
#' data("larynx", package="KMsurv")
#' tne1 <- tne(Surv(time, delta) ~ stage, data=larynx)
#' comp(tne1, scores=c(4, 2, 3, 1))
#' attr(tne1, "tft")
#' ### see effect of FH test
#' data("alloauto", package="KMsurv")
#' tne1 <- tne(Surv(time, delta) ~ type, data=alloauto)
#' comp(tne1, p=c(0, 1), q=c(1, 1))
#'
#'
#' @rdname comp
#' @aliases comp.tne
#' @method comp tne
#' @export
comp.tne <- function(x,
                     p=1,
                     q=1,
                     scores=seq.int(attr(x, "ncg"))){
    stopifnot(inherits(x, "tne"))
    stopifnot(attr(x, "ncg") >= 2)
    stopifnot(all(c(p, q) >= 0))
    stopifnot(all(c(p, q) <= 1))
    stopifnot(length(p)==length(q))
    ## number of F-H tests
    fh1 <- length(p)
    ## WEIGHTS
    data.table::setkey(x, t)
    wt1 <- data.table::data.table(matrix(rep(1, x[, length(unique(t))] * (5 + fh1)),
                                         ncol=5 + fh1))
    FHn <- paste("FH_p=", p, "_q=", q, sep="")
    n1 <- c("lrt", "n", "sqrtN", "S1", "S2", FHn)
    data.table::setnames(wt1, n1)
    ## Gehan-Breslow generalized Wilcoxon, weight = n
    data.table::set(wt1, j=2L, value=x[, max(n), by=t][, V1])
    ## Tarone-Ware, weight = sqrt(n)
    data.table::set(wt1, j=3L, value=wt1[, sqrt(.SD), .SDcols=2])
    ## Peto-Peto, weight = S(t) = modified estimator of survival function
    data.table::set(wt1, j=4L, value=x[, cumprod((1 - (sum(e) / (max(n) + 1)))), by=t][, V1])
    ## modified Peto-Peto (by Andersen), weight = S(t)n / n+1
    data.table::set(wt1, j=5L, value=wt1[, S1] * x[, max(n) / (max(n) + 1), by=t][, V1])
    ## Fleming-Harrington
    S3 <- x[, cumprod((max(n) - sum(e)) / max(n)), by=t][, V1]
    ## weight of first 1st element is 1 as depends on [i-1]
    S3 <- c(1, S3[seq.int(length(S3) - 1L)])
    ## Fleming-Harrington
    ## easier to read as a loop here
    for (j1 in seq.int(fh1)){
        data.table::set(wt1, j=5L + j1,
                        value=S3^p[j1] * ((1 - S3)^q[j1]))
    }
    ## hold results for log-rank tests
    res <- data.table::data.table(matrix(0, nrow=5 + fh1, ncol=8))
    data.table::setnames(res, c("W", "Q", "Var", "Z", "pNorm", "chiSq", "df", "pChisq"))
    data.table::set(res, j=1L, value=n1)
    ## events minus predicted
    if (!"pred" %in% names(attributes(x))) predict(x)
    eMP1 <- x[, sum(e) - (max(ncg) * sum(e) / max(n)), by=list(t, cg)]
    rbindlist(split(eMP1, eMP1[, cg])), by="t"), alllow.cartesian=TRUE)
    merge(eMP1[, sum(V1), by=t], eMP1[, list(cg,t)], by="t", allow.cartesian=T)
    x1 <- copy(x[, list(t, cg)])
    x1[, "eMP" := attr(x, "pred")[, "eMP"]]
    t1 <- x[, sort(unique(t))]
    y1 <- x[, unique(cg)]
    ## get no. at risk for each unique time and covariate group
    n1 <- lapply(t1, FUN=function(t1) (setkey(x1[t==t1, sum(eMP), by=cg][y1], cg)))
    apply(x1, 1, print)
    n1 <- t(sapply(n1, function(x) x[, V1]))

    y1 <- x[, unique(cg)]
    x1[, sum(eMP), by=t]
    ## covariance
    if (!"COV" %in% names(attributes(x))) COV(x)
    cov1 <- attr(x, "COV")
    ## number of covariate groups
    ncg1 <- attr(x, "ncg")
### 2 groups only:
    if (ncg1==2){
### log-rank family
        ## make observed - expected for group 1
        eMP1 <- eMP1[, .SD, .SDcols=(ncol(eMP1) - 1)]
        ## log-rank, weight = 1
        res[1, "Q" := sum(eMP1)]
        res[1, "Var" := sum(cov1)]
        ## weight = n
        res[2, "Q" := sum(eMP1 * wt1[, n])]
        res[2, "Var" := sum(cov1 * wt1[, n^2])]
        res[3, "Q" := sum(eMP1 * wt1[, sqrtN])]
        res[3, "Var" := sum(cov1 * wt1[, n])]
        res[4, "Q" := sum(eMP1 * wt1[, S1])]
        res[4, "Var" := sum(cov1 * wt1[, S1^2])]
        res[5, "Q" := sum(eMP1 * wt1[, S2])]
        res[5, "Var" := sum(cov1 * wt1[, S2^2])]
        for (i in seq.int(fh1)){
            res[5L + i, "Q" := sum(eMP1 * wt1[, .SD, .SDcols=5L + i])]
            res[5L + i, "Var" := sum(cov1 * wt1[, .SD, .SDcols=5L + i]^2)]
        }
### supremum tests
### aka Renyi statistics
### (analagous to 2-sample Kolmogorov-Smirnov test)
        res1 <- data.table::data.table(matrix(0, nrow=5 + fh1, ncol=5))
        data.table::setnames(res1, c("W", "maxAbsZ", "Var", "Q", "pSupBr"))
        data.table::set(res1, j=1L, value=c("1", "n", "sqrt(n)", "S1", "S2", FHn))
        ## log-rank weights
        res1[1, "maxAbsZ" := max(abs(cumsum(eMP1)))]
        res1[2, "maxAbsZ" := max(abs(cumsum(eMP1 * wt1[, n])))]
        res1[3, "maxAbsZ" := max(abs(cumsum(eMP1 * wt1[, sqrtN])))]
        res1[4, "maxAbsZ" := max(abs(cumsum(eMP1 * wt1[, S1])))]
        res1[5, "maxAbsZ" := max(abs(cumsum(eMP1 * wt1[, S2])))]
        for (i in seq.int(fh1)){
            res1[5L + i, "maxAbsZ" :=  max(abs(cumsum(eMP1 * wt1[, .SD, .SDcols=5L + i])))]
        }
        ## results
        res1[, Var := res1[, Var]]
        res1[, "Q" := maxAbsZ / sqrt(Var)]
        res1[, "pSupBr" := sapply(Q, probSupBr)]
        data.table::setattr(res1, "class", c("sup", class(res1)))
    }
### >2 groups
    if (ncg1 > 2){
        df1 <- seq.int(ncg1 - 1L)
        eMP1 <- eMP1[, .SD, .SDcols=seq.int(ncg1 + 1L, ncol(eMP1))]
        ## hold results
        res2 <- data.table::data.table(matrix(0, nrow=5 + fh1, ncol=4))
        data.table::setnames(res2, c("W", "chiSq", "df", "pChisq"))
        data.table::set(res2, j=1L, value=n1)
        ## we save results below as these are also used by
        ## tests for trend
        ## log-rank, weight = 1
        eMP1cs <- colSums(eMP1)
        cov1rs <- rowSums(cov1, dims=2)
        res2[1, "chiSq" := eMP1cs[df1] %*% solve(cov1rs[df1, df1]) %*% eMP1cs[df1]]
        eMP1n <- colSums(sweep(eMP1, MARGIN=1, STATS=wt1[, n], FUN="*"))
        ## do this in two steps the first time as easier to read
        cov1n <- sweep(cov1, MARGIN=3, STATS=wt1[, n^2], FUN="*")
        cov1n <- rowSums(cov1n, dims=2)
        res2[2, "chiSq" := eMP1n[df1] %*% solve(cov1n[df1, df1]) %*% eMP1n[df1]]
        eMP1sn <- colSums(sweep(eMP1, MARGIN=1, STATS=wt1[, sqrtN], FUN="*"))
        cov1sn <- rowSums(sweep(cov1, MARGIN=3, STATS=wt1[, n], FUN="*"), dims=2)
        res2[3, "chiSq" := eMP1sn[df1] %*% solve(cov1sn[df1, df1]) %*% eMP1sn[df1]]
        eMP1S1 <- colSums(sweep(eMP1, MARGIN=1, STATS=wt1[, S1], FUN="*"))
        cov1S1 <- rowSums(sweep(cov1, MARGIN=3, STATS=wt1[, S1^2], FUN="*"), dims=2)
        res2[4, "chiSq" := eMP1S1[df1] %*% solve(cov1S1[df1, df1]) %*% eMP1S1[df1]]
        eMP1S2 <- colSums(sweep(eMP1, MARGIN=1, STATS=wt1[, S2], FUN="*"))
        cov1S2 <- rowSums(sweep(cov1, MARGIN=3, STATS=wt1[, S2^2], FUN="*"), dims=2)
        res2[5, "chiSq" := eMP1S2[df1] %*% solve(cov1S2[df1, df1]) %*% eMP1S2[df1]]
        ## Fleming-Harrington
        for (i in seq.int(fh1)){
            eMP1FH <- colSums(sweep(eMP1,
                                    MARGIN=1,
                                    STATS=wt1[, unlist(.SD, use.names=FALSE), .SDcols=5L + i],
                                    FUN="*"))
            cov1FH <- rowSums(sweep(cov1,
                                    MARGIN=3,
                                    STATS=wt1[, unlist(.SD, use.names=FALSE), .SDcols=5L + i],
                                    FUN="*"),
                              dims=2)
            res2[5 + i, "chiSq" := eMP1FH[df1] %*% solve(cov1FH[df1, df1]) %*% eMP1FH[df1]]
        }
        ## results
        res2[, "df" := attr(x, "ncg") - 1L]
        res2[, "pChisq" := 1 - stats::pchisq(chiSq, df)]
        data.table::setattr(res2, "class", c("lrt", class(res2)))
### Tests for trend 
        if (is.null(scores)) scores <- seq.int(ncg1)
        ## scores - all combinations
        sAC1 <- as.matrix(expand.grid(scores, scores))
        ## scores - product of all combinations
        scoProd1 <- apply(sAC1, 1, prod)
        ## Log-rank
        res[1, "Q" := sum(eMP1cs * scores)]
        res[1, "Var" := sum(cov1rs * scoProd1)]
        res[2, "Q" := sum(eMP1n * scores)]
        res[2, "Var" := sum(cov1n * scoProd1)]
        res[3, "Q" := sum(eMP1sn * scores)]
        res[3, "Var" := sum(cov1sn * scoProd1)]
        res[4, "Q" := sum(eMP1S1 * scores)]
        res[4, "Var" := sum(cov1S1 * scoProd1)]
        res[5, "Q" := sum(eMP1S2 * scores)]
        res[5, "Var" := sum(cov1S2 * scoProd1)]
        ## Fleming-Harrington
        for (i in seq.int(fh1)){
            eMP1FH <- colSums(sweep(eMP1,
                                    MARGIN=1,
                                    STATS=wt1[, unlist(.SD, use.names=FALSE), .SDcols=5L + i],
                                    FUN="*"))
            cov1FH <- rowSums(sweep(cov1,
                                    MARGIN=3,
                                    STATS=wt1[, unlist(.SD, use.names=FALSE), .SDcols=5L + i],
                                    FUN="*"),
                              dims=2)
            res[5 + i, "Q" := sum(eMP1FH * scores)]
            res[5 + i, "Var" := sum(cov1FH * scoProd1)]
        }
    }
    ## results
    res[, "Z" := Q / sqrt(Var)]
    res[, "pNorm" := 2 * (1 - stats::pnorm(abs(Z)))]
    res[, "chiSq" := Q^2 / Var]
    res[, "df" := 1]
    res[, "pChisq" := 1 - stats::pchisq(chiSq, df)]
    data.table::setattr(res, "class", c("lrt", class(res)))
    data.table::setattr(x, "lrw", wt1)
    if (ncg1==2){
        data.table::setattr(x, "lrt", res)
        data.table::setattr(x, "sup", res1)
    } else {
        data.table::setattr(x, "lrt", res2)
        data.table::setattr(x, "tft", res)
    }
    return(attr(x, "lrt"))
}
### Helper functions
## Probability of supremum of absolute value of
## standard Brownian motion process B(t)
## 1e4 is good enough for all practical purposes for k
probSupBr <- function(x){
    k <- c(0L, seq.int(1e4))
    1 - (4 / pi) * (sum(((( - 1)^k) / (2 * k + 1)) * exp(-(((pi^2) * (2 * k + 1)^2) / (8 * x^2)))))
}


          P_1       P_2      eMP_1      eMP_2
 1: 2.1680672 3.8319328  2.1680672 -2.1680672
 2: 0.4174757 0.5825243 -0.5825243  0.5825243
 3: 0.8571429 1.1428571  0.8571429 -0.8571429
 4: 0.8988764 1.1011236 -0.1011236  0.1011236
 5: 0.9113924 1.0886076 -1.0886076  1.0886076
 6: 0.4520548 0.5479452 -0.5479452  0.5479452
 7: 0.4696970 0.5303030  0.4696970 -0.4696970
 8: 0.0000000 0.0000000  0.0000000  0.0000000
 9: 0.9090909 1.0909091 -1.0909091  1.0909091
10: 0.4489796 0.5510204 -0.5510204  0.5510204
11: 0.4444444 0.5555556 -0.5555556  0.5555556
12: 0.4500000 0.5500000 -0.5500000  0.5500000
13: 0.0000000 0.0000000  0.0000000  0.0000000
14: 0.0000000 0.0000000  0.0000000  0.0000000
15: 0.0000000 0.0000000  0.0000000  0.0000000
16: 0.8800000 1.1200000 -0.1200000  0.1200000
17: 0.4347826 0.5652174 -0.5652174  0.5652174
18: 0.4500000 0.5500000 -0.5500000  0.5500000
19: 0.0000000 0.0000000  0.0000000  0.0000000
20: 0.0000000 0.0000000  0.0000000  0.0000000
21: 0.0000000 0.0000000  0.0000000  0.0000000
22: 0.0000000 0.0000000  0.0000000  0.0000000
23: 0.4444444 0.5555556 -0.5555556  0.5555556
24: 0.0000000 0.0000000  0.0000000  0.0000000
25: 0.0000000 0.0000000  0.0000000  0.0000000
26: 0.4000000 0.6000000 -0.6000000  0.6000000
27: 0.0000000 0.0000000  0.0000000  0.0000000
28: 0.0000000 0.0000000  0.0000000  0.0000000
