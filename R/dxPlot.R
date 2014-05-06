##' @name dxPlot
##' @title Diagnostic plots for coxph model
##' @rdname dxPlot
##' @export dxPlot
dxPlot <- function(x, ...){
    UseMethod("dxPlot")
}
##' @rdname dxPlot
##' @aliases dxPlot.coxph
##' @method dxPlot coxph
##' @S3method dxPlot coxph
##'
##' @include multiPlot.R
##'
##' @param x An object of class \code{coxph}
##' @param ... Additional arguments (not implemented)
##' @param ties  Method of handling ties when refitting model. Must be
##' one of \code{breslow} or \code{efron}
##' @param defCont Definition of continuous variable.
##' Variables with more than \emph{n} unique values will
##' be split into quantiles to facilitate graphs.
##' (This does not apply to factor variables)
##' @param noQuantiles No. of quantiles into which to
##' split continuous variables
##' @param noPerPage Number of plots per page
##' @param width Width of screen (display device) in pixels.
##' Set to \code{NULL} for default plot size
##' @param height Height of screen (display device) in pixels.
##' Set to \code{NULL} for default plot size
##'
##' @details
##' The \bold{Cox-Snell} residuals are used to assess the fit of a
##' proportional hazards model.
##' \cr
##' The residuals are generated from the (non time-dependent)
##' covariates \eqn{Z}, a matrix
##' with one row per observation (total \eqn{n}) and additional indicators
##' of time \eqn{T} and status \eqn{\delta}{D}.
##' estimated coefficients \eqn{b}, where \eqn{b} is a vector of
##' length \eqn{p} (no. of predictors) :
##' \deqn{ r_j = \hat{H}_0(T_j) \exp ( \sum_{k=1}^p Z_{jk}b_k ), \quad j=1,...,n }{
##'  r[j] = H0(T[j]) exp ( SUM Z[jk] b[k] ), j=1,...,n }
##' Here \eqn{\hat{H}_0}{H0} is Breslows estimator of the `baseline'
##' hazard (i.e. all coefficients are zero).
##' \cr \cr
##' If the coefficients are close to their true values, then \eqn{r_j}{r[j]}
##' should follow a unit-exponential distribution, i.e. \eqn{H_0(t)=t}{H0(t)=t}.
##' \cr \cr
##' To verify this, we calculate the Nelson-Aalen estimator of the
##' cumulative hazard rate of the \eqn{r_j}{r[j]}s. A plot of this
##' estimator against \eqn{r_j}{r[j]} should be a straight line
##' through the origin with a slope of \eqn{1}.
##' \cr \cr
##' The \bold{martingale} residual is used to help determine the best functional
##' form of a covariate in a \code{coxph} model.
##' The Cox model assumes that the hazard function satisfies:
##' \deqn{\lambda_{i}(t) = \lambda_0(t) \exp(X_i\beta)}{
##'  h[i](t) = h[0](t) exp (X[i]Beta)}
##' That is, for a continuous variable, an unit increase in the variable produces the same
##' change in risk across the value of the variable. (E.g. an increase of age of 5 years leads to the same
##' change in hazard no matter what the increase is from or to).
##' \cr \cr
##' To verify this is the case, a null model is fitted
##' (i.e no coefficients, similar to intercept-only model in linear regression).
##' Martingale residuals are calcuated for this. Plots of these residuals against the values of each of the
##' predictors in the model are shown.
##' \cr \cr
##' If the correct model for covariate \eqn{j}
##' is based on a smooth function \eqn{f()}, i.e.
##' \eqn{ \exp(f(X_j)\beta_j)}{exp(f(X[j])Beta[j])}
##' then the following should hold:
##' \deqn{E(M_i|X_{ij}=x_j) \approx cf(x_j)}{
##'  E( M[i]| X[i,j]=x[j] ) = cf(x[j])}
##' Where \eqn{M} is the martingale residual and \eqn{c} depends
##' on the amount of censoring and is
##' roughly independent of \eqn{x_j}{x[j]}.
##' \cr \cr
##' A \code{lowess} smoothed line is added to the plot.
##' This should be approximately linear if the
##' assumption of proportional hazards is met.
##' If the plot shows a threshold, a discretised version of the covariate
##' may be preferable.
##' \cr \cr
##' The assumption of \bold{proportional hazards} can be checked in a
##' number of ways.
##' These methods work by stratifying a covariate \eqn{G} into
##' \eqn{k} disjoint strata. A stratified \code{coxph} model is fitted to these
##' strata and one is selected as a reference.
##' \cr
##' The \bold{cumulative hazard} \eqn{\hat{H}_g(t)}{H[g](t)} is
##' plotted for each stratum \eqn{g}.
##' These should be a constant multiple of the reference stratum
##' \eqn{\hat{H}_1(t)}{H[1](t)} over time.
##' \cr \cr
##' Another way to compare these is to plot the
##' \bold{differences in log cumulative hazard}, that is:
##' \deqn{\log \hat{H}_g(t) - \log \hat{H}_1(t), \quad g = 2 , ..., k}{
##'  log(H[1](t)) - log(H[g](t)) for g = 2 to k}
##' Each curve should be horizontal and constant over time.
##' Curves above zero indicate
##' increased hazard in the stratum \eqn{g} vs the reference at that time.
##' \cr \cr
##' Finally \bold{Andersen plots} show
##' \deqn{\log \hat{H}_1(t) vs \log \hat{H}_g(t), \quad g = 2,...,k}{
##'  log(H[g](t)) vs log(H[1](t)) , g = 2 to k}
##' If proportional hazards are present, these should be straight lines through
##' the origin. If the curve is convex this shows that
##' \eqn{ \hat{H}_g(t) \div \hat{H}_1(t)}{
##' H[g](t) / H[1](t)} is an increasing function of \eqn{t}{time}.
##' Thus if convex, the hazard rate in \eqn{g} is increasing vs the reference.
##'
##' @return Graphs as described above. Plotted with base graphics.
##'
##' @note
##' Caution - for plots to verify proportional hazards: the variance of the curves is
##' not constant over time.
##'
##' @references Examples are from
##' Klein J, Moeschberger M 2003
##' \emph{Survival Analysis}, 2nd edition.
##' New York: Springer.
##' Example 11.1 - 11.7, pg 355-66
##' @references Last example is from:
##' Therneau T, Grambsch P 2000.
##' \emph{Modeling Survival Data}, 1st edition.
##' New York: Springer.
##' Section 5.1.2, pg 91.
##' @references
##' Andersen PK, Borgan O, Gill R, Keiding N 1982.
##' Linear Nonparametric Tests for Comparison of Counting Processes, with
##' Applications to Censored Survival Data, Correspondent Paper.
##' \emph{International Statistical Review} \bold{50}(3):219--44.
##' \href{http://www.jstor.org/stable/1402489}{JSTOR}
##' @examples
##' data(bmt, package="KMsurv")
##' bmt <- within(bmt, {
##' z1 <- z1 -28
##' z2 <- z2- 28
##' z3 <- z1*z2
##' z4 <- as.double( group== 2 )
##' z5 <- as.double( group== 3 )
##' z6 <- z8
##' z7 <- (z7 / 30) - 9
##' z8 <- z10
##' })
##' c1 <- coxph(Surv(t2,d3) ~ z1+z2+z3+z4+z5+z6+z7+z8,
##'             method="breslow",
##'             data=bmt)
##' dxPlot(c1)
##' data(alloauto, package="KMsurv")
##' c1 <- coxph(Surv(time,delta) ~ factor(type),
##'             method='breslow',
##'             data=alloauto)
##' dxPlot(c1)
##' c1 <- coxph(formula = Surv(time, status == 2) ~ age + log(bili), data=pbc)
##' dxPlot(c1)
dxPlot.coxph <- function(x, ...,
                         ties="breslow",
                         defCont=2,
                         noQuantiles=5,
                         noPerPage=2,
                         height=NULL,
                         width=NULL
                         ){
    stopifnot(inherits(x, "coxph"))
### set up plotting params
    .multiPlot(noPerPage=noPerPage, width=width, height=height)
### no. events (take from 2nd column of Surv object)
    e1 <- model.frame(x)[, 1][, "status"]
###----------------------------------------
### Cox-Snell residuals
    coxSnellRes <- e1 - survival:::residuals.coxph(x, type="martingale")
### intercept-only model
    s1 <- survfit(Surv(coxSnellRes, e1) ~ 1)
### hazard for intercept-only model
    H0 <- cumsum(s1$n.event / s1$n.risk)
    main1 <- "Cox-Snell residuals vs cumulative hazard of residuals \n Should follow line at 45 degrees if well fit"
    main2 <- "Complete model"
    plot(s1$time, H0, type='s',
         xlab="Cox-Snell residual",
         ylab="Nelson-Aalen estimate of cumulative hazard of residuals",
         main=main2,
         col='blue')
    abline(0, 1, col='red', lty=2)
    graphics::mtext(main1, line=0.3, outer=TRUE)
### counter to check whether to add mtext
    count1 <- 1
###----------------------------------------
### Cox-Snell stratified by each co-variate
    for (i in 2:length(model.frame(x))){
        n1 <- names(model.frame(x))[i]
        mf1 <- model.frame(x)[ ,i]
### split if continous (default in quantiles)
        if ( length(unique(mf1)) > defCont && !(is.factor(mf1)) ){
            mf1 <- Hmisc::cut2(mf1, g=noQuantiles)
        }
        mf1 <- as.factor(mf1)
        ns1 <- levels(mf1)
        s1 <- survfit(Surv(coxSnellRes, e1) ~ mf1)
        main2 <- paste0("Stratified on ", n1)
### is mf1 has been divided into factor
### plot one line per level of factor
        plot(s1, fun="cumhaz",
             xlab="Cox-Snell residual",
             ylab="Nelson-Aalen estimate of cumulative hazard of residuals",
             main=main2,
             col=rainbow(length(ns1), start=0.5, end=1)
             )
### add legend
        graphics::legend("bottomright", ns1,
                         lty=1,
                         title=paste0("Levels of ", n1),
                         col = rainbow(length(ns1), start=0.5, end=1)
                         )
        graphics::abline(0, 1, lty=2)
        count1 <- count1+1
        if (count1 > noPerPage){
            count1 <- 1
            graphics::mtext(main1, line=0.3, outer=TRUE)
        }
    }
### ----------------------------------------
### martingale residuals
    main1 <- "Covariate vs martingale residuals (from model without this covariate)
If not linear, suggests discretised version of covariate preferable \n (note no plots for binary covariates)"
    n1 <- names(x$coefficients)
###
    for (i in 2:length(model.frame(x))){
### skip binary covariates
        if ( length(unique(model.frame(x)[, i])) <=2 ) next
        f1 <- deparse(stats::formula(x))
### refit as intercept only
        f2 <- stats::as.formula(".~1")
        c2 <- stats::update(x, formula=f2)
### get residuals (hidden method in survival)
        r1 <- survival:::residuals.coxph(c2, type="martingale")
### get names of the coefficients from model.frame
### note excluding Surv
        mf1 <- model.frame(x)[, i]
        n2 <- names(model.frame(x))[i]
        plot(mf1, r1,
             xlab=n2,
             ylab="Martingale residual = excess events (deaths)",
             main=paste0("covariate ", n2)
             )
### add smoothing line
### use iter=0 to prevent removal of potential outliers
        graphics::lines(stats::lowess(mf1, r1, iter=0), col="red")
        count1 <- count1+1
        if (count1 > noPerPage){
            count1 <- 1
            graphics::mtext(main1, line=0.3, outer=TRUE)
        }
    }
###----------------------------------------
### cumulative hazard
### one curve per coefficient
    main1 <- "Cumulative hazard by time, stratfied by covariate \n If hazards proportional then curves should be constant multiples of baseline"
    for (i in 2:length(model.frame(x))){
        n1 <- names(model.frame(x))[i]
        mf1 <- model.frame(x)[, i]
        if (length(unique(mf1)) > defCont){
            mf1 <- Hmisc::cut2(mf1, g=noQuantiles)
        }
### put original time + events into data.frame
        df1 <- data.frame(t1=unclass(model.frame(x)[[1]])[, 1],
                          e1=unclass(model.frame(x)[[1]])[, 2],
                          mf1=mf1)
        df1$mf1 <- as.factor(df1$mf1)
        ns1 <- levels(df1$mf1)
### ### refit as stratified model
        f1 <- stats::as.formula("Surv(t1, e1) ~ strata(mf1)")
        c2 <- coxph(f1, data=df1)
        s1 <- survfit(c2, type="aalen")
        plot(s1, fun='cumhaz', log=TRUE,
             ylab="Cumulative hazard",
             xlab="Time",
             main=paste0("Covariate ", n1),
             col=rainbow(length(ns1), start=0.5, end=1)
             )
### add legend
        graphics::legend("bottomright", ns1,
                         title = paste0("Levels of: ", n1),
                         lty=1,
                         col = rainbow(length(ns1), start=0.5, end=1)
                         )
        count1 <- count1+1
        if (count1 > noPerPage){
            count1 <- 1
            graphics::mtext(main1, line=0.3, outer=TRUE)
        }
    }
###----------------------------------------
### differences in log-hazard curves
### note repeats code above which is inefficient
### but this is easier to read...
    main1 <- "Difference in log cumulative hazard rates by covariate \n (should be constant over time)
If >0 shows survival advantage for reference stratum"
### store merged frames and names
    m1l <- ns1l <- vector("list", (ncol(model.frame(x))-1))
    for (i in 2:length(model.frame(x))){
        n1 <- names(model.frame(x))[i]
        mf1 <- model.frame(x)[, i]
        if (length(unique(mf1)) > defCont){
            mf1 <- Hmisc::cut2(mf1, g=noQuantiles)
        }
### put original time + events into data.frame
        df1 <- data.frame(t1=unclass(model.frame(x)[[1]])[, 1],
                          e1=unclass(model.frame(x)[[1]])[, 2],
                          mf1=mf1)
        df1$mf1 <- as.factor(df1$mf1)
        ns1 <- levels(df1$mf1)
### ### refit as stratified model
        f1 <- stats::as.formula("Surv(t1, e1) ~ strata(mf1)")
        c2 <- coxph(f1, data=df1)
        s1 <- survfit(c2, type="aalen")
### get baseline hazard (= -log survival)
        s2 <- data.frame(t=s1$time, -log(s1$surv))
### names of strata
        ns1 <- paste0(n1, "=", ns1 )
        s2$strat <- unlist(sapply(1:length(s1$strata), function(i) rep(ns1[i], s1$strata[i]) ) )
### split by strata
        s3 <- by(s2[, -3], s2[, 3], FUN=identity)
### rename 2nd column to avoid duplicate column names
### needed if merging >3 data.frames
        for(j in seq_along(s3)) {
            colnames(s3[[j]])[2] <- ns1[j]
        }
### merge list based on time
        m1 <- Reduce(function(...)
                     merge(..., by="t",
                           all=TRUE), s3)
### replace NAs in first row with zeros
        m1[1, is.na(m1[1, ])] <- 0
### last observation carried forward
        m1 <- zoo::na.locf(m1)
### get differences between reference and other strata
        diffs1 <- sapply(3:ncol(m1), function (i) -log(m1[, 2]) +log(m1[, i]) )
### one line per difference
        matplot(m1$t, diffs1,
                type="s",
                lty=1,
                xlab="Time",
                ylab="Difference in log cumulative hazard rates",
                col=rainbow(length(ns1), start=0.7, end=1),
                main=paste0("Covariate ", n1))
### add legend
        l1 <- sapply(1:ncol(diffs1),
                     function(i) paste0(ns1[i+1], "  -  ", ns1[1])
                     )
        tit1 <- "Stratum - reference stratum"
        graphics::legend("bottomright", legend=l1,
                         lty=1,
                         title=tit1,
                         cex=0.75,
                         col = rainbow(length(ns1), start=0.7, end=1)
                         )
        graphics::abline(0, 0, lty=2)
        count1 <- count1+1
        if (count1 > noPerPage){
            count1 <- 1
            graphics::mtext(main1, line=0.3, outer=TRUE)
        }
### store result for Anderson plot
        m1l[[i-1]] <- m1
        ns1l[[i-1]] <- ns1
    }
###----------------------------------------
### Andersen plots
    main1 <- "Log cumulative hazard rates vs reference stratum \n (should be linear plot through origin)
If convex (towards top left) shows ratio of hazards is increasing over time"
    for (i in 2:length(model.frame(x))){
        n1 <- names(model.frame(x))[i]
        m1 <- m1l[[i-1]]
        ns1 <- ns1l[[i-1]]
        matplot(m1[, 2], m1[, -c(1,2)],
                lty=1,
                xlim=c(0,1),
                xlab=paste0("Cumulative hazard (reference): \n", ns1[1]),
                ylim=c(0,1),
                ylab="Cumulative hazard (see legend)",
                type="s",
                col=rainbow(length(ns1)-1, start=0.7, end=1),
                main=paste0("Covariate ", n1)
                )
        abline(0, 1, lty=2)
        graphics::legend("topleft", legend=ns1[-1],
                         lty=1,
                         col=rainbow(length(ns1[-1]), start=0.7, end=1)
                         )
        count1 <- count1+1
        if (count1 > noPerPage){
            count1 <- 1
            graphics::mtext(main1, line=0.3, outer=TRUE)
        }
     }
}
