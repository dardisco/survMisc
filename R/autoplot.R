##' @name autoplot
##' @export autoplot.survfit
##' @aliases autoplot.survfit
##' @method autoplot survfit
##' @description Uses \code{ggplot2} to plot survival curves (Kaplan-Meier plot)
##' @title Generate a ggplot for \code{survfit} object
##' @param object An object of class \code{survfit}
##' @param ... Additional arguments
##' @param alpha Transparency to use for confidence intervals or bands
##' @param shape Shape of marks to indicate censored onservations.
##' \cr Default is 3 which gives vertical ticks.
##' \cr Use 10 for circular marks.
##' @param xlab Label for x axis
##' @param ylab Label for y axis
##' @param title Title for graph
##' @param legendLabs Legend labels. These can be used to replace the names
##' of the strata from the fit. Should be given in the same order as those strata.
##' @param CI Include confidence intervals (plotted as lines).
##' \cr These are taken from the \code{survfit} object.
##' @param bands Include confidence bands (plotted as ribbons i.e.
##' filled and joined by diagonal lines rather than steps)
##' @param pval Add \eqn{p} value from log-rank test to lower left of plot
##' @param plotTable Add table below plot.
##' This contains one row for each time increment for each strata,
##' giving the number at risk at that time.
##' @param divideTime Divide time by this number to give time increments.
##' A larger number may be preferable with large data sets to avoid
##' crowding.
##' @param returnTable Return \code{data.frame} used for table.
##' @return A \code{ggplot} and optionally a \code{data.frame} as above.
##' @author Chris Dardis. Based on existing work by
##' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
##' @examples
##' data(kidney, package="KMsurv")
##' s1 <- survfit(Surv(time=time, event=delta) ~ type, data=kidney)
##' autoplot(s1)
##' autoplot(s1, CI=TRUE, pval=TRUE, plotTable=TRUE, divideTime=5,
##'  legendLabs=c("surgical", "percutaneous"),
##'  title="Time to infection following catheter placement \n
##'    by type of catheter, for dialysis patients")
##' s1 <- survfit(Surv(time=time, event=delta) ~ 1, data=kidney)
##' autoplot(s1)
##' data(rectum.dat, package="km.ci")
##' s1 <- survfit(Surv(time, status) ~ 1, data=rectum.dat)
##' ### change confidence intervals to log Equal-Precision confidence bands
##' km.ci::km.ci(s1, method="logep")
##' autoplot(s1, bands=TRUE)
autoplot.survfit <- function(object, ...,
                             alpha=0.5,
                             shape=3,
                             xlab="Time",
                             ylab="Survival",
                             title="Marks show times with censoring",
                             legendLabs=NULL,
                             CI=FALSE,
                             bands=FALSE,
                             pval=FALSE,
                             plotTable=FALSE,
                             divideTime=1,
                             returnTable=FALSE){
    stopifnot(inherits(object, "survfit"))
    if(!is.null(legendLabs)) stopifnot(
        length(legendLabs)==length(object$strata))
### generate data to plot
### declare variables (for R CMD check)
### st1 is vector for strata identification
    st1 <- stNames <- surv <- n.risk <- n.censor <- upper <- lower <- NULL
### change names for strata to legendLabs if required
    if(is.null(legendLabs)){
        stNames <- names(object$strata)
    } else {
        stNames <- legendLabs
    }
### add vector for one strata according to number of rows of strata
    st1 <- unlist(sapply( 1:length(object$strata),
                         function (i) rep(stNames[i], object$strata[i]) ))
### if only one strata (intercept only model)
    if (is.null(object$strata)) st1 <- as.factor(rep(1, length(object$time)))
### create data.table with data from survfit
### add column for strata
### (using data.table here as avoids duplication when adding rows later)
    dt1 <- data.table(time=object$time,
                      n.risk=object$n.risk,
                      n.event=object$n.event,
                      n.censor = object$n.censor,
                      surv=object$surv,
                      upper=object$upper,
                      lower=object$lower,
                      strata=factor(st1))
### add additional 2 rows for strata
### for time=0 and first time before any events
    for(i in 1:length(unique(st1))){
        dt2 <- dt1[dt1$strata==unique(st1)[i]]
        dt3 <- data.table(time=c(0, dt2$time[1]),
                          n.risk=rep(dt2$n.risk[1], 2),
                          n.event=c(0, 0),
                          n.censor=c(0, 0),
                          surv=c(1, 1),
                          upper=c(1, 1),
                          lower=c(1, 1),
                          rep(dt2$strata[i], 2)
                          )
        dt1 <- rbindlist(list(dt1,dt3))
    }
### order by strata and survival estimate
    dt1 <- dt1[order(dt1$strata, (1-dt1$surv)), ]
### plot single lines only
    if (!CI && !bands){
        g1 <- ggplot(data=dt1, aes(group=strata, colour=strata)) +
         geom_step(aes(x=time, y=surv), direction="hv") +
### add lines to show times where subjects censored
         geom_point(data=subset(dt1, n.censor>=1),
                   aes(x=time, y=surv), shape=shape) +
### use palette Dark2 for prominent shades
### (suitable for colorblind)
         scale_colour_brewer(type="qual", palette="Dark2",
                             guide=guide_legend(keywidth=3, keyheight=3)) +
         scale_x_continuous(xlab) +
         scale_y_continuous(ylab) +
         ggtitle(title) +
         theme(legend.text=element_text(size=15),
               legend.title=element_text(size=15)
               )
    }
    if (CI && !bands){
### use palette Dark2 for prominent shades
        g1 <- ggplot(data=dt1, aes(colour=strata, group=strata)) +
         geom_step(aes(x=time, y=surv), direction="hv", fill=strata) +
         geom_step(aes(x=time, y=upper),
                   direction="hv", linetype=10, alpha=alpha) +
         geom_step(aes(x=time, y=lower),
                   direction="hv", linetype=10, alpha=alpha) +
         geom_point(data=subset(dt1, n.censor>=1),
                    aes(x=time, y=surv), shape=shape) +
         scale_colour_brewer(type="qual", palette="Dark2",
                             guide=guide_legend(keywidth=3, keyheight=3)) +
         scale_x_continuous(xlab) +
         scale_y_continuous(ylab) +
         ggtitle(title) +
         theme(legend.text=element_text(size=15),
               legend.title=element_text(size=15)
               )
    }
    if (bands){
### use palette Set2 for lighter shades as large fill area
        g1 <- ggplot(data=dt1, aes(colour=strata, group=strata)) +
         geom_step(aes(x=time, y=surv), direction="hv") +
         geom_ribbon(aes(x=time, ymax=upper, ymin=lower, fill=strata),
                     alpha=alpha) +
         geom_point(data=subset(dt1, n.censor>=1),
                    aes(x=time, y=surv), shape=shape) +
         scale_fill_brewer(type="qual", palette="Set2") +
         scale_colour_brewer(type="qual", palette="Dark2",
                             guide=guide_legend(keywidth=3, keyheight=3)) +
         scale_x_continuous(xlab) +
         scale_y_continuous(ylab) +
         ggtitle(title) +
         theme(legend.text=element_text(size=15),
               legend.title=element_text(size=15)
               )
    }
### p value for log-rank test
    if(pval) {
        sd1 <- survival::survdiff(eval(object$call$formula),
                                  data=eval(object$call$data))
        p1 <- stats::pchisq(sd1$chisq,
                            length(sd1$n) - 1,
                            lower.tail=FALSE)
        p1txt <- ifelse(p1 < 0.0001,
                        "p < 0.0001",
                        paste("Log-rank test \n p =", signif(p1, 3))
                        )
        g1 <- g1 + annotate("text",
                            x = 0.1 * max(dt1$time),
                            y = 0.2,
                            label = p1txt)
    }
### times to show on table
    times1 <- seq(0, max(object$time), by=divideTime)
### make table for plot
### if intercept-only model
    if (is.null(object$strata)){
    df1 <- data.frame(
        strata = as.factor(rep(1, length(times1))),
        time = summary(object, times = times1, extend = TRUE)$time,
        n.risk = summary(object, times = times1, extend = TRUE)$n.risk
        )
    } else {
        df1 <- data.frame(
            strata = summary(object, times=times1, extend=TRUE)$strata,
            time = summary(object, times=times1, extend=TRUE)$time,
            n.risk = summary(object, times=times1, extend=TRUE)$n.risk
            )
        if(!is.null(legendLabs))
### change names of strata to legend labels
            df1$strata <- factor(df1$strata, labels=legendLabs)
    }
    if(plotTable) {
### table grob (graphical object)
        tg1 <- ggplot(df1, aes(x=time, y=strata,
                               label = format(n.risk, nsmall = 0))) +
                geom_text(size = 3.5) +
                scale_y_discrete(breaks=as.character(levels(df1$strata)),
                                 labels=levels(df1$strata)) +
                scale_x_continuous(limits=c(0, max(object$time)),
                                   breaks=times1) +
                ggtitle("Number at risk by time") +
                theme_grey() +
                theme(
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    plot.title = element_text(size = rel(0.75)),
                    axis.title.x = element_blank(),
### change background to white
### panel.background = element_blank(),
                    axis.title.y = element_blank()
                    )
    }
### plot without table
    if(!plotTable) print(g1)
### plot with table
    if(plotTable){
### split g1 into components
        g2 <- ggplot_gtable(ggplot_build(g1))
### get legend items
        leg <- which(sapply(g2$grobs, function(x) x$name) == "guide-box")
        legend <- g2$grobs[[leg]]
        grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                                 legend,
                                 tg1 + theme(legend.position="none"),
                                 ncol=2,
                                 widths=c(0.85, 0.15),
                                 heights=c(2, 0.75)
                                 ))
    }
    if(returnTable) return(df1)
}
###
