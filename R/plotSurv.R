##' @name plot.Surv
##' @method plot Surv
##' @export
##' @title Plot Survival object
##' @description
##' Plots an object of class \code{Surv}.
##' \cr
##' Different methods apply to different types of \code{Surv} objects.
##' @export
##' @param x A \code{Surv} object
##' @param l Length of arrow. Length is \code{l/nrow(x)}
##' @param ... Additional arguments. These are passed to \code{graphics::arrows} when
##' drawing right- or left-censored observations.
##' @return A graph (base graphics).
##' The type of graph depends on the \code{type} of the \code{Surv} object.
##' This is given by \code{attr(s, which="type")} :
##' \item{counting}{Lines with an arrow pointing right if right censored}
##' \item{right}{Lines with an arrow pointing right if right censored}
##' \item{left}{Lines with an arrow pointing left if left censored}
##' \item{interval}{If censored:
##'  \itemize{
##'    \item Lines with an arrow pointing right if right censored.
##'    \item Lines with an arrow pointing left if left censored.
##'  }
##' If not censored:
##'  \itemize{
##'    \item Lines if observations of more than one time point
##'    \item Points if observation of one time only (i.e. start and end times are the same)
##'  }
##' }
##'
##' @keywords plot
##' 
##' @examples
##' df0 <- data.frame(t1=c(0, 2, 4, 6, NA, NA, 12, 14),
##'                   t2=c(NA, NA, 4, 6, 8, 10, 16, 18))
##' s5 <- Surv(df0$t1, df0$t2, type="interval2")
##' plot(s5)
plot.Surv <- function(x, l=3, ...){
    if(!class(x)=="Surv") stop("Only applies to class 'Surv'")
    ## if(!attr(x,which="type")=="right") warning("Applies to right censored data")
    type1 <- attr(x,which="type")
    nr1 <- nrow(x)
    ## largest observed time (used for x axis on plot)
    max1 <- ifelse(type1=="interval"|type1=="counting",
                   max(x[,1:2][!is.na(x[,1:2])]),
                   max(x[,1]))
    ## add jitter to ends of x axis
    j1 <- stats::runif(1, max1/50, 2 * max1 / 50)
    maxX <- max1+j1
    minX <- 0-j1
    ## blank plot with appropriate axes
    graphics::plot(x=0, type="n",
         xlim=c(minX, maxX),
         ylim=c(0, nrow(x)),
         xlab="time",
         ylab="observation")
    graphics::grid()
    ## make lines
    Seg <- function(i, x0=0, x1=x[i, 1]) graphics::segments(x0=x0, y0=i, x1=x1, y1=i)
    ## make arrows, default point to right as this is more common
    ## xpd allows arrows to be on the edge of plot margins
    Arr <- function(i, code=2, x0=0, x1=x[i,1]){
        graphics::arrows(x0=x0, y0=i,
                         x1=x1, y1=i,
                         length=(l/nr1),
                         angle=67.5,
                         code=code,
                         xpd=TRUE,
                         ...)
    }
    ## for right censored data
    plotR <- function(x){
        for (i in 1:nrow(x)){
            if(x[i,"status"]==0){
                Seg(i)
            } else {
                Arr(i)
            }}}
    ## for left censored data
    plotL <- function(x){
            for (i in 1:nrow(x)){
                if(x[i,"status"]==0){
                    Arr(i,code=1)
                } else {
                    Seg(i)
                }}}
    ## for counting format data
    plotC <- function(x){
        for (i in 1:nrow(x)){
            if(x[i,"status"]==0){
                Seg(i,x0=x[i,1],x1=x[i,2])
            } else {
                Arr(i,x0=x[i,1],x1=x[i,2])
            }}}
    ## for interval format data
    plotI <- function(x){
        for (i in 1:nrow(x)){
            ## switch based on type of observation
            ## (use point if exact observation)
            switch(x[[i,3]]+1,
                   Arr(i, x0=x[i, 1], x1=max1),
                   graphics::points(x[i,1], i, pch=21, bg=1),
                   Arr(i, code=1, x0=0, x1=x[i,1]),
                   Seg(i, x0=x[i,1], x1=x[i,2]))
            }}
    ## plot based on Surv type
    switch(type1,
           right = plotR(x),
           left = plotL(x),
           counting = plotC(x),
           interval = plotI(x))
    ## add title
    int1 <- ifelse(type1=="interval",
                   "\n point = exactly observed event time",
                   "")
    text1 <- " censored survival data\n arrow = censored observation"
    tit1 <- paste(type1, text1, int1, sep="")
    graphics::title(tit1)
}
