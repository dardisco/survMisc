##' @name gamTerms
##' @export
##' @title Individual terms of a Generalized Additive
##' or Cox Proportional Hazards Model
##' 
##' @description
##' Returns the individual terms of a \code{gam} or \code{coxph} object,
##' along with the standard errors, in a way useful for plotting.
##' 
##' @param fit The result of a Generalized Additive Model (\code{gam})
##' or a Cox proportional hazards model (\code{coxph})
##' @param se If \code{TRUE}, also return the standard errors
##' @param link If \code{TRUE}, then the individual terms are centered so that
##' the average of the inverse-link of the data,
##' i.e., the data on the original scale, has mean zero
##' @param weights A vector of case weights.
##' \cr
##' If not supplied (the default), the data is centered so that the
##' weighted mean is zero.
##' @param data A \code{data.frame} in which to evaluate the model.
##' \cr
##' If missing, \code{eval(fit$call$data)} is used.
##' 
##' @return A \code{list} with one element per term.
##' \cr
##' Each element is a \code{matrix}
##' whose columns are \code{x}, \code{y}, and (optionally) \code{se(y)}.
##' \cr
##' There is one row per unique \code{x} value, and the matrix is sorted by these values.
##' (This makes it easy to plot the results).
##' \cr
##' The first element of the list, \code{constant},
##' contains an overall mean for the decomposition.
##' 
##' @seealso
##' \code{\link{air}}
##' \cr
##' \code{?gam::gam}
##' \cr
##' \code{?gam::plot.gam}
##' 
##' @examples
##' data(air, package="survMisc")
##' gfit <- gam::gam(ozone ~ gam::s(temperature) + gam::s(wind), data=air)
##' temp <- gamTerms(gfit)
##' identical(names(temp), c("constant", "temperature", "wind"))
##' ### air has 111 rows, but only 28 unique wind speeds:
##' dim(temp$wind)
##' ### plot the fit versus square root of wind speed
##' yy <- cbind(temp$wind[, 2],
##'             temp$wind[, 2] - 1.96 * temp$wind[, 3],
##'             temp$wind[, 2] + 1.96 * temp$wind[, 3])
##' ### Adding the constant makes this a plot of
##' ### actual y (ozone) at the mean temp
##' yy <- yy + temp$constant
##' graphics::matplot(sqrt(temp$wind[, 1]), yy, lty=c(1, 2, 2),
##'                  type='l', col=1, xaxt='n', xlab='Wind Speed', ylab='Ozone')
##' temp <- seq(3, 19, 2)
##' graphics::axis(1, sqrt(temp), format(temp))
##' 
##' @author Terry Therneau, Dirk Larson. Updated/adapted from S-plus by Chris Dardis.
gamTerms <- function(fit,
                     se=TRUE,
                     link=FALSE,
                     weights,
                     data) {
### reconstruct the data, without transformations
### get precictors (right hand side of formula)
    rhs1 <- all.vars(delete.response(terms(formula(fit))))
    keep <- match(c("", "formula", "data", "subset", "na.action"),
                  names(fit$call),
                  nomatch=0)
### get original dataset
### data <- get(deparse(fit$call$data), env=parent.frame())
    if(missing(data)) data <- eval(fit$call$data)
### Now get the terms, and match them
    fit$na.action <- NULL
### predicted values
    pred1 <- predict(fit, type="terms", se.fit = se)
### vnames <- attr(terms(fit$formula), "term.labels")
    vnames <- rhs1
### We need to remove any offset from the vnames list
    terms2 <- terms(as.formula(formula(fit)))
    if (length(attr(terms2, 'offset')) > 0)  {
        j <- attr(terms2, 'offset') - attr(terms2, 'response')
        vnames <- vnames[-j]
    }
    tname1 <- attr(terms2, "term.labels")
###
    if(se){
        tfit <- pred1$fit
    } else {
        tfit <- pred1
    }
###
    if(is.matrix(tfit)){
        tname2 <- dimnames(tfit)[[2]]
    } else {
        tfit <- as.matrix(tfit)
        tname2 <- names(fit$assign)
        tname2 <- tname2[tname2 != "(Intercept)"]
    }
###
    if(nrow(data) > nrow(tfit)) {
        keep <- match(row.names(data), dimnames(tfit)[[1]])
        data <- data[!is.na(keep),  , drop=FALSE]
    }
###----------------------------------------
### Find the mean of each column, and subtract it out
### S doesn't quite get this right
### The extras get added back into the constant
### If link=T, then do it on the y scale rather than the linear predictor
    if (!missing(weights) && length(weights) != nrow(tfit))
        stop("Wrong length for weight vector")
    if (link) {
### the name of the link
        temp <- as.vector(fit$family[2])
### declare variable
        glm.links <- NULL
        indx <- match(temp, unlist(glm.links[1, ]))
        if (is.na(indx)) stop("Link function not recognized")
###
        ilink <- glm.links['inverse', indx][[1]]
        lnk   <- glm.links['link', indx][[1]]
###
        if (missing(weights)) {
            tmean <- lnk(colMeans(ilink(tfit)))
        } else {
            tmean <- lnk(colSums(weights * ilink(tfit))/sum(weights))
        }
    } else {
        if (missing(weights)){
            tmean <- colMeans(tfit)
        } else {
            tmean <- colSums(weights*tfit) / sum(weights)
        }
    }
        tfit <- tfit - rep(tmean, each=nrow(tfit))
###
### Now walk through the columns one by one, and construct
### a data frame for each one
    nterm <- length(tname2)
    outlist <- list(constant = attr(pred1$fit, "constant") + sum(tmean) )
    for(i in 1:nterm) {
        k <- match(tname2[i], tname1)
        xx1 <- data[[vnames[k]]]
        xx2 <- sort(unique(xx1))
        keep <- match(xx2, xx1)
        if(se) {
            if(is.matrix(pred1$se.fit)){
                 zz <- data.frame(x = xx2, y = tfit[keep, i],
                                  se = pred1$se.fit[keep, i])
                 } else {
                     zz <- data.frame(x = xx2, y = tfit[keep, i],
                                      se = pred1$se.fit[keep])
                 }
        } else {
            zz <- data.frame(x = xx2, y = tfit[keep, i])
        }
        outlist[[vnames[k]]] <- zz
    }
###
    return(outlist)
}
