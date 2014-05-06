##' @name gamTerms
##' @export
##' @title Individual terms of a Generalized Additive
##' or Cox Proportional
##' Hazard Model
##' @description
##' Returns the individual terms of a \code{gam} or \code{coxph} object,
##' along with the standard errors, in a way useful for plotting.
##' @param fit The result of a Generalized Additive Model (\code{gam})
##' or a Cox proportional hazards model (\code{coxph})
##' @param se.fit if \code{TRUE}, also return the standard errors
##' @param link if \code{TRUE}, then the individual terms are centered so that
##' the average of the inverse-link of the data,
##' i.e., the data on the original scale, has mean zero
##' @param weights a vector of case weights, data is centered so that the
##' weighted mean is zero
##' @return A list with one element per term. Each element is a matrix
##' whose columns are x, y, and (optionally) se(y),
##' with one row per unique x value, sorted by the first column.
##' (This makes it easy to plot the results).
##' The first element of the list, \code{constant},
##' contains an overall mean for the decomposition
##' @seealso \code{gam}, \code{plot.gam}
##' @examples
##' data(air)
##' gfit <- gam::gam(ozone ~ gam::s(temperature) + gam::s(wind), data=air)
##' temp <- gamTerms(gfit)
##' identical( names(temp), c("constant", "temperature", "wind") )
##' ### air has 111 rows, but only 28 unique wind speeds:
##' dim(temp$wind)
##' ### plot the fit versus square root of wind speed
##' yy <- cbind(temp$wind[, 2],
##'             temp$wind[, 2] - 1.96*temp$wind[, 3],
##'             temp$wind[, 2] + 1.96*temp$wind[, 3])
##' ### Adding the constant makes this a plot of actual y (ozone)
##' ### at the mean temp
##' yy <- yy + temp$constant
##' graphics::matplot(sqrt(temp$wind[, 1]), yy, lty=c(1, 2, 2),
##' type='l', col=1, xaxt='n', xlab='Wind Speed', ylab='Ozone')
##' temp <- seq(3, 19, 2)
##' graphics::axis(1, sqrt(temp), format(temp))
##' @author Terry Therneau, Dirk Larson, updated from S-plus by Chris Dardis
##' @keywords plot
gamTerms <- function(fit,
                     se.fit=TRUE,
                     link=FALSE,
                     weights) {
### reconstruct the data, without transformations
    Terms <- all.vars(delete.response(terms(formula(fit))))
    keep <- match(c("", "formula", "data", "subset", "na.action"),
                  names(fit$call), nomatch = 0)
### get original dataset
    data <- eval(fit$call$data)
### Now get the terms, and match them
    fit$na.action <- NULL
### predicted values
    termp <- predict(fit, type = "terms", se.fit = se.fit)
### vnames <- attr(terms(fit$formula), "term.labels")
    vnames <- Terms
### We need to remove any offset from the vnames list
    Terms2 <- terms(as.formula(formula(fit)))
    if (length(attr(Terms2, 'offset')) >0)  {
        j <- attr(Terms2, 'offset') - attr(Terms2, 'response')
        vnames <- vnames[-j]
        }
    tname1 <- attr(Terms2, "term.labels")
###
    if(se.fit){
        tfit <- termp$fit
    } else {
        tfit <- termp
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
        data <- data[!is.na(keep),  , drop = F]
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
    outlist <- list(constant = attr(termp$fit, "constant") + sum(tmean) )
    for(i in 1:nterm) {
        k <- match(tname2[i], tname1)
        xx1 <- data[[vnames[k]]]
        xx2 <- sort(unique(xx1))
        keep <- match(xx2, xx1)
        if(se.fit) {
            if(is.matrix(termp$se.fit)){
                 zz <- data.frame(x = xx2, y = tfit[keep, i],
                                  se = termp$se.fit[keep, i])
                 } else {
                     zz <- data.frame(x = xx2, y = tfit[keep, i],
                                      se = termp$se.fit[keep])
                 }
        } else {
            zz <- data.frame(x = xx2, y = tfit[keep, i])
        }
        outlist[[vnames[k]]] <- zz
    }
###
    return(outlist)
}
