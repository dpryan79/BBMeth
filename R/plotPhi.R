#' @rdname plotPhi
#' @export
plotPhi <- function(BBMF, minPhi=NA, maxPhi=NA, minMu=0.05, maxMu=0.95) {
    if(!require(ggplot2))
        stop("Sorry, this function requires ggplot2, which either isn't
              installed or can't be loaded.")
    if(class(BBMF) != "BBMethFit")
        stop("BBMF must be a BBMethFit object")
    if(minMu < 0 | maxMu > 1) {
        message("minMu and maxMu should be between 0 and 1 (inclusive), since
                 the inverse-logit is used. These thresholds will be logit
                 transformed prior to use.")
        minMu = logit(minMu)
        maxMu = logit(maxMu)
    }

    #Guess bounds if none are provided
    if(is.na(minPhi) & is.na(maxPhi)) {
        IDX <- which(
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) >= minMu &
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) <= maxMu
        )
        minPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.01, 0.99))
        maxPhi = minPhi[2]
        minPhi = minPhi[1]
    } else if(is.na(minPhi)) {
        IDX <- which(
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) >= minMu &
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) <= maxMu
        )
        minPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.01))
    } else if(is.na(maxPhi)) {
        IDX <- which(
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) >= minMu &
            c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) <= maxMu
        )
        maxPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.99))
    }

    #Filter out extreme values, namely those with too low/high methylation (phi
    #estimates are bad there) or too low/high phi (again, largely meaningless
    #fits)
    IDX <- which(
        c(BBMF@mmPhi %*% t(est.phi(BBMF))) >= minPhi &
        c(BBMF@mmPhi %*% t(est.phi(BBMF))) <= maxPhi &
        c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) >= minMu &
        c(BBMF@mmMu %*% t(invlogit(est.mu(BBMF)))) <= maxMu
    )

    d <- data.frame(
        mus=invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))[IDX]),
        phis=(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX])
    )
    g <- ggplot(d, aes(x=mus, y=phis))
    g <- g + stat_density2d(aes(fill=..density..), geom="tile", contour=F)

    if(length(BBMF@prior) > 0) { #We've received a fit from optim
        prior = BBMF@prior
        d2 <- data.frame(x=c(seq(from=minMu, to=maxMu, by=0.001)))
        d2$middle = prior[1]
        d2$upper = d2$middle+prior[2]
        d2$lower = d2$middle-prior[2]
        #Clip to the specified bounds
        d2$upper[d2$upper > maxPhi] = maxPhi
        d2$middle[d2$middle > maxPhi] = maxPhi
        d2$lower[d2$lower > maxPhi] = maxPhi
        d2$upper[d2$upper < minPhi] = minPhi
        d2$middle[d2$middle < minPhi] = minPhi
        d2$lower[d2$lower < minPhi] = minPhi
        g <- g + geom_line(aes(x=x, y=middle), colour="red", show_guide=F, d2)
        g <- g + geom_line(aes(x=x, y=upper), linetype="dotted", colour="red", show_guide=F, d2)
        g <- g + geom_line(aes(x=x, y=lower), linetype="dotted", colour="red", show_guide=F, d2)
    }
    g <- g + xlim(minMu, maxMu) + ylim(minPhi, maxPhi)
    g <- g + xlab(expression(mu)) +
        ylab(expression(paste("logit(",phi,")",sep="")))
    g <- g + scale_fill_gradient(name="Density")
    suppressWarnings(plot(g))
}
