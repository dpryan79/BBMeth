fitGaussian <- function(phi) {
    gaussian <- function(x, phi) {
        mu = x[1]
        SD = exp(x[2])
        ll = -sum(dnorm(phi, mean=mu, sd=SD, log=T))
        return(ll)
    }
    res <- optim(c(mean(phi, na.rm=T), sd(phi, na.rm=T)), gaussian, phi=phi, hessian=F, control=list(maxit=1000))
    return(res)
}

#' @rdname estimatePrior
#' @export
estimatePrior <- function(BBMF, minPhi=NA, maxPhi=NA, minMu=0.2, maxMu=0.8, do.plot=F) {
    if(class(BBMF) != "BBMethFit")
        stop("BBMF must be a BBMethFit object.")
    if(minMu < 0 | maxMu > 1) {
        message("minMu and maxMu should be between 0 and 1 (inclusive), since the inverse-logit is used. These thresholds will be logit transformed prior to use.")
        minMu = logit(minMu)
        maxMu = logit(maxMu)
    }

    #Guess appropriate bounds if none are provieded
    if(is.na(minPhi) && is.na(maxPhi)) {
        IDX = which(invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) >= minMu &
            invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) <= maxMu)
        minPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.05, 0.95))
        maxPhi=minPhi[2]
        minPhi=minPhi[1]
    } else if(is.na(minPhi)) {
        IDX = which(invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) >= minMu &
            invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) <= maxMu)
        minPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.05))
    } else if(is.na(maxPhi)) {
        IDX = which(invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) >= minMu &
            invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) <= maxMu)
        maxPhi = quantile(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX], c(0.95))
    }

    #Filter as desired
    IDX = which(
        c(BBMF@mmPhi %*% t(est.phi(BBMF))) >= minPhi &
        c(BBMF@mmPhi %*% t(est.phi(BBMF))) <= maxPhi &
        invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) >= minMu &
        invlogit(c(BBMF@mmMu %*% t(est.mu(BBMF)))) <= maxMu
    )
    res = fitGaussian(c(BBMF@mmPhi %*% t(est.phi(BBMF)))[IDX])

    res$par[2] = exp(res$par[2])

    if(res$par[2] > 1)
        message("The estimated prior standard deviation seems high. Double-check if this is correct with plotPhi() and manually change minPhi and/or maxPhi if not.")
    if(do.plot) {
        plotPhi(BBMF, minPhi=minPhi, maxPhi=maxPhi, minMu=minMu, maxMu=maxMu)
    }
    BBMF@prior = res$par

    return(BBMF)
}

#Compute the log-likelihood given a prior
loglik.shrunken_phi <- function(x, M, Cov, prior, mmMu, mmPhi, lower_threshold, mu_dist) {
    upper_threshold=1-lower_threshold
    mu = x[1:ncol(mmMu)]
    phi = x[(ncol(mmMu)+1):length(x)]
    #Determine a and b
    mu = invlogit(mmMu %*% mu)
    phi = mmPhi %*% phi
    il_phi = invlogit(phi)

    #If mu is 0 then 'a' can become 0!!!
    if(any(mu<lower_threshold))
        mu[mu<lower_threshold] = lower_threshold
    #If mu is 1 then 'b' becomes 0!!!
    if(any(mu>upper_threshold))
        mu[mu>upper_threshold] = upper_threshold

    #If phi is 0, then 'a' and 'b' become infinite and the log-likelihood is also infinite!!!
    if(any(il_phi < lower_threshold))
        il_phi[il_phi<lower_threshold] = lower_threshold
    #If phi is 1 then 'a' can be infinite!!!
    if(any(il_phi > upper_threshold))
        il_phi[il_phi>upper_threshold] = upper_threshold

    tmp = (1/(il_phi))-1
    a = mu*tmp
    #If mu and phi are close enough to 0, then 'a' will be NaN rather than 1
    b = tmp*(1-mu) #tmp-a
    #Determine the prior log likelihood
    prior_phi = dnorm(phi, mean=prior[1], sd=prior[2], log=T) #This is on the logit scale
    prior_mu = log(f(mu))
    #Unshrunken current log likelihood
    vals = lchoose(Cov,M)-lbeta(b,a)+lbeta(Cov-M+b, M+a)
    if(any(abs(lbeta(b,a)) == Inf) | any(is.na(vals))) {
        message("Warning, a log-likelihood is +/-Inf or NaN! Diagnostic information follows:")
        idx <- which(abs(lbeta(b,a)) == Inf | is.na(vals))
        message("A")
        message(paste(a[idx], collapse=", "))
        message("B")
        message(paste(b[idx], collapse=", "))
    }

    ll = -sum(vals + prior_phi + prior_mu) #the summed -log-likelihood
    return(ll)
}

#BBM and BBMF describe a single location
fit.shrunken_phi <- function(BBM, BBMF, lower_threshold, mu_dist) {
    mmMu = BBMF@mmMu
    mmPhi = BBMF@mmPhi
    out <- list(
        mu = c(rep(NA, ncol(mmMu))),
        phi = c(rep(NA, ncol(mmPhi))),
        ll = NA,
        convergence = 100, #Denoting not full rank
        SE.mu = c(rep(NA, ncol(mmMu))),
        SE.phi = c(rep(NA, ncol(mmPhi)))
    )

    #If the model matrix won't end up being full rank, then abort here
    if(BBMF$convergence == 100)
        return(out)

    M = c(assays(BBM)[["M"]])
    Cov = c(assays(BBM)[["Cov"]])
    prior = BBMF@prior

    #Remove samples with 0 coverage, we'll still be full rank
    RMV <- which(Cov==0)
    if(length(RMV) > 0) {
        M <- M[-RMV]
        Cov <- Cov[-RMV]
        mmMu <- mmMu[-RMV,,drop=F]
        mmPhi <- mmPhi[-RMV,,drop=F]
    }

    #Fit
    est.mu = est.mu(BBMF)
    est.phi = est.phi(BBMF)
    fit <- optim(c(est.mu, est.phi), loglik.shrunken_phi, M=M, Cov=Cov, mmMu=mmMu, mmPhi=mmPhi, prior=prior, mu_dist=mu_dist,
        lower_threshold=1e-16, hessian=T,
        lower=logit(lower_threshold), upper=logit(1-lower_threshold),
        control=list(maxit=5000))

    #If a mu estimate changed by more than ~0.2, then it's possible that the phi gradient pulled things in the wrong direction
    est.mus = invlogit(mmMu %*% t(est.mu)) #est.mu is a row vector!
    mus = invlogit(mmMu %*% fit$par[1:ncol(mmMu)])
    if(any(abs(mus-est.mus) > 0.5)) {
        est.phi = fit$par[(1+ncol(mmMu)):length(fit$par)]
        fit <- optim(c(est.mu, est.phi), loglik.shrunken_phi, M=M, Cov=Cov, mmMu=mmMu, mmPhi=mmPhi, prior=prior, mu_dist=mu_dist,
            lower_threshold=1e-16, hessian=T,
            lower=logit(lower_threshold), upper=logit(1-lower_threshold),
            control=list(maxit=5000))
    }

    #Fill the output list
    out[["mu"]] = fit$par[1:ncol(mmMu)]
    out[["phi"]] = fit$par[(1+ncol(mmMu)):length(fit$par)]
    out[["ll"]] = -fit[["value"]]
    out[["convergence"]] <- fit[["convergence"]]

    #Compute the SE from the hessian
    if(fit$convergence == 0) {
        SE = tryCatch({sqrt(diag(solve(fit$hessian)))}, error=function(err){NA})
        if(!any(is.na(SE))) { #We need the hessian to be invertible
            out[["SE.mu"]] = SE[1:ncol(mmMu)]
            out[["SE.phi"]] = SE[(ncol(mmMu)+1):length(SE)]
        } else {
            out[["convergence"]] = 2
        }
    }

    return(out)
}

fit.shrunken.wrapper <- function(BBM, BBMF, lower_threshold, mu_dist) {
    if(nrow(BBM) > 1) {
        f = c(1:nrow(BBM))
        return(mapply(fit.shrunken_phi, split(BBM, f), split(BBMF, f), MoreArgs=list(lower_threshold=lower_threshold, mu_dist=mu_dist), SIMPLIFY=F))
    } else {
        return(fit.shrunken_phi(BBM, BBMF, lower_threshold=lower_threshold, mu_dist=mu_dist))
    }
}
 
#' @rdname shrinkPhi
#' @export
shrinkPhi <- function(BBMF, BBM, presplit=T, lower_threshold=1e-6, BPPARAM=MulticoreParam(workers=1)) {
    if(class(BBMF) != "BBMethFit")
        stop("BBMF must be a BBMethFit object.")
    if(class(BBM) != "BBMeth")
        stop("BBM must be a BBMeth object.")
    if(nrow(BBM) != length(BBMF))
        stop("BBM isn't compatible with BBMF. There must be one row in BBM for
            every entry in BBMF (i.e., nrow(BBM) == length(BBMF))")
    if(length(BBMF@prior) != 2) {
        BBMF=estimatePrior(BBMF)
    }

    #Estimate the background mu distribution
    d = density(c(invlogit(unique(BBMF@mmMu) %*% t(BBMF@est.mu))), from=0, to=1)
    fun = approxfun(d$x, d$y, yleft=0, yright=0)
    rm(d)

    #Be verbose
    BPPARAM$verbose=T

    if(BPPARAM$workers>1 & presplit==T) {
        f = c(rep(1:BPPARAM$workers, each=ceiling(nrow(BBM)/BPPARAM$workers), length.out=nrow(BBM)))
    } else {
        if(BPPARAM$workers>1)
            message("N.B., it can take a while for all of the workers to start running. Please be patient.")
        f = c(1:nrow(BBM))
    }

    res <- bpmapply(fit.shrunken.wrapper, split(BBM, f), split(BBMF, f), MoreArgs=list(lower_threshold=lower_threshold, mu_dist=fun), SIMPLIFY=F, BPPARAM=BPPARAM)

    if(length(res) != nrow(BBM)) #More than one worker
        res = unlist(res, recursive=F)

    #Add/replace values
    mcols(BBMF)$ll = c(sapply(res, function(x) return(x$ll)))
    mcols(BBMF)$convergence = c(sapply(res, function(x) return(x$convergence)))
    BBMF@mu = matrix(unlist(lapply(res, function(x) return(x$mu))), ncol=ncol(BBMF@mmMu), byrow=T)
    BBMF@phi = matrix(unlist(lapply(res, function(x) return(x$phi))), ncol=ncol(BBMF@mmPhi), byrow=T)
    BBMF@se.mu = matrix(unlist(lapply(res, function(x) return(x$SE.mu))), ncol=ncol(BBMF@mmMu), byrow=T)
    BBMF@se.phi = matrix(unlist(lapply(res, function(x) return(x$SE.phi))), ncol=ncol(BBMF@mmPhi), byrow=T)

    return(BBMF)
}
