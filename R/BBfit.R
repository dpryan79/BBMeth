#Assign each sample to different groups
split_samples <- function(mm) {
    mm2 <- unique(mm) #nrow(mm2) == number of values to estimate
    IDX <- apply(mm, 1, function(x) which(apply(mm2, 1, function(y) all(y-x==0))))
    return(IDX)
}
    
#Method of moments estimation of mu
mu.method.moments <- function(M, Cov, f, inv_mm) {
    X <- M/Cov
    if(any(X < 1e-16))
        X[X < 1e-16] = 1e-16
    if(any(X > 1-1e-16))
        X [X > 1-1e-16] = 1-1e-16
    X <- matrix(X)

    #X is a n_sample x 1 matrix
    #The coefficients are C in mm %*% C = X
    #So, multiply by the pseudoinverse of mm
    return(c(inv_mm %*% logit(X)))
}

#Method of moments estimation of phi
#Just regress the normalized residuals
#phi only makes sense for sites with Cov>2!
phi.method.moments <- function(M, Cov, f, inv_mm, est.mu, mmMu) {
    #Get the mu estimates (a column matrix)
    mus = c(invlogit(mmMu %*% matrix(est.mu)))
    if(any(mus < 1e-16))
        mus[mus < 1e-16] = 1e-16
    if(any(mus > 1-1e-16))
        mus[mus > 1-1e-16] = 1-1e-16
#    residuals = ((M/Cov-mus)^2)/(nrow(mmMu)-ncol(mmMu)) #degrees of freedom normalized residuals
    vars = (M-mus*Cov)^2 #Variance
    one_if_binom = vars/(Cov*mus*(1-mus)) #If phi is 0 or mu is 1 or 0, then this will be ~1
    phis = one_if_binom/(Cov-1)
    if(any(is.na(phis) | phis==Inf))
        phis[is.na(phis) | phis==Inf] <- 1e-16
    if(any(phis<1e-16))
        phis[phis<1e-16] <- 1e-16
    if(any(phis>1-1e-16))
        phis[phis>1-1e-16] <- 1-1e-16
    
    #As for mu.methods.moments
#    return(c(inv_mm %*% logit(residuals)))
    return(c(inv_mm %*% logit(phis)))
}

#Log-likelihood of the beta-binomial function given a model matrix (mm) and parameters
loglik.betabin <- function(mu, phi, M, Cov, lower_threshold = 1e-16) {
    upper_threshold = 1-lower_threshold

    #If mu is ~0 then 'a' becomes 0
    if(any(mu<lower_threshold))
        mu[mu<lower_threshold] = lower_threshold
    #If mu is ~1 then 'b' becomes 0
    if(any(mu>upper_threshold))
        mu[mu>upper_threshold] = upper_threshold
    #If phi is ~1 then 'a' becomes 0
    if(any(phi > upper_threshold))
        phi[phi>upper_threshold] = upper_threshold
    #If phi is 0, then 'a' becomes Inf and 'b' becomes NaN
    if(any(phi<lower_threshold))
        phi[phi<lower_threshold] = lower_threshold

    tmp = (1/(phi))-1
    a = (mu)*tmp
    #If mu and phi are close enough to 0, then 'a' will be NaN rather than 1
    if(any(is.na(a)))
        a[is.na[a]] = 1
    b = tmp-a
    vals <- lchoose(Cov,M)-lbeta(b,a)+lbeta(Cov-M+b, M+a)
    if(any(vals %in% c(Inf, -Inf, NaN))) {
        IDX <- which(vals %in% c(Inf, -Inf, NaN))
        print(c("mu", mu[IDX]))
        print(c("phi", phi[IDX]))
        print(c("M", M[IDX]))
        print(c("Cov", Cov[IDX]))
    }
    return(-sum(vals))
}

#This is just a wrapper for the above function, since optim passes the estimates as a single vector
loglik.betabin2 <- function(x, M, Cov, mm, mmPhi) {
    end <- ncol(mm)
    mu <- x[1:end]
    phi <- x[(end+1):length(x)]
    mu <- invlogit(mm%*%mu)
    phi <- invlogit(mmPhi%*%phi)
    rv <- loglik.betabin(mu, phi, M, Cov)
    if(rv %in% c(-Inf, Inf, NaN)) {
        message("returning NaN")
        print(c("mu", mu))
        print(c("phi", phi))
        rv <- NaN
    }
    return(rv)
}

fit.BB <- function(BBM, mm, mmPhi, groupsMu, groupsPhi, inv_mmMu, inv_mmPhi, lower_threshold, momentOnly, refit=50, verbose=F) {
    M = c(assays(BBM)[["M"]])
    Cov = c(assays(BBM)[["Cov"]])
    out <- list(
        mu = c(rep(NA, ncol(mm))),
        phi = c(rep(NA, ncol(mmPhi))),
        convergence = 100, #Denoting no longer full rank!
        ll = NA,
        mm = mm,
        mmPhi = mmPhi
    )
    #Remove any samples with 0 coverage
    fullrank = ncol(mm) #Assuming the model matrix is already full-rank...
    fullrankPhi = ncol(mmPhi) #Assuming the model matrix is already full-rank...
    RMV = which(Cov==0)
    if(length(RMV)) {
        M = M[-RMV]
        Cov = Cov[-RMV]
        mm = mm[-RMV,,drop=F]
        mmPhi = mmPhi[-RMV,,drop=F]
        if(fullrank > qr(mm)$rank) return(out)
        if(fullrankPhi > qr(mmPhi)$rank) return(out)
        groupsMu = split_samples(mm)
        groupsPhi = split_samples(mmPhi)
        s = svd(mm)
        inv_mmMu = s$v %*% diag(1/s$d) %*% t(s$u)
        s = svd(mmPhi)
        inv_mmPhi = s$v %*% diag(1/s$d) %*% t(s$u)
    }

    #Get starting estimates
    #Need to switch to just using qr.solve()
    init.mu = mu.method.moments(M, Cov, groupsMu, inv_mmMu)
    init.phi = phi.method.moments(M, Cov, groupsPhi, inv_mmPhi, init.mu, mm)

    #If we're just producing method of moments approximations, then return here
    if(momentOnly) {
        out[["mu"]] <- init.mu
        out[["phi"]] <- init.phi
        out[["convergence"]] <- 0
        return(out)
    }

    #Full beta-binomial fit
    fit <- optim(c(init.mu, init.phi), loglik.betabin2, M=M, Cov=Cov, mm=mm, mmPhi=mmPhi,
        lower=logit(lower_threshold), upper=logit(1-lower_threshold),
        hessian=F, control=list(maxit=1000))
    if(fit$convergence != 0 & refit>0) { #Did we converge?
        for(i in c(1:refit)) { #Try starting at the last position
            fit <- optim(fit$par, loglik.betabin2, M=M, Cov=Cov, mm=mm, mmPhi=mmPhi,
                lower=logit(lower_threshold), upper=logit(1-lower_threshold),
                hessian=F, control=list(maxit=1000))
            if(fit$convergence==0) break
        }
    }

    #If the fit changes mu too much, then it's likely that the phi gradient pulled a coefficient in the wrong direction
    est.mus = invlogit(mm %*% init.mu)
    mus = invlogit(mm %*% fit$par[1:ncol(mm)])
    if(any(abs(mus-est.mus) > 0.5)) {
        est.phi = fit$par[(1+ncol(mm)):length(fit$par)]
        fit <- optim(c(init.mu, est.phi), loglik.betabin2, M=M, Cov=Cov, mm=mm, mmPhi=mmPhi,
            lower=logit(lower_threshold), upper=logit(1-lower_threshold),
            hessian=F, control=list(maxit=1000))
        if(fit$convergence != 0 & refit>0) { #Did we converge?
            for(i in c(1:refit)) { #Try starting at the last position
                fit <- optim(fit$par, loglik.betabin2, M=M, Cov=Cov, mm=mm, mmPhi=mmPhi,
                    lower=logit(lower_threshold), upper=logit(1-lower_threshold),
                    hessian=F, control=list(maxit=1000))
                if(fit$convergence==0) break
            }
        }
    }
    mus = invlogit(mm %*% fit$par[1:ncol(mm)])
    if(any(abs(mus-est.mus) > 0.5))
        fit$convergence=3 #for debugging

    #Fill the output list
    out[["mu"]] <- fit$par[1:ncol(mm)]
    out[["phi"]] <- fit$par[(ncol(mm)+1):length(fit$par)]
    out[["mu"]] <- as.vector(out[["mu"]]) #pretty things up a bit
    out[["phi"]] <- as.vector(out[["phi"]])
    out[["convergence"]] <- fit[["convergence"]]

    return(out)
}

#A wrapper around fit.BB, such that we can either initially split everything by row or just send a chunk to a worker and have it work on that.
fit.BB.wrapper <- function(BBM, mmMu, mmPhi, groupsMu, groupsPhi, inv_mmMu, inv_mmPhi, lower_threshold, momentOnly, refit=50) {
    if(nrow(BBM) > 1) {
        f = c(1:nrow(BBM))
        BBM2 = split(BBM, f)
        return(mapply(fit.BB, BBM2, MoreArgs=list(mm=mmMu, mmPhi=mmPhi, groupsMu=groupsMu, groupsPhi=groupsPhi, inv_mmMu=inv_mmMu, inv_mmPhi=inv_mmPhi,
            momentOnly,lower_threshold=lower_threshold, refit=refit, verbose=T), SIMPLIFY=F))
    } else {
        return(fit.BB(BBM, mmMu, mmPhi=mmPhi, groupsMu=groupsMu, groupsPhi=groupsPhi, inv_mmMu=inv_mmMu, inv_mmPhi=inv_mmPhi,
            momentOnly,lower_threshold=lower_threshold, refit=refit))
    }
}

#Fit a BBMeth object with a beta-binomial model, ideally using multiple cores
#Return a BBMethFit object
#' @rdname BBfit
#' @export
BBfit <- function(BBM, mmMu, mmPhi=NA, sameMM=T, refit=50, presplit=T, lower_threshold=1e-6, BPPARAM=MulticoreParam(workers=1), momentOnly=T) {
    if(class(BBM) != "BBMeth")
        stop("BBM must be a BBMeth object")
    if(!is.matrix(mmMu))
        stop("mmMu must be a matrix")
    if(!is.matrix(mmPhi)) {
        if(sameMM) {
            mmPhi = mmMu
        } else {
            mmPhi = matrix(c(rep(1, nrow(mmMu))), ncol=1)
            colnames(mmPhi) ="(Intercept)"
        }
    }

    #Be verbose
    BPPARAM$verbose=T

    #Ensure that we're of full rank
    if(min(dim(mmMu)) > qr(mmMu)$rank)
        stop("mmMu is not of full rank, it's likely that some of the columns are linear combinations of each other.")
    if(min(dim(mmPhi)) > qr(mmPhi)$rank)
        stop("mmPhi is not of full rank, it's likely that some of the columns are linear combinations of each other.")
    if(ncol(BBM) < qr(mmMu)$rank + qr(mmPhi)$rank)
        stop("No residual degrees of freedom!")

    #We for estimating the starting parameters, we need to invert the model
    #matrix, and determine the group membership of the various samples. It's
    #most efficient to only do that once.
    groupsMu = split_samples(mmMu)
    groupsPhi = split_samples(mmPhi)
    s = svd(mmMu)
    inv_mmMu = s$v %*% diag(1/s$d) %*% t(s$u)
    s = svd(mmPhi)
    inv_mmPhi = s$v %*% diag(1/s$d) %*% t(s$u)

    #Split the input for each worker as desired...
    if(presplit==T & BPPARAM$workers>1) {
        f = c(rep(1:BPPARAM$workers, each=ceiling(nrow(BBM)/BPPARAM$workers), length.out=nrow(BBM)))
    } else {
        if(BPPARAM$workers>1)
            message("N.B., it can take a while for all of the workers to start running. Please be patient.")
        f = c(1:nrow(BBM))
    }

    res <- bpmapply(fit.BB.wrapper, split(BBM, f), MoreArgs=list(mm=mmMu,
        mmPhi=mmPhi, refit=refit, groupsMu=groupsMu, groupsPhi=groupsPhi,
        inv_mmMu=inv_mmMu, inv_mmPhi=inv_mmPhi, lower_threshold=lower_threshold,
        momentOnly=momentOnly),
        SIMPLIFY=F, BPPARAM=BPPARAM, USE.NAMES=T)

    if(length(res) != nrow(BBM)) #Single row input or as many workers as rows
        res <- unlist(res, recursive=F)

    gr = rowData(BBM)
    mcols(gr)$convergence = c(sapply(res, function(x) return(x$convergence)))
    out = new("BBMethFit",
        seqnames = seqnames(gr),
        ranges = ranges(gr),
        strand = strand(gr),
        elementMetadata = elementMetadata(gr),
        seqinfo = seqinfo(gr),
        est.mu = matrix(unlist(lapply(res, function(x) return(x$mu))), ncol=ncol(mmMu), byrow=T),
        est.phi = matrix(unlist(lapply(res, function(x) return(x$phi))), ncol=ncol(mmPhi), byrow=T),
        mmMu = mmMu,
        mmPhi = mmPhi
    )
    colnames(out@est.mu) = colnames(mmMu)
    colnames(out@est.phi) = colnames(mmPhi)

    return(out)
}
