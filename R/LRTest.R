LRTcoefs <- function(full, reduced) {
    #Determine which columns to look at
    idxMu = which(!(colnames(full@mmMu) %in% colnames(reduced@mmMu)))
    idxPhi = which(!(colnames(full@mmPhi) %in% colnames(reduced@mmPhi)))
    if(length(c(idxMu,idxPhi)) == 0)
        stop("The models used by the two objects seem to be identical")

    out = matrix(NA, ncol=2, nrow=length(full))
    Numerator = rep(0, length(full))
    Denominator = rep(0, length(full))

    if(length(idxMu) > 0) {
        rs = rowSums(full@mmMu)
        idx2 = which(apply(full@mmMu[,idxMu,drop=F], 1, function(x) all(x == 1)))
        final_idx = which(rs == min(rs[idx2]))
        if(length(final_idx) == 0) {
            if(length(idxMu) == 1) {
                message("It seems that the coefficient isn't estimable. Returning NAs for it")
            } else {
                message("It seems that the coefficients aren't estimable. Returning NAs for them.")
            }
        } else {
            full_value = invlogit(mu(full) %*% t(full@mmMu[final_idx[1],,drop=F]))
            reduced_value = invlogit(mu(reduced) %*% t(reduced@mmMu[final_idx[1],,drop=F]))
            out[,1] = c(full_value - reduced_value)
        }
        Numerator = Numerator + c(rowSums(full@mu[,idxMu,drop=F]))
        Denominator = Denominator + c(rowSums(full@se.mu[,idxMu,drop=F]))
    }
    if(length(idxPhi) > 0) {
        rs = rowSums(full@mmPhi)
        idx2 = which(apply(full@mmPhi[,idxPhi,drop=F], 1, function(x) all(x == 1)))
        final_idx = which(rs == min(rs[idx2]))
        if(length(final_idx) == 0) {
            if(length(idxPhi) == 1) {
                message("It seems that the phi coefficient isn't estimable. Returning NAs for it")
            } else {
                message("It seems that the phi coefficients aren't estimable. Returning NAs for them.")
            }
        } else {
            full_value = invlogit(phi(full) %*% t(full@mmPhi[final_idx[1],,drop=F]))
            reduced_value = invlogit(phi(reduced) %*% t(reduced@mmPhi[final_idx[1],,drop=F]))
            out[,2] = c(full_value - reduced_value)
        }
        Numerator = Numerator + c(rowSums(full@phi[,idxPhi,drop=F]))
        Denominator = Denominator + c(rowSums(full@se.phi[,idxPhi,drop=F]))
    }
    W = Numerator/Denominator
    W[is.nan(W)] = NA
    out = cbind(out, W)

    #Ensure that rows where either model didn't converge are NA
    IDX = which(!(full$convergence %in% c(0,2)) | !(reduced$convergence %in% c(0,2)))
    out[IDX,] = NA

    return(out)
}

#' @rdname LRTest
#' @export
LRTest <- function(full, reduced, method="BH", coef=NA, type="mu") {
    if(class(full) != "BBMethFit")
        stop("The full model must be a BBMethFit object")
    if(class(reduced) != "BBMethFit")
        stop("The reduced model must be a BBMethFit object")
    if(length(full) != length(reduced))
        stop("The objects have different numbers of rows")

    #Likelihood ratio
    d = 2*(full$ll - reduced$ll)

    #Determine the degrees of freedom
    dof = (qr(full@mmMu)$rank+qr(full@mmPhi)$rank) - (qr(reduced@mmMu)$rank+qr(reduced@mmPhi)$rank)
    if(dof == 0)
        stop("No residual degrees of freedom!")

    #For cases where at least one of the groups didn't converge, set results to NA
#    RMV = which(!full$convergence %in% c(0,2) | !reduced$convergence %in% c(0,2))
    RMV = which(full$convergence != 0 | reduced$convergence != 0)

    pvals = 1 - pchisq(d, df=dof)
    pvals[RMV] = NA
    padj = p.adjust(pvals, method=method)

    #Make the final GRanges object
    changes = LRTcoefs(full, reduced)
    out = GRanges(seqnames=seqnames(full), ranges=ranges(full), strand=strand(full))
    mcols(out)$deltaMu = changes[,1]
    mcols(out)$deltaPhi = changes[,2]
    mcols(out)$test_statistic = d
    mcols(out)$W = changes[,3]
    mcols(out)$pval = pvals
    mcols(out)$padj = padj
    mcols(out)[RMV,] = NA #Ignore rows that didn't converge
    return(out)
}
