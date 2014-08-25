#' @rdname variogram
#' @export, BPPARAM=BPPARAM
#This is technically the semivariogram
variogram <- function(res, maxGap=100, BPPARAM=MulticoreParam(workers=1), span=0.5, degree=2, minCount=10) {
    if(!require(bumphunter))
        stop("Sorry, the bumphunter package can't be loaded")
    if(maxGap < 2)
        stop("maxGap must be at least 2, otherwise you're only comparing a point
            with itself")
    end(res) = start(res)
    RMV = which(is.na(res$pval))
    res = res[-RMV]
#    res$Z = qnorm(res$pval, lower.tail=F) #make significant p-values positive
    cl = clusterMaker(seqnames(res), start(res), assumeSorted=F, maxGap=maxGap)
    res = split(res, cl)
    #Precompute distance matrices
    ms = bplapply(res, function(x) {
        pos=start(x)-min(start(x))+1
        out=as.matrix(dist(pos))
        out[lower.tri(out)]=NA
        return(out)
    }, BPPARAM=BPPARAM)
    f2 <- function(res, m, i=2) {
        nr = nrow(m)
        IDX = which(m==i)
        if(length(IDX)) { #set a minimum
            query=IDX%%nr #the row
            subject=ceiling(IDX/nr) #the column
            return(res$W[query]-res$W[subject])
        }
    }
    f <- function(i, res, ms) {
        diffs = unlist(mapply(f2, res, ms, MoreArgs=list(i=i)))
        if(length(diffs)>minCount)
            return(median((diffs^2)/0.91))
        return(NA)
    }
    maxh = max(unlist(bplapply(res, function(x) max(start(x)-min(start(x))), BPPARAM=BPPARAM)))
    h = 2:maxh
    vars = unlist(bplapply(h, f, res=res, ms=ms, BPPARAM=BPPARAM))

    #Format things for printing and regression
    d = data.frame(h=h, vars=vars)
    RMV = which(is.na(d$vars)) #Remove NAs
    if(length(RMV))
        d = d[-RMV,,drop=F]
    d$est = loess(vars~h, data=d, span=span, degree=degree)$fitted
    with(d, plot(vars~h, xlab="h", ylab=expression(paste(gamma, "(h)"))))
    return(d)
}

makeNull <- function(BBMF1, BBMF2, BBM, samples=NA, maxGap=100, nSites=10000, BPPARAM=MulticoreParam(workers=1), makeVariogram=T) {
    BBMF1 = BBMF1[c(1:nSites)]
    BBMF2 = BBMF2[c(1:nSites)]
    BBM = BBM[c(1:nSites)]
    prior1 = BBMF1@prior
    prior2 = BBMF2@prior
    if(all(is.na(samples))) {
        samples = sample(1:nrow(BBMF1@mmMu), nrow(BBMF1@mmMu))
        message(sprintf("The new model matrix row order is: %s", paste(samples, collapse=",")))
    }

    message("Performing initial method of moments reestimation...")
    BBMF1 = BBfit(BBM, BBMF1@mmMu[samples,], mmPhi=BBMF1@mmPhi[samples,], BPPARAM=BPPARAM, momentOnly=T)
    BBMF2 = BBfit(BBM, BBMF2@mmMu[samples,], mmPhi=BBMF2@mmPhi[samples,], BPPARAM=BPPARAM, momentOnly=T)
    BBMF1@prior = prior1
    BBMF2@prior = prior2
    message("Incorporating previous prior...")
    BBMF1 = shrinkPhi(BBMF1, BBM, BPPARAM=BPPARAM)
    BBMF2 = shrinkPhi(BBMF2, BBM, BPPARAM=BPPARAM)

    resNull = LRTest(BBMF1, BBMF2)
    if(makeVariogram)
        variogram(resNull, maxGap=maxGap, BPPARAM=BPPARAM)

    return(resNull)
}

#Should shrink mu as well
#Use a dataset where there actually is a difference!!!
#make a graph of how the new group are made of the old groups
