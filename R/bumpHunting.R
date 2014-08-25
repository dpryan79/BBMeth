#' @rdname correlogram
#' @export
correlogram <- function(res, min.dist=2, max.dist=500, do.plot=T) {
    #Remove NAs
    RMV = which(is.na(res$pval))
    res = res[-RMV]

    if(min.dist<2)
        min.dist=2

    vals = qnorm(res$pval, lower.tail=F)
    vars = c(rep(NA, max.dist-min.dist))
    for(i in min.dist:max.dist) {
        #Find positions a given distance away
        res2 = GRanges(seqnames=seqnames(res), ranges=IRanges(start=start(res)+i, width=width(res)))
        start(res2) = start(res)+i
        m = findOverlaps(res, res2)
        if(length(m) > 0) {
            vars[i-min.dist+1] = cor(res$deltaMu[queryHits(m)],res$deltaMu[subjectHits(m)])
        }
    }
    d = data.frame(x=c(min.dist:max.dist), y=vars)
    #fit
    fit = summary(nls(y ~ A*exp(-x/B), data=d, start=list(A=max(d$y), B=100)))$coefficients
    cutoff=qexp(0.9, rate=1/(fit[2,1]-fit[2,2]), lower.tail=T)
    print(fit)
    d2 = data.frame(x=d$x, y=fit[1,1]*exp(-d$x/fit[2,1]))
    if(do.plot) {
        if(!require(ggplot2)) {
            message("Can't load the ggplot2 packge, so no graph will be produced")
        } else {
            g = ggplot(d, aes(x=x, y=y)) + geom_point()
            g = g + xlab("Distance") + ylab("Correlation")
            g = g + geom_line(data=d2, colour="red")
            #Draw a vertical line at the clustering cutoff
            g = g + geom_vline(xintercept=cutoff, colour="red", linetype="dotted")
            print(g)
        }
    }
    return(cutoff)
}

#' @rdname bumphunting
#' @export
#' This is just a convenience wrapper around bumphunter
#' mu.threshold and phi.threshold refer to thresholds to call something interesting
#' The thresholds may also be assymetric (e.g., mu.threshold=c(-0.2, 0.1))
bumpHunting <- function(res, maxGap=NA, mu.threshold=0.1, phi.threshold=NA) {
    if(!require(bumphunter))
        stop("Sorry, the bumphunter package couldn't be loaded")
    if(is.na(mu.threshold) & is.na(phi.threshold))
        stop("Either mu.threshold or phi.threshold must be defined!")
    if(is.na(maxGap))
        maxGap = correlogram(res, do.plot=F)

    #Determine which type of bumps we're looking for
    col=2
    cutoff=phi.threshold
    if(any(is.na(phi.threshold))) {
        col=1
        cutoff=mu.threshold
    }

    #Remove NAs
    RMV = which(is.na(res$pval))
    res = res[-RMV]

    #Make some clusters
    cl = clusterMaker(as.character(seqnames(res)), start(res), assumeSorted=F, maxGap=maxGap)

    #Segment the clusters and make a nice table
    tab = regionFinder(mcols(res)[,col], as.character(seqnames(res)), start(res), cl, cutoff=cutoff)

    #Ignore regions with fewer than 3 datapoints
    RMV = which(tab$L < 3)
    tab = tab[-RMV,,drop=F]
    return(tab)
}

#Given a model matrix with unique rows, return which sample can be used as an example of which coefficient
mm2samples <- function(mm) {
    rs = rowSums(mm)
    cols = c(rep(NA, ncol(mm)))
    found = rep(0, ncol(mm))
    for(i in unique(rs)) {
        for(j in which(rs == i)) {
            for(k in which(mm[j,] == 1)) {
                if(found[k] == 0) {
                    found[k] = 1
                    cols[k] = j
                    break #Stop at the first column
                }
            }
        }
    }
    return(colnames(mm)[cols])
}

#' @rdname plotBump
#' @export
plotBump <- function(BBM, BBMF, bump, type="Mu", extend=1000, colour=1, shape=NA, show.coefs=T) {
    if(!require(ggplot2))
        stop("ggplot2 could not be loaded")
    if(class(BBM) != "BBMeth")
        stop("The first argument must be a BBMeth object.")
    if(class(BBMF) != "BBMethFit")
        stop("The second argument must be a BBMethFit object.")

    #Extend the range of the bump by 25% in each direction
    gr = GRanges(seqname=Rle(bump$chr), ranges=IRanges(start=bump$start, end=bump$end))
    seqlevels(gr) = seqlevels(BBMF)
    gr2 = gr
    if(start(gr)-extend < 1) {
        start(gr2) = 1
    } else {
        start(gr2) = start(gr)-extend
    }
    end(gr2) = end(gr) + extend

    #Find overlaps
    IDX = queryHits(findOverlaps(rowData(BBM), gr2))

    if(type=="Mu") {
        fit.vals = c(unique(BBMF@mmMu) %*% t(BBMF@mu[IDX,,drop=F]))
        mm = unique(BBMF@mmMu)
    } else {
        fit.vals = c(unique(BBMF@mmPhi) %*% t(BBMF@phi[IDX,,drop=F]))
        mm = unique(BBMF@mmPhi)
    }
    fit.vals = invlogit(fit.vals)

    #Create a dataframe
    d = data.frame(
        x = rep(start(rowData(BBM)[IDX]), each=ncol(BBM)),
        est.vals = c(t(counts(BBM[IDX])/counts(BBM[IDX], type="Cov"))),
        sample = rep(row.names(colData(BBM)), length.out=length(fit.vals))
    )
    d2 = data.frame(x=c(start(gr), start(gr), end(gr), end(gr)), y=c(0,1,1,0))

    #Add the coloring factor for the est.vals
    if(!is.numeric(colour))
        colour = which(colour %in% colnames(colData(BBM)))
    d$colour = factor(rep(colData(BBM)[,colour], length.out=length(IDX)))

    #Tweak the shape for est.vals, if needed
    if(!is.na(shape)) {
        if(!is.numeric(shape))
            shape = which(colnames(colData(BBM)) %in% shape)
        d$shape = rep(colData(BBM)[,shape], length.out=length(IDX))
    }

    #The estimated values
    d3 = data.frame(x=rep(start(rowData(BBM)[IDX]), each=ncol(mm)),
        y = fit.vals)
    #Get the coefficient association of each unique row in the model matrix
    d3$colour = factor(rep(mm2samples(mm), length.out=length(d3$y)))
    d3$colour = relevel(factor(d3$colour), colnames(mm)[1])

    #Make the graph
    g = ggplot(d, aes(x=x, y=est.vals))
    if(!is.na(shape)) {
        g = g + geom_point(aes(shape=shape, colour=colour))
        g = g + scale_shape(guide=guide_legend(title=colnames(colData(BBM))[shape]))
    } else {
        g = g + geom_point(aes(colour=colour))
    }
    g = g + scale_colour_discrete(guide=guide_legend(title=colnames(colData(BBM))[colour]))
    g = g + geom_polygon(aes(y=y), d2, alpha=0.2, stat="identity")
    if(show.coefs == T) {
        g = g + geom_line(aes(y=y, linetype=colour), d3)
        g = g + scale_linetype(guide=guide_legend(title="Coefficient"))
    }
    g = g + ylim(0,1) + xlim(start(gr2), end(gr2))
    g = g + ylab("Methylation %")
    g = g + xlab(seqnames(gr2)[1])
    print(g)
}
