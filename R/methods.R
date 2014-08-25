#' Accessors for the 'est.phi' slot in a BBMethFit object.
#'
#' The est.phi slot holds non-shrunken over-dispersion parameter. Each row
#' represents a single observation (generally a single C or CpG) and one column
#' for each column in the phi model matrix.
#'
#' @usage
#' \S4method{est.phi}{BBMethFit}(object)
#'
#' \S4method{est.phi}{BBMethFit, matrix}(object)<-value
#'
#' @docType methods
#' @name est.phi
#' @rdname est.phi
#' @aliases est.phi, est.phi,BBMethFit-method est.phi<-,BBMethFit,matrix-method
#'
#' @param object a \code{BBMethFit} object.
#' @param value a real (floating point) matrix
#' @author Devon Ryan
#'
est.phi.BBMethFit <- function(object) return(object@est.phi)

#' @rdname est.phi
#' @export
setMethod("est.phi", signature(object="BBMethFit"), est.phi.BBMethFit)

#' @name est.phi
#' @rdname est.phi
#' @exportMethod "est.phi<-"
setReplaceMethod("est.phi", signature(object="BBMethFit", value="matrix"),
  function( object, value ) {
   object@est.phi <- value
   object
})

#phi
phi.BBMethFit <- function(object) return(object@phi)
setMethod("phi", signature(object="BBMethFit"), phi.BBMethFit)
setReplaceMethod("phi", signature(object="BBMethFit", value="matrix"),
  function( object, value ) {
   object@phi <- value
   object
})

#est.mu
est.mu.BBMethFit <- function(object) return(object@est.mu)
setMethod("est.mu", signature(object="BBMethFit"), est.mu.BBMethFit)
setReplaceMethod("est.mu", signature(object="BBMethFit", value="matrix"),
  function( object, value ) {
   object@est.mu <- value
   object
})

#Mu
mu.BBMethFit <- function(object) return(object@mu)
setMethod("mu", signature(object="BBMethFit"), mu.BBMethFit)
setReplaceMethod("mu", signature(object="BBMethFit", value="matrix"),
  function( object, value ) {
   object@mu <- value
   object
})

#' @rdname coverage
#' @export
setMethod(coverage, signature(x="BBMeth"), function(x) data.frame(All=colSums(assays(x)[["Cov"]])/nrow(x)))
setMethod(coverage, signature(x="BBMethList"), function(x) as.data.frame(do.call(cbind, lapply(x, function(y) colSums(assays(y)[["Cov"]])/nrow(y)))))

#' @rdname split
#' @rdname unlist
#' @export
setMethod(split, signature(x="BBMeth", f="ANY"),
    function(x,f,drop) {
        x = as(x, "SummarizedExperiment")
        out = lapply(split(x, f, drop=drop), function(x) as(x, "BBMeth"))
        out = as(out, "BBMethList")
    }
)
setMethod(split, signature(x="BBMeth", f="missing"),
    function(x,drop) {
        f = as.factor(seqnames(rowData(x)))
        x = as(x, "SummarizedExperiment")
        out = lapply(split(x, f, drop=drop), function(x) as(x, "BBMeth"))
        out = as(out, "BBMethList")
    }
)
setMethod(unlist, signature(x="BBMethList"),
    function(x, recursive, use.names) {
        x = as(x, "list")
        return(do.call(rbind, x))
    }
)
setMethod(split, signature(x="BBMethFit", f="ANY"),
    function(x,f) {
        if(length(f) != length(x)) {
            f2 = rep(f, each=ceiling(length(x)/length(f)), length.out=length(x))
            fl = split(c(1:length(x)), f2)
        } else{
            fl = split(c(1:length(x)), f)
        }
        out = lapply(fl, function(f3) x[f3])
        return(out)
    }
)
        

#' @rdname [
#' @export
setMethod(`[`, signature(x="BBMeth", i="numeric"),
    function(x,i) {
        out <- SummarizedExperiment(assays=SimpleList(M=assays(x)[[1]][i,,drop=F], Cov=assays(x)[[2]][i,,drop=F]),
            rowData=rowData(x)[i], colData=colData(x), exptData=exptData(x))
        out = as(out, "BBMeth")
        return(out)
    }
)
setMethod(`[`, signature(x="BBMeth", i="logical"),
    function(x,i) {
        i2 = which(i==T)
        out <- SummarizedExperiment(assays=SimpleList(M=assays(x)[[1]][i2,,drop=F], Cov=assays(x)[[2]][i2,,drop=F]),
            rowData=rowData(x)[i2], colData=colData(x), exptData=exptData(x))
        out = as(out, "BBMeth")
        return(out)
    }
)
setMethod(`[`, signature(x="BBMethFit", i="numeric"),
    function(x,i) {
        out = new("BBMethFit",
            seqnames = seqnames(x)[i],
            ranges = ranges(x)[i],
            strand = strand(x)[i],
            elementMetadata = elementMetadata(x)[i,,drop=F],
            seqinfo = seqinfo(x),
            est.mu = est.mu(x)[i,,drop=F],
            est.phi = est.phi(x)[i,,drop=F],
            mmMu = x@mmMu,
            mmPhi = x@mmPhi,
            prior = x@prior
        )
        if(nrow(phi(x)) > 0) {
            phi(out) = phi(x)[i,,drop=F]
            mu(out) = mu(x)[i,,drop=F]
        }
        return(out)
    }
)
setMethod(`[`, signature(x="BBMethFit", i="logical"),
    function(x,i) {
        i2 = which(i==T)
        out = new("BBMethFit",
            seqnames = seqnames(x)[i2],
            ranges = ranges(x)[i2],
            strand = strand(x)[i2],
            elementMetadata = elementMetadata(x)[i2,,drop=F],
            seqinfo = seqinfo(x),
            est.mu = est.mu(x)[i2,,drop=F],
            est.phi = est.phi(x)[i2,,drop=F],
            mmMu = x@mmMu,
            mmPhi = x@mmPhi,
            prior = x@prior
        )
        if(nrow(phi(x)) > 0) {
            phi(out) = phi(x)[i,,drop=F]
            mu(out) = mu(x)[i,,drop=F]
        }
        return(out)
    }
)

#' @rdname head
#' @export
setMethod(head, signature(x="BBMeth"), 
    function(x,n=6L) {
        l = length(x)
        if(n>=l) {
            return(x[c(1:n)])
        } else {
            return(x[c(1:l)])
        }
    }
)

#' @rdname tail
#' @export
setMethod(tail, signature(x="BBMeth"), 
    function(x, n=6L) {
        l = length(x)
        if(n>=l) {
            return(x)
        } else {
            return(x[c((l-n+1):l)])
        }
    }
)

#' @rdname window
setMethod(window, signature(x="BBMeth"), 
    function(x, start=NA, end=NA, width=NA) {
        if(sum(is.na(c(start,end,width))) >= 2)
            stop("At least 2 of 'start', 'end' and 'width' must be defined.")
        if(!is.na(start)) { #(start&end) | (start&width)
            if(is.na(end))
                end=start+width-1
        } else { #width&end
            start=end-width+1
        }
        l = length(x)
        if(start>x | end>x)
            stop("Invalid sequence coordinates.
  Please make sure the supplied 'start', 'end' and 'width' arguments
  are defining a region that is within the limits of the sequence.")
        return(x[c(start:end)])
    }
)

#' @rdname counts
setMethod(counts, signature(object="BBMeth"),
    function(object, type="M") {
        if(type != "M" & type != "Cov")
            stop("'type' must be either 'M' or 'Cov'")
        return(assays(object)[[type]])
    }
)
