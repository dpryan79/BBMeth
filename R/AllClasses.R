#' @rdname BBMeth
#' @export
setClass("BBMeth",
    contains = "SummarizedExperiment")

setValidity("BBMeth", function(object) {
    if(!("M" %in% names(assays(object)))) 
        return("The 'M' slot must contain a matrix of methylation counts!")
    if(!("Cov" %in% names(assays(object)))) 
        return("The 'Cov' slot must contain a matrix of coverage counts!")
    if(!is.numeric(assays(object)[["M"]]))
        return("The methylation counts must be numeric!")
    if(!is.integer(assays(object)[["M"]]))
        return("The methylation counts must be integers!")
    if(min(assays(object)[["M"]])<0)
        return("The methylation counts must non-negative!")
    if(!is.numeric(assays(object)[["Cov"]]))
        return("The coverage counts must be numeric!")
    if(!is.integer(assays(object)[["Cov"]]))
        return("The coverage counts must be integers!")
    if(min(assays(object)[["Cov"]])<0)
        return("The coverage counts must non-negative!")
    if(dim(assays(object)[["M"]]) != dim(assays(object)[["Cov"]])) 
        return("The methylation and coverage matrices must have the same dimensions!")
    if(length(rowData(object)) != nrow(assays(object)[["M"]]))
        return("The length of the coordinates must match the number of rows in the count matrices!")
    if(ncol(colData(object)) < 1)
        return("colData must not be empty!")
    TRUE
})

#' BBMeth object and constructor
#'
#' The \code{BBMeth} class is a subclass of \code{SummarizedExperiment}. The
#' assay data for this class consists of two matrices of non-negative integers
#' named "M" and "Cov". "M" holds methylation counts and "Cov" holds coverage
#' counts per position. Position information is held in a GRanges object as the
#' rowData.
#'
#' @param M a matrix of non-negative integer methylation counts
#' @param Cov a matrix of non-negative integer coverage counts
#' @param gr A \code{GRanges-class} denoting the position associated with each
#' row of the count matrices.
#' @param groups a \code{data.frame} or \code{DataFrame} containing the group
#' membership of each sample. Rows are samples and columns parameters
#' @param sampleTable a data.frame with at least 3 columns. The first contains
#' the sample names, one per row. The second column contains the filename
#' associated with each sample. The file should be in bedGraph format with the
#' 5th column containing a methylated C counts and the 6th column containing
#' unmethylated C counts (i.e., the count of Ts). The third and subsequent
#' columns denote the group membership of each samples (e.g., WT vs. Mutant,
#' Treated vs Control, Male vs Female, etc.). Any factor or covariate that you
#' wish to use in your statistical model must be present here.
#' @param BPPARAM A \code{BiocParallelParam} instance, which defaults to
#' \code{MulticoreParam(workers=1)} if not specified. See \code{bplapply} for details.
#'
#' @return A BBMeth object.
#'
#' @aliases BBMeth, BBMeth-class
#'
#' @docType class
#'
#' @rdname BBMeth
#' @export
BBMeth <- function(M, Cov, gr, groups) {
    if(!is.matrix(M))
        stop("M must be a matrix!")
    if(!is.matrix(Cov))
        stop("Cov must be a matrix!")
    if(is.data.frame(groups)) {
        groups <- as(groups,"DataFrame")
    } else if(!(class(groups) == "DataFrame")) {
        stop("groups must be a data.frame!")
    }

    BBM <- SummarizedExperiment(assays=SimpleList(M=M, Cov=Cov), rowData=gr, colData=groups)
    BBM <- as(BBM, "BBMeth")
    colnames(BBM) <- colnames(M)
    return(BBM)
}

#' @rdname BBMeth
#' @export
BBMethFromBedGraphs <- function(sampleTable, BPPARAM=MulticoreParam(workers=1)) {
    if(ncol(sampleTable) < 3) {
        stop("sampleTable must have at least 3 columns (samples names, file names, group affiliations)")
    }

    message("Importing bedGraphs...")
    files <- bplapply(as.character(sampleTable[,2]), function(x) rtracklayer::import.bedGraph(x), BPPARAM=BPPARAM)

    #Take care of any mismatch between covered contigs (this occurs frequently when _random or _patch contigs
    snames = character()
    for(i in c(1:length(files)))
        snames <- union(snames,as.character(seqlevels(files[[i]])))
    for(i in c(1:length(files)))
        seqlevels(files[[i]]) <- snames

    #Create the final GRanges object
    message("Merging ranges...")
    gr = files[[1]]
    for(i in c(2:length(files))) 
        gr = unique(c(gr, files[[i]]))
    mcols(gr) <- NULL

    #Make the count matrices
    message("Combining counts...")
    M = matrix(c(0), nrow=length(gr), ncol=length(files))
    Cov = matrix(c(0), nrow=length(gr), ncol=length(files))
    for(i in seq(along=files)) {
        m <- findOverlaps(gr, files[[i]])
        q <- queryHits(m)
        s <- subjectHits(m)
        M[q, i] <- mcols(files[[i]])[s,2]
        Cov[q, i] <- mcols(files[[i]])[s,2]+mcols(files[[i]])[s,3]
    }
    colnames(M) <- as.character(sampleTable[,1])
    colnames(Cov) <- as.character(sampleTable[,1])

    #Make the groups
    groups = sampleTable[,-c(1,2), drop=FALSE]

    return(BBMeth(M=M, Cov=Cov, gr=gr, groups=groups))
}

#' @rdname BBMethFit
#' @export
setClass("BBMethFit",
    contains="GRanges",
    representation(
        est.mu="matrix",
        est.phi="matrix",
        mu="matrix",
        phi="matrix",
        se.mu="matrix",
        se.phi="matrix",
        mmMu="matrix",
        mmPhi="matrix",
        prior="vector"
    )
)

#' @rdname BBMethList
setClass("BBMethList",
    prototype = prototype(elementType = "BBMeth"),
    contains = "list")
