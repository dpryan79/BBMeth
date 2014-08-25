\docType{class}
\name{BBMeth-class}
\alias{BBMeth}
\alias{BBMeth-class}
\alias{BBMethFromBedGraphs}
\title{BBMeth object and constructors}
\usage{
BBMeth(M, Cov, gr, groups)

BBMethFromBedGraphs(sampleTable, BPPARAM=bpparam())
}
\arguments{
  \item{M}{a matrix of non-negative integer methylation counts.}

  \item{Cov}{a matrix of non-negative integer coverage counts.}

  \item{gr}{a \code{GRanges} object holding the position information for each
  row of the aforementioned matrices.}

  \item{groups}{a \code{data.frame} or \code{DataFrame} denoting the group/treatment
  affiliation of each column of \code{M} and \code{Cov}. This must contain at least
  one column and have as many rows as the count matrices have columns.}

  \item{sampleTable}{a \code{data.frame} with three or more columns. The first
  column holds sample names and should either be a character vector or otherwise
  coercible to one. The second column is the file name of a bedGraph file
  holding methylation metrics for that sample. The 5th column of each bedGraph
  file should contain the count of methylated Cs at a given position and the 6th
  column the count of unmethylated Cs (i.e., Ts). The 4th column of the bedGraph
  files are ignored. The third and subsequent columns of sampleTable should hold
  group affiliations for each sample (e.g., WT vs. Mutant, Treated vs Control,
  Male vs Female, etc.). Any factor or covariate that you wish to use in your
  statistical model must be present here.}

  \item{BPPARAM}{a \code{BiocParallelParam} instance, which defaults to
  \code{MulticoreParam(workers=1)}. See \code{?bplapply} for details.}
}
\value{
A BBMeth object.
}
\description{The \code{BBMeth} class is a subclass of \code{SummarizedExperiment}
intended to hold the original methylation and coverage counts, position
information and sample affiliation of a bisulfite-sequencing experiment. The assay data for this class consists of
two matrices of non-negative integers named "M" and "Cov". "M" holds methylation
counts and "Cov" holds coverage counts per position. Position information is
held in a GRanges object as the rowData. Sample affiliation with treatment
groups and covariates is held in the colData. Methods useful for working with a
\code{SummarizedExperiment} will also work with a \code{BBMeth} object.

The constructor functions can take either a set of matrices with an associated
GRanges object and a data.frame or a data.frame listing bedGraph filenames, as
described above.
}
\section{Warning}{
Loading multiple whole-genome bedGraph files can require quite a bit of memory.
If you receive an error of the form "negative widths are not allowed" when using
\code{BBMethFromBedGragphs} then it's likely that you've simply run out of
memory. Since a \code{BBMeth} object is just a subclass of \code{SummarizedExperiment},
you can merge multiple of them with \code{cbind}.
}
\examples{
????
}
