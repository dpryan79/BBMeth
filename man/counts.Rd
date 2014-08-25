\name{counts}
\alias{counts}
\title{Extract matrix of counts from a BBMeth object.}
\description{This function accepts a BBMeth object and optionally a "type" and
 returns the matrix of raw counts, with rows as positions and columns as
 samples.}
\usage{
  counts(BBM, type="M")
}
\arguments{
  \item{BBM}{
    An BBMeth object.
  }
  \item{type}{
    Either "M" or "Cov". "M" (the default), will produce the count matrix of
    methylated Cs. "Cov" will produce the coverage matrix.
  }
}
\value{
  A numeric matrix with rows as positions and columns as samples.
}
\examples{
  #Filter a BBMeth object, keeping rows with a minimum average 4x depth
  IDX <- which(rowSums(counts(BBM, type=Cov)) >= 4*ncol(BBM))
  BBM.filtered <- BBM[IDX]
}
