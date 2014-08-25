\docType{class}
\name{BBMethList-class}
\alias{BBMethList}
\alias{class:BBMethList}
\alias{BBMethList-class}
\title{BBMethList objects}
\description{The \code{BBMethList} class is a specific instance of a \code{list}
where each element of the list is a \code{BBMeth} object. You can use the
standard list-based functions on this class, such as \code{lapply}.}
\section{Coercion}{
  \describe{
    \item{}{
      \code{as(from, "BBMethList")}: Supported \code{from} types
      include:
      \describe{
        \item{"list"}{Each element of \code{from} is coerced to
          a \code{BBMeth}.}
      }
    }
  }
}
\section{"Splitting and combining"}{
  In the code snippets below, \code{BBM} is a \code{BBMeth} object and \code{BBML} a \code{BBMethList}.

  \describe{
    \item{}{
      \code{split(BBM, f)}:
      Split \code{BBM} by \code{f}. Since it's often most useful to split \code{BBM} by
      chromosome/contig, ommitting \code{f} will do this automatically.
    }
    \item{}{
      \code{BBM <- unlist(BBML)}:
      Unlist a \code{BBMethList} to create a \code{BBMeth} object. Note that the
      \code{seqLevels} of the \code{rowData} of each element must be compatible.
    }
  }
}
\seealso{
  \code{\link{coverage}}
}
\author{Devon P. Ryan}
