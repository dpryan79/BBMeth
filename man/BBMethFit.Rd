\docType{class}
\name{BBMethFit-class}
\alias{BBMethFit}
\alias{BBMeth-class}
\title{BBMethFit object}
\description{The \code{BBMethFit} class is an extension of a GRanges object with
  the addition of a number of variably-sized matrices. This class hold the
  results of fitting a BBMeth object with a beta-binomial distribution.
}
\details{The \code{BBMethFit} object holds the results of fitting methylation
data with a beta-binomial model. This includes the estimated mu and phi (raw and
shrunken) parameters as well as the log-likelihood and the model matrices. The
matrices holding the mu, phi, and shrunken_phi estimates will have one column
for each of the columns in the corresponding model matrix and one row for each
position. Position information is held as a GRanges object, which the BBMethFit
class extends. Consequently, any method for accessing position information from
a GRanges object will also work on a BBMethFit object.

The \code{convergence} codes are as follows: 0, converged properly; 1, hit the
maximum number of iterations; 10, degeneracy of the Nelder-Mead simplex; 100,
did not converge. All but 100 are directly output by \code{\link{optim}}.

The log-likelihood values are actually negative, such that values closer to 0
are larger.
}
\section{Slots}{
  \describe{
    \item{@mu}{A matrix holding the estimated mean methylation values from each
      group. Values are normally logit transformed (see \code{\link{BBfit}}).
    }
    \item{@phi}{A matrix holding the estimated over-dispersion parameters, as for
      mu above.
    }
    \item{@shrunken_phi}{A matrix holding the shrunken estimates of phi, if those
      have been calculated.
    }
    \item{@mmMu}{The model.matrix used for fitting mu.}
    \item{@mmPhi}{The model.matrix used for fitting phi.}
  }
}
\section{Accessors}{
  In addition to the standard \code{GRanges} accessors (e.g., \code{\link{length}},
\code{seqnames()}, etc.), the following are available:

  \describe{
    \item{mu(obj)}{The matrix of mu estimates.}
    \item{phi(obj)}{The matrix of phi estimates.}
    \item{shrunken_phi(obj)}{The matrix of shrunken phi estimates (if present).}
    \item{ll(obj)}{The log-likelihood values for each position. If the fit did not
 converge for this position then the value is meaningless.}
    \item{convergence(obj)}{A vector of convergence flags. See details below.}
  }

  These may also be used for replacement (e.g., mu(obj) <- values).
}
\section{Subsetting}{
  In the examples below, \code{obj} is a BBMethFit object.

  \describe{
    \item{obj[i]}{A BBMethFit object consisting of the ith position (or vector of positions).}
    \item{head(obj, n=6L)}{The first "n" elements of obj as a BBMethFit object.}
    \item{tail(obj, n=6L)}{The last "n" elements of obj as a BBMethFit object.}
    \item{window(obj, start=NA, end=NA, width=NA)}{Extract the subset of a BBMethFit object, where 2 of "start", "end", and "width" must be defined. "width" specifies the width of the window to extract.}
  }
}
\section{See also}{
  inter-range-methods
}
