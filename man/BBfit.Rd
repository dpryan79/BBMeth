\name{BBfit}
\alias{BBfit}
\alias{fit.BB}
\title{Fit a beta-binomial model to a BBMeth object, returning a BBMethFit object.}
\description{This function accepts a BBMeth object and one or two model matrices
 and fits a beta-binomial model accordingly. The returned object contains the
 fit mu and phi values, which are logit transformed.}
\usage{
  BBfit(BBM, mmMu, mmPhi=NA, sameMM=TRUE,
    refit=50, presplit=T,
    BPPARAM=MulticoreParam(workers=1))
}
\arguments{
  \item{BBM}{
    An BBMeth object.
  }
  \item{mmMu}{
    The model matrix for use with fitting mu, the expected methylation value.
  }
  \item{mmPhi}{
    A model matrix use for fitting phi, the overdispersion parameter. If this is
    omitted it will be the same as that used for mmMu. However, see the sameMM
    option below.
  }
  \item{sameMM}{
    If this is set to FALSE and mmPhi is omitted, then a model matrix consisting
    of a single column of ones will be used for fitting phi. This has the effect
    of setting phi to a constant, regardless of how many groups are present. If
    mmPhi is not omitted, then this option has no effect.
  }
  \item{refit}{
    An integer >= 0, specifying how many times to try continue trying to refit
    when convergence does't occur.
  }
  \item{presplit}{
    If True and multiple cores/CPUs/computers are being used, the dataset will
    be split (more or less) evenly by the number of worker cores/CPUs/computers.
    Otherwise, the dataset is split by position and sent to workers as they have
    free resources. If the cost of communicating short messages between workers
    is very low or the time required per fit is very erratic then it may be
    beneficial to set presplit to F.
  }
  \item{BPPARAM}{
     A "BiocParallelParam" instance. See \code{?bplapply} for details.
  }
}
\details{
  The mean (mu) and overdispersion (phi) of each position is determined by
  fitting with a beta-binomial distribution. The model matrix used to fit these
  two parameters needn't be the same, which facilitates looking for changes in
  variability separately from mean methylation. The returned fits are normally
  still logit transformed, which is convenient for downstream processing.
}
\value{
  BBfit returns a \code{\link{BBMethFit}} object.
}
\examples{
  ????
  Examples go here.
}
