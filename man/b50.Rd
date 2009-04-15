\name{b50, p50 and b64}

\encoding{latin1}

\alias{b50}
\alias{p50}
\alias{b64}

\docType{data}

\title{Simulated Data sets which Illustrate the Usage of the Package geoRglm}

\usage{
data(b50)
data(p50)
data(b64)
}

\description{The simulated data sets \code{b50} and \code{p50} are used
  in several examples throughout the package
  documentation. The simulated data set \code{b64} was used in Diggle,
  Ribeiro and Christensen (2003), and Christensen and Ribeiro (2002).
}

\format{
  The object is a list
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the simulated data from the binomial-logit normal
      model for \code{b50} and \code{b64}, and the simulated data from the Poisson-log normal
      model for \code{p50}.  }
    \item{\code{units.m}}{\eqn{n}-dimensional vector giving the number
      of trials for \code{b50}  and \code{b64}, and the observation times
      for \code{p50}.  }
    \item{\code{cov.model}}{the correlation model used for the simulation.  }
    \item{\code{nugget}}{the values of the \code{nugget} parameter used for the simulation.  }
    \item{\code{cov.pars}}{the covariance parameters used for the simulation.  }
    \item{\code{kappa}}{the value of the smoothness parameter \eqn{\kappa}{kappa} used for the simulation.  }
  }
}

\source{Simulated data sets.}


\references{
Peter J. Diggle, Paulo J. Ribeiro and Ole F. Christensen (2003). An
introduction to model-based geostatistics. In : \emph{Spatial statistics and computational methods} (ed. J. Møller),
Springer Verlag, 43-86.

Christensen, O. F. and Ribeiro Jr, P. J. (2002). geoRglm - a package for generalised linear spatial models. \emph{R News}, 2(2), 26-28.
  
Further information about \bold{geoRglm} can be found at:\cr
\url{http://gbi.agrsci.dk/~ofch/geoRglm}.
}

\keyword{datasets}
