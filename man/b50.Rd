\name{b50 and y50}

\alias{b50}
\alias{y50}

\non_function{b50, y50}

\title{Simulated Data sets which Illustrate the Usage of the Package geoRglm}

\usage{
data(b50)
data(y50)
}

\description{The simulated data sets \code{b50} and \code{y50} are used
  in several examples throughout the package
  documentation.
}

\format{
  The object is a list
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the simulated data from the binomial-logit normal
      model for \code{b50}, and the simulated data from the Poisson-log normal
      model for \code{y50}.  }
    \item{\code{units.m}}{\eqn{n}-dimensional vector giving the number
      of trials for \code{b50} and the observation times
      for \code{y50}.  }
    \item{\code{cov.model}}{the correlation model used for the simulation.  }
    \item{\code{nugget}}{the values of the \code{nugget} parameter used for the simulation.  }
    \item{\code{cov.pars}}{the covariance parameters used for the simulation.  }
    \item{\code{kappa}}{the value of the smoothness parameter \eqn{\kappa}{kappa} used for the simulation.  }
  }
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\keyword{datasets}
