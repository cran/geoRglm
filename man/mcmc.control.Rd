\name{mcmc.control}

\alias{mcmc.control}

\title{ Defines options for the MCMC-algorithm }

\description{
  This auxiliary function defines options for the MCMC-algorithm used by
  \code{\link{pois.krige.bayes}}, \code{\link{binom.krige.bayes}} and \code{\link{krige.conv}}.
}

\usage{
mcmc.control(S.scale, Htrunc, S.start, burn.in, thin, n.iter, phi.start,
             phi.scale)
}


\arguments{
\item{S.scale}{proposal variance for the update of \eqn{S} in the algorithm.  }
    \item{Htrunc}{truncation parameter for the MCMC-algorithm. Either a number or a vector of the same dimension as the data. 
     Default is Htrunc = 2*data + 5.  }
    \item{S.start}{starting value for \code{S} in the MCMC-algorithm. Default value exist.  }
    \item{burn.in}{length of the burn-in period.  Default is 0.  }
    \item{thin}{thinning: only every \code{thin} scan is stored. Default is \code{thin=10}.  }
    \item{n.iter}{number of iterations performed. Default is
      \code{n.iter=1000*thin}.  }  
    \item{phi.start}{starting value for \eqn{\phi}{phi}. Default is the median of \code{prior$phi.discrete}. }
    \item{phi.scale}{proposal variance for the update of \eqn{\phi}{phi}
in the algorithm.  }
}

\value{
  A list with processed arguments to be passed to the main function.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}}, \code{\link{binom.krige.bayes}}
and \code{\link{krige.conv}}.  }

\keyword{spatial}

