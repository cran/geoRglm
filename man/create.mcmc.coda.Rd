\name{create.mcmc.coda}

\alias{create.mcmc.coda}


\title{Create an mcmc object for the CODA package}
\description{
  This function creates an mcmc object for the CODA package for
  output from the functions \code{glsm.mcmc}, \code{binom.krige.bayes} and \code{pois.krige.bayes}. The
  functions in CODA can then be used to investigate convergence and
  mixing of the MCMC algorithm.
}

\usage{
create.mcmc.coda(x, mcmc.input)
}

\arguments{
  \item{x}{an output object from \code{glsm.mcmc}, \code{binom.krige.bayes} or \code{pois.krige.bayes}.  }
  \item{mcmc.input}{input parameters for the MCMC algorithm. It can take an output from \code{mcmc.control} or
    a list with elements. Only \code{thin} and \code{burn.in} are
    used, and both have default values (\code{thin=10}, \code{burn.in=0}).   }
}


\value{
  An object of class \code{mcmc} to be used for CODA.
}

\references{ 
Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}


\examples{
## se example in help file for glsm.mcmc
}

\keyword{spatial}













