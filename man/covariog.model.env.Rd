\name{covariog.model.env}

\alias{covariog.model.env}
\alias{pois.log.grf}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Envelope for Empirical Covariogram for the Poisson-log normal model}
\description{
  Computes envelope for empirical covariogram by simulating
  data for given model parameters. This function is for the Poisson-log
  normal model.  
}

\usage{ 
covariog.model.env(geodata, coords = geodata$coords, units.m = "default",
     obj.covariog, model.pars, nsim = 500, prob = c(0.025, 0.975),
     messages.screen = TRUE)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing element \code{coords} 
    as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data set.
    If not provided the argument
    \code{coords} must be given instead. 
    The list may also contain an argument \code{units.m} as described below.
  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing
    Euclidean coordinates of the \eqn{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{units.m}{ \eqn{n}-dimensional vector of observation times for the data. By default (\code{units.m = "default"}), 
  it takes \code{geodata$units.m} in case this  exist and else the value 1 for every observation.
  }
  \item{obj.covariog}{an object of the class \code{"covariogram"}, typically
    an output of the function \code{\link{covariog}}.  }
  \item{model.pars}{a list with model specification and parameter
    values. 
    The required components of the list are:
    \itemize{
      \item \code{beta}, the mean parameter. Default is \code{beta=0}.
      \item \code{cov.model}, the covariance model. Default is \code{"exponential"}.
      \item \code{cov.pars}, the covariance parameters
      \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}.
      \item \code{kappa}, the extra covariance parameters for some of the
      covariance models. Default is \code{kappa = 0.5}.
      \item \code{nugget}, the relative nugget variance. Default is \code{nugget = 0}.
    }
  }
  \item{nsim}{number of simulations used to compute the envelope. Default is \code{nsim = 500}.  }
  \item{prob}{the quantiles used for constructing the envelopes. Default is 2.5\% and 97.5\%  }
  \item{messages.screen}{logical. If \code{TRUE}, the default, status messages
    are printed while the function is running.  }
}

\details{  
  The envelope is computed assuming a Poisson-log normal model.
  Simulated values are generated at the data locations, given the model parameters.
  The empirical covariogram is computed for each simulation
  using the same binning as for the original covariogram of the data.
  The envelope is computed by taking, at each lag, the quantile-values of the covariograms for the simulated data.
}

\value{ 
  An object of the \code{\link{class}} \code{"covariogram.envelope"} which is a
  list with the components:
  \item{u}{a vector with distances. }
  \item{v.lower}{a vector with the upper-quantile covariogram values for each
    distance in \code{u}.  }
  \item{v.upper}{a vector with the lower-quantile covariogram values for each
    distance in \code{u}.  }
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{covariog}} for covariogram calculation and \code{\link{plot.covariogram}} for graphical output.  }

\examples{
\testonly{library(geoR)}
if(is.R()) data(y50)
covar <- covariog(y50, uvec = c(1:10))
parmval <- list(cov.model = "exponential", cov.pars = c(0.8,0.1), 
                beta = 1)
class(parmval) <- "covariomodel"
konvol <- covariog.model.env(y50, obj.covariog = covar, 
                             model.pars = parmval)
plot(covar, envelope.obj = konvol)
lines(parmval, max.dist = 10, lty = 1)
}


\keyword{spatial}
