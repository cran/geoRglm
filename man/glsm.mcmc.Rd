\name{glsm.mcmc}

\alias{glsm.mcmc}


\title{Conditional Simulation for a generalised linear spatial model}
\description{
  This function performs conditional simulation (by MCMC)
  in a generalised linear spatial model
  for fixed parameters.
}

\usage{
glsm.mcmc(geodata, coords = geodata$coords, data = geodata$data, 
         units.m = "default",  model, mcmc.input, messages)
}

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \bold{geoR} data set.
    If not provided the arguments \code{coords} and \code{data} must be provided instead.  
    The list may also contain an argument \code{units.m} as described
    below.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  } 
  \item{units.m}{ \eqn{n}-dimensional vector of observation times for the data. By default (\code{units.m = "default"}), 
    it takes \code{geodata$units.m} in case this exist and else a vector of 1's.  }
  \item{model}{defines the model components. Either an object of
    \code{class} \code{likGLSM}; typically output from
    \code{likfit.glsm}, or a list containing
    the arguments : 
    \cr
    \itemize{
      \item{trend}{specifies the trend (covariate) values at the data
	locations. See documentation of \code{\link[geoR]{trend.spatial}}
	for further details.
	Default is \code{trend = "cte"}.   }
      \item{beta}{numerical value of the mean (vector) parameter.  }
      \item{cov.model}{string indicating the name of the model for the
	correlation function. Further details  in the
	documentation for \code{\link[geoR]{cov.spatial}}.  }
      \item{cov.pars}{a vector with the 2 covariance parameters \eqn{\sigma^2}{sigma^2},
	and \eqn{\phi}{phi} for the underlying Gaussian field.  }
      \item{kappa}{additional smoothness parameter required by the following correlation
	functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} and
	\code{"gneiting.matern"}.   }
      \item{nugget}{the value of the nugget parameter
	\eqn{\tau^2}{tau^2} for the underlying Gaussian field. Default is
	\code{nugget = 0}.  }
      \item{aniso.pars}{parameters for geometric anisotropy
	correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
	a two elements vector with values for the anisotropy parameters
	must be provided. Anisotropy correction consists of a
	transformation of the data and prediction coordinates performed
	by the function \code{\link[geoR]{coords.aniso}}.  }
      \item{family}{equal to either \code{"poisson"} or
	\code{"binomial"}  }
      \item{link}{equal to either \code{"canonical"} (default), \code{"log"},
	\code{"boxcox"} or \code{"logit"}. For \code{"canonical"} then
	in general
	the canonical link function is used (\code{"log"} for the
	Poisson distribution and \code{"logit"} for the binomial distribution), but when
	\code{lambda} is also specified then the Box-Cox class is used
	(a mis-use of the terminology "canonical", really).  }
      \item{lambda}{numeric value of the Box-Cox transformation parameter.
	The value \eqn{\lambda = 1}{lambda = 1} corresponds to
	no transformation and the default value \eqn{\lambda = 0}{lambda = 0} corresponds to
	the log-transformation.  Only used when \code{family = "poisson"} }
    } 
  }
  \item{mcmc.input}{input parameter for the MCMC algorithm. It can take an output from \code{mcmc.control} or
    a list with elements as for the arguments in
    \code{\link{mcmc.control}}. See documentation for
    \code{\link{mcmc.control}}.\cr ATTENTION: the argument \code{S.scale} is necessary while
    all the others have default values.   }
  \item{messages}{logical. Indicates whether or not status messages are printed on the screen (or other output device)
    while the function is running.  }
}


\details{
  For simulating the conditional distribution of \eqn{S} given \eqn{y}, the Langevin-Hastings algorithm 
with the parametrisation in Papaspilliopoulus, Roberts and Skold (2003)
is used. This algorithm is a Metropolis-Hastings algorithm, where the
proposal distribution uses gradient information from the
log-posterior distribution. 

The proposal variance (called \code{S.scale}; see \code{\link{mcmc.control}})
for the algorithm needs to be scaled
such that approximately 60 percent of the proposals are accepted. We
also recommend that the user to studies plots of the autocorrelations.
 
}

\value{
  A list with the following components:
  
  \item{simulations}{an \eqn{n \times n.sim}{n x n.sim} matrix with
    \eqn{n.sim} being the number of MCMC simulations. containing \eqn{S_i}. Each column 
    corresponds to a conditional simulation of the conditional distribution of \eqn{S_i} at the data locations. 
    }
  \item{acc.rate}{matrix with acceptance rates from MCMC.
    Only returned when no prediction locations are given.  }
  \item{model}{Information about the model parameters, link function and
    error distribution family used.  }
  \item{geodata}{Information about the data.  }
  \item{call}{the function call.  }
}

\references{
  O. Papaspiliopoulus and G. O. Roberts and M. Skold
  (2003). Non-centered parameterizations for hierarchical models and
  data augmentation. \emph{Bayesian statistics 7} (eds. J. M. Bernardo,
  S. Bayarri, J. O. Berger, A. P. Dawid, D. Heckerman, A. F. M. Smith
  and M. West), Oxford University Press, 307-326.

  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://gbi.agrsci.dk/~ofch/geoRglm}.
}

\author{
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{binom.krige}} for prediction with fixed parameters in the
  Binomial-normal model, \code{\link{pois.krige}} for prediction with fixed parameters in the
  Poisson normal model.
}

\examples{
\dontshow{library(geoR)}
if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) set.seed(1234)
data(b50)
test <- glsm.mcmc(b50, model = list(family="binomial",
             cov.pars = c(1,1), beta = c(1,0), trend =~ rnorm(50),
             cov.model="spherical", nugget=0.3),
          mcmc.input = mcmc.control(S.scale = 0.2, thin = 1))
## visulalising the MCMC output using the coda package
test.coda <- create.mcmc.coda(test, mcmc.input = list(thin = 1))
\dontrun{
plot(test.coda)
autocorr.plot(test.coda) 
}
\dontshow{
plot(test.coda,ask=FALSE)
autocorr.plot(test.coda,ask=FALSE)
}
}

\keyword{spatial}













