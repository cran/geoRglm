\name{pois.likfit.sim}

\alias{pois.likfit.sim}
\alias{maxim.aux1}
\alias{lik.sim}
\alias{lik.sim.boxcox}
\alias{intgr.lik.sim} 
\alias{intgr.lik.sim.boxcox}


%- Also NEED an `\alias' for EACH other topic documented here.

\title{Monte Carlo Maximum Likelihood Estimation for the Poisson normal model}

\description{This function performs Monte Carlo maximum likelihood for the spatial Poisson-log normal model, based 
on a Monte Carlo sample from the conditional distribution.
}

\usage{
pois.likfit.sim(mcmc.obj, trend = "cte", cov.model = "exponential", 
       kappa = 0.5, ini.phi, fix.nugget.rel = FALSE, nugget.rel = 0, 
       aniso.pars = NULL, fix.lambda = TRUE, lambda = NULL, 
       lik.type = "standard", limits = likfit.limits(), 
       messages.screen = TRUE, ...)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{mcmc.obj}{object with the Monte Carlo simulations and corresponding approximating density. 
    This object should be an output from the function
    \code{\link{prepare.lik.sim}}.  }
  \item{trend}{specifies the covariate values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details. Default is  \code{trend.d="cte"}.   }
  \item{cov.model}{a string specifying the model for the correlation
    function. For further details see documentation for \code{\link[geoR]{cov.spatial}}.  }
  \item{kappa}{additional smoothness parameter required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"gneiting.matern"} and \code{"cauchy"}.  }
  \item{ini.phi}{initial value for the covariance parameter \eqn{\phi}{phi}.  }
  \item{fix.nugget.rel}{logical, saying whether the parameter
    \eqn{\tau_R^2}{tau_R^2} (relative nugget) should be regarded as fixed
    (\code{fix.nugget.rel = TRUE}) or should be
    estimated (\code{fix.nugget.rel = FALSE}). Default is \code{fix.nugget.rel = FALSE}.  }
  \item{nugget.rel}{value of the relative nugget parameter.
    Regarded as a fixed value if \code{fix.nugget.rel = TRUE}, otherwise
    as the initial value for the minimization algorithm. Default is \code{nugget.rel = 0}.  }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = NULL} the correction will be the same as for the generated sample in \code{mcmc.obj}.
    Otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link[geoR]{coords.aniso}}.  }
  \item{fix.lambda}{logical, indicating whether the Box-Cox transformation parameter
    \eqn{\lambda}{lambda} should be regarded as fixed
    (\code{fix.lambda = TRUE}) or should be be estimated (\code{fix.lambda = FALSE}). Default is  \code{fix.lambda = TRUE}.  }
  \item{lambda}{value of parameter \eqn{\lambda}{lambda} in the Box-Cox class of link functions.
    Regarded as a fixed value if \code{fix.lambda = TRUE}, otherwise as the initial value for the
    minimization algorithm. Default is \code{lambda = NULL}, in which case the used link function will be the same as for the 
    generated sample in \code{mcmc.obj}.  }
  \item{lik.type}{type of likelihood used, \code{"standard"} or \code{"integrated"}. Using the latter option, the integrated 
    likelihood is maximised. Here the integrated likelihood incorporates uncertainty on \eqn{\beta}{beta} and 
    \eqn{\sigma^2}{sigma^2} corresponding to uniform priors. Default is \code{"standard"}.  }
  \item{limits}{values defining lower and upper limits for the model parameters used in the numerical minimization.
    The auxiliary function \code{\link[geoR]{likfit.limits()}} is used to set the limits.  }
  \item{messages.screen}{logical. Indicates whether status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{additional parameters to be passed to the minimization
    function. Typically  arguments of the type \code{control()} which controls the
    behavior of the minimization algorithm. For further details, see the documentation 
    for the minimization function \code{\link{optim}}.  }
}

\details{
This function estimates the parameters in the Poisson normal
model, using a Monte Carlo approximation to the likelihood. Further details can be found in Christensen (2001).

Lower and upper limits for parameter values can be specified using the function \code{likfit.limits()}.
For example, including \code{limits = likfit.limits(phi=c(0, 10))} in the function call
will change the limits for the parameter \eqn{\phi}{phi}.
Default values are used if the argument \code{limits} is not provided.
     
Only when the \code{mcmc.obj} object contains an object \code{mu} giving the intensity, is it possible to use other link
functions than the link function used for the generated sample in \code{mcmc.obj}
     
We strongly recommend that the user does not 
provide self-made input objects for \code{mcmc.obj}, but only uses objects created by \code{\link{prepare.lik.sim}}. In case the user
really wants to create his own objects, he should study the source code very carefully to understand how it works.
}

\value{
  A list with the following components:
  \item{cov.model}{a string with the name of the correlation
    function.  }
  \item{beta}{estimate of the parameter \eqn{\beta}{beta}.
    This can be a scalar or vector depending on the
    covariates (trend) specified in the model.  }
  \item{cov.pars}{a vector with the estimates of the parameters
    \eqn{\sigma^2}{sigma^2} 
    and \eqn{\phi}{phi}, respectively.  }
  \item{nugget.rel}{value of the relative nugget parameter \eqn{\tau_R^2}{tau_R^2}.
    This is an estimate if \code{fix.nugget.rel = FALSE}, and otherwise a given fixed value.  }
  \item{kappa}{value of the smoothness parameter. Valid only when
    the correlation function is one of: \code{"matern"},
    \code{"powered.exponential"}, \code{"cauchy"}
    or \code{"gneiting.matern"}.  }
  \item{lambda}{values of the parameter for the Box-Cox class of link functions. A fixed value if
    \code{fix.lambda = TRUE}, otherwise the estimated value.  }
  \item{aniso.pars}{values of the anisotropy parameters used.  }
  \item{loglik}{the value of the maximized likelihood.  }
  \item{call}{the function call.  }
}

\references{
  Christensen, O. F. (2001). Monte Carlo maximum likelihood in model-based geostatistics.
  Research report R-01-2017, Department of Mathematical Sciences, Aalborg University (submitted for journal publication).
   
Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{prepare.lik.sim}} on how to prepare the object \code{mcmc.obj}, \code{\link{pois.krige.bayes}} for 
Bayesian inference in the Poisson normal model, and \code{\link[geoR]{likfit}} for 
parameter estimation in the Gaussian spatial model.
}

\examples{
if(is.R()) data(y50)
\testonly{library(geoR)}
\dontrun{
mcmc.4 <- mcmc.control(S.scale = 0.000035, n.iter = 1000)
prior.4 <- prior.glm.control(phi.prior = "fixed", phi = 150)
bayes.4 <- pois.krige.bayes(y50, prior = prior.4, mcmc.input = mcmc.4)     
mcmcobj.4 <- prepare.lik.sim(bayes.4, y50)   
lik.4 <- pois.likfit.sim(mcmcobj.4, ini.phi = 10, fix.nugget.rel = TRUE) 
lik.4.spherical.nugget <- pois.likfit.sim(mcmcobj.4, ini.phi = 5.59, 
                              cov.model = "spherical", nugget.rel =
0.385)}
\testonly{
y9 <- list(data=y50$data[1:9],coords=y50$coords[1:9,])
class(y9) <- "geodata"
mcmc.4 <- mcmc.control(S.scale = 0.000035, n.iter = 100)
prior.4 <- prior.glm.control(phi.prior = "fixed", phi = 150)
bayes.4 <- pois.krige.bayes(y9, prior = prior.4, mcmc.input = mcmc.4)     
mcmcobj.4 <- prepare.lik.sim(bayes.4, y9)   
lik.4 <- pois.likfit.sim(mcmcobj.4, ini.phi = 10, fix.nugget.rel = TRUE) 
}
}

\keyword{spatial}
