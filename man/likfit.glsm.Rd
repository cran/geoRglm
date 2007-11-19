\name{likfit.glsm}

\alias{likfit.glsm}
\alias{.func.val}
\alias{.NewtonRhapson.step}
\alias{.maxim.aux1}
\alias{.lik.sim}
\alias{.lik.sim.boxcox}


\title{Monte Carlo Maximum Likelihood Estimation in a Generalised Linear
  Spatial Model}

\description{This function performs Monte Carlo maximum likelihood in a
  generalised linear spatial model, based 
  on a Monte Carlo sample from the conditional distribution.
}

\usage{
likfit.glsm(mcmc.obj, trend = mcmc.obj$trend, cov.model = "matern", 
       kappa = 0.5, ini.phi, fix.nugget.rel = FALSE, nugget.rel = 0, 
       aniso.pars = NULL, fix.lambda = TRUE, lambda = NULL, 
       limits = pars.limits(), messages, \dots)
}

\arguments{
  \item{mcmc.obj}{object with the Monte Carlo simulations and corresponding approximating density. 
    This object should be an output from the function
    \code{\link{prepare.likfit.glsm}}.  }
  \item{trend}{specifies the covariate values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details. Default is that the trend is the same as in the
    \code{mcmc.obj} object.   }
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
    as the initial value for the maximization algorithm. Default is \code{nugget.rel = 0}.  }
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
  \item{limits}{values defining lower and upper limits for the model parameters used in the numerical minimization.
    The auxiliary function \code{\link[geoR]{pars.limits()}} is used to
    set the limits.  }
  \item{messages}{logical. Indicates whether status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{additional parameters to be passed to the optimisation
    function. Typically  arguments of the type \code{control()} which controls the
    behavior of the optimization algorithm. For further details, see the documentation 
    for the minimization function \code{\link{optim}}.  }
}

\details{
This function estimates the parameters in the Poisson/Binomial normal
model, using a Monte Carlo approximation to the likelihood. Further details can be found in Christensen (2004).

Parameter estimation is done numerically using the \R function
\code{\link{optim}} with box-constraints, i.e. \code{method="L-BFGS-B"}.

Lower and upper limits for parameter values can be specified using the function \code{pars.limits()}.
For example, including \code{limits = pars.limits(phi=c(0.01, 10))} in the function call
will specify the limits for the parameter \eqn{\phi}{phi}.
Default values are used if the argument \code{limits} is not provided.
     
Only when the \code{mcmc.obj} object contains an object \code{mu} giving the intensity, is it possible to use other link
functions than the link function used for the generated sample in
\code{mcmc.obj}

We strongly recommend that the user does not 
provide self-made input objects for \code{mcmc.obj}, but only uses objects created by \code{\link{prepare.likfit.glsm}}. In case the user
really wants to create his own objects, he should study the source code
very carefully to understand how it works.

Summary and print methods for summarising and printing the output also exist. 
}

\value{
  A list with the following components:
  \item{family}{the error distribution (Poisson or Binomial).}
   \item{link}{the name of the link function. }
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
  \item{trend}{the trend  }
  \item{parameters.summary}{a data-frame with all model parameters, their
    status (estimated or fixed) and values.  }
  \item{loglik}{the value of the maximized likelihood.  }
  \item{npars}{number of estimated parameters.  }
  \item{info.minimisation}{results returned by the minimisation function.}
  \item{call}{the function call.  }
}

\references{
  Christensen, O. F. (2004). Monte Carlo maximum likelihood in model-based geostatistics.
  Journal of computational and graphical statistics \bold{13} 702-718.
   
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.daimi.au.dk/~olefc/geoRglm}.
}

\author{
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{prepare.likfit.glsm}} on how to prepare the object \code{mcmc.obj}, \code{\link{glsm.mcmc}} for 
MCMC simulation in generalised linear spatial model, and \code{summary.likGLSM} for
summarising the output. See also \code{\link[geoR]{likfit}} for 
parameter estimation in the Gaussian spatial model.
}

\examples{
data(p50)
\dontshow{library(geoR)}
\dontrun{
mcmc.5 <- mcmc.control(S.scale = 0.6, thin=20, n.iter=50000, burn.in=1000)
model.5 <- list(cov.pars=c(0.6, 0.01), beta=1, family="poisson")
outmcmc.5 <- glsm.mcmc(p50, model= model.5, mcmc.input = mcmc.5)     
mcmcobj.5 <- prepare.likfit.glsm(outmcmc.5)   
lik.5 <- likfit.glsm(mcmcobj.5, ini.phi = 0.1, fix.nugget.rel = TRUE)
print(lik.5)
summary(lik.5)
lik.5.sph.nugget <- likfit.glsm(mcmcobj.5, ini.phi = 1, 
                           cov.model = "spherical", nugget.rel = 0.385)
print(lik.5.sph.nugget)
summary(lik.5.sph.nugget)
}
\dontshow{
y9 <- list(data=p50$data[1:9],coords=p50$coords[1:9,])
class(y9) <- "geodata"
mcmc.5 <- mcmc.control(S.scale = 0.6, n.iter = 100)
model.5 <- list(cov.pars=c(0.1, 1), beta=1.1, family="poisson")
outmcmc.5 <- glsm.mcmc(y9, model= model.5, mcmc.input = mcmc.5)     
mcmcobj.5 <- prepare.likfit.glsm(outmcmc.5)   
lik.5 <- likfit.glsm(mcmcobj.5, ini.phi = 0.4, fix.nugget.rel = TRUE)
print(lik.5)
summary(lik.5)
data(b50)
b9 <- list(data=b50$data[1:9],coords=b50$coords[1:9,])
trend5 <- 1:9-5
model5 <- list(cov.pars=c(0.1, 1), beta=c(1.2,0),family="binomial",trend=~trend5)
outmcmc5 <- glsm.mcmc(b9, model= model5, mcmc.input = mcmc.5)  
mcmcobj5 <- prepare.likfit.glsm(outmcmc5)   
lik5 <- likfit.glsm(mcmcobj5, trend=model5$trend, ini.phi = 1,
                    fix.nugget.rel = TRUE)
print(lik5)
summary(lik5)
}
}

\keyword{spatial}
