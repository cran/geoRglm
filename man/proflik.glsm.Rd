\name{proflik.glsm}


\alias{proflik.glsm}


\title{Computes Profile Likelihood for generalised linear spatial models}

\description{Computes two dimenaional profile likelihood for the
  parameters
  (phi, nugget.rel) for a model previously
  derived using the function \code{likfit.glsm}.
}

\usage{
proflik.glsm(mcmc.obj, obj.likfit.glsm, phi.values, nugget.rel.values,
             messages, \dots)
}



\arguments{
  \item{mcmc.obj}{object with the Monte Carlo simulations and corresponding approximating density. 
    This object should be an output from the function
    \code{\link{prepare.likfit.glsm}}.  }
  \item{obj.likfit.glsm}{ Output file from
    \code{\link{likfit.glsm}}.
  }
  \item{phi.values}{set of values of the parameter \code{phi} for which the
    profile likelihood will be computed.  }
  \item{nugget.rel.values}{set of values of the relative nugget parameter
    for which the profile likelihood will be computed.
    Only used if \code{obj.likfit.glsm} was created with the option \code{fix.nugget = FALSE}.
  }
  \item{messages}{logical. Indicates whether messages should be printed on the screen (or output device)
    while the function is running. Note that for this function
    additional messages can be obtained by setting the global option
    \code{verbose=TRUE}  }
  \item{\dots}{additional parameters to be passed to the maximization
    function. Typically  arguments of the type \code{control()} which controls the
    behavior of the minimization algorithm. For further details, see the documentation 
    for the minimization function \code{\link{optim}}.  }
}


\value{
  An object of the class `"proflik"' which is  a list. The element
  contains values of the pair of parameters and the corresponding value of the profile likelihood.
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.au.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\lik{likfit.geoRglm}} for the parameter estimation, and
  \code{\link[geoR]{proflik}} for the profile likelihood in the Gaussian spatial model.
}

\examples{
data(p50)
\dontshow{library(geoR)}
\dontrun{
mcmc.5 <- mcmc.control(S.scale = 0.6, thin=20, n.iter=50000, burn.in=1000)
model.5 <- list(cov.pars=c(0.6, 0.01), beta=1, family="poisson")
outmcmc.5 <- glsm.mcmc(p50, model= model.5, mcmc.input = mcmc.5)     
mcmcobj.5 <- prepare.likfit.glsm(outmcmc.5)   
lik.5.sph.nugget <- likfit.glsm(mcmcobj.5, ini.phi = 1, 
             cov.model = "spherical", nugget.rel = 0.385)
pr.lik.5.sph.nugget <- proflik.glsm(mcmcobj.5, lik.5.sph.nugget,
         phi.values = seq(0.5,5,l=10), nugget.rel.values=seq(0.5,5,l=10))
plot(pr.lik.5.sph.nugget)
}
\dontshow{
y9 <- list(data=p50$data[1:9],coords=p50$coords[1:9,])
class(y9) <- "geodata"
mcmc.5 <- mcmc.control(S.scale = 0.6, n.iter = 100)
model.5 <- list(cov.pars=c(0.1, 1), beta=1.1, family="poisson")
outmcmc.5 <- glsm.mcmc(y9, model= model.5, mcmc.input = mcmc.5)     
mcmcobj.5 <- prepare.likfit.glsm(outmcmc.5)   
lik.5.sph.nugget <- likfit.glsm(mcmcobj.5, ini.phi = 1, 
                           cov.model = "spherical", nugget.rel = 0.385)
pr.lik.5.sph.nugget <- proflik.glsm(mcmcobj.5, lik.5.sph.nugget,
         phi.values = seq(0.5,5,l=3), nugget.rel.values=seq(0.5,5,l=3))
plot(pr.lik.5.sph.nugget)
}
}

\keyword{spatial}
