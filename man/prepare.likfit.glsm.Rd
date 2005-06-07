\name{prepare.likfit.glsm}


\alias{prepare.likfit.glsm}


\title{Prepare for Monte Carlo MLE}

\description{This function takes an output object from \code{\link{glsm.mcmc}}, and the corresponding data object of class 
  \code{geodata}, and prepares the 
  variables needed for the Monte Carlo maximum likelihood function
  \code{\link{likfit.glsm}}.  
}

\usage{
prepare.likfit.glsm(mcmc.output, use.intensity = FALSE)
}


\arguments{
  \item{mcmc.output}{an output file from the function \code{\link{glsm.mcmc}}.  }
  \item{use.intensity}{logical. If \code{use.intensity = TRUE} then the integration variable in the Monte Carlo approximation 
    will not be \eqn{S} but the intensity \eqn{g_{\lambda}^{-1}(S)}{g_{lambda}^{-1}(S)}. The latter makes it possible to
    use other link functions in  \code{\link{likfit.glsm}} than the one used in \code{mcmc.output}. 
    Default is \code{use.intensity = FALSE}.  }
}

\value{
  An object containing the sample and the approximating density to be used in \code{\link{likfit.glsm}}.
}

\author{
  Ole F. Christensen \email{olefc@birc.au.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{likfit.glsm}} for how to use the output from this function, and \code{\link{glsm.mcmc}} for 
  generating the object \code{mcmc.output}.
}

\examples{
\dontrun{
data(p50)
mcmc.4 <- mcmc.control(S.scale = 0.000035, n.iter = 1000)
kr4 <- list(family="poisson", cov.pars=c(1,1), beta=0)
condsim.4 <- glsm.mcmc(p50, mcmc.input = mcmc.4, model = kr4)     
mcmcobj.4 <- prepare.likfit.glsm(condsim.4)   
lik.4 <- likfit.glsm(mcmcobj.4, ini.phi = 10, fix.nugget.rel = TRUE) 
lik.4.spherical.nugget <- likfit.glsm(mcmcobj.4, ini.phi = 5.59, 
                              cov.model = "spherical", nugget.rel = 0.385)
}
}


\keyword{spatial}
