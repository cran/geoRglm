\name{prepare.lik.sim}

\alias{prepare.lik.sim}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Prepare for Monte Carlo MLE}

\description{This function takes an output object from \code{\link{pois.krige.bayes}}, and  the corresponding data object of class 
  \code{geodata}, and prepares the needed 
  variables for the Monte Carlo maximum likelihood function
  \code{\link{pois.likfit.sim}}.  
}

\usage{
prepare.lik.sim(bayes.output, geodata, use.intensity = FALSE)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{bayes.output}{an output file from the function \code{\link{pois.krige.bayes}}.  }
  \item{geodata}{the \bold{geoR} data set used for creating \code{bayes.output}.  }
  \item{use.intensity}{logical. If \code{use.intensity = TRUE} then the integration variable in the Monte Carlo approximation 
    will not be \eqn{S} but the intensity \eqn{g_{\lambda}^{-1}(S)}{g_{lambda}^{-1}(S)}. The latter makes it possible to
    use other link functions in  \code{\link{pois.likfit.sim}} than the one used in \code{bayes.output}. 
    Default is \code{use.intensity = FALSE}.  }
}

\value{
  An object containing the sample and the approximating density to be used in \code{\link{pois.likfit.sim}}.
}

\references{  
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.likfit.sim}} for how to use the output from this function, and \code{\link{pois.krige.bayes}} for 
  generating the object \code{bayes.output}.
}

\examples{
\dontrun{
if(is.R()) data(y50)
mcmc.4 <- mcmc.control(S.scale = 0.000035, n.iter = 1000)
prior.4 <- prior.glm.control(phi.prior = "fixed", phi = 150)
bayes.4 <- pois.krige.bayes(y50, prior = prior.4, mcmc.input = mcmc.4)     
mcmcobj.4 <- prepare.lik.sim(bayes.4, y50)   
lik.4 <- pois.likfit.sim(mcmcobj.4, ini.phi = 10, fix.nugget.rel = TRUE) 
lik.4.spherical.nugget <- pois.likfit.sim(mcmcobj.4, ini.phi = 5.59, 
                              cov.model = "spherical", nugget.rel = 0.385)
}
}


\keyword{spatial}
