\name{output.glm.control}

\alias{output.glm.control}

\title{ Defines output options }

\description{
  This auxiliary function defines output options for
  \code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.
}

\usage{
output.glm.control(sim.posterior, sim.predict, keep.mcmc.sim, quantile,
                   threshold, inference, messages.screen)
}


\arguments{
  \item{sim.posterior}{logical. Indicates
    whether or not the MCMC-sample from the posterior distribution of
    the parameters should be returned.  Default is \code{sim.posterior = TRUE}.  }
  \item{sim.predict}{logical. Defines whether simulations are
    drawn from the predictive distribution. Only valid if prediction
    locations are provided in the argument \code{locations}. Default is \code{sim.predict = FALSE}.   }
  \item{keep.mcmc.sim}{logical. Indicates
    whether or not the MCMC-sample from the posterior distribution of
    the signal should be returned. Here the signal
     is the inverse link function of the random effect
     \eqn{g^{-1}(S)}{g^{-1}(S)}, which equals
    \eqn{\exp(S)}{exp(S)} for the Poison model with \eqn{\log}{log}-link,
    \eqn{(\lambda S +1)^{\lambda}}{(lambda S + 1)^{lambda}} for the Poison model with Box-Cox-link
    or \eqn{\exp(S)/(1+\exp(S))}{exp(S)/(1+exp(S))} for the binomial model with logit-link. Default is
    \code{keep.mcmc.sim = TRUE}.  }
  \item{quantile}{indicates whether quantiles of the
    simulations from the
    predictive distributions are computed and returned.
    If a vector with numbers in the interval
    \eqn{[0,1]}{[0,1]} is provided, the output includes the object
    \code{quantile}, which contains values of corresponding
    estimated quantiles. For
    example, if \code{quantile = c(0.25,0.50,0.75)} the function
    returns the quartiles of the distributions at each of the
    prediction locations.
    If \code{quantile = TRUE} the default values \code{c(0.025, 0.5,
      0.975)}, are assumed. 
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{threshold}{one or more values with threshold values can be given to this argument.
    If so, an object called
    \code{probability} is included in the output. This object contains, for each
    prediction location, the probability that the variable is less
    than or equal to the threshold value given in the argument. Default is that 
    nothing is computed.  } 
  \item{inference}{logical. Indicates whether or not inference
    (summary of the parameters, and prediction) is performed or
    not. Default is \code{inference=TRUE}.
    Setting \code{inference=FALSE} is useful in an initial stage, when 
    tuning the MCMC algorithm (choosing \code{S.scale} and \code{phi.scale} appropriate). 
  }
  \item{messages.screen}{logical. Indicates whether or not status messages are printed on the screen (or other output device)
    while the function is running.  }
}


\value{
  A list with processed arguments to be passed to the main function.
}


\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.  }

\keyword{spatial}

