\name{output.glm.control}

\alias{output.glm.control}

\title{ Defines output options }

\description{
  This auxiliary function defines output options for
  \code{\link{pois.krige.bayes}}, \code{\link{binom.krige.bayes}} and \code{\link{glsm.krige}} (not all arguments have an effect for the later).
}

\usage{
output.glm.control(sim.posterior, sim.predict, keep.mcmc.sim, quantile,
                   threshold, inference, messages)
}


\arguments{
  \item{sim.posterior}{logical. Indicates
    whether or not the MCMC-sample from the posterior distribution of
    the parameters should be returned for \code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}. Default is \code{sim.posterior = TRUE}.  }
  \item{sim.predict}{logical. Defines whether simulations are
    drawn from the predictive distribution. Only valid if prediction
    locations are provided in the argument \code{locations}. Default is \code{sim.predict = FALSE}.   }
  \item{keep.mcmc.sim}{logical. Indicates
    whether or not the MCMC-sample from the posterior distribution of
    the signal should be returned for \code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}. Here the signal
     is the inverse link function of the random effect
     \eqn{g^{-1}(S)}{g^{-1}(S)}, which equals
    \eqn{\exp(S)}{exp(S)} for the Poison model with \eqn{\log}{log}-link,
    \eqn{(\lambda S +1)^{\lambda}}{(lambda S + 1)^{lambda}} for the Poison model with Box-Cox-link
    and \eqn{\exp(S)/(1+\exp(S))}{exp(S)/(1+exp(S))} for the binomial model with logit-link. Default is
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
    The default is \code{quantile = TRUE} where the values \code{c(0.025, 0.5, 0.975)} are used.
    If \code{quantile = FALSE} no quantiles are computed (and hence
    neither median nor uncertainty is included in the output).
    Only used if prediction locations are provided in the argument
    \code{locations}. Not implemented for \code{\link{glsm.krige}}. }
  \item{threshold}{one or more values with threshold values can be given to this argument.
    If so, an object called
    \code{probability} is included in the output. This object contains, for each
    prediction location, the probability that the variable is less
    than or equal to the threshold value given in the argument. Default is that 
    nothing is computed. Not implemented for \code{\link{glsm.krige}}. } 
  \item{inference}{logical. Indicates whether or not inference
    (summary of the parameters, and prediction) is performed or
    not. Default is \code{inference=TRUE}.
    Setting \code{inference=FALSE} is useful in an initial stage, when 
    tuning the MCMC algorithm (choosing \code{S.scale} and \code{phi.scale} appropriate). Not implemented for \code{\link{glsm.krige}}.
  }
  \item{messages}{logical. Indicates whether or not status messages are printed on the screen (or other output device)
    while the function is running.  }
}


\value{
  A list with processed arguments to be passed to the main function.
}

\author{
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}}, \code{\link{binom.krige.bayes}} and \code{\link{binom.krige.bayes}}.  }

\keyword{spatial}

