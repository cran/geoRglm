\name{prior.glm.control}

\alias{prior.glm.control}
\alias{prior.glm.check.aux}

\title{ Defines prior options }

\description{
  This auxiliary function defines prior options for
  \code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.
}

\usage{
prior.glm.control(beta.prior = c("flat", "normal", "fixed"),
              beta = NULL, beta.var.std = NULL,
              sigmasq.prior = c("uniform", "sc.inv.chisq", "reciprocal", "fixed"),
              sigmasq = NULL, df.sigmasq = NULL,
              phi.prior = c("uniform", "exponential","fixed",
                            "squared.reciprocal", "reciprocal"),
              phi = NULL, phi.discrete = NULL,
              tausq.rel = 0)
}


\arguments{
  \item{beta.prior}{prior distribution for the mean (vector) parameter \eqn{\beta}{beta}. The options are \code{"flat"} (default),  
    \code{"normal"} or \code{"fixed"}.  }
  \item{beta}{hyper-parameter for the prior distribution of the mean (vector) parameter \eqn{\beta}{beta}. 
    Only used if \code{beta.prior = "normal"} or \code{beta.prior = "fixed"}. For the latter \code{beta} defines the value of
    the known mean.  }
  \item{beta.var.std}{(co)variance hyper-parameters for the prior of the mean
    (vector) parameter \eqn{\beta}{beta}. Only used if \code{beta.prior = "normal"}.  }
  \item{sigmasq.prior}{prior distribution for the parameter \eqn{\sigma^2}{sigma^2}. The options are \code{"uniform"} (default), 
    \code{"sc.inv.chisq"}, \code{"reciprocal"} (gives improper posterior), or \code{"fixed"}.  }
  \item{sigmasq}{fixed value of the parameter \eqn{\sigma^2}{sigma^2} when
    \code{sigmasq.prior = "fixed"}. Parameter \eqn{S^2_{\sigma}}{} in the scaled inverse-\eqn{\chi^2}{chi^2} prior distribution for \eqn{\sigma^2}{sigma^2}.  }
  \item{df.sigmasq}{parameter \eqn{n_{\sigma}}{} in the scaled
    inverse-\eqn{\chi^2}{chi^2} prior distribution for
    \eqn{\sigma^2}{sigma^2}.  }
  \item{phi.prior}{prior distribution for the range parameter \eqn{\phi}{phi}. 
    Options are: \code{"uniform"} (\eqn{\propto 1}),
    \code{"exponential"} (\eqn{\exp(-\nu * \phi)}{exp(- nu * phi)}),
    \code{"fixed"} (known value of  \eqn{\phi}{phi}),
    \code{"squared.reciprocal"} (\eqn{1/\phi^2}{1/phi^2}),
    \code{"reciprocal"} (\eqn{1/\phi}{1/phi}).
    Alternativelly, a user defined discrete distribution can be
    specified by providing a vector of probabilities. These
    probabilities corrresponds to a prior distribution with support \code{phi.discrete}. \cr     
    If the \code{"fixed"} the argument \eqn{\phi}{phi}
    should be provided and it is regarded as fixed when performing predictions.

  }
  \item{phi}{fixed value of the parameter \eqn{\phi}{phi} when
    \code{phi.prior = "fixed"}. Mean of the prior distribution when \code{phi.prior = "exponential"}.   }
  \item{phi.discrete}{support points for the discretisation of the
    prior for the parameter \eqn{\phi}{phi}.  }
  \item{tausq.rel}{the value of the relative nugget parameter
    \eqn{\tau_R^2}{tau_R^2}. Default is \code{tausq.rel = 0}.  }
}


\value{
  A list with processed arguments to be passed to the main function.
}


\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.  }

\keyword{spatial}

