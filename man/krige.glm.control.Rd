\name{krige.gkm.control}


\alias{krige.glm.control}
\alias{krige.glm.check.aux}



\title{Defines options and model for prediction}
\description{
  This auxiliary function defines options and model for
  \code{\link{pois.krige}} and \code{\link{binom.krige}}. 
}

\usage{
krige.glm.control(type.krige = "sk", trend.d = "cte", trend.l = "cte",
                  obj.model = NULL, beta, cov.model, cov.pars, kappa,
                  nugget, micro.scale, dist.epsilon = 1e-10, 
                  aniso.pars, lambda)
}

\arguments{
  \item{type.krige}{type of prediction to be performed (minimal mean
    square error prediction). Options are
    \code{"sk"} and \code{"ok"} corresponding to prediction with fixed
    parameters (\code{type.krige = "sk"}), which is the default, or prediction with a uniform
    prior on \eqn{\beta}{beta} (\code{type.krige = "ok"}). 
    Prediction using a model with covariates can be done by specifying the
    covariate model using the arguments \code{trend.d} and
    \code{trend.l}.  }
  \item{trend.d}{specifies the trend (covariate) values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details.
    Default is \code{trend.d = "cte"}.   }
  \item{trend.l}{specifies the trend (covariate) values at prediction
    locations. It must be of the same type as for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{obj.model}{a list with the model parameters.  }
  \item{beta}{numerical value of the mean (vector) parameter.
    Only used if \code{type.krige="sk"}.  }
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
  \item{micro.scale}{micro-scale variance. If specified, the
    nugget is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}.
    This has effect on prediction where the ``signal'' part of \eqn{S}
    (without the measurement error part of the nugget) is predicted. The
    default is \code{micro.scale = nugget}.  }
  \item{dist.epsilon}{a numeric value. Locations which are separated
    by a distance less than this value
    are considered co-located.  }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link[geoR]{coords.aniso}}.  }
  \item{lambda}{numeric value of the Box-Cox transformation parameter
    for \code{pois.krige}.
    The value \eqn{\lambda = 1}{lambda = 1} corresponds to
    no transformation and \eqn{\lambda = 0}{lambda = 0} corresponds to
    the log-transformation. 
    Prediction results are back-transformed and
    returned is the same scale as for the original data.  }
}



\value{
  A list with processed arguments to be passed to the main
  function. 
}


Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige}} and \code{\link{binom.krige}}.  }


\keyword{spatial}













