\name{model.glm.control}

\alias{model.glm.control}

\title{ Defines model options }

\description{
  This auxiliary function defines model options for
  \code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.
}

\usage{
model.glm.control(trend.d = "cte", trend.l = "cte", cov.model = "matern",
                  kappa = 0.5, aniso.pars = NULL, lambda = 0)
}


\arguments{
   \item{trend.d}{specifies the covariate values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details. Default is  \code{trend.d="cte"}.   }
  \item{trend.l}{specifies the covariate values at prediction
    locations. It must be of the same type as for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{cov.model}{string indicating the name of the model for the
    correlation function. Further details  in the
    documentation for \code{\link[geoR]{cov.spatial}}.  }
  \item{kappa}{additional smoothness parameter required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"gneiting.matern"} and \code{"cauchy"}.   }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = NULL} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link[geoR]{coords.aniso}}.  }
  \item{lambda}{parameter controlling the link function for
    \code{pois.krige.bayes}. \code{lambda = 0} is the default, which corresponds to a log link. 
    Other values gives a link function from the Box-Cox class.  }  
}


\value{
  A list with processed arguments to be passed to the main function.
}

\references{ 
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.au.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}} and \code{\link{binom.krige.bayes}}.  }

\keyword{spatial}

