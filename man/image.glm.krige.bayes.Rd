\name{image.glm.krige.bayes}

\alias{image.glm.krige.bayes}
\alias{persp.glm.krige.bayes}

\title{Plots Results of the Predictive Distribution}

\description{
  This function produces an image or perspective plot of a selected
  element
  of the predictive distribution
  returned by the functions \code{\link{binom.krige.bayes}} and \code{\link{pois.krige.bayes}}.
}

\usage{
\method{image}{glm.krige.bayes}(x, locations, borders,
                  values.to.plot=c("median", "uncertainty",
                     "quantiles", "probabilities", "simulation"),
                  number.col, coords.data, xlim, ylim, x.leg, y.leg,
                  \dots)



\method{persp}{glm.krige.bayes}(x, locations, borders,
                  values.to.plot=c("median", "uncertainty",
                     "quantiles", "probabilities", "simulation"),
                  number.col, \dots)
}

\arguments{
  \item{x}{an object of the class \code{glm.krige.bayes}, typically an
    output of the functions \code{\link{binom.krige.bayes}} or \code{\link{pois.krige.bayes}}.  }
  \item{locations}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates of the prediction locations, which should define a
    regular grid
    in order to be plotted by \code{\link{image}} or
    \code{\link{persp}}. By default does not need to be provided and
    evaluates the attribute "prediction.locations" from the input object.  }
  \item{borders}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates defining the borders of a region inside the grid defined
    by \code{locations}. Elements in the
    argument \code{values} are assigned to locations internal to the
    borders  and \code{NA}'s to the external ones.  }
  \item{values.to.plot}{select the element of the predictive
    distribution to be plotted. See DETAILS below.}
  \item{number.col}{Specifies the number of the column to be plotted.
    Only used if previous argument is set to one of \code{"quantiles"}, 
    \code{"probabilities"}  or \code{"simulation"}.  }
  \item{coords.data}{optional. If an \eqn{n \times 2}{n x 2}
    matrix with the data coordinates is
    provided, points indicating the data locations are included in the
    plot.  }
  \item{xlim, ylim} {the ranges to be encompassed by the x and y axes.  }
  \item{x.leg, y.leg}{limits for the legend in the horizontal and vertical directions.  }
  \item{\dots}{extra arguments to be passed to the plotting function
    \code{\link{image}} or \code{\link{persp}}.  }
  
}

\details{
  The functions \code{\link{binom.krige.bayes}} and \code{\link{pois.krige.bayes}} return
  summaries and other results about the predictive distributions.
  The argument \code{values.to.plot} specifies which result will be
  plotted. It can be passed to the function in two different forms:
  \itemize{
    \item a vector with the object containing the values to be plotted, or
    \item one of the following options:
    \code{"median"}, 
    \code{"uncertainty"},
    \code{"quantiles"},
    \code{"probability"} or
    \code{"simulation"}.  }
  For the last three options, if the results are stored in matrices,
  a column number must be provided using the argument \code{number.col}.

  The documentation for the functions \code{\link{binom.krige.bayes}} and
  \code{\link{pois.krige.bayes}} provide
  further details about these options.

}

\value{
  An \code{\link{image}} or \code{\link{persp}} plot is produced on the
  current graphics device. No values are returned. 
}


\references{
 Further information about \bold{geoRglm} can be found at:\cr
 \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}


\seealso{\code{\link[geoR]{image.krige.bayes}} for plotting output from \code{\link[geoR]{krige.bayes}}}

\examples{
#See examples in the documentation for the functions binom.krige.bayes and pois.krige.bayes .
}

\keyword{spatial}

