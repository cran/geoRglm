\name{hist.glm.krige.bayes}

\alias{hist.glm.krige.bayes}


\title{Plots Sample from Posterior Distributions}

\description{
  Plots histograms and/or density estimation with samples
  from the posterior distribution of the model parameters for output
  from the functions \code{\link{binom.krige.bayes}} and \code{\link{pois.krige.bayes}}
}

\usage{
\method{hist}{glm.krige.bayes}(x, pars, density.est = TRUE, histogram = TRUE, \dots)
}



\arguments{
  \item{x}{an object of the class \code{glm.krige.bayes}, with an output from
    the funtions \code{\link{binom.krige.bayes}} or \code{\link{pois.krige.bayes}}. }
  \item{pars}{a vector with the names of one or more of the model parameters.
    Defaults to all model parameters.  }
  \item{density.est}{logical indication whether a line with the density 
    estimation should be added to the plot.  }
  \item{histogram}{logical indicating whether the histogram is included in the 
    plot.  }
  \item{\dots}{further arguments for the plotting functions and or for the 
    density estimation.  }
}

\value{
  Produces a plot in the currently graphics device.\cr
  Returns an \code{\link{invisible}} list with the components:
  
  \item{histogram}{with the output of the function \code{\link{hist}} for 
    each parameter}
  \item{density.estimation}{with the output of the function 
    \code{\link{density}} for each parameter}
  
}


\author{
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{
  \code{\link{binom.krige.bayes}}, \code{\link{pois.krige.bayes}}, \code{\link[geoR]{hist.krige.bayes}}.
}

\examples{
## See documentation for binom.krige.bayes and pois.krige.bayes
}

\keyword{spatial}
\keyword{dplot}
