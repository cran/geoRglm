\name{rongelap}

\alias{rongelap}

\non_function{rongelap}

\title{Radionuclide Concentrations on Rongelap Island}

\usage{data(rongelap)}

\description{ 
  This data-set was used by Diggle, Tawn and Moyeed (1998) to illustrate the methods
  discussed in the paper. The radionuclide concentration data set consists 
  of measurements of \eqn{\gamma}{gamma}-ray counts at \eqn{157} locations. 
}
  
\format{
  The object is a list
  with the following components:
  \describe{
    \item{\code{coords}}{the coordinates of data locations.  }
    \item{\code{data}}{the data.  }
    \item{\code{units.m}}{\eqn{n}-dimensional vector of observation-times for the data.  }
    \item{\code{borders}}{a matrix with the coordinates defining the coastline on
    Rongelap Island.  }
  }
  
}

\source{
For further details on the radionuclide concentration data, see Diggle, Harper and Simon (1997) and Diggle {\em et al.} (1998).
}

\references{
  Diggle, P. J., Harper, L. and Simon, S. L. (1997). Geostatistical analysis of residual contamination from nuclea
  testing. In: \emph{Statistics for the environment 3: pollution assesment and control} (eds. V. Barnet and K. F. Turkmann), Wiley, 
  Chichester, 89-107.
  
  Diggle, P. J., Tawn, J. A. and Moyeed, R. A. (1998). Model-based
  geostatistics (with Discussion). \emph{Applied Statistics},
  \bold{47}, 299--350.
}

\keyword{datasets}
\keyword{spatial}
