\name{pois.krige.bayes}

\alias{pois.krige.bayes}
\alias{mcmc.bayes.pois.log}
\alias{mcmc.bayes.pois.boxcox}
\alias{mcmc.bayes.conj.pois.log}
\alias{mcmc.bayes.conj.pois.boxcox}
\alias{krige.bayes.extnd}
\alias{cut0.calc.mixed.gauss}
\alias{cut0.calc.mixed.t}
\alias{phi.number}
\alias{multgauss}
\alias{BCinv.aux}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Bayesian Posterior Simulation and Prediction for the Poisson
  Spatial model}
\description{This function performs posterior simulation (by MCMC) and
  spatial prediction in the Poisson spatial model (with link function
  from the Box-Cox class).
}

\usage{
pois.krige.bayes(geodata, coords = geodata$coords, data = geodata$data, 
         units.m = "default", locations = "no", 
         model, prior, mcmc.input, output)}


\arguments{
  \item{geodata}{a list containing elements \code{coords} and \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \bold{geoR} data set.
    If not provided the arguments \code{coords} and \code{data} must be given instead. 
    The list may also contain an argument \code{units.m} as described below. }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{units.m}{\eqn{n}-dimensional vector of observation times for the data. By default (\code{units.m = "default"}),
  it takes \code{geodata$units.m} in case this exist and else a vector of 1's.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data frame with the 2-D
    coordinates of the \eqn{N} prediction locations.  }
  \item{model}{a list defining the components of the model. It can take an output from \code{model.glm.control} or
    a list with elements as for the arguments in \code{model.glm.control}.
    See documentation for \code{model.glm.control}.
    All arguments have default values }
  \item{prior}{specification of priors for the model parameters.
    It can take an output from \code{prior.glm.control} or
    a list with elements as for the arguments in \code{prior.glm.control}. See documentation for
   \code{prior.glm.control}. ATTENTION: When \code{phi.prior = "fixed"} then \code{phi} must be specified, and when 
    \code{phi.prior} is not \code{"fixed"} then \code{phi.discrete} must be specified. All other parameters have default values.  }
  \item{mcmc.input}{input parameter for the MCMC algorithm. It can take an output from \code{mcmc.control} or
    a list with elements as for the arguments in
    \code{mcmc.control}. See documentation for
    \code{mcmc.control}.\cr ATTENTION: the argument \code{S.scale} must be specified,
    the argument \code{phi.start} must specified when 
    \code{prior$phi} is not \code{"fixed"}, while all the others have default values.  }
  \item{output}{parameters for controlling the output. It can take an output from \code{output.glm.control} or
    a list with elements as for the arguments in \code{output.glm.control}.
    See documentation for \code{output.glm.control}.  }
}

\details{

 \code{pois.krige.bayes} is a function for Bayesian geostatistical 
  analysis in the Poisson spatial model.
  
  It can be used for an analysis with fixed values of the parameters. However, the
  function \code{pois.log.krige} may be preferable in this case.
  
  The Bayesian algorithm is using a discretized version of the prior
  distribution for the parameter \eqn{\phi}{phi}. This means that the prior for \eqn{\phi}{phi} is always a proper prior.
  
For simulating from the posterior distribution of \eqn{S} given \eqn{y}, 
we use a Langevin-Hastings type algorithm. 
This algorithm is a Metropolis-Hastings algorithm, where the
proposal distribution uses gradient information from the posterior. The algorithm is described below.
For shortness of presentation, we present only the MCMC algorithm for
the case where \eqn{\beta}{beta} follows a uniform prior and the link
function is the canonical log-link.

When \eqn{\beta}{beta} follows a uniform prior and the prior for \eqn{\sigma^2}{sigma^2} is a scaled inverse-\eqn{\chi^2}{chi^2}
distribution, 
the marginalised improper density of \eqn{S} is
\deqn{
f_I(s) \propto |D^T V^{-1}D|^{1/2}|V|^{-1/2}\{n_{\sigma}S^2_{\sigma} 
+ s^T (V^{-1}-V^{-1}D(D^T V^{-1}D)^{-1}D^T V^{-1})s \}^{-(n-p+n_{\sigma})/2},
}{
f_I(s) propto |D^T*V^{-1}*D|^{1/2}*|V|^{-1/2}*{n.sigma*S2.sigma
+ s^T*(V^{-1}-V^{-1}*D*(D^T*V^{-1}*D)^{-1}*D^T*V^{-1})*s }^{-(n-p+n.sigma)/2},
}
where \eqn{V} is the correlation matrix of \eqn{S}. The uniform prior for \eqn{\sigma^2}{sigma^2} corresponds 
to \eqn{S^2_{\sigma}=0}{S2.sigma=0} and 
\eqn{n_{\sigma}=-2}{n.sigma=-2}, and the reciprocal prior for \eqn{\sigma^2}{sigma^2} corresponds to \eqn{S^2_{\sigma}=0}{S2.sigma=0} 
and \eqn{n_{\sigma}=0}{n.sigma=0}.

We use the reparametrisation \eqn{S = Q\Gamma}{S = Q*Gamma}, where \eqn{Q} is the Cholesky factorisation of \eqn{V} so that \eqn{V=QQ^T}.
Posterior simulations of \eqn{S} are obtained by transforming
MCMC simulations from the conditional distribution of \eqn{\Gamma}{Gamma} given
\eqn{Y=y}.

The log posterior density of
 \eqn{\Gamma}{Gamma} given \eqn{Y=y} is
\deqn{
\log f(\gamma|y) = const(y) - \frac{1}{2} \gamma^T (I_n - V^{-1/2}D(D^T V^{-1}D)^{-1}D^T V^{-1/2})\gamma 
 +\sum_{i=1}^n y_i s_i - \exp(s_i),
}{
log f(gamma|y) = const(y) - 1/2 * gamma^T*(I_n - V^{-1/2}*D*(D^T*V^{-1}*D)^{-1}*D^T*V^{-1/2})*gamma 
 +\sum_i y_i*s_i -exp(s_i),
}
where
\eqn{(s_1,\ldots,s_n)^T= Q \gamma}{(s_1,...,s_n)^T = Q*gamma}.

For the truncated Langevin-Hastings update we use a truncated form of the gradient (truncating by \eqn{H_i}) of the log target density,
\deqn{
\nabla(\gamma )^{trunc}= - (I_n - Q^{-1}D(D^T V^{-1}D)^{-1}D^T(Q^{-1})^T)\gamma + Q^T\left\{y_i-\exp(s_i)\wedge H_i \right\}_{i=1}^n .
}{
nabla(gamma )^{trunc}= - (I_n - Q^{-1}*D*(D^T*V^{-1}*D)^{-1}*D^T*(Q^{-1})^T)*gamma + Q^T*{y_i-min{\exp(s_i),H_i} }_{i=1}^n .
}

The proposal \eqn{\gamma'}{gamma'} follows a multivariate normal distribution with mean vector
\eqn{\xi(\gamma)=\gamma+(h/2)\nabla(\gamma)^{trunc}}{xi(gamma)=gamma+(h/2)*nabla(gamma)^{trunc}} and covariance matrix \eqn{hI},
where \eqn{h>0} is a user-specified ``proposal variance'' (called
\code{S.scale}; see \code{mcmc.control}).

When \code{phi.prior} is not \code{"fixed"}, we update the parameter \eqn{\phi}{phi} by a random walk Metropolis step.
Here \code{mcmc.input$phi.scale} (see \code{mcmc.control}) is the proposal variance, which needs to be sufficiently large so that 
that algorithm easily can move between the discrete values in \code{prior$phi.discrete} (see \code{prior.glm.control}).

\bold{CONTROL FUNCTIONS}

  The function call includes auxiliary control functions which allows
  the user to specify and/or change the specification of 1) model
  components
  (using \code{model.glm.control}), 2) prior
  distributions (using \code{prior.glm.control}), 3) options for the
  MCMC algorithm (using \code{mcmc.control}), and 4) options for the
  output (using \code{output.glm.control}).
  Default values are available in most of the cases.
  The arguments for the control functions are described in their
  respective help files.
  
For the prediction part of the algorithm, we use the median of the 
predictive distribution as the predictor and 1/4 of the length of the 95 percent predictive interval as a measure of the prediction 
uncertainty. Below we describe the procedure for doing calculating these quantities.
 
First we perform a Bayesian Gaussian prediction with the given priors on \eqn{\beta}{beta}
and \eqn{\sigma^2}{sigma^2} on each of the simulated \eqn{S}{S}-``datasets'' from the 
posterior distribution (and in case \eqn{\phi}{phi} is not fixed, for each sampled \eqn{\phi}{phi} value).
This Gaussian prediction is done by calling the 
function \code{krige.bayes.extnd}, which 
is an extension of \code{\link[geoR]{krige.bayes}} 
allowing for more than one ``data set''. 

For calculating the probability above a threshold for the predictive distribution given the data, we first calculate this 
probability
for each of the simulated \eqn{S}{S}-``datasets''.
This is done using the fact that the predictive distribution 
for each of the simulated \eqn{S}{S}-``datasets'' is a multivariate \eqn{t}-distribution.
Afterwards the probability above a threshold is calculated by taking the empirical mean of these conditional probabilities.

Now the median and the 2.5 percent and 97.5 percent quantiles can be calculated by an iterative procedure, where first a guess
of the value is made, and second this guess is checked by calculating the probability of being less than this value. In case
the guess is too low, it is adjusted upwards, and vise verse.
}

\value{A list with the following components:
  
  \item{posterior}{A list with results for the posterior distribution of the
    model parameters and the random effects at the data locations. The components are:\cr  
  
  \itemize{
    \item{beta}{summary of posterior distribution
      for the parameter \eqn{\beta}{beta}.  }
    \item{sigmasq}{summary of the posterior distribution
      for the parameter \eqn{\sigma^2}{sigma^2}.  }
    \item{phi}{summary of the posterior distribution
      of the parameter \eqn{\phi}{phi}.  }
    \item{simulations}{sample from the posterior distribution of \eqn{\exp(S)/(1+\exp(S))}{exp(S)/(1+exp(S))} at the data locations. 
      Returned only if \code{keep.mcmc.sim = TRUE}.  }
    \cr
    \cr
  }
  }
  \item{predictive}{A list with results for the predictive distribution at the
    prediction locations (if provided). The
    components are:
    \cr
    \itemize{    
      \item{simulations}{a numerical matrix. Each column contains a simulation
        from the predictive distribution. Returned
        only if \code{sim.predict = TRUE}.  }
      \item{median}{a vector with the estimated median at the prediction
        locations.  }
      \item{uncertainty}{a vector with the estimated uncertainty at the
        prediction locations, defined as the length of the \eqn{95\%} prediction interval divided by 4. }
      \item{quantile}{A matrix or vector with quantile estimators.  }
      \item{probability}{A matrix or vector with probabilities above a threshold. Returned
        only if the argument \code{threshold} is used.  }
      \cr
      \cr
    }
  }
  \item{model}{model components used as defined by \code{model.glm.control}.  }
  \item{prior}{priors used for the model parameters.  }
  \item{mcmc.input}{input parameters used for the MCMC algorithm.  }
  \item{.Random.seed}{system random seed before running the function.
    Allows reproduction of results. If
    the \code{\link{.Random.seed}} is set to this value and the function is run
    again, it will produce exactly the same results.  }
  \item{call}{the function call.  }
  
}

\references{
Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.log.krige}} for prediction with fixed parameters in the
  Poisson-log normal model.  \code{\link[geoR]{krige.bayes}} for
  Bayesian prediction in the Gaussian spatial model. 
}

\examples{
if(is.R()) data(y50)
\testonly{library(geoR)}
## MCMC with fixed phi
prior.5 <- prior.glm.control(phi.prior = "fixed", phi = 0.1)
mcmc.5 <- mcmc.control(S.scale = 0.01, thin = 1)
test.5 <- pois.krige.bayes(y50, prior = prior.5, mcmc.input = mcmc.5)
\dontrun{
## Now chose S.scale (Acc-rate=0.60 is preferable).
mcmc.5.new <- mcmc.control(S.scale = 0.13, thin = 100)
test.5.new <- pois.krige.bayes(y50, locations = t(cbind(c(2.5,3), c(-6050,-3270))), 
               prior = prior.5, mcmc.input = mcmc.5.new, 
               output = list(threshold = 10, quantile = c(0.49999, 0.99)))
## MCMC with random phi.
## Note here that we can start with the S.scale from above.
mcmc.6 <- mcmc.control(S.scale = 0.13, n.iter = 2000, thin = 100, 
                       phi.scale = 0.01)
prior.6 <- prior.glm.control(phi.discrete = seq(0.02, 1, 0.02))
test.6 <- pois.krige.bayes(y50, prior = prior.6, mcmc.input = mcmc.6)
## Acc-rate=0.60 , acc-rate-phi = 0.25-0.30  are preferable
mcmc.6.new <- mcmc.control(S.scale=0.08, n.iter = 400000, thin = 200,
                       burn.in = 5000, phi.scale = 0.12, phi.start = 0.5)
prior.6 <- prior.glm.control(phi.prior = "uniform", 
               phi.discrete = seq(0.02, 1, 0.02))
test.6.new <- pois.krige.bayes(y50, 
                 locations = t(cbind(c(2.5,3.5), c(-60,-37))), 
                 prior = prior.6, mcmc.input = mcmc.6.new)
}
\testonly{
y9 <- grf(grid = expand.grid(x = seq(1, 3, l = 3), y = seq(1, 3, l = 3)), cov.pars = c(0.1, 0.2))
y9$data <- rpois(9, lambda = exp(y9$data+1))
mcmc5 <- mcmc.control(S.scale = 0.075, thin = 1, n.iter=10)
test5 <- pois.krige.bayes(y9, locations =
t(cbind(c(2.5,3),c(-6050,-3270))), prior = list(phi.prior = "fixed", phi= 0.1), mcmc.input = mcmc5,
output = list(threshold = 10, quantile = c(0.05,0.99))) 
prior6 <- prior.glm.control(phi.prior = "uniform", phi.discrete = c(0.2,0.25), tausq.rel = 0.05) 
mcmc6 <- list(S.scale = 0.075, thin = 1, n.iter=10, burn.in = 1, phi.scale = 1.2, phi.start = 0.5) 
test6 <- pois.krige.bayes(y9, locations =t(cbind(c(2.5,3.5),c(-60,-37))), prior = prior6, mcmc.input = mcmc6)
run2.y9 <- pois.krige.bayes(y9, prior = list(phi.prior = "fixed", phi= 0.1), mcmc.input = mcmc5,
model = list(trend.d ="1st", trend.l ="1st", aniso.pars =c(1,2), kappa
=1), locations = t(cbind(c(2.5,3.5),c(-60,-37))) ,
output = output.glm.control(sim.predict = TRUE)) }
}
\keyword{spatial}
