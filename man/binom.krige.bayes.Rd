\name{binom.krige.bayes}

\alias{binom.krige.bayes}
\alias{.mcmc.bayes.binom.logit}
\alias{.mcmc.bayes.conj.binom.logit}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Bayesian Posterior Simulation and Prediction for the Binomial
  Spatial model}
\description{This function performs posterior simulation (by MCMC) and spatial prediction in the binomial-logit spatial
  model.
}

\usage{
binom.krige.bayes(geodata, coords = geodata$coords, data = geodata$data, 
         units.m = "default", locations = "no", 
         model, prior, mcmc.input, output)
}


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
  \item{units.m}{\eqn{n}-dimensional vector giving the number of trials for the data. 
    By default (\code{units.m = "default"}),
    it takes \code{geodata$units.m} in case this exist and else a vector of 1's.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data frame, or a
    list with the 2-D
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

 \code{binom.krige.bayes} is a function for Bayesian geostatistical 
  inference in the binomial-logit spatial model.
  
  The Bayesian algorithm is using a discretized version of the prior
  distribution for the parameter \eqn{\phi}{phi}. This means that the prior for \eqn{\phi}{phi} is always a proper prior.
  
For simulating from the posterior distribution of \eqn{S} given \eqn{y}, 
we use a Langevin-Hastings type algorithm. 
This algorithm is a Metropolis-Hastings algorithm, where the
proposal distribution uses gradient information from the posterior. The algorithm is described below.
For shortness of presentation, we only present the MCMC algorithm for the case where \eqn{\beta}{beta} follows a uniform prior.

When \eqn{\beta}{beta} follows a uniform prior and the prior for \eqn{\sigma^2}{sigma^2} is a scaled inverse-\eqn{\chi^2}{chi^2} distribution, 
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
 +\sum_{i=1}^n y_i s_i - \log(1+\exp(s_i)),
}{
log f(gamma|y) = const(y) - 1/2 * gamma^T*(I_n - V^{-1/2}*D*(D^T*V^{-1}*D)^{-1}*D^T*V^{-1/2})*gamma 
 +\sum_i y_i*s_i -log(1+exp(s_i)),
}
where
\eqn{(s_1,\ldots,s_n)^T= Q \gamma}{(s_1,...,s_n)^T = Q*gamma}.

For the Langevin-Hastings update we use the gradient of the log target density,
\deqn{
\nabla(\gamma )^{trunc}= - (I_n - Q^{-1}D(D^T V^{-1}D)^{-1}D^T(Q^{-1})^T)\gamma + Q^T\left\{y_i-\exp(s_i)/(1+\exp(s_i))\right\}_{i=1}^n .
}{
nabla(gamma )^{trunc}= - (I_n - Q^{-1}*D*(D^T*V^{-1}*D)^{-1}*D^T*(Q^{-1})^T)*gamma + Q^T*{y_i-\exp(s_i)/(1+exp(s_i))}_{i=1}^n .
}

The proposal \eqn{\gamma'}{gamma'} follows a multivariate normal distribution with mean vector
\eqn{\xi(\gamma)=\gamma+(h/2)\nabla(\gamma)^{trunc}}{xi(gamma)=gamma+(h/2)*nabla(gamma)^{trunc}} and covariance matrix \eqn{hI},
where \eqn{h>0} is a user-specified ``proposal variance'' (called
\code{S.scale}; see \code{mcmc.control}).

When \code{phi.prior} is not \code{"fixed"}, we update the parameter \eqn{\phi}{phi} by a random walk Metropolis step.
Here \code{mcmc.input$phi.scale} (see \code{mcmc.control}) is the proposal variance, which needs to be sufficiently large so that 
the algorithm easily can move between the discrete values in \code{prior$phi.discrete} (see \code{prior.glm.control}).

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
  
In the prediction part of the function we want to predict
\eqn{\exp(S^*)/(1+\exp(S^*))}{exp(S^*)/(1+exp(S^*))} at locations of interest.
For the prediction part of the algorithm, we use the median of the 
predictive distribution as the predictor and 1/4 of the length of the 95 percent predictive interval as a measure of the prediction 
uncertainty. Below we describe the procedure for calculating these quantities.
 
First we perform a Bayesian Gaussian prediction with the given priors on \eqn{\beta}{beta}
and \eqn{\sigma^2}{sigma^2} on each of the simulated \eqn{S}{S}-``datasets'' from the 
posterior distribution (and in case \eqn{\phi}{phi} is not fixed, for each sampled \eqn{\phi}{phi} value).
This Gaussian prediction is done by calling the function \code{krige.bayes.extnd}, which 
is an extension of \code{\link[geoR]{krige.bayes}} 
allowing for more than one ``data set''.

For calculating the probability below a threshold for the predictive distribution given the data, we first calculate this 
probability
for each of the simulated \eqn{S}{S}-``datasets''.
This is done using the fact that the predictive distribution 
for each of the simulated \eqn{S}{S}-``datasets'' is a multivariate \eqn{t}-distribution.
Afterwards the probability below a threshold is calculated by taking the empirical mean of these conditional probabilities.

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
    \item{acc.rate}{The acceptance rates.  }
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
      \item{probability}{A matrix or vector with probabilities below a threshold. Returned
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
 \url{http://www.daimi.au.dk/~olefc/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.au.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{binom.krige}} for prediction with fixed parameters in the
  binomial-logit normal model, \code{\link{pois.krige.bayes}} for Bayesian prediction in the
  Poisson normal model, \code{\link[geoR]{krige.bayes}} for
  Bayesian prediction in the Gaussian spatial model. 
}

\examples{
data(b50)
\dontshow{library(geoR)}
if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) set.seed(1234)
\dontrun{
mcmc.10 <- mcmc.control(S.scale = 0.09, n.iter = 1000, phi.scale = 0.2, 
              phi.start = 4.5)
prior.10 <- prior.glm.control(phi.discrete = seq(0.2,5,0.2))
test.10 <- binom.krige.bayes(b50, locations=t(cbind(c(2.5,3.5),c(-1,3.5),c(2.5,1.5),c(-1,1.5))),
              prior = prior.10, mcmc.input = mcmc.10)
image(test.10)
persp(test.10)
}
\dontshow{
b9 <- grf(grid = expand.grid(x = seq(1, 7, l = 3), y = seq(1, 7, l = 3)), cov.pars = c(0.1, 2))
b9$units.m <- rep(4,9)
b9$data <- rbinom(9, size = rep(4,9), prob = exp(b9$data)/(1+exp(b9$data)))
##
b9.pr <- prior.glm.control(phi.discrete = c(2.5, 3))
b9.mcmc <- mcmc.control(S.scale = 0.1, phi.scale = 1.4,n.iter=10,thin=1) 
grid <- expand.grid(x = c(1, 7), y = c(1, 6))
run.b9 <- binom.krige.bayes(b9, locations = grid, prior = b9.pr, mcmc.input = b9.mcmc) 
image.kriging(locations = grid, values = run.b9$pred$median)
pred.grid <- expand.grid(x = seq(0.0125, 0.9875, l = 2), y = seq(0.0125, 0.9875, l = 2))
prior.b9 <- prior.glm.control(beta.prior = "normal", beta = 0, beta.var.std = 1, phi.prior = "exponential", phi = 2.5, phi.discrete = c(2.5,3), sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5) 
mcmc.b9 <- mcmc.control(S.scale = 0.05, phi.scale = 1.5, thin = 1, burn.in = 1, n.iter=10)
run.b9 <- binom.krige.bayes(b9, prior = prior.b9, mcmc.input = mcmc.b9, locations = pred.grid, output = list(sim.predict = TRUE))
run.mean <- apply(run.b9$pre$sim,1,mean)
rum.med.ctrl <- apply(run.b9$pre$sim,1,median)
A <- apply(run.b9$pre$sim,2,mean)
####
prior2.b9 <- prior.glm.control(beta.prior = "normal", beta = c(0,0,0), beta.var.std = matrix(cbind(c(1,0.25,0),c(0.25,0.5,0.125),c(0,0.125,0.5)),3,3), phi.prior = "exponential", phi = 2.5, phi.discrete = c(2.5,3), sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5) 
mcmc2.b9 <- mcmc.control(S.scale = 0.02, phi.scale = 1.5, thin = 1, n.iter=10)
run2.b9 <- binom.krige.bayes(b9, prior = prior2.b9, mcmc.input = mcmc2.b9, model = list(trend.d ="1st", trend.l ="1st"), locations = pred.grid, output = output.glm.control(sim.predict = TRUE))
##
prior3.b9 <- prior.glm.control(beta.prior = "normal", beta = 0,
beta.var.std = 100, phi.prior = "exponential", phi = 2.5, phi.discrete=c(2.5,3),
sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5,tausq.rel=0.1) 
mcmc3.b9 <- mcmc.control(S.scale = 0.0007, phi.scale = 1.5, thin = 1,
n.iter=10)
run3.b9 <- binom.krige.bayes(b9, prior = prior3.b9, mcmc.input =
mcmc3.b9, locations = pred.grid, output = output.glm.control(sim.predict
= TRUE), model = model.glm.control(aniso.pars = c(1,2), kappa =1))
##
prior4.b9 <- prior.glm.control(beta.prior = "flat", phi.prior = "exponential", phi = 2.5, phi.discrete = c(2.5,3), sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5) 
mcmc4.b9 <- mcmc.control(S.scale = 0.2, phi.scale = 1.5, thin = 1, n.iter=10)
run4.b9 <- binom.krige.bayes(b9, prior = prior4.b9, mcmc.input = mcmc4.b9, locations = pred.grid, output = list(sim.predict = TRUE))
mean(run3.b9$pos$phi$s)
mean(run4.b9$pos$phi$s)
AA <- apply(run.b9$pre$sim,2,mean)
####
prior.b9 <- prior.glm.control(beta.prior = "normal", beta = 0, beta.var.std = 1, phi.prior = "fixed", phi = 2.5, sigmasq.prior = "sc.inv.chisq", df.sigmasq = 50, sigmasq = 0.5) 
mcmc.b9 <- mcmc.control(S.scale = 0.05, phi.scale = 1.5, thin = 1, n.iter=10)
run3.b9 <- binom.krige.bayes(b9, prior = prior.b9, mcmc.input = mcmc.b9, locations = pred.grid)
prior4.b9 <- prior.glm.control(beta.prior = "flat", phi.prior = "fixed", phi = 2, sigmasq.prior = "sc.inv.chisq", df.sigmasq = 50, sigmasq = 0.5) 
mcmc4.b9 <- mcmc.control(S.scale = 0.2, phi.scale = 0.015, thin = 1, n.iter =10)
run4.b9 <- binom.krige.bayes(b9, prior = prior4.b9, mcmc.input = mcmc4.b9, locations = pred.grid, output = list(sim.predict = TRUE))
##
prior3.b9 <- prior.glm.control(beta.prior = "normal", beta = 0, beta.var.std = 100, phi.prior = "fixed", phi = 2, sigmasq.prior = "fixed", sigmasq = 0.5)
mcmc3.b9 <- list(S.scale = 0.0007, thin = 1, n.iter =10)
run3.b9 <- binom.krige.bayes(b9, prior = prior3.b9, mcmc.input = mcmc3.b9, locations = pred.grid)
prior4.b9 <- prior.glm.control(beta.prior = "flat", phi.prior = "fixed", phi = 2, sigmasq.prior = "fixed", sigmasq = 0.5)
mcmc4.b9 <- mcmc.control(S.scale = 0.2, thin = 1, n.iter =10)
run4.b9 <- binom.krige.bayes(b9, prior = prior4.b9, mcmc.input =
mcmc4.b9, locations = pred.grid, output = output.glm.control(sim.predict
= TRUE))
}
}

\keyword{spatial}
