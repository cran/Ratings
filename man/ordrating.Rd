\name{ordrating}
\alias{ordrating}


\title{ MCMC for Unidimensional Ordinal IRT Model }


\description{
  This function generates a sample from the posterior distribution of a
  unidimensional ordinal item response theory (IRT) model, with Gaussian
  priors on the ability and item parameters.  The user
  supplies data and priors, and a sample from the posterior
  distribution is returned as an mcmc object, which can be
  subsequently analyzed with functions provided in the coda
  package.
}

\usage{
ordrating(Y, beta.constraint = NULL, theta.neg.index = NULL, 
	     theta.pos.index = NULL, vinva = 0.2, vinvb = 0.2, ma = 0, 
	     mb = 1, theta.start = NULL, gamma.start = NULL, 
	     burnin = 1000, mcmc = 10000, thin = 1, tune = 1, 
	     verbose = 0, seed = NA)
}

\arguments{
  \item{Y}{ Matrix of data to be analyzed. Entries must be integers from 
            \eqn{1, \ldots, C}{1, ..., C} or \code{NA} where \eqn{C}{C} 
            is the number of ordinal categories. Items are on the rows and 
            subjects are on the columns. }

  \item{beta.constraint}{ Possible inequality constraint on all 
    \eqn{\beta}{beta} parameters. \code{beta.constraint = NULL} (default) 
    implies no constraint, \code{beta.constraint > 0} implies all 
    \eqn{\beta}{beta} values are constrained to be positive, and 
    \code{beta.constraint < 0} implies all \eqn{\beta}{beta} values are 
    constrained to be negative.  }

  \item{theta.neg.index}{ An index of \eqn{\theta}{theta} that is constrained 
   to be negative. To constrain the \eqn{\theta}{\theta} for the subject whose 
   data is in the jth column of \code{Y} to be negative one would set 
   \code{theta.neg.index = j}. }

  \item{theta.pos.index}{ An index of \eqn{\theta}{theta} that is constrained 
   to be positive. To constrain the \eqn{\theta}{\theta} for the subject whose 
   data is in the jth column of \code{Y} to be positive one would set 
   \code{theta.pos.index = j}. }

  \item{vinva}{ The prior inverse variance for each \eqn{\alpha_r}{alpha_r}. }

  \item{vinvb}{ The prior inverse variance for each \eqn{\beta_r}{beta_r}. }

  \item{ma}{ The prior mean for each \eqn{\alpha_r}{alpha_r}. }

  \item{mb}{ The prior mean for each \eqn{\beta_r}{beta_r}. }

  \item{theta.start}{ Starting values for \eqn{\theta}{theta}. If \code{NULL} 
  starting values are calculated based on the sample means of \code{Y}. }

  \item{gamma.start}{ Starting values for \eqn{\gamma}{gamma}. If \code{NULL} starting values are chosen to be equally spaced between 1 and the number of ordinal categories \eqn{C}{C}.  }

  \item{burnin}{ The number of burn-in iterations for the sampler. }

  \item{mcmc}{ The number of MCMC iterations for the sampler. }

  \item{thin}{ The thinning interval used in the simulation.  The number of
    Gibbs iterations must be divisible by this value. }

  \item{tune}{ The scalar tuning parameter for the Metropolis-Hastings
    sampling. Must be strictly positive.}

  \item{verbose}{ A switch which determines whether or not the progress of
      the sampler is printed to the screen.   If \code{verbose} is greater
      than 0 then every \code{verbose}th iteration will be printed to the
      screen. }

  \item{seed}{ The seed for the random number generator.  If NA, the Mersenne
    Twister generator is used with default seed 12345; if an integer is 
    passed it is used to seed the Mersenne twister.  The user can also
    pass a list of length two to use the L'Ecuyer random number generator,
    which is suitable for parallel computation.  The first element of the
    list is the L'Ecuyer seed, which is a vector of length six or NA (if NA 
    a default seed of \code{rep(12345,6)} is used).  The second element of 
    list is a positive substream number. See the MCMCpack 
    specification for more details. }
}

\details{
Let \eqn{r = 1, \ldots, R}{r = 1, ..., R} index items, 
    \eqn{p = 1, \ldots, P}{p = 1, ..., P} index subjects, and 
    \eqn{c = 1, \ldots, C}{c = 1, ..., C} index ordinal rating categories. 

\eqn{Y}{Y} is an \eqn{R \times P}{R x P} matrix with elements in 
\eqn{\{1, \ldots, C, NA \}}{{1, ..., C, NA}}. \code{NA} denotes missing data 
that are assumed to be missing at random (MAR). 

The distribution of \eqn{Y}{Y} is governed by a \eqn{R \times P}{R x P} matrix 
of latent variables \eqn{Y^*}{Y*} and a series of cutpoints 
\eqn{\gamma}{gamma}. 

More specifically, 

\deqn{y^*_{rp} = \alpha_r + \beta_r \theta_p + \epsilon_{rp}}{y*_{rp} = alpha_r + beta_r theta_p + epsilon_{rp}} 

where 
\eqn{\epsilon_{rp} \sim N(0, 1)}{epsilon_{rp} ~iid N(0, 1)}.
It is assumed that \eqn{y_{rp} = c}{y_{rp = c}} if and only if 

\deqn{y^*_{rp} \in (\gamma_{c-1}, \gamma_{c}]}{y*_{rp} in (gamma_{c-1}, gamma_c]}. 

The priors used for this model are that each \eqn{\alpha_r}{alpha_r} is iid Gaussian, each \eqn{\beta_r}{beta_r} is iid Gaussian, each \eqn{\theta_p}{theta_p} is standard normal, and the elements of \eqn{\gamma}{\gamma} are improper uniform with all parameters assumed to be a priori independent.


\code{ordrating} simulates from the posterior distribution using
  a Metropolis-Hastings within Gibbs sampling algorithm. The algorithm
  employed is based on work by Cowles (1996).  Note that the first
  element  \eqn{\gamma_1}{gamma_1} is normalized to zero, and thus not 
  returned in the mcmc object.

\code{ordrating} fits a model that is a special case of the model fit by \code{MCMCordfactanal} in the \code{MCMCpack} package. The primary differences are the types of identifying constraints employed, the dimensionality of the \eqn{\theta}{theta} and \eqn{\beta}{beta} parameters, and the computational speed. Because \code{ordrating} fits a narrower class of models it can be optimized for speed much more effectively. 

}


\value{
   An mcmc object that contains the posterior sample.  This 
   object can be summarized by functions provided by the \code{coda} package.
}

\references{

  Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation 
 and Interpretation of Online Ratings Data with Model-based Figures." 
 The American Statistician.   

  M. K. Cowles. 1996. ``Accelerating Monte Carlo Markov Chain Convergence for
  Cumulative-link Generalized Linear Models.'' \emph{Statistics and Computing.}
  6: 101-110.
     
  Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling.'' 
  Springer: New York.
  
    Kevin M. Quinn. 2004. ``Bayesian Factor Analysis for Mixed Ordinal and
  Continuous Responses.'' \emph{Political Analysis}. 12: 338-353.

  Shawn Treier and Simon Jackman. 2003. ``Democracy as a Latent Variable.'' 
  Paper presented at the Midwest Political Science Association Annual Meeting.
     
   Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2002.
   \emph{Output Analysis and Diagnostics for MCMC (CODA)}.
   \url{http://www-fis.iarc.fr/coda/}.
 }


\seealso{\code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
         \code{\link[MCMCpack]{MCMCordfactanal}}}

\examples{
\dontrun{
## Mondo Times example from Ho & Quinn (nd).

data(Mondo)

ord.out <- ordrating(Mondo, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

plot(ord.out)
summary(ord.out) 




## subsetting the Mondo data to include only raters who rated 5 or more 
## outlets 

Mondo.sub <- Mondo[apply(!is.na(Mondo), 1, sum) >= 5, ]
## also getting rid of outlets that are not rated now
Mondo.sub <- Mondo.sub[,apply(is.na(Mondo.sub), 2, mean) != 1] 

ord.out <- ordrating(Mondo.sub, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

plot(ord.out)
summary(ord.out) 

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
