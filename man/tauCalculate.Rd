\name{tauCalculate}
\alias{tauCalculate}

\title{ Posterior Predictive Probabilities from ordrating }

\description{
Calculates the quantity labeled \eqn{\tau_{pc}}{tau_{pc}} in Ho and Quinn (n.d.).  \eqn{\tau_{pc}}{tau_{pc}} can be thought of as the probability that a randomly chosen rater (from the set of observed raters) will give product \eqn{p}{p} a rating of \eqn{c}{c} given the observed data. 
}

\usage{
tauCalculate(out, ndraws = 500)
}

\arguments{
  \item{out}{ An output object from the \code{ordrating} function. }
  \item{ndraws}{ The number of Monte Carlo draws used to calculate the posterior predictive probabilities. Must be less than or equal to the number of rows in \code{out}. }
}

\details{
Section 3.2 of Ho and Quinn (n.d.) provides full details for how 
\eqn{\tau_{pc}}{tau_{pc}} is calculated. The basic ideas are the following. 

\code{tauCalculate} takes the MCMC output from \code{ordrating} and calculates the sample average (over all row units in \eqn{Y}{Y}) of the posterior predictive probability of a particular column unit in \eqn{Y}{Y} being rated as \eqn{c}{c}. This is done for all column units and ratings categories.   
}

\value{
A \eqn{P \times C}{P x C} matrix where \eqn{P}{P} is the number of columns in \eqn{Y}{Y} and \eqn{C}{C} is the number of ordinal rating categories. Each row of this matrix gives the posterior predictive probability that a randomly chosen rater, from the set of observed raters, will give product \eqn{p}{p} a rating of \eqn{c}{c}. 
}


\references{ 
  Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation 
 and Interpretation of Online Ratings Data with Model-based Figures." 
 The American Statistician.   
}


\seealso{ \code{\link[Ratings]{ordrating}},
          \code{\link[Ratings]{barplotModelBased}},
          \code{\link[Ratings]{starplotModelBased}}}


\examples{
\dontrun{
## Mondo Times example from Ho & Quinn (nd).
## may have to increase stack limit to run this example on some machines

data(Mondo)

ord.out <- ordrating(Mondo, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

tau <- tauCalculate(ord.out, 500)




## subsetting the Mondo data to include only raters who rated 5 or more 
## outlets (should avoid any stacksize problems)

Mondo.sub <- Mondo[apply(!is.na(Mondo), 1, sum) >= 5, ]
## also getting rid of outlets that are not rated now
Mondo.sub <- Mondo.sub[,apply(is.na(Mondo.sub), 2, mean) != 1] 

ord.out <- ordrating(Mondo.sub, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

tau <- tauCalculate(ord.out, 500)


}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

