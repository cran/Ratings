\name{barplotModelBased}
\alias{barplotModelBased}

\title{ Model-based Barplots of Ratings Data }

\description{
Creates the model-based barplot described in Ho and Quinn (n.d.). 
}

\usage{
barplotModelBased(tau.mat, barcol = "darkgray", top.limit = 0.5, 
                  scale.factor = 1, ...)
}

\arguments{
  \item{tau.mat}{A \eqn{P \times C}{P x C} matrix where \eqn{P}{P} is the number of products being rated and \eqn{C}{C} is the number of rating categories. The \eqn{p}{p}th row, \eqn{c}{c}th column of \code{tau.mat} should give the probability that product \eqn{p}{p} is equal to category \eqn{c}{c}. }
  \item{barcol}{ The color for the bars in the barplot. }
  \item{top.limit}{ The maximum value of the \eqn{y}{y} axis. Should be between 0 and 1. }
  \item{scale.factor}{ Scaling factor that adjusts the scaling of the \eqn{y} axis. Any positive value is allowed. }
  \item{\dots}{ Other arguments to plotting functions. }
}


\seealso{ \code{\link[Ratings]{ordrating}},
          \code{\link[Ratings]{tauCalculate}},
          \code{\link[Ratings]{starplotModelBased}}}


\references{ 
Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation and Interpretation of Online Ratings Data with Model-based Figures." The American Statistician.  
}   	     
 
\examples{
\dontrun{
## Mondo Times example from Ho & Quinn (nd).
## may have to increase stack limit to run this example on some machines

data(Mondo)

ord.out <- ordrating(Mondo, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

## get rating probabilities
tau <- tauCalculate(ord.out, 500)

## just the labeled outlets
tau.sub <- tau[-grep("thetaOutlet", rownames(tau)),]

## clean up names
rownames(tau.sub) <-  gsub("theta", "", rownames(tau.sub))

barplotModelBased(tau.sub)

## more informative labels
colnames(tau.sub) <- c("Awful", "Poor", "Average", "Very Good", "Great")

barplotModelBased(tau.sub)





## subsetting the Mondo data to include only raters who rated 5 or more 
## outlets (should avoid any stacksize problems)

Mondo.sub <- Mondo[apply(!is.na(Mondo), 1, sum) >= 5, ]
## also getting rid of outlets that are not rated now
Mondo.sub <- Mondo.sub[,apply(is.na(Mondo.sub), 2, mean) != 1] 

ord.out <- ordrating(Mondo.sub, beta.constraint=1, tune=.035, 
                     ma=1, mb=-5, vinva=1, vinvb=0.05,
		     gamma.start=c(-300, 0, 1.5, 3.0, 4.5, 300),
		     thin=20, burnin=20000, mcmc=100000, verbose=1000)

## get rating probabilities
tau <- tauCalculate(ord.out, 500)

## just the labeled outlets
tau.sub <- tau[-grep("thetaOutlet", rownames(tau)),]

## clean up names
rownames(tau.sub) <-  gsub("theta", "", rownames(tau.sub))

barplotModelBased(tau.sub)

## more informative labels
colnames(tau.sub) <- c("Awful", "Poor", "Average", "Very Good", "Great")

barplotModelBased(tau.sub)


}
}


\keyword{ hplot }

