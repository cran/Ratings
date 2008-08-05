\name{starplotModelBased}
\alias{starplotModelBased}

\title{ Model-based Starplots of Ratings Data }

\description{
Creates the model-based starplot described in Ho and Quinn (n.d.). 
}

\usage{
starplotModelBased(tau.mat, colvec = NULL, starsize = 0.2, 
                   interpolation.level = 200, ...)
}

\arguments{
  \item{tau.mat}{
A \eqn{P \times C}{P x C} matrix where \eqn{P}{P} is the number of products being rated and \eqn{C}{C} is the number of rating categories. The \eqn{p}{p}th row, \eqn{c}{c}th column of \code{tau.mat} should give the probability that product \eqn{p}{p} is equal to category \eqn{c}{c}.
}
  \item{colvec}{ Vector of rgb colors used for the plot. }
  \item{starsize}{ Positive scalar that regulates the size of the stars in the plot.  }
  \item{interpolation.level}{ Regulates the smoothness of the color scale. \code{interpolation.level = 1} simply uses the colors in \code{colvec}. Values of \code{interpolation.level > 1} provide increasing amounts of interpolation between the values in \code{colvec}. Larger values of \code{interpolation.leve} produce smoother transitions between colors. }
  \item{\dots}{ Other arguments to plotting functions. }
}


\seealso{ \code{\link[Ratings]{ordrating}},
          \code{\link[Ratings]{tauCalculate}},
          \code{\link[Ratings]{barplotModelBased}}}


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

starplotModelBased(tau.sub)

## more informative labels
colnames(tau.sub) <- c("Awful", "Poor", "Average", "Very Good", "Great")

starplotModelBased(tau.sub)





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

starplotModelBased(tau.sub)

## more informative labels
colnames(tau.sub) <- c("Awful", "Poor", "Average", "Very Good", "Great")

starplotModelBased(tau.sub)

## a different color scheme
mycol <- rgb(red=c(255, 243, 231, 219, 207, 159, 82.5, 30, 0, 0, 255),
              green=c(255.0, 250.1, 245.2, 240.3, 235.4, 215.8, 175, 144,
                0, 0, 69),
              blue=c(255, 255, 255, 255, 255, 255, 255, 255, 238, 183, 0),
              maxColorValue=255
              )

starplotModelBased(tau.sub, colvec=mycol)

}
}


\keyword{ hplot }



