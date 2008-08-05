\name{barplotObserved}
\alias{barplotObserved}

\title{ Observed Data Barplots of Ratings Data }

\description{
Creates the observed data barplot discussed in Ho and Quinn (n.d.).
}

\usage{
barplotObserved(Y, catnames, barcol = "darkgray", top.limit = 0.5, 
                scale.factor = 1, ...)
}

\arguments{
  \item{Y}{ Matrix of data to be analyzed. Entries must be integers from 
            \eqn{1, \ldots, C}{1, ..., C} or \code{NA} where \eqn{C}{C} 
            is the number of ordinal categories. Raters are on the rows and 
            products are on the columns.}
  \item{catnames}{ A character vector of names for the rating categories. Should have \eqn{C}{C} elements. First element should give the name for a rating of 1, second element the names for a rating of 2, etc.}
  \item{barcol}{ The color for the bars in the barplot. }
  \item{top.limit}{ The maximum value of the \eqn{y}{y} axis. Should be between 0 and 1. }
  \item{scale.factor}{ Scaling factor that adjusts the scaling of the \eqn{y} axis. Any positive value is allowed. }
  \item{\dots}{ Other arguments to plotting functions. }
}

\seealso{ \code{\link[Ratings]{ordrating}},
          \code{\link[Ratings]{barplotModelBased}},
          \code{\link[Ratings]{starplotObserved}}}


\references{ 
Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation and Interpretation of Online Ratings Data with Model-based Figures." The American Statistician.  
}   	     
 
\examples{
\dontrun{
## Mondo Times example from Ho & Quinn (nd).

data(Mondo)

## just the first 12 cols
Mondo.sub <- Mondo[,1:12]

barplotObserved(Mondo.sub, top.limit=1, scale.factor=.7,
                  catnames = c("Awful", "Poor", "Average", 
                               "Very Good", "Great"))

}
}


\keyword{ hplot }

