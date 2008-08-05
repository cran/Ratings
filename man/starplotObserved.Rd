\name{starplotObserved}
\alias{starplotObserved}

\title{ Observed Data Starplots of Ratings Data }

\description{
Creates the observed data starplot discussed in Ho and Quinn (n.d.).
}

\usage{
starplotObserved(Y, catnames, starcol = "darkgray", starsize = 0.2, ...)
}

\arguments{
  \item{Y}{ Matrix of data to be analyzed. Entries must be integers from 
            \eqn{1, \ldots, C}{1, ..., C} or \code{NA} where \eqn{C}{C} 
            is the number of ordinal categories. Raters are on the rows and 
            products are on the columns.}
  \item{catnames}{ A character vector of names for the rating categories. Should have \eqn{C}{C} elements. First element should give the name for a rating of 1, second element the names for a rating of 2, etc.}
  \item{starcol}{ Color of stars. }
  \item{starsize}{  Positive scalar that regulates the size of the stars in the plot. }
  \item{\dots}{ Other arguments to plotting functions. }
}

\seealso{ \code{\link[Ratings]{ordrating}},
          \code{\link[Ratings]{starplotModelBased}},
          \code{\link[Ratings]{barplotObserved}}}


\references{ 
Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation and Interpretation of Online Ratings Data with Model-based Figures." The American Statistician.  
}   	     
 
\examples{
\dontrun{
## Mondo Times example from Ho & Quinn (nd).

data(Mondo)

## just the first 12 cols
Mondo.sub <- Mondo[,1:12]

starplotObserved(Mondo.sub, starcol="green",
                  catnames = c("Awful", "Poor", "Average", 
                               "Very Good", "Great"))

}
}


\keyword{ hplot }


