\name{Mondo}
\alias{Mondo}
\docType{data}
\title{Mondo Times News Media Content Quality Data}

\description{
This data matrix contains ratings of the content quality of news media as collected by Mondo Times. 
}

\usage{
data(Mondo)
}

\format{
Rows of this data matrix correspond to raters and columns to news outlets. The entry in the ith row and jth columns corresponds to the rating given by rater i to news outlet j. These ratings are on a five point scale: 1 = "Awful", 2 = "Poor", 3 = "Average", 4 = "Very Good", and 5 = "Great". NA values indicate that the rater in question did not rate the news outlet in question. 

The first 12 columns of the data matrix have the full names of the news outlets so that the analysis in Ho and Quinn (forthcoming) can be replicated. The remaining news outlets are not named due to a confidentiality agreement with Mondo Times. 
}

\source{
  Mondo Times \url{http://www.mondotimes.com/} 
}
\references{
  Ho, Daniel E. and Kevin M. Quinn. forthcoming. "Improving the Presentation and Interpretation of Online Ratings Data with Model-based Figures." The American Statistician. 
}
\examples{
data(Mondo)
}
\keyword{datasets}
