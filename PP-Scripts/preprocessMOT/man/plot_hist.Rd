\name{plot_hist}
\alias{plot_hist}
\title{
Plotting the histogram from the data matrix
}
\description{
Plotting the histogram using the data matrix and medain centered data matrix
}
\usage{
plot_hist(data, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  data matrix
}
  \item{file}{
   pdf file for histogram file
}
}
\details{
Plotting the histogram of the data and the median centered data
}
\references{
Anguraj Sadanandam, et all (2016), MOTassigner
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
data =matrix(rnorm(100),10,10)
file="hist.pdf"
plot_hist(data, file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
