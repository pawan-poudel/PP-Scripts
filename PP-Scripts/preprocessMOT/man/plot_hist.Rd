\name{plot_hist}
\alias{plot_hist}
\title{
Plotting a histogram of the data matrix
}
\description{
Plotting a histogram of the original data matrix and median centered data matrix.
}
\usage{
plot_hist(data, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  data matrix.
}
  \item{file}{
   PDF filename for histogram plot generated.
}
}
\details{
Plotting a histogram of the original data matrix and median centered data matrix.
}
\references{
Anguraj Sadanandam, et al. (2016), MOTassigner
}
\author{
Pawan Poudel paawanpoudel@gmail.com
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
data = matrix(rnorm(100),10,10)
file = "hist.pdf"
plot_hist(data, file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
