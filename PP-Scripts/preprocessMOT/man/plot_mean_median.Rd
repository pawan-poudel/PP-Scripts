\name{plot_mean_median}
\alias{plot_mean_median}
\title{
Plot the mean and median values for the data.
}
\description{
Plot the mean and median values for the data.
}
\usage{plot_mean_median(data_1, file)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_1}{
  data matrix for plotting the mean and median.
}
  \item{file}{
  PDF filename for the plot generated.
}
}
\details{
Plot the mean and median values for the data.
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
data_1 = matrix(rnorm(100),10,10)
file = "plot_file.pdf"
plot_mean_median(data_1, file)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mean }
\keyword{ ~median }% __ONLY ONE__ keyword per line
