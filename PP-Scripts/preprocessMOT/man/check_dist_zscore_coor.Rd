\name{check_dist_zscore_coor}
\alias{check_dist_zscore_coor}
\title{
Check the data distribution, zscore, batch effects and repeating samples in the data matrix
}
\description{
This function takes the filenames of the data matrix and then performs quality control using following parameters: data distribution (using zscores, mean and median plots), batch effects (using PCA plots) and  sample samples correlation (removing repeating samples). 
}
\usage{check_dist_zscore_coor(file)}
\arguments{
  \item{file}{
	Absolute path of the file (containing data matrix)
}
}

\details{
Check the data distribution, zscore, batch effects and repeating samples in the data matrix
}

\references{Anguraj Sadanandam, et all (2016), MOTassigner}
\author{
Pawan Poudel paawanpoudel@gmail.com
}
\examples{
	# Just the filename containing data matrix is required
	file="/Users/ppoudel/test/data.matrix.txt"
	check_dist_zscore_coor(file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~normalisation }
\keyword{ ~correlation }% __ONLY ONE__ keyword per line
