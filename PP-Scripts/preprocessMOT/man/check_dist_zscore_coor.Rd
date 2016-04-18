\name{check_dist_zscore_coor}
\alias{check_dist_zscore_coor}
\title{
Check the Data Distribution, Z-Score, Batch Effects and Repeating Samples in the Data Matrix
}
\description{
This function takes the filename of the data matrix and then performs quality control using following parameters: data distribution (using z-scores, mean and median plots), batch effects (using PCA plots) and  duplicate samples (sample-sample correlation). 
}
\usage{check_dist_zscore_coor(file, output_dir)}
\arguments{
  \item{file}{
	absolute path of the file (containing the data matrix).
}
  \item{output_dir}{
	absolute path of the destination directory.
}
}

\details{
Check the data distribution, z-score, batch effects and repeating samples in the data matrix.
}

\references{Anguraj Sadanandam, et al. (2016), MOTassigner}
\author{
Pawan Poudel paawanpoudel@gmail.com
}
\examples{
	# Just the filename containing data matrix is required
	file = "/Users/ppoudel/test/data.matrix.txt"
	output_dir = "/Users/ppoudel/test/"
	check_dist_zscore_coor(file, output_dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~normalisation }
\keyword{ ~correlation }% __ONLY ONE__ keyword per line
