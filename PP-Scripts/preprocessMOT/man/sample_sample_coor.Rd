\name{sample_sample_coor}
\alias{sample_sample_coor}
\title{
	Sample-sample correlation
}
\description{
	Checks correlation between samples and removes samples with correlation score > 0.99, in order to remove duplicate samples.
}
\usage{sample_sample_coor(data, file1, filename_coor)}
\arguments{
  \item{data}{
  data matrix to check for sample-sample correlation.
}
  \item{file1}{
  output file to write the names of matching samples.
}
  \item{filename_coor}{
  final correlation file.
}
}
\details{
Checks and removes samples with correlation score > 0.99.
}
\references{Anguraj Sadanandam, et al. (2016), MOTassigner}
\author{
Pawan Poudel paawanpoudel@gmail.com
}

\examples{
x = matrix(rnorm(100),10,10)
rownames(x) = paste0("Genes", c(1:10))
colnames(x) = paste0("Sample", c(1:10))
sample_sample_coor(data = x, file1 = "matching samples.txt", filename_coor = "final unique samples.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
