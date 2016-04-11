\name{correlationAssigner}
\alias{correlationAssigner}
\title{
	Correlation assinger
}
\description{
This function uses the PAM centroid scores (of each subtypes) to classify samples into one of the subtypes based on the similarity in the gene expression pattern.
}
\usage{
correlationAssigner(filename, signatures, outputdir, sd = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
	filename of the gene expression matrix that you want to classify
}
  \item{signatures}{
	filename of the pam centroids
}
  \item{outputdir}{
	final directory that you want to put your output
}
  \item{sd}{
	standard deviation threshold required, default=0
}

}
\details{
This script is used to classify the validation data set using the PAM centroids scores
}

\references{
MOTassigner, 2016
}
\author{
Pawan Poudel
}
\note{
This script is a part of MOTcorrelationassigner
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

  filename="/Users/ppoudel/test/correlation_assigner/organs/endo/Endometrium.txt"
  signatures="/Users/ppoudel/test/correlation_assigner/centroids/PAM-centroids.txt"
  outputdir="/Users/ppoudel/test/correlation_assigner/output"
  
  correlationAssigner(filename=filename, signatures=signatures, outputdir=outputdir, sd = 0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~correlationassigner }
\keyword{ ~correlation }% __ONLY ONE__ keyword per line
