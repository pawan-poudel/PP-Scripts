\name{run_CorrelationAssigner}
\alias{run_CorrelationAssigner}

\title{
Correlation assigner
}
\description{
This is the main function that takes the input and then runs the other functions
}
\usage{
run_CorrelationAssigner(project_folder, signatures, sd = 0, outputdir, sig_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{project_folder}{
	location where the data is present. Please make sure that data from each organ are place in seperate folder named by organs
}
  \item{signatures}{
 pam centroids file
}
  \item{sd}{
 standard deviation cutoff [default is set to 0]
}
  \item{outputdir}{
	location of the folder where you would like to put your output
}
  \item{sig_file}{
	additional message or informtaion that you would like to append to the output files
}
}
\details{
This is a part of MOTassigner pipeline
}
\references{
MOTassigner, 2016
}
\author{
Pawan Poudel
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

  project_folder="/Users/ppoudel/test/correlation_assigner/organs"
  signatures="/Users/ppoudel/test/correlation_assigner/centroids/PAM-centroids.txt"
  outputdir="/Users/ppoudel/test/correlation_assigner/output"
  sig_file="Pancan-test"
  
  run_CorrelationAssigner(project_folder=project_folder, signatures=signatures, sd = 0, outputdir=outputdir, sig_file=sig_file)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~runCorrelation }
\keyword{ ~correlationassigner }% __ONLY ONE__ keyword per line
