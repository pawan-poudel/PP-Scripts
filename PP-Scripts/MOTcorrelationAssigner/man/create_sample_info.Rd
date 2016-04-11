\name{create_sample_info}
\alias{create_sample_info}
\title{
Create the tab seperated sample info file with the samples and organ information
}
\description{
This script is a part of MOTcorrelationassigner pipeline. The script uses foldername to assign the organ information thus please segregate the samples using the organ name
}
\usage{
create_sample_info(project_folder, outputdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{project_folder}{
  location of the folder where the gene expression matrix lies. Please name the folder based on the organ information
}
  \item{outputdir}{
  location of the output directory
}
}
\details{
This script is a part of MOTassinger pipeline.
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
  
  project_folder="/Users/ppoudel/test/correlation_assigner/organs/"
  outputdir="/Users/ppoudel/test/correlation_assigner/output/"
  
  create_sample_info(project_folder, outputdir)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~createsampleinfo }
\keyword{ ~sampleinfo }% __ONLY ONE__ keyword per line
