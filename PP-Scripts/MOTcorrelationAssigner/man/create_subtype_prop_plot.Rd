\name{create_subtype_prop_plot}
\alias{create_subtype_prop_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This script create the proportion plot after classifying the data set to one of the MOTassigner subtypes
}
\description{
This function creates a subtype proportion plot for each organs
}
\usage{
create_subtype_prop_plot(file, output_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
	filename of the annotation file that you want to plot
}
  \item{output_dir}{
	final location to put the output
}

}
\details{
This script is a part of MOTassigner pipeline
}

\references{
MOTassinger 2016
}
\author{
Pawan Poudel
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
  file="/Users/ppoudel/test/correlation_assigner/annodir/anno.txt"
  output_dir="/Users/ppoudel/test/correlation_assigner/output/"
  
  create_subtype_prop_plot(file=file, output_dir=output_dir)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
