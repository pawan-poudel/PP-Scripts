\name{impute_screenExpr}
\alias{impute_screenExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Removes the repeating genes/probes and impute the missing values (if required)
}
\description{
This script removes the non variable genes using the standard deviation criteria (default=0)
}
\usage{
impute_screenExpr(exprFiles, sdCutoffs, outdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprFiles}{
    name of the expression file
}
  \item{sdCutoffs}{
  standard deviation cutoff [default=0]
}
  \item{outdir}{
  location of output directory
}
}
\details{
This script is a part of MOTassigner pipeline
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

  exprFiles="/Users/ppoudel/test/correlation_assigner/organs/endo/Endometrium.txt"
  sdCutoffs=0
  outdir="/Users/ppoudel/test/correlation_assigner/output/"
  
  impute_screenExpr(exprFiles=exprFiles, sdCutoffs==0, outdir=outdir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~impute }
\keyword{ ~screenExpr }% __ONLY ONE__ keyword per line
