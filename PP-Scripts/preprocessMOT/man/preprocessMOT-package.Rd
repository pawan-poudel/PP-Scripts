\name{preprocessMOT-package}
\alias{preprocessMOT-package}
\alias{preprocessMOT}
\docType{package}
\title{
This packages is a part of MOTassigner paper, which performs QC of the dataset used in the paper.
}
\description{
	
This package takes the raw cel files (preferably Hgu133plus2), and then creates the final processed QC matrix. It filters the samples based on (1 -/+0.05) NUSE critera, check batch effects, plots data distribution and performs sample sample correlations to remove the repeating samples. Additionally, QC can be performed on any data matrix with the genes in the first column followed by the samples.
}
\details{
\tabular{ll}{
Package: \tab preprocessMOT\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-03-03\cr
License: \tab GPL=>2 \cr
}
	QC, check batch effects, sample sample correlation
}
\author{
Pawan Poudel
Maintainer: Pawan Poudel <paawanpoudel@gmail.com>
}
\references{
	This package is the part of MOTassigner paper, which performs QC of the dataset used in the paper.
}

\keyword{preprocess }
