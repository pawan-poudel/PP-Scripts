\name{makePcaPlot}
\alias{makePcaPlot}
\title{
	Create the PCA plot using the ggplot
}
\description{
Create the PCA plot using the ggplot
}
\usage{makePcaPlot(x, sample_info, plot_file, title = "")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	data matrix
}
  \item{sample_info}{
	should contain 2 columns- sample name(first), group info (second)
}
  \item{plot_file}{
  	name of the pdf file for the PCA plot

}
  \item{title}{
  	title for the plot
	 }
}
\details{
Create the PCA plot
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
x=matrix(rnorm(100),10,10)
rownames(x)=paste0("Genes", c(1:10))
colnames(x)=paste0("Sample", c(1:10))
sample_info=data.frame(colnames(x), rep(1, ncol(x)))
makePcaPlot( x=x, sample_info=sample_info,plot_file="test.pdf",title="test-plot")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~pca }% __ONLY ONE__ keyword per line
