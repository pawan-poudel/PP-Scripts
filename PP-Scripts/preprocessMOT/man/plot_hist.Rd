\name{plot_hist}
\alias{plot_hist}
\title{

Plotting the histogram from the data matrix
}
\description{
Plotting the histogram using the data matrix and medain centered data matrix
}
\usage{
plot_hist(data, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  data matrix
}
  \item{file}{
   pdf file for histogram file
}
}
\details{

Plotting the histogram of the data and the median centered data

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, file) 
{
    data = apply(data[, -1], 1, as.numeric)
    rowMed = apply(data, 1, median, na.rm = TRUE)
    data1 = data - rowMed
    pdf(file)
    par(mfrow = c(1, 2))
    hist(data, breaks = 100, xlab = expression("gene"), main = "Distribution of genes", 
        freq = FALSE, cex.lab = 1)
    lines(density(data), col = "red", lwd = 2)
    box()
    hist(data1, breaks = 100, xlab = expression("genes"), main = "Distribution of genes after median centering", 
        freq = FALSE, cex.lab = 1)
    lines(density(data1), col = "red", lwd = 2)
    box()
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
