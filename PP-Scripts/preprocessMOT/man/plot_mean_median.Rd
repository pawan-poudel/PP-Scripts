\name{plot_mean_median}
\alias{plot_mean_median}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Plot the mean and median values for the data
}
\description{

Plot the mean and median values for the data
}
\usage{
	
plot_mean_median(data_1, file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_1}{
  data matrix for plotting the mean and median
}
  \item{file}{
  pdf file for the plot
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function (data_1, file) 
{
    data = apply(data_1[, -1], 1, as.numeric)
    means <- rowMeans(data)
    med <- apply(data, 1, median)
    stds <- apply(data, 1, sd)
    indsd <- rep(1, length(stds))
    ind <- 1
    indsd[stds < ind] <- 0
    pdf(file)
    par(mfrow = c(1, 1))
    par(mar = c(4, 3, 2, 2), oma = c(0, 0, 0, 2), mgp = c(0, 
        0.5, 0), xaxt = "s")
    plot(1:length(means), means, ylim = c(min(c(stds, means, 
        med)), max(c(stds, means, med)) + 1), type = "b", xlab = "", 
        ylab = "", main = "Average and variance profiles", col = 1, 
        axes = F)
    points(1:length(stds), stds, type = "b", col = "green")
    text(1:length(stds), stds, rownames(data), col = indsd, cex = 0.8)
    points(1:length(med), med, type = "b", col = "red")
    box()
    axis(2, col = "black")
    axis(1, at = seq(1, nrow(data), ceiling(nrow(data)/25)), 
        labels = seq(1, nrow(data), ceiling(nrow(data)/25))[1:length(seq(1, 
            nrow(data), ceiling(nrow(data)/25)))], las = 2)
    mtext("genes", side = 1, col = "black", line = 2.2, cex = 0.8, 
        las = 1)
    legend("topright", bty = "n", c("means", "median", "standard deviations"), 
        col = 1:3, text.col = 1:3, pch = 1)
    dev.off()
    z_score = apply(data[, -1], 1, standardise_data)
    zscore_plot_f = gsub(".pdf", "zscore.pdf", file)
    pdf(zscore_plot_f)
    par(mfrow = c(1, 2))
    hist(data, breaks = 100, xlab = expression("gene"), main = "Distribution of genes", 
        freq = FALSE, cex.lab = 1)
    lines(density(data), col = "red", lwd = 2)
    box()
    hist(z_score, breaks = 100, xlab = expression("genes"), main = "Zscore", 
        freq = FALSE, cex.lab = 1)
    lines(density(z_score), col = "red", lwd = 2)
    box()
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mean }
\keyword{ ~median }% __ONLY ONE__ keyword per line
