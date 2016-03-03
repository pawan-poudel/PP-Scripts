\name{makePcaPlot}
\alias{makePcaPlot}
\title{
	Create the PCA plot using the ggplot
}
\description{
Create the PCA plot using the ggplot
}
\usage{
makePcaPlot(x, sample_info, plot_file, title = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	data matrix
}
  \item{sample_info}{
	should contain 2 columns- sample name(first), group info (second)
}
  \item{title}{
	title for the plot
}
  \item{plot_file}{
	 name of the pdf file for the PCA plot
}
}
\details{
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
function (x, sample_info, title = "", plot_file) 
{
    sample_info <- sample_info[order(sample_info[, 2]), ]
    data.1 <- subset(x, select = as.character(sample_info[, 1]))
    data.1 <- apply(data.1, 2, as.numeric)
    ind <- rowVars(data.1) < .Machine$double.eps
    dx <- data.1[!ind, ]
    data <- t(apply(dx, 1, scale))
    mydata <- t(data)
    rownames(mydata) <- colnames(data.1)
    mydata.pca <- prcomp(mydata, retx = TRUE, center = TRUE, 
        scale. = TRUE)
    scores <- mydata.pca$x
    cex = 3
    group_name <-  sample_info[, 2]
    df <- data.frame(rownames(scores), scores[, c(1, 2)], group_name)
    colnames(df) <- c("samples", "PCA1", "PCA2", "Groups")
    p <- ggplot(df, aes(PCA1, PCA2, colour = Groups)) + geom_point(size = I(2)) + 
        labs(x = "PCA-1", y = "PCA-2", color = "Groups") + scale_fill_manual(values = as.factor(df$Groups))
    p <- p + labs(title = "PCA plot") + scale_colour_Publication() + 
        theme_Publication()
    score_file <- gsub(".pdf", "pca.scores.txt", plot_file)
    write.table(scores, score_file, sep = "\t")
    ggsave(p, file = plot_file)
    return(df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~pca }
\keyword{ ~pca }% __ONLY ONE__ keyword per line
