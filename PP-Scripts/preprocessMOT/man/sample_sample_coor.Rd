\name{sample_sample_coor}
\alias{sample_sample_coor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Sample sample correlation
}
\description{
	Checks and removes the samples with the correlation score >0.99
}
\usage{
sample_sample_coor(data, file1, filename_coor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{file1}{
%%     ~~Describe \code{file1} here~~
}
  \item{filename_coor}{
%%     ~~Describe \code{filename_coor} here~~
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
function (data, file1, filename_coor) 
{
    dx = apply(data[, -1], 2, as.numeric)
    my_genes = data[, 1]
    M <- cor(dx)
    new_M = M[, ] > 0.99
    similar_samples_data_frame = which(new_M, arr.ind = TRUE)
    similar_samples_index = similar_samples_data_frame[which(similar_samples_data_frame[, 
        1] != similar_samples_data_frame[, 2]), ]
    log_f = paste0(file1, "_R_logs.txt")
    if (length(similar_samples_index) == 0) {
        write.table(data, filename_coor, sep = "\t", row.names = FALSE, 
            quote = FALSE)
        cat("This doesnot have repeated sampels", file = log_f, 
            sep = "\n", append = TRUE)
    }
    else {
        p = unique(similar_samples_index[, 1])
        q = unique(similar_samples_index[, 2])
        i = intersect(p, q)
        cat("Removing the repeating samples \n")
        cat("Following are the very similar samples ", similar_samples_index, 
            sep = "\n")
        cat("Following samples have been removed ", as.vector(similar_samples_index), 
            "\n")
        cat(colnames(data)[i], file = log_f, sep = "\n")
        cat(as.vector(similar_samples_index), file = log_f, append = TRUE)
        corr_score <- gsub(".txt", "correlation_scores_samples.txt", 
            file1)
        data.Unique = data[, -i]
        write.table(data.Unique, filename_coor, sep = "\t", row.names = FALSE, 
            quote = FALSE)
        write.table(similar_samples_data_frame, corr_score, sep = "\t", 
            row.names = FALSE, quote = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
