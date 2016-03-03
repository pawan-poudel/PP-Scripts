\name{impute_missing_values}
\alias{impute_missing_values}
\title{
   Imputing the missing values
}
\description{
    Imputing the missing values
}
\usage{
impute_missing_values(pancan12_genomeMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pancan12_genomeMatrix}{
  	data matrix with the genes in the first column followed by samples
}

}
\details{

Imputing the missing values if required

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
function (pancan12_genomeMatrix) 
{
    rownames(pancan12_genomeMatrix) = pancan12_genomeMatrix[, 
        1]
    pancan12_genomeMatrix_1 = apply(pancan12_genomeMatrix[, -1], 
        2, as.numeric)
    rownames(pancan12_genomeMatrix_1) = pancan12_genomeMatrix[, 
        1]
    pancan12_genomeMatrix.imputed = impute.knn(as.matrix(pancan12_genomeMatrix_1))
    genes = rownames(pancan12_genomeMatrix.imputed$data)
    rownames(pancan12_genomeMatrix.imputed$data) = NULL
    data_to_write = cbind(genes, pancan12_genomeMatrix.imputed$data)
    rownames(data_to_write) = NULL
    colnames(data_to_write)[1] = "Genes"
   
    return(data_to_write)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~impute }
\keyword{ ~impute }% __ONLY ONE__ keyword per line
