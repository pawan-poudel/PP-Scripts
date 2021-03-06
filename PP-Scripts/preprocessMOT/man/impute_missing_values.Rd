\name{impute_missing_values}
\alias{impute_missing_values}
\title{
   Impute missing values
}
\description{
    Impute missing values.
}
\usage{
impute_missing_values(pancan12_genomeMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pancan12_genomeMatrix}{
  	a genes by samples data matrix, containing the gene names in the first column.
}
}
\details{
Imputing missing values, if required.
}

\references{
Hastie T, Tibshirani R, Narasimhan B and Chu G. impute: impute: Imputation for microarray data. R package version 1.44.0.

Anguraj Sadanandam, et al. (2016), MOTassigner
}
\author{
Pawan Poudel paawanpoudel@gmail.com
}

\examples{
  x = matrix(rnorm(6),2,3)
  pancan12_genomeMatrix = data.frame(c("Gene1","Gene2"), x)
  imputed_data = impute_missing_values(pancan12_genomeMatrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~impute }
\keyword{ ~impute }% __ONLY ONE__ keyword per line
