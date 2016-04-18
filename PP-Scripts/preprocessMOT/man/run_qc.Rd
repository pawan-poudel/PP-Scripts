\name{run_qc}
\alias{run_qc}
\title{
    Run the quality control pipeline.
}
\description{
Reads the Affymetrix .CEL files, checks the NUSE score, annotates the probes using hgu133plus2.db and then calls other functions to perform additional QC.
}
\usage{run_qc(folder_location, final_dir)}
\arguments{
  \item{folder_location}{
  absolute path (folder name) where the .CEL files lie.
}
  \item{final_dir}{
  final path for output of the quality control files.
}
}
\details{

This function calls reads the Affymetrix CEL files and then check for the NUSE score, annotates the probes and the calls other function to perform complete QC.

}
\references{
Anguraj Sadanandam, et al. (2016), MOTassigner
}
\author{
Pawan Poudel paawanpoudel@gmail.com
}
\examples{
folder_location = "/Users/ppoudel/DATA/External/Collections/Microarrays/Affymetrics/Raw/MOT-5/PANCREAS_GSE151471/"
final_dir = "/Users/ppoudel/test/"
run_qc(folder_location, final_dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
