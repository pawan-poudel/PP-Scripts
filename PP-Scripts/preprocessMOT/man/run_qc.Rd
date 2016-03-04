\name{run_qc}
\alias{run_qc}
\title{
    Run the quality control 
}
\description{
Reads the affymetrix CEL files, and then checks for the NUSE score, annotates the probes using hgu133plus2.db and then calls other function to perform additional QC.
}
\usage{run_qc(folder_location, final_dir)}
\arguments{
  \item{folder_location}{
  Absolute path (just input the folder name) where the CEL files lies
}
  \item{final_dir}{

  final path that you would like to put the QCed files
}
}
\details{

This function calls reads the affymetrix CEL files and then check for the NUSE score, annotates the probes and the calls other function to perform complete QC.

}
\references{
Anguraj Sadanandam, et all (2016), MOTassigner
}
\author{
Pawan Poudel
}
\examples{
folder_location="/Users/ppoudel/DATA/External/Collections/Microarrays/Affymetrics/Raw/MOT-5/PANCREAS_GSE151471/"
final_dir="/Users/ppoudel/test/"
run_qc(folder_location,final_dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
