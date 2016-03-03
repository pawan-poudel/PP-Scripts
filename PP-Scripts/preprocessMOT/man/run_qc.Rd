\name{run_qc}
\alias{run_qc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Run the quality control 
}
\description{
    Reads the affymetrix CEL files, and then checks for the NUSE score, annotates the probes using hgu133plus2.db and then calls other function to perform additional QC.
}
\usage{
run_qc(folder_location, final_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder_location}{
  absolute path (just input the folder name) where the CEL files lies
}
  \item{final_dir}{

  ffinal path that you would like to put the QCed files
}
}
\details{

This function calls reads the affymetrix CEL files and then check for the NUSE score, annotates the probes and the calls other function to perform complete QC.

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
function (folder_location, final_dir) 
{
    setwd(folder_location)
  #get the gse identifier
  f=basename(getwd())
  fd <- paste0(final_dir,"/",f)
  dir.create(fd, showWarnings = FALSE)
  
  fd1 <- paste0(fd,"/QC/")
  dir.create(fd1, showWarnings=FALSE)
  
  # reading the affymetrics object
  data.raw=ReadAffy()
  
  PLMset<- fitPLM(data.raw,normalize=FALSE)
  nuse_score=NUSE(PLMset, type="stats")
  
  i=which(nuse_score[1,]< 1.05 & nuse_score[1,]>.95)
  samples_to_keep=names(i)
  
  pf <- paste0(fd1, Sys.Date(),"_",f, "_", "NUSE.pdf")
  nuse_f <- paste0(fd1, Sys.Date(),"_",f, "_", "NUSE.txt")
  
  # plotting the nuse socre
  pdf(pf)
  NUSE(PLMset)
  dev.off()
  # writing the nuse information to the file
  write.table(t(nuse_score),nuse_f, sep="\t")
  
  # cleaning the data and performing the RMA
  data.raw.clean=data.raw[,samples_to_keep]
  data_rma=rma(data.raw.clean) #the expression values are in log2 scale
  data_raw_exp=exprs(data.raw.clean)
  
  data_exp=exprs(data_rma)
  
  probes_n=rownames(data_exp)
  probes_data_to_write=data.frame(probes_n, data_exp)
  
  sample_info=data.frame(samples_to_keep, paste0("Organ-",f,rep(1, length(samples_to_keep))))
  
  inf=paste0(f,"Batch effect before normalisaton and removing NUSE more than +/- 1.05 " )
  
  batch_plot_f=paste0(fd1, "/", Sys.Date(),"_",f, "_", "RMA_NORMALISED_GENES_EXPRESSION_HGU_PLUS2_Arrays-Batch-effect-after-normalisation.pdf")
  
  makePcaPlot(x = data_exp, sample_info=sample_info, title = inf, plot_file = batch_plot_f)
  
  
  x=unlist(as.list(hgu133plus2SYMBOL))
  match_values=match(rownames(data_exp),names(x))
  matching_index=which(!is.na(match_values))
  
  xx=x
  names(xx)=NULL
  rownames(data_exp)=xx[match_values[matching_index]]
  
  
  fp=paste0(fd1, "/", Sys.Date(),"_",f, "_", "RMA_NORMALISED_GENE_EXPRESSION_HGU_PLUS2_Arrays_Probes.txt")
  write.table(data_exp, fp, sep="\t", row.names=FALSE, quote = FALSE)
  
  
  data_exp_final <- cbind(rownames(data_exp), data_exp)
  colnames(data_exp_final)[1] <- "Genes"
  rownames(data_exp_final) <- NULL
  
  # writing the RMA normalised data  
  ff=paste0(fd1, "/", Sys.Date(),"_",f, "_", "RMA_NORMALISED_GENE_EXPRESSION_HGU_PLUS2_Arrays.txt")
  write.table(data_exp_final, ff, sep="\t", row.names=FALSE, quote = FALSE)
  
  # performing QC part 2
  check_dist_zscore_coor(ff)
  
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
