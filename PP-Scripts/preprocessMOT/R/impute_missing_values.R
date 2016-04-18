impute_missing_values <-
function(pancan12_genomeMatrix){
  
  pancan12_genomeMatrix_1=apply(pancan12_genomeMatrix[,-1], 2, as.numeric)
  genes=pancan12_genomeMatrix[,1]
  pancan12_genomeMatrix.imputed=impute.knn(as.matrix(pancan12_genomeMatrix_1))
  rownames(pancan12_genomeMatrix.imputed$data)=NULL
  data_to_write=cbind(genes,pancan12_genomeMatrix.imputed$data)
  colnames(data_to_write)[1]="Genes"

  #write.table(data_to_write, file, sep='\t')
  return(data_to_write)
  
  
}
