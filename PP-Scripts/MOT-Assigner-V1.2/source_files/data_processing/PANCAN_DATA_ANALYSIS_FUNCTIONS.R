#/*
# PANCAN_DATA_ANALYSIS_FUNCTIONS.R
#  This file contains the functions that are used for the data analysis for the pancan 12 or orther datasets
#  Following functions are present in the modules
#
#   1: impute_missing_values
#   2: plot_batch_effects
#   3: plot_hist
#   4: plot_mean_median
#   5: sample_sample_coor
#   6: screen_genes
#   7: merge_file
#   8: create_sample_info
#   9: run_combat

# Created by Pawan Poudel on 06/02/2015.
#
#*/

# loading the libraries

cat('LOADING THE REQUIRED LIBRARIES\n')
library('sva')
require('impute')
require('lattice')
require('pheatmap')
require("matrixStats")
require("reshape")
require("ggplot2")



#source("/Users/ppoudel/Dropbox (ICR)/SysPreMed/Pawan/CANCER-PIPE/scripts/")

#-----------------------------------------------------------------------------------------------------
#FUNCTIONS STARTS HERE
#
#-------------------------------------------------------------------------------------------------------


theme_Publication <- function(base_size=16, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.text.x = element_text(angle=60), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.5, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#------------------------------
# SPLITTING THE DATA MATRIX ON THE BASIS OF ORGAN INFO
#----------------------------------

split_matrix<- function(data, sample_info){
    
    
    
    
    #get the organ information from the sample info file
    organs=as.character(unique(sample_info[,3]))
    
    for (l in 1:length(organs)){
        
        organ_ids=as.character(data.frame(sample_info[sample_info$Organs ==organs[l],])[,1]) # getting the description section
        
        m=match(colnames(data), organ_ids)
        w=which(!is.na(m))
        
        organ_gene_expression=subset(data, select=colnames(data)[w])
        
        organ_gene_expression=cbind(data[,1], organ_gene_expression)
        colnames(organ_gene_expression)[1]="gene"
        
        filename=paste0(organs[l], "/", "CCLE_CELL_LINES_", organs[l], ".txt")
        
        dir.create(organs[l])
        
        write.table(organ_gene_expression, filename, quote=FALSE, sep="\t", row.names=FALSE)
        
        cat("Done creating the file for organ:", organs[l])
        
        
    }

}


#
makePcaPlot<- function(x, sample_info, title = "", plot_file) {
  
  #arranging the data
  sample_info <- sample_info[order(sample_info[,2]), ]


  data.1 <- subset(x, select=as.character(sample_info[,1]))
  data.1 <- apply(data.1, 2, as.numeric)
  
  ind <- rowVars(data.1) < .Machine$double.eps
  dx <- data.1[!ind,]
  
  
  data <- t(apply(dx, 1, scale))
  mydata <- t(data)
  rownames(mydata) <- colnames(data.1)
  mydata.pca <- prcomp(mydata, scale.=TRUE)
  
  #percent <- round((((mydata.pca$sdev)^2 / sum(mydata.pca$sdev^2))*100)[1:2])
  #loadings <- mydata.pca$rotation
  
  #rownames(loadings) <- colnames(mydata)
  scores <- mydata.pca$x
  cex = 3
  
  group_name <- gsub("\\..+", "", sample_info[,2])
  
  df <- data.frame(rownames(scores), scores[,c(1,2)], group_name)
  colnames(df) <- c("samples", "PCA1", "PCA2", "Groups")
  
  p <- ggplot(df, aes(PCA1, PCA2, colour = Groups)) + geom_point(size = I(2)) +
    labs(x="PCA-1", y="PCA-2", color="Groups")+
    scale_fill_manual(values = as.factor(df$Groups))
  
  p <- p + labs(title = "PCA plot") +scale_colour_Publication()+ theme_Publication()
  
  
  ggsave(p, file=plot_file)
  return(df)
}


#--------------------------------
#GET THE SUBTYPE INFO FROM CORRELATION MAPPED FILE
#----------------------------------

produce_pie_chart_silhoutte <- function(barplot_file){
  
  data_table <- read.table(barplot_file, sep="\t", header=TRUE)
   
   data_table$membership <- paste0("Subtype-", data_table$membership)
  counts <- table(data_table$membership)
  counts_prop <- counts/sum(counts)*100
  
  df.1 <-round (counts_prop , digits = 2 )
  
  df.2 <- melt(df.1)
  
  colnames(df.2) <- c("subtypes", "value")
  g <- ggplot(df.2, aes(x=subtypes, y=value, fill=subtypes, label=value) )+scale_fill_brewer(palette = "Set2")
  g <- g + geom_bar(stat="identity") + geom_text() +labs(title="Percentage of samples in each subtypes")+ coord_polar()+scale_colour_Publication()+ theme_Publication()+scale_y_continuous(breaks=NULL)
  
  pie_file <- gsub(".txt", "positive_silhoutte_piechart.pdf",barplot_file)
 
  
  # now plotting a mosaic plot
  ggsave(g, file=pie_file)
 
}

#--------------------------------
#GET THE SUBTYPE INFO FROM CORRELATION MAPPED FILE
#----------------------------------

annotate_sample_subtype=function(correlation_classification){
    
    # getting rid for the first element
    X_anatomical_origin = gsub("^[[:alnum:]]+_", "", correlation_classification$samples)
    membership=gsub("\\.score", "", correlation_classification$labels)
    membership=gsub("X", "", membership)
    
    combined_data=data.frame(correlation_classification$samples,X_anatomical_origin, membership)
    colnames(combined_data)=c("sample","X_anatomical_origin", "membership")
    
    return(combined_data)
}



#-------------------------------------------------------------------
# EXPLORATORY DATA ANALYSIS PLOTS: HISTOGRAM
#-------------------------------------------------------------------------------
plot_hist=function(data, file){
    
    data=apply(data[,-1], 1, as.numeric)
    rowMed=apply(data, 1, median, na.rm=TRUE)
    data1=data-rowMed
    
    pdf(paste0(file,"_density_plot.pdf"))
    par(mfrow=c(1,2))
    hist(data,breaks=100, xlab=expression("gene"), main="Distribution of genes", freq=FALSE, cex.lab=1)
    lines(density(data), col="red",lwd=2)
    box()
    hist(data1,breaks=100, xlab=expression("genes"), main="Distribution of genes after median centering", freq=FALSE, cex.lab=1)
    lines(density(data1), col="red",lwd=2)
    box()
    dev.off()
}


# mapping the hiseq genes to the affymetrix genes
#---------------------------------------------------------------------------------
# This function maps the hiseq to the annotation present in affymetrix genes.
# Input:list containing gene names
# Output:mapped list of hgnc and affy ids




map_hiseq_to_affy=function(genes){
    
    ensembl=useMart("ensembl")
    ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
    
    mapped_mRNA = getBM(c("hgnc_symbol","affy_hg_u133_plus_2"), filters="hgnc_symbol",
    values=genes, mart=ensembl)
    i=which(mapped_mRNA[,2]== "")
    mapped_genes=mapped_mRNA[-i,]
    mapped_genes=mapped_genes[which(!is.na(mapped_genes[,2])),]
    return(unique(mapped_genes[,1]))
    
}


# Merging the files
#----------------------------------------------------------------------------------------
# This merge the list of files from multiple organs or multiple data sets
# Input: list of filenames that you would like to merge,
#   the fuction assumes that the files are tab seperated expression matrix with the gene names in the first colums
# Output:Tab seperated merged (on gene names) expression values with the gene in the first columns


merge_file = function(file_list){
    
    
    final_df=data.frame()
   
    genes_1=read.delim2(file_list[1], sep="\t", row.names=NULL, header=TRUE)
    index=which(!is.na(genes_1[,1]))
    
    final_df=data.frame(genes_1[index,])
    
    
    colnames(final_df)[1]="Genes"
    
    if(length(file_list)>=2){
      
      for (f in 2:length(file_list))
        {
        
        data_1=read.delim2(file_list[f], sep="\t", row.names=NULL, header=TRUE)
        colnames(data_1)[1]="Genes"
        
        final_df=merge(final_df, data_1, by.x="Genes", by.y="Genes")
        
      }
      
      return(final_df)
      
    }else{

      return(final_df)
    
    }
}


cat("Here")

#combat
#---------------------------------------------------------------------------------

run_combat=function(df, sample_f){
    
    colnames(sample_f)=c("samples", "batch")
    batch=sample_f$batch
    mod = model.matrix(~1, data=sample_f)
    df_1=apply(df[,-1], 2, as.numeric)
    data_b_c_1=ComBat(dat=df_1, batch=batch, mod=mod)
    
    #median centering the data for normalisation
    rowMed=  apply(data_b_c_1,1,median)
    
    df_f=data_b_c_1 - rowMed
    dbc=data.frame(df[,1],df_f )
    colnames(dbc)[1]="Genes"
    return(dbc)
        
}


# Silhoutte width plot
#---------------------------------------------------------------------------------

require(RColorBrewer)

plot_silhoutte <- function(fil) {
  
  sil_plot <- gsub(".txt", ".pdf", fil)
  sil_txt <- gsub(".txt", ".silhoutte_summary.txt", fil)
  
  df <- read.delim2(fil, sep="\t", header = TRUE)
  
  
  x <- ncol(df)
  colnames(df)[x] <- "sil_width"
  df <- df[order(df$membership, df$sil_width), ]
  df$name <- factor(rownames(df), levels = rownames(df))
  
  
  df$membership <- as.factor(df$membership)
  df$sil_width <- as.numeric(df$sil_width)
  
  mapping <- aes_string(x = "name", y = "sil_width", fill="membership")
  ave <- tapply(df$sil_width, df$membership, mean)
  n <- tapply(df$membership, df$membership, length)
  sil.sum <- data.frame(cluster = names(ave), size = n, ave.sil.width = round(ave, 2))
  print(sil.sum)
  
  #colours <- c("#2121D9", "#9999FF", "#D92121", "#21D921", "#FFFF4D", "#FF9326")
  
  p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
    labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", "\n Average silhouette width: ", round(mean(df$sil_width),  2))) + ggplot2::ylim(c(NA, 1))
  
  
  q <-   p + scale_fill_brewer(palette = "Set2")+ scale_colour_Publication() + theme_Publication()
  
  
  ggsave(q, file=sil_plot)
  
  write.table(sil.sum, sil_txt, sep = "\t")
  
}


#performing the silhoutte
#---------------------------------------------------------------------------------

do_silhoutte=function(data_file, consensus_file,silhoutte_dir, info, k )
{
    #naming the result files
    result_file=paste0(silhoutte_dir, "/", info,"_", k,"_silhoutte_result.txt")
    silhoutte_selected_data_file=paste0(silhoutte_dir, "/",  info,"_", k,"_positive_silhoutte_samples.txt")
   
    
    nmf_data=read.delim(data_file, sep="\t", header=TRUE)
    consensus_data=read.delim(consensus_file, sep="\t", header=TRUE)
    col=dim(nmf_data)[2]
    row=dim(nmf_data)[1]
    
    #ordering the consensus data
    consensus_data=consensus_data[order(consensus_data[,2]),]
    
    consensus_data_subset=data.frame(samples=consensus_data[,1], lables=consensus_data[,2]) 
    num=length(unique(consensus_data[,2]))
    
    
    m=match(consensus_data_subset$samples, colnames(nmf_data))
    w=which(!is.na(m))
    
    subtype=as.factor(consensus_data_subset$lables[w])  
    nmf_data_subset=subset(nmf_data, select = as.character(consensus_data_subset$samples[w]))
    my_genes=nmf_data[,1]
    
    nmf_data_subset=apply(nmf_data_subset, 2, as.numeric)
    rownames(nmf_data_subset)=my_genes
    
    si=silhouette(as.numeric(subtype), dist(t(nmf_data_subset), "euclidean"))
    
    si_1=matrix(as.numeric(si), nrow=(col-1), ncol=3)
    si_1=cbind(consensus_data, si_1)
    
    # now writing the silhoutte results
    write.table(si_1, result_file, sep="\t",row.names=FALSE, quote = FALSE)
    
    # now selecting the samples with the positive silhoutte score    
    w=which(si_1[,5] >=0)
    #silhoutte selected samples
    silhoutte_selected_samples=si_1[w,1]
    
    silhoutte_selected_data=subset(nmf_data, select=silhoutte_selected_samples)
    silhoutte_selected_data_to_write=cbind(nmf_data[,1],silhoutte_selected_data)
    
    write.table(silhoutte_selected_data_to_write, silhoutte_selected_data_file, sep="\t", quote = FALSE, row.names=FALSE)
    cat(silhoutte_selected_data_file, "====== file")
    
    return(result_file)

}

