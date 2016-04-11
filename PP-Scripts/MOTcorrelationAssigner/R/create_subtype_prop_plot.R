
create_subtype_prop_plot <-
function(file, output_dir){
    
    df <- read.delim( file, sep="\t", header = TRUE)
    
    colnames(df) <- c("samples", "organs","membership")
    
    df[,3] <- gsub("\\.score", "",df[,3])
    df[,3] <- gsub("\\.", "\\-",df[,3])
    
    
    # getting the color information for each subtypes
    colours <- c( "#fdbb84", "#e34a33", "#67000d", "#08519c", "#addd8e", "#df65b0")
    subtype <- c("Basal", "Classical", "Inflammatory", "Stem-like", "Stem-PPAR", "TA")
    
    
    my_col <- data.frame(subtype, colours)
    
    # matching the color information with the 6 subtypes and assigning colours accordingly
    
    m1 <- match( my_col[,1], df[,3])
    w1 <- which(!is.na(m1))
    
    sub.col <- my_col[w1,2]
    
    # creating the bar plot
    p <- ggplot(df, aes(x= organs, fill=membership))+ geom_bar(position="fill",aes(fill = membership))
    p <- p + scale_fill_manual(values = as.character(sub.col))+ theme_Publication()+ labs(title="Proportions of samples in each subtype",x="Organs", y = "Proportions")
    
    
    file_1 <- gsub(".+\\/", "", file)
    prop_pdf_file <- paste0(output_dir, "/", file_1)
    prop_pdf_file <- gsub(".txt",".pdf",prop_pdf_file )
    
    ggsave(p, filename=prop_pdf_file)
    
    
}

