
# Input parameters --------------------------------------------------------





# Illumina run quality control --------------------------------------------




# Sequence data quality control -------------------------------------------


library(ShortRead)

runs <- dir("data/", pattern="run_")

for (i in seq(along=runs)){
  path <- paste0("data/",runs[i]) # CHANGE ME to the directory containing your demultiplexed forward-read fastq files
  
  #path <- paste0("data/",runs[i])
  
  #Plot number of reads
  dat <- as.data.frame(countLines(dirPath=path, pattern=".fastq")) %>%
    rownames_to_column()  %>%
    `colnames<-`(c("Sample", "Reads")) %>%
    filter(str_detect(Sample,"R1"))
  
  #Plot pooling
  
  gg.pooling <- ggplot(data=dat, aes(x=Sample,y=Reads),stat="identity") + 
    geom_bar(aes(fill=Reads),stat="identity")  + 
    scale_fill_viridis(name = "Reads", begin=0.1) + 
    theme(axis.text.x = element_text(angle=90, hjust=1), plot.title=element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5))+ 
    geom_hline(aes(yintercept = mean(Reads)))  +
    xlab("sample name")+
    ylab("Number of reads") + 
    labs(title= paste0("Pooling for : ", runs[i]), subtitle = paste0("Total Reads: ", sum(dat$Reads), " Average reads: ",  sprintf("%.0f",mean(dat$Reads))," Standard deviation: ", sprintf("%.0f",sd(dat$Reads)))) +
    coord_flip()
  
  plot(gg.pooling)
}


# Read quality and lengths ------------------------------------------------

runs <- dir("data/", pattern="run_")
readcounts <- vector("list", length=length(runs))

for (i in seq(along=runs)){
  path <- paste0("data/",runs[i],"/trimmed" )# CHANGE ME to the directory containing your demultiplexed forward-read fastq files
  
  filtFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
  filtRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))
  p1 <- plotQualityProfile(filtFs, aggregate = TRUE) + ggtitle(paste0(runs[i]," Forward Reads")) 
  p2 <- plotQualityProfile(filtRs, aggregate = TRUE) + ggtitle(paste0(runs[i]," Reverse Reads"))
  
  #output plots
  dir.create("output/figures/")
  pdf(paste0("output/figures/",runs[i],"_prefilt_quality.pdf"), width = 11, height = 8 , paper="a4r")
  plot(p1+p2)
  dev.off()
  
  #Get lengths
  readcounts[[i]] <- cbind(width(readFastq(file.path(path, fastqFs))), width(readFastq(file.path(path, fastqRs))))
  
}
