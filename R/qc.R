
# Input parameters --------------------------------------------------------

#Take these as an input of bash script?

param <- read.csv("parameters.csv")


#parameter inputs - 

# Trim primers ------------------------------------------------------------

# For this part of script - want to take input parameters from the param file
# Copy fasta files from sequencing folder to temp, do trimming

#Loop over runs - Maxlength set to remove any untrimmed reads
runs <- dir("data/", pattern="run_")

for (i in seq(along=runs)){
  path <- paste0("data/",runs[i]) # CHANGE ME to the directory containing your demultiplexed forward-read fastq files
  
  #Trim forward primers - Set a maxlength to remove all those that werent trimemd
  fastqFs <- sort(list.files(path, pattern="_R1", full.names = TRUE))
  fastqRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))
  
  bbtools_trim(install="bin/bbmap", fwd=fastqFs,rev=fastqRs, primers=c("CCTACGGGNGGCWGCAG","GACTACHVGGGTATCTAATCC"), copyundefined=TRUE, outpath="trimmed",ktrim="l", ordered=TRUE,mink=FALSE, hdist=2, overwrite=TRUE, samelength=TRUE, forcetrimright = 300, maxlength = 285)
}

#output diagnostics plots from trimming

#output trimmed files into temp folder


# Filter samples ----------------------------------------------------------

# input parameters from param file

# input files from trimmed temp folder
runs <- dir("data/", pattern="run_")
filtered_out <- vector("list", length=length(runs))

for (i in seq(along=runs)){
  path <- paste0("data/",runs[i],"/trimmed/") # CHANGE ME to the directory containing your demultiplexed forward-read fastq files
  filtpath <- file.path(path, "filtered") # Filtered forward files go into the path/filtered/ subdirectory
  dir.create(filtpath)
  fastqFs <- sort(list.files(path, pattern="R1_001.*"))
  fastqRs <- sort(list.files(path, pattern="R2_001.*"))
  
  if(length(fastqFs) != length(fastqRs)) stop(paste0("Forward and reverse files for ",runs[i]," do not match."))
  
  filtered_out[[i]] <- (filterAndTrim(fwd=file.path(path, fastqFs), filt=file.path(filtpath, fastqFs),
                                      rev=file.path(path, fastqRs), filt.rev=file.path(filtpath, fastqRs),
                                      maxEE=c(2,3),truncQ = 0,truncLen=c(280,200), maxN = 0,  rm.phix=TRUE, compress=TRUE, verbose=TRUE))
  
  # post filtering plot
  filtFs <- sort(list.files(filtpath, pattern="R1_001.*", full.names = TRUE))
  filtRs <- sort(list.files(filtpath, pattern="R2_001.*", full.names = TRUE))
  p1 <- plotQualityProfile(filtFs, aggregate = TRUE) + ggtitle(paste0(runs[i]," Forward Reads")) 
  p2 <- plotQualityProfile(filtRs, aggregate = TRUE) + ggtitle(paste0(runs[i]," Reverse Reads"))
  
  #output plots
  dir.create("output/figures/")
  pdf(paste0("output/figures/",runs[i],"_postfilt_quality.pdf"), width = 11, height = 8 , paper="a4r")
  plot(p1+p2)
  dev.off()
  
  #Get lengths post filter
  readcounts[[i]] <- cbind(width(readFastq(file.path(path, filtFs))), width(readFastq(file.path(path, filtRs))))
}


print(filtered_out)


# output sequences in temp folder

# remove trimmed 

# DADA2 error model -------------------------------------------------------


runs <- dir("data/", pattern="run_")
set.seed(100)

for (i in seq(along=runs)){
  path <- paste0("data/",runs[i],"/trimmed/" )# CHANGE ME to the directory containing your demultiplexed forward-read fastq files
  filtpath <- file.path(path, "filtered") # Filtered forward files go into the path/filtered/ subdirectory
  
  filtFs <- list.files(filtpath, pattern="R1_001.*", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="R2_001.*", full.names = TRUE)
  
  # Learn error rates from samples
  # nread tells the function how many reads to use in error learning, this can be increased for more accuracy at the expense of runtime
  
  errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
  
  ##Print error plots to see how well the algorithm modelled the errors in the different runs
  print(plotErrors(errF, nominalQ=TRUE)+ ggtitle(paste0(runs[i]," Forward Reads")))
  print(plotErrors(errR, nominalQ=TRUE)+ ggtitle(paste0(runs[i]," Reverse Reads")))
  
  #Error inference and merger of reads - Using pseudo pooling for increased sensitivity
  
  dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
  dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")
  
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
  
  # Construct sequence table
  
  seqtab <- makeSequenceTable(mergers)
  
  saveRDS(seqtab, paste0(path,"/seqtab.rds")) # CHANGE ME to where you want sequence table saved
}


# Remove Chimeras ---------------------------------------------------------


for (i in seq(along=runs)){
  path <- paste0("data/",runs[i],"/trimmed/" )
  seqs <- list.files(path, pattern="seqtab.rds", full.names = TRUE)
  
  assign(paste("st", i, sep = ""),readRDS(seqs))
  stlist <- append(stlist, paste("st", i, sep = ""), after=length(seqs))
}

st.all <- mergeSequenceTables(st1, st2, st3)

#Test collapsed
st.all <- collapseNoMismatch(st.all, minOverlap = 20, orderBy = "abundance",
                             vec = TRUE, verbose = TRUE)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=FALSE, verbose=TRUE)

#Check output of chimera removal
print(paste(sum(seqtab.nochim)/sum(st.all),"of the abundance remaining after chimera removal"))

#Check complexity
hist(seqComplexity(seqtab.nochim), 100)


#Look at seqlengths
plot(table(nchar(getSequences(seqtab.nochim))))

dir.create("output/rds/")
saveRDS(seqtab.nochim, "output/rds/seqtab_final_CUT.rds") # CHANGE ME to where you want sequence table saved



# Assign taxonomy ---------------------------------------------------------

seqtab.nochim <- readRDS("output/rds/seqtab_final_CUT.rds")

# Assign Kingdom:Genus taxonomy using RDP classifier
tax <- assignTaxonomy(seqtab.nochim, "reference/silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot=60, outputBootstraps=FALSE)
colnames(tax) <- c("Root", "Phylum", "Class", "Order", "Family", "Genus")

##add species to taxtable using exact matching
tax_plus <- addSpecies(tax, "reference/silva_species_assignment_v132.fa.gz", allowMultiple=TRUE)


##join genus and species name in species rank column - need to make this part of addspecies
sptrue <- !is.na(tax_plus[,7])
tax_plus[sptrue,7] <- paste(tax_plus[sptrue,6],tax_plus[sptrue,7], sep=" ")

tax_plus <- propagate_tax(tax_plus,from="Phylum")

#add Genus_SPP
#for(col in seq(7,ncol(tax_plus))) { 
# propagate <- is.na(tax_plus[,col]) & !is.na(tax_plus[,col-1])
#  tax_plus[propagate,col:ncol(tax_plus)] <-  "spp."
#}

#Check Output
taxa.print <- tax_plus # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Write taxonomy table to disk
saveRDS(tax_plus, "output/rds/tax_RDP_final_CUT.rds") 


# Create phylogeny --------------------------------------------------------

# if create phylo = true

seqtab.nochim <- readRDS("output/rds/seqtab_final_CUT.rds")

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# Write taxonomy table to disk
saveRDS(fitGTR, "output/rds/phytree.rds") 


# Clean up ----------------------------------------------------------------

#delete all temporary files

