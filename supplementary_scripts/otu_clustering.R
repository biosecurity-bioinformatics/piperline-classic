# Load packages
library(tidyverse)
library(taxreturn)
library(seqateurs)
library(Biostrings)

# Read in dataset 
dataset <- readxl::read_excel("clustering_example.xlsx") #Set to your data file

# Extract the sequence column for clustering and turn it into a DNAStringSet object
seqs <- DNAStringSet(dataset$sequence)

# cluster the sequences at 97% similarity
clustered <- cluster_otus(seqs, similarity = 0.97, cores = 1)

# Join the clustering results to the original dataset - then move the columns to be after sequences
dataset2 <- dataset %>%
  left_join(clustered, by="sequence") %>%
  relocate(cluster, cluster_size, .after = "sequence")

# Write out the updated dataset as a csv
write_csv(dataset2, "clustering_result.csv")  #Set to desired output
