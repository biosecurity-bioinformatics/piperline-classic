# Load necessary libraries. use install.packages("packagename") if they need to be installed first.
library(tidyverse)
library(readxl)

# Read in sheet 1 and select just the columns we need
sheet1 <- read_excel("summary_filtered_Jan_Feb_2023.xlsx")%>%
  dplyr::select(sequence, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species, starts_with("L6LGD") )

# Read in sheet 2 and select just the columns we need
sheet2 <- read_excel("summary_filtered_Feb_March_2023.xlsx")%>% 
  dplyr::select(sequence, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species, starts_with("L69TG") )

merged <- sheet1 %>% 
  bind_rows(sheet2) %>% #Join both sheets together
  pivot_longer( # Pivot to long format
    cols = -c("sequence", "Root", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    names_to = "sample",
    values_to = "reads"
    ) %>% 
  filter(!is.na(reads)) %>% # Filter out any entries with no reads
  pivot_wider( # Pivot back to wide format
    names_from = "sample",
    values_from="reads",
    values_fill = list(reads=0) # Fill any missing values with zero
    ) 

# Check all sequences are present in merged samples
starting_seqs <- length(unique(c(sheet1$sequence, sheet2$sequence)))
merged_seqs <- length(unique(merged$sequence))

all.equal(starting_seqs, merged_seqs) # TRUE - all sequences are present in merged output

# Check all columns (samples) are present in output
starting_colnames <- unique(c(colnames(sheet1), colnames(sheet2)))
merged_colnames <- colnames(merged)
all.equal(starting_colnames, merged_colnames) # TRUE - all columns (samples) are present in merged output

# Write out the merged datasheet
write_csv(merged, "summary_filtered_Jan_March_2023.csv")

