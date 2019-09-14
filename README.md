# piperline
Metabarcoding pipeline template

This workflow is divided into 3 sections:

01_QC - This takes an input of run_locations in the sequencing group folder - an option for Re-Dulmultiplex also. What this does is:
for loop length (run_locations given)
copy the data to a temp folder, run a number of QC steps and output
1- Run level illumina style QC
2- Read level fastqc style
3- Sample read summary to choose DADA2 filtering parameters. - this could be a seperate shiny app to display the filterign options similar to dada2, but also interactive visualisation of expected error filtering
This then outputs all these in a summary folder or HTML, and deletes all temp files. You then copy the outputs back to your local PC for viewing in a web browser.

02_Analyse - This script also takes an input of run_locations in the sequencing group folder. And takes a number of parameters including:
DADA2 Parameters:
    --trimLeft              How far to trim on the left (default = 0)
    --maxN                  (default = 0)
    --maxEE                 (default = Inf)
    --truncLenF             (default = 0)
    --truncLenR             (default = 0)
    --truncQ                (default = 2)
Taxonomic assignment parameters: MEthod, and database and other necessary parameters. See https://github.com/jgolob/maliampi for a great options of this. These options should be decided by looking at the QC outputs of the first script
The output of this is the final seqtab, and qc tracker file. This then gets copies back to your local PC for visualsiation in the final shiny app, and further analysis in R/phyloseq

03_visualise - This wraps a shiny package similar to ranacapa https://gauravsk.shinyapps.io/ranacapa/ 
This takes an input of: QC files from the first script, the Seqtab and QC tracker from the second script, and a sample metadata file and displays some basic visualisations. This is the "first pass" of your analysis, and anything more complex or for publication can be further edited in R

The main pages on this shiny app will be:
QC
- Illumina run level QC outputs - similar to what you see on basespace - this is an important qualtiy control measure
- Read level QC - fastqc output, reads per sample, index switching measurements
- Results of read tracking through pipeline - Reads lost in filtering, mergign etc, number of ASV's out
- Results of Taxonomic assignment - method used, database used & Version, number of ASV's assigned to each taxonomic rank visualised.

Results:
- Table visualsiation of Entire seqtab + a column for assignment confidence, Enable filtering similar to excel
- Interactive filtering? - View replicates, View outliers, View prevalence of features
-  Explore sequencing depth - taxon rarefication species accumulation chart,
- Krona chart of results
- Heatmap of results - allow interactive agglomeration to each taxonomic rank and re-plot (can this be done in shinyphyloseq)
- Barchart of results - allow viewing of raw reads and Relative abundance transformed
- For both of the above - Could have a drag and drop table underneath to arrange samples - also a subset to certain samples option, and perhaps have an export plot button to render a ggplot pdf and download
- 
