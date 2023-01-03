library(targets)
library(tarchetypes)
source("R/functions.R")
options(tidyverse.quiet = TRUE)

# Load packages -----------------------------------------------------------
tar_option_set(packages = c(
  "devtools",
  "ggplot2",
  "gridExtra",
  "data.table",
  "tidyverse", 
  "stringdist",
  "patchwork",
  "vegan",
  "seqinr",
  "patchwork",
  "stringi",
  "phangorn",
  "magrittr",
  "galah",
  "phyloseq",
  "DECIPHER",
  "Biostrings",
  "ShortRead",
  "ggtree",
  "savR",
  "dada2",
  "ngsReports",
  "taxreturn",
  "seqateurs",
  "speedyseq"
), workspace_on_error = TRUE)


# Targets pipeline
list(
  # Track files
  tar_file(samdf_file, "sample_data/Sample_info.csv"),
  tar_file(params_file, "sample_data/loci_params.csv"),
  
  # Load sample tracking files
  tar_target(samdf, read_csv(samdf_file)),
  tar_target(params,read_csv(params_file)),
  
  # Create directories
  tar_target(create_dirs,  step_validate_folders(getwd())),
  
  # Look for sequencing reads
  tar_files(
    fastq_paths,
    purrr::map(list.dirs("data", recursive=FALSE),
                          list.files, pattern="_R1_", full.names = TRUE) %>%
      unlist() 
  ),

  # Check sequencing reads match sample sheet & Create temporary samdf file
  tar_file(
    make_temp_samdf1,
    {
      step_check_files(samdf, fastq_paths) %>%
      write_csv("temp/temp_samdf1.csv")
      return("temp/temp_samdf1.csv")
      }
  ),
  
  # Track temporary samdf
  tar_target(temp_samdf1, read_csv(make_temp_samdf1)),

  tar_group_by(temp_samdf1_grouped, temp_samdf1, fcid),
  
## Per-Sample QC -----------------------------------------------------------
#  tar_target(sample_qc, 
#             temp_samdf1 %>% mutate(sample_qc = purrr::pmap(list(sample_id, fcid),
#              .f = ~step_sample_qc(sample_id = ..1, fcid=..2, multithread=FALSE, quiet=TRUE)
#              )),
#            pattern = map(temp_samdf1), iteration = "vector"),
#  
#  tar_group_by(sample_qc_grouped, sample_qc, fcid),
#  
#  tar_target(multiqc, 
#             sample_qc_grouped %>%
#               group_by(fcid) %>%
#               nest() %>%
#               mutate(multi_qc = purrr::map(fcid, step_multiqc, quiet=FALSE)),
#             pattern = map(sample_qc_grouped), iteration = "vector"),
#  

# Sequencing run QC -------------------------------------------------------
  tar_target(seq_qc,
             temp_samdf1_grouped %>%
                group_by(fcid) %>%
                nest() %>%
                mutate(seq_qc = purrr::map(fcid, step_seq_qc)),
             pattern = map(temp_samdf1_grouped), iteration = "vector"),
  
  tar_target(switching_qc,
             temp_samdf1_grouped %>%
               group_by(fcid) %>%
               nest() %>%
               mutate(switching_qc = purrr::map(fcid, step_switching_calc, multithread=FALSE, quiet=TRUE)),
             pattern = map(temp_samdf1_grouped), iteration = "vector"),
             

# Demultiplex and trim primers --------------------------------------------
  tar_target(primer_trim,
            {
            process <- temp_samdf1 %>%
             mutate(primer_trim = purrr::pmap(dplyr::select(., sample_id, for_primer_seq, rev_primer_seq, pcr_primers, fcid),
                                     .f = ~step_primer_trim(sample_id = ..1, for_primer_seq=..2, rev_primer_seq=..3, pcr_primers = ..4,
                                                             input_dir = paste0("data/",..5), output_dir =  paste0("data/",..5,"/01_trimmed"),
                                                             qc_dir=paste0("output/logs/",..5), quiet = FALSE))) 
            #write_csv(process, "temp/primer_trim_res.csv") # This is overwriting - need a merge function
            # Return filepaths for tracking
            outF <- process$primer_trim %>%  
              bind_rows()%>%
              pull(fwd_out) 
            outR <- process$primer_trim %>%
              bind_rows()%>%
              pull(rev_out) 
            # Check for empty files
            outF <- outF[file.exists(outF)]
            outR <- outR[file.exists(outR)]
            # Return list of completed paths
            return(c(outF,outR))
            },
    pattern = map(temp_samdf1), format="file", iteration = "vector"),

 # tar_file(primer_trim_results, "temp/primer_trim_res.csv"),
  # Check sequencing reads match sample sheet & Create temporary samdf file
  tar_file(
    make_temp_samdf2,
    {
      temp_samdf1 %>%
      select(-where(is.list)) %>%
      step_demux_samdf() %>%
      #left_join(params %>% dplyr::select(pcr_primers, ref_db)) %>%
      #step_add_params(params) %>%
      step_check_files(primer_trim) %>%
      write_csv("temp/temp_samdf2.csv")
      return("temp/temp_samdf2.csv")
    }
  ),

  # Track temporary samdf
  tar_target(temp_samdf2, read_csv(make_temp_samdf2)),
  
# Filter reads ------------------------------------------------------------
  tar_target(read_filter,
             {
             process <- temp_samdf2 %>%
             #  dplyr::slice(1:5)%>%
             mutate(read_filter = purrr::pmap(dplyr::select(., sample_id, fcid),
                   .f = ~step_read_filter(sample_id = ..1,
                   input_dir = paste0("data/",..2,"/01_trimmed/"), output_dir = paste0("data/",..2,"/02_filtered"),
                   maxEE = 1, truncLen = 150, rm.lowcomplex = 0,
                   quiet = FALSE)))
             #write_csv(process, "temp/filter_res.csv")
             # Return filepaths for tracking
             outF <- process$read_filter %>%  
               bind_rows()%>%
               pull(fwd_out) 
             outR <- process$read_filter %>%
               bind_rows()%>%
               pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed paths
             return(c(outF,outR))
             },                  
             pattern = map(temp_samdf2), format="file", iteration = "vector"),

  tar_file(
    make_temp_samdf3,
    {
      temp_samdf2 %>%
        select(-where(is.list)) %>%
        step_check_files(read_filter) %>%
        write_csv("temp/temp_samdf3.csv")
      return("temp/temp_samdf3.csv")
    }
  ),
  
  # Track temporary samdf
  tar_target(temp_samdf3, read_csv(make_temp_samdf3)),

## Pre-filtering quality plots ---------------------------------------------

# Sample a random set of samples for read quality plotting
tar_target(read_samples,{
  group_sizes <- temp_samdf3 %>%
    group_by(fcid) %>%
    group_size()
  if(all(group_sizes > 10)){
    n_samples <- 10
  } else {
    n_samples = min(group_sizes)
  }
  out <- temp_samdf3 %>%
    group_by(fcid) %>%
    slice_sample(n=n_samples) 
}),

tar_target(prefilt_qualplots,
          read_samples %>%
           mutate(prefilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/01_trimmed/"), truncLen=NULL, quiet = FALSE, n = 10000)
            )),
            pattern = map(read_samples), iteration = "vector"),

## Write out prefilt qualplots
tar_target(write_prefilt_qualplots, {
  prefilt_qualplots %>% 
    group_by(fcid) %>%
    nest() %>%
    purrr::pwalk(list(fcid, data),
      .f = ~{
      pdf(file=paste0("output/logs/",..1,"/prefilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
      print(..2$prefilt_qualplots)
      try(dev.off(), silent=TRUE)
      })
  out <- paste0("output/logs/",unique(prefilt_qualplots$fcid),"/prefilt_qualplots.pdf")
  return(out)
  }, format="file", iteration = "vector"),

## Post-filtering quality plots --------------------------------------------
tar_target(postfilt_qualplots,
           read_samples %>%
            mutate(postfilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/02_filtered/"), truncLen=NULL, quiet = FALSE, n = 10000)
           )),
            pattern = map(read_samples), iteration = "vector"),
 
# Write out postfilt qualplots
tar_target(write_postfilt_qualplots, {
  postfilt_qualplots %>% 
    group_by(fcid) %>%
    nest() %>%
    purrr::pwalk(list(fcid, data),
                 .f = ~{
                   pdf(file=paste0("output/logs/",..1,"/postfilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
                   print(..2$postfilt_qualplots)
                   try(dev.off(), silent=TRUE)
                 })
  out <- paste0("output/logs/",unique(postfilt_qualplots$fcid),"/postfilt_qualplots.pdf")
  return(out)
}, format="file", iteration = "vector"),

# Infer sequence variants with DADA2 --------------------------------------

  # Create new samdf from output of previous step? Edit check function to take the 
  tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid),
 
  # How to make it just redo one of the dadas if only one runs filtered file changed changed?
  tar_target(dada,{
             process <- temp_samdf3_grouped %>%
             group_by(fcid) %>%
             nest() %>%
             mutate(dada2 = purrr::map(fcid,
                                        .f = ~step_dada2(fcid = .x,
                                                         input_dir = paste0("data/",.x,"/02_filtered"),
                                                         output = paste0("output/rds/",.x,"_seqtab.rds"),
                                                         qc_dir = paste0("output/logs/",.x),
                                                         quiet = FALSE)
             )) %>%
               unnest(data)
             #write_csv(process, "temp/filter_res.csv")
             out <- paste0("output/rds/",unique(temp_samdf3_grouped$fcid),"_seqtab.rds")
             return(out)
             },
             pattern = map(temp_samdf3_grouped), format="file", iteration = "vector"),
  
##  Merge infered variants from each run and subset to target loci ---------
 tar_target(subset_seqtab, {
            process <- temp_samdf3 %>%
             ungroup() %>%
             group_by(pcr_primers) %>%
             nest() %>%
             mutate(subset_seqtab = purrr::map(pcr_primers, 
                   .f = ~{
                   #seqtabs <- list.files("output/rds/", pattern="seqtab.rds", full.names = TRUE)
                   if(length(dada) > 1){
                   st.all <- mergeSequenceTables(tables=dada)
                   } else if(length(dada) == 1) {
                   st.all <- readRDS(dada)
                   }
                   st.all <- st.all[str_detect(rownames(st.all), .x),]
                   st.all <- st.all[,colSums(st.all) > 0]
                   saveRDS(st.all, paste0("output/rds/",.x,"_seqtab.rds"))
                   return(st.all)
              }))
            out <- paste0("output/rds/",unique(process$pcr_primers),"_seqtab.rds")
            return(out)
            }, format="file", iteration = "vector"),
 

##  Filter ASV's per locus -------------------------------------------------
  tar_target(filtered_seqtab, {
           process <- temp_samdf3 %>%
            dplyr::select(-one_of("exp_length", "phmm", "coding", "genetic_code"))%>%
            left_join(params %>% dplyr::select(pcr_primers, exp_length, phmm, coding, genetic_code)) %>%
            group_by(pcr_primers, exp_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq) %>%
            nest() %>%
            mutate(subset_seqtab = purrr::map(pcr_primers, ~{
                 readRDS(subset_seqtab[str_detect(subset_seqtab, .x)])
            })) %>%
            ungroup()%>%
            mutate(filtered_seqtab = purrr::pmap(dplyr::select(.,pcr_primers, subset_seqtab, exp_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq),
                .f = ~step_filter_asvs(
                seqtab = ..2,
                output = paste0("output/rds/",..1,"_seqtab.cleaned.rds"),
                qc_dir = "output/logs/",
                min_length = ..3-10,
                max_length = ..3+10,
                phmm = ..4,
                check_frame = ..5,
                genetic_code = ..6,
                primers = c(..7, ..8),
                multithread = FALSE, 
                quiet = FALSE)
          )) %>%
          unnest(data)
         out <- paste0("output/rds/",unique(process$pcr_primers),"_seqtab.cleaned.rds")
         return(out)
  }, format="file", iteration = "vector"),
  

## Merge all loci into final seqtab ----------------------------------------
  tar_target(merged_seqtab, {
    process <- temp_samdf3 %>%
      nest(data=everything()) %>%
      mutate(final_seqtab = purrr::map(data, ~{
        #seqtabs <- fs::dir_ls("output/rds/", glob="*.cleaned.rds")
        seqtabs <- filtered_seqtab
        seqtabs <- seqtabs[seqtabs %>%
                             purrr::map_lgl(function(y){
                               any(str_detect(y, unique(.x$pcr_primers)))
                             })]
        if(length(seqtabs) > 1){
          st.all <- mergeSequenceTables(tables=seqtabs)
        } else if(length(seqtabs) == 1) {
          st.all <- readRDS(seqtabs)
        }
        saveRDS(st.all, "output/rds/seqtab_final.rds")
        return(TRUE)
      })) %>%
      unnest(data)
    out <- paste0("output/rds/seqtab_final.rds")
    return(out)
   }, format="file", iteration = "vector"),
  
  # Create grouped seqtab
  #tar_group_by(merged_seqtab_grouped, merged_seqtab, target_gene),

# Assign taxonomy ---------------------------------------------------------
 # # Track the taxonomy files - could join these similar to the way im joining the eqtab again
 # tar_file(blast_db,
 #          params  %>%
 #              pull(blast_db) %>%
 #              unique(),
 # ),
 # tar_file(ref_db,
 #          params  %>%
 #            pull(ref_db) %>%
 #            unique(),
 # ),

# Shouldnt need to add the parames earlier, should just be able to join here

## IDTAXA -------------------------------------------------------------------
 tar_target(tax_idtaxa,{ 
   process <- temp_samdf3 %>%
     dplyr::select(-one_of("target_gene", "ref_db"))%>%
     left_join(params %>% dplyr::select(pcr_primers, target_gene, ref_db)) %>%
     tidyr::separate_rows(ref_db, sep=";") %>%
     group_by(target_gene, pcr_primers, ref_db) %>%
     nest()  %>% 
     mutate(filtered_seqtab = purrr::map(pcr_primers, ~{
       readRDS(filtered_seqtab[str_detect(filtered_seqtab, .x)])
     }))  %>%
     mutate(idtaxa = purrr::pmap(list(target_gene, filtered_seqtab, ref_db),
                                 .f = ~step_assign_taxonomy(
                                   seqtab = ..2,
                                   database = ..3,
                                   ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                                   output = paste0("output/rds/",basename(..3) %>% str_remove("\\..*$"),"_taxtab.rds"),
                                   qc_dir = "output/logs/",
                                   threshold = 60,
                                   multithread = FALSE, 
                                   quiet = FALSE)
     ))%>%
     unnest(data)
   out <- paste0("output/rds/",unique(basename(process$ref_db)  %>% str_remove("\\..*$")),"_taxtab.rds")
   return(out)
   }, format="file", iteration = "vector"),

## BLAST -------------------------------------------------------------------
 tar_target(tax_blast,
            {
           process <- temp_samdf3 %>%
             dplyr::select(-one_of("target_gene", "blast_db"))%>%
             left_join(params %>% dplyr::select(pcr_primers, target_gene, blast_db)) %>%
             tidyr::separate_rows(blast_db, sep=";") %>%
             group_by(target_gene, pcr_primers, blast_db) %>%
             nest() %>% 
             mutate(filtered_seqtab = purrr::map(pcr_primers, ~{
               readRDS(filtered_seqtab[str_detect(filtered_seqtab, .x)])
             }))  %>%
             mutate(blast = purrr::pmap(list(target_gene, filtered_seqtab, blast_db),
                                        .f = ~step_blast_tophit(
                                          seqtab = ..2,
                                          database = ..3,
                                          ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                          output = paste0("output/rds/",basename(..3) %>% str_remove("\\..*$"),"_blast.rds"),
                                          qc_dir = "output/logs/",
                                          identity = 97,  
                                          coverage=95,
                                          evalue=1e06,
                                          max_target_seqs=5,
                                          max_hsp=5, 
                                          multithread = FALSE, 
                                          quiet = FALSE)
           ))
          out <- paste0("output/rds/",unique(basename(process$blast_db)  %>% str_remove("\\..*$")),"_blast.rds")
          return(out)
        }, format="file", iteration = "vector"),
 
## Aggregate taxonomic assignment methods-----------------------------------------------
  tar_target(joint_tax,
             {
             process <- temp_samdf3 %>%
               dplyr::select(-one_of("target_gene"))%>%
               left_join(params %>% dplyr::select(pcr_primers, target_gene, ref_db, blast_db)) %>%
               group_by(target_gene) %>%
               nest() %>% 
               mutate(idtaxa = purrr::map(data, ~{
                 ref_dbs <- .x %>%
                   tidyr::separate_rows(ref_db, sep=";") %>%
                   pull(ref_db) %>%
                   unique() %>%
                   basename() %>% 
                   str_remove("\\..*$")
                 taxtabs <- tax_idtaxa[str_detect(tax_idtaxa, ref_dbs)] %>%
                   purrr::map(readRDS)
                if(length(taxtabs) == 1){
                  out <- taxtabs[[1]]
                } else if(length(taxtabs) == 2){
                  out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                } else if(length(taxtabs) == 3){
                  temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                  out <- coalesce_tax(temptax, taxtabs[[3]])
                } 
              }),
               blast = purrr::map(data, ~{
                 blast_dbs <- .x %>%
                   tidyr::separate_rows(blast_db, sep=";") %>%
                   pull(blast_db) %>%
                   unique() %>%
                   basename() %>% 
                   str_remove("\\..*$")
                 taxtabs <- tax_blast[str_detect(tax_blast, blast_dbs)]%>%
                  purrr::map(readRDS)
                 if(length(taxtabs) == 1){
                   out <- taxtabs[[1]]
                 } else if(length(taxtabs) == 2){
                   out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                 } else if(length(taxtabs) == 3){
                   temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                   out <- coalesce_tax(temptax, taxtabs[[3]])
                 } 
               })) %>%
               mutate(joint_tax = purrr::pmap(list(target_gene, idtaxa, blast),
                                              .f = ~step_join_tax_blast(
                                                tax = ..2,
                                                blast_spp = ..3,
                                                ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                                output = paste0("output/rds/",..1,"_taxblast.rds"),
                                                propagate_tax = TRUE)
               )) %>%
               unnest(data)
             out <- paste0("output/rds/",unique(process$target_gene),"_taxblast.rds")
             return(out)
            }, format="file", iteration = "vector"),
  
## Merge taxonomy tables  -----------------------------------------------------
 tar_target(merged_tax,
            {
            process <- temp_samdf3 %>%
            dplyr::select(-one_of("target_gene"))%>%
            left_join(params %>% dplyr::select(pcr_primers, target_gene)) %>%
            nest(data=everything())%>%
            mutate(merged_tax = purrr::map(data,
                                          .f = ~{
                                            taxtabs <- joint_tax
                                            taxtabs <- taxtabs[taxtabs %>%
                                                                 purrr::map_lgl(function(y){
                                                                   any(str_detect(y, unique(.x$target_gene)))
                                                                 })]
                                            tax <- taxtabs %>%
                                              purrr::map(readRDS) %>%
                                              bind_rows() %>%
                                              as.matrix()
                                            saveRDS(tax, "output/rds/final_tax.rds")
                                            return(tax)
                                          })) %>%
          unnest(data)            
          out <- "output/rds/final_tax.rds"
          return(out)
          }, format="file", iteration = "vector"),

# Create phyloseq object --------------------------------------------------
  tar_file(seqtab_path, merged_seqtab),
  tar_file(taxtab_path, merged_tax),
  tar_target(samdf_path, make_temp_samdf2),
  
  tar_target(ps,{
    process <- step_phyloseq(
       seqtab_path = seqtab_path,
       taxtab_path = taxtab_path,
       samdf_path = samdf_path,
       seqs_path=NULL,
       phy_path=NULL)
    out <- "output/rds/ps.rds"
    saveRDS(process, out)
    return(out)
  }, format="file", iteration = "vector"),

## Output unfiltered results -----------------------------------------------
  tar_target(ps_summary, {
    out <- ps %>%
     readRDS()%>%
     step_output_summary(out_dir="output/results/unfiltered", type="unfiltered")
    return(out)
  }, format="file", iteration = "vector"),
  
## Filter phyloseq ---------------------------------------------------------
# Subset taxa in phyloseq object
tar_target(ps0,
           ps %>%
           readRDS() %>%
           subset_taxa(
             Phylum == "Arthropoda"
           ) %>%
          filter_taxa(function(x) mean(x) > 0, TRUE)
),

# Filter low abundance samples
tar_target(ps1,
           {
           process <- ps0 %>%
           step_filter_phyloseq(
             min_reads=1000, 
             quiet=FALSE
           )
           out <- "output/rds/ps_filtered.rds"
           saveRDS(process, out)
           return(out)
           }, format="file", iteration = "vector"),

## Output filtered results -------------------------------------------------
tar_target(ps_filt_summary, {
  out <- ps1 %>%
    readRDS()%>%
    step_output_summary(out_dir="output/results/filtered", type="filtered")
  return(out)
}, format="file", iteration = "vector"),

## Output imappests format -----------------------------------------------
tar_target(ps_imap_output, {
           out <- ps1 %>%
             readRDS()%>%
           step_output_imap(out_dir="output/results/final")
  return(out)
}, format="file", iteration = "vector")
)