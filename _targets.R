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
    fastq_path,
    purrr::map(list.dirs("data", recursive=FALSE),
                          list.files, pattern="_R1_", full.names = TRUE) %>%
      unlist() 
  ),

  tar_target(temp_samdf1, step_check_files(samdf, fastq_path)),
  
  tar_group_by(temp_samdf1_grouped, temp_samdf1, fcid),

# Sequencing QC -------------------------------------------------------
  tar_target(seq_qc,
             temp_samdf1_grouped %>%
                dplyr::group_by(fcid) %>%
                tidyr::nest() %>%
                dplyr::mutate(seq_qc = purrr::map(fcid, step_seq_qc, quiet=FALSE, write_all=FALSE)),
             pattern = map(temp_samdf1_grouped), iteration = "vector"),
  
  tar_target(switching_qc,
             temp_samdf1_grouped %>%
               dplyr::group_by(fcid) %>%
               tidyr::nest() %>%
               dplyr::mutate(switching_qc = purrr::map(fcid, step_switching_calc, multithread=FALSE, quiet=TRUE)),
             pattern = map(temp_samdf1_grouped), iteration = "vector"),
             

# Demultiplex and trim primers --------------------------------------------
tar_target(primer_trim,
           {
             temp_samdf1 %>%
               dplyr::mutate(primer_trim = purrr::pmap(dplyr::select(., sample_id, for_primer_seq, rev_primer_seq, pcr_primers, fcid),
                                                .f = ~step_primer_trim(sample_id = ..1, for_primer_seq=..2, rev_primer_seq=..3, pcr_primers = ..4,
                                                                       input_dir = paste0("data/",..5), output_dir =  paste0("data/",..5,"/01_trimmed"),
                                                                       qc_dir=paste0("output/logs/",..5), quiet = FALSE)))%>%
               dplyr::select(sample_id, sample_name, fcid, primer_trim)
           },
           pattern = map(temp_samdf1), iteration = "vector"),

# Return filepath for tracking
tar_target(primer_trim_path,
           {
             outF <- primer_trim$primer_trim %>%  
               dplyr::bind_rows()%>%
               dplyr::pull(fwd_out) 
             outR <- primer_trim$primer_trim %>%
               dplyr::bind_rows()%>%
               dplyr::pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed path
             return(c(outF,outR))
             
           },
           pattern = map(primer_trim), format="file", iteration = "vector"),

## Make temporary samdf
 tar_target(temp_samdf2, {
   temp_samdf1 %>%
                    dplyr::select(-where(is.list)) %>%
                    step_demux_samdf() %>%
                    step_check_files(primer_trim_path)
   }),

## Filter reads ------------------------------------------------------------
  tar_target(read_filter,
             {
             temp_samdf2 %>%
                dplyr::select(sample_id, sample_name, pcr_primers, fcid) %>%
                dplyr::left_join(params %>% 
                          dplyr::select(target_gene, pcr_primers, read_min_length, read_max_length, read_max_ee, 
                                        read_trunc_length, read_trim_left, read_trim_right), by = "pcr_primers") %>%
                 dplyr::mutate(read_filter = purrr::pmap(dplyr::select(., sample_id, fcid, read_min_length, read_max_length, read_max_ee, 
                                                                read_trunc_length, read_trim_left, read_trim_right),
                   .f = ~step_filter_reads(
                    sample_id = ..1,
                    input_dir = paste0("data/",..2,"/01_trimmed/"),
                    output_dir = paste0("data/",..2,"/02_filtered"),
                    min_length = ..3,
                    max_length = ..4,
                    max_ee = ..5,
                    trunc_length = ..6,
                    trim_left = ..7,
                    trim_right = ..8,
                    rm.lowcomplex = 0,
                    quiet = FALSE)))%>%
                 dplyr::select(sample_id, sample_name, fcid, read_filter)
             },                  
             pattern = map(temp_samdf2), iteration = "vector"),

# Return filepath for tracking
tar_target(read_filter_path,
           {
             outF <- read_filter$read_filter %>%  
               dplyr::bind_rows()%>%
               dplyr::pull(fwd_out) 
             outR <- read_filter$read_filter %>%
               dplyr::bind_rows()%>%
               dplyr::pull(rev_out) 
             # Check for empty files
             outF <- outF[file.exists(outF)]
             outR <- outR[file.exists(outR)]
             # Return list of completed path
             return(c(outF,outR))
             
           },
           pattern = map(read_filter), format="file", iteration = "vector"),

# Make temporary samdf
 tar_target(temp_samdf3, {
   temp_samdf2 %>%
           dplyr::select(-where(is.list)) %>%
           step_check_files(read_filter_path) 
 }),


## Pre-filtering quality plots ---------------------------------------------

# Sample a random set of samples for read quality plotting - Sample the different pcr_primers too!
tar_target(read_samples,{
  group_sizes <- temp_samdf3 %>%
    dplyr::group_by(fcid, pcr_primers) %>%
    dplyr::group_size()
  if(all(group_sizes > 10)){
    n_samples <- 10
  } else {
    n_samples = min(group_sizes)
  }
  out <- temp_samdf3 %>%
    dplyr::group_by(fcid, pcr_primers) %>%
    dplyr::slice_sample(n=n_samples) 
}),

tar_target(prefilt_qualplots,
          read_samples %>%
           dplyr::mutate(prefilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/01_trimmed/"), truncLen=NULL, quiet = FALSE, n = 10000)
            )),
            pattern = map(read_samples), iteration = "vector"),

## Write out prefilt qualplots
tar_target(write_prefilt_qualplots, {
  prefilt_qualplots %>% 
    dplyr::group_by(fcid) %>%
    tidyr::nest() %>%
    purrr::pwalk(list(fcid, data),
      .f = ~{
      pdf(file=paste0("output/logs/",..1,"/", ..1,"_prefilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
      print(..2$prefilt_qualplots)
      try(dev.off(), silent=TRUE)
      })
  out <- paste0("output/logs/",unique(prefilt_qualplots$fcid),"/", unique(prefilt_qualplots$fcid),"_prefilt_qualplots.pdf")
  return(out)
  }, format="file", iteration = "vector"),

## Post-filtering quality plots --------------------------------------------
tar_target(postfilt_qualplots,
           read_samples %>%
            dplyr::mutate(postfilt_qualplots = purrr::pmap(list(sample_id, fcid),
                 .f = ~plot_read_quals(sample_id = ..1,
                 input_dir = paste0("data/",..2,"/02_filtered/"), truncLen=NULL, quiet = FALSE, n = 10000)
           )),
            pattern = map(read_samples), iteration = "vector"),
 
# Write out postfilt qualplots
tar_target(write_postfilt_qualplots, {
  postfilt_qualplots %>% 
    dplyr::group_by(fcid) %>%
    tidyr::nest() %>%
    purrr::pwalk(list(fcid, data),
                 .f = ~{
                   pdf(file=paste0("output/logs/",..1,"/", ..1, "_postfilt_qualplots.pdf"), width = 11, height = 8 , paper="a4r")
                   print(..2$postfilt_qualplots)
                   try(dev.off(), silent=TRUE)
                 })
  out <- paste0("output/logs/",unique(postfilt_qualplots$fcid),"/", unique(postfilt_qualplots$fcid),"_postfilt_qualplots.pdf")
  return(out)
}, format="file", iteration = "vector"),

## Infer sequence variants with DADA2 --------------------------------------

  # Group temporary samdf by fcid
  tar_group_by(temp_samdf3_grouped, temp_samdf3, fcid),
 
  # How to make it just redo one of the dadas if only one runs filtered file changed changed?
  tar_target(dada,{
             temp_samdf3_grouped %>%
             dplyr::group_by(fcid) %>%
             tidyr::nest() %>%
             dplyr::mutate(dada2 = purrr::map(fcid,
                                        .f = ~step_dada2(fcid = .x,
                                                         input_dir = paste0("data/",.x,"/02_filtered"),
                                                         output = paste0("output/rds/",.x,"_seqtab.rds"),
                                                         qc_dir = paste0("output/logs/",.x),
                                                         quiet = FALSE,
                                                         write_all = FALSE)
             ))
             },
             pattern = map(temp_samdf3_grouped), iteration = "vector"),


# Return filepath for tracking
tar_target(dada_path,
           {
             return(paste0("output/rds/",unique(dada$fcid),"_seqtab.rds"))
           },
           pattern = map(dada), format="file", iteration = "vector"),

  
##  Merge infered variants from each run and subset to target loci ---------
 tar_target(subset_seqtab, {
            process <- temp_samdf3 %>%
             dplyr::ungroup() %>%
             dplyr::group_by(pcr_primers) %>%
             tidyr::nest() %>%
             dplyr::mutate(subset_seqtab = purrr::map(pcr_primers, 
                   .f = ~{
                   #seqtabs <- list.files("output/rds/", pattern="seqtab.rds", full.names = TRUE)
                   if(length(dada_path) > 1){
                   st.all <- mergeSequenceTables(tables=dada_path)
                   } else if(length(dada_path) == 1) {
                   st.all <- readRDS(dada_path)
                   }
                   st.all <- st.all[str_detect(rownames(st.all), .x),]
                   st.all <- st.all[,colSums(st.all) > 0]
                   saveRDS(st.all, paste0("output/rds/",.x,"_seqtab.rds"))
                   out <- rowSums(st.all) %>%
                     tibble::enframe(name="fq", value="subset_seqtab_reads")
                   return(out)
              })) %>% 
              tidyr::unnest(data, subset_seqtab) %>%
              dplyr::select(sample_id, sample_name, fcid, subset_seqtab_reads)  %>%
              dplyr::mutate(path = paste0("output/rds/",unique(pcr_primers),"_seqtab.rds"))
            }, iteration = "vector"),


# Return filepath for tracking
tar_target(subset_seqtab_path,
           {
             return(unique(subset_seqtab$path))
           }, format="file"),


#  Filter ASV's per locus -------------------------------------------------
 tar_target(filtered_seqtab, {
          temp_samdf3 %>%
           dplyr::select(-one_of("asv_min_length", "asv_max_length", "phmm", "coding", "genetic_code"))%>%
           dplyr::left_join(params %>% dplyr::select(pcr_primers, asv_min_length, asv_max_length, phmm, coding, genetic_code)) %>%
           dplyr::group_by(pcr_primers, asv_min_length, asv_max_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq) %>%
           tidyr::nest() %>%
           dplyr::mutate(subset_seqtab = purrr::map(pcr_primers, ~{
                readRDS(subset_seqtab_path[stringr::str_detect(subset_seqtab_path, .x)])
           })) %>%
           dplyr::ungroup()%>%
           dplyr::mutate(filtered_seqtab = purrr::pmap(dplyr::select(.,pcr_primers, subset_seqtab, asv_min_length, asv_max_length, phmm, coding, genetic_code, for_primer_seq, rev_primer_seq),
               .f = ~step_filter_asvs(
               seqtab = ..2,
               output = paste0("output/rds/",..1,"_seqtab.cleaned.rds"),
               qc_dir = "output/logs/",
               min_length = ..3,
               max_length = ..4,
               phmm = ..5,
               check_frame = ..6,
               genetic_code = ..7,
               primers = c(..8, ..9),
               multithread = FALSE, 
               quiet = FALSE)
         )) %>% 
         unnest_wider(filtered_seqtab) %>%
         dplyr::mutate(filtered_asvs = purrr::map(filtered_asvs, ~{
           .x %>%
             dplyr::select(-sample_id) # remove sample_id from the nested column to avoid duplicate name columns
         })) %>%
         tidyr::unnest(c(data, filtered_asvs))%>%
         dplyr::select(sample_id, sample_name, fcid, reads_starting, reads_chimerafilt, pcr_primers, reads_lengthfilt,
                       reads_phmmfilt, reads_framefilt, reads_final, plot)%>%
         dplyr::mutate(path = paste0("output/rds/",pcr_primers,"_seqtab.cleaned.rds"))
 }, iteration = "vector"),
 
# Return filepath for tracking
tar_target(filtered_seqtab_path,
           {
             return(unique(filtered_seqtab$path))
           }, format="file"),

# Write out seqtab filtering plots
tar_target(write_seqtab_qualplots, {
  pdf("output/logs/ASV_cleanup_summary.pdf", width = 11, height = 8 , paper="a4r")
    print(unique(filtered_seqtab$plot))
  try(dev.off(), silent=TRUE)
  out <- "output/logs/ASV_cleanup_summary.pdf"
  return(out)
}, format="file", iteration = "vector"),


## Merge all loci into final seqtab ----------------------------------------
  tar_target(merged_seqtab_path, {
    process <- temp_samdf3 %>%
      tidyr::nest(data=everything()) %>%
      dplyr::mutate(final_seqtab = purrr::map(data, ~{
        seqtabs <- filtered_seqtab_path
        seqtabs <- seqtabs[seqtabs %>%
                             purrr::map_lgl(function(y){
                               any(stringr::str_detect(y, unique(.x$pcr_primers)))
                             })]
        # Remove empty seqtabs
        empty <- seqtabs[purrr::map_lgl(seqtabs, function(z){
          !ncol(readRDS(z)) > 0
        })]
        if(length(empty) > 0 ){
          warning(paste0("No sequences in: ", empty, " skipping"))
        }
        if(length(seqtabs[!seqtabs %in% empty]) > 1){
          st.all <- dada2::mergeSequenceTables(tables=seqtabs[!seqtabs %in% empty])
        } else if(length(seqtabs[!seqtabs %in% empty]) == 1) {
          st.all <- readRDS(seqtabs[!seqtabs %in% empty])
        }
        saveRDS(st.all, "output/rds/seqtab_final.rds")
        return(TRUE)
      })) %>%
      tidyr::unnest(data)
    out <- paste0("output/rds/seqtab_final.rds")
    return(out)
   }, format="file", iteration = "vector"),

# Assign taxonomy ---------------------------------------------------------
  tar_file(ref_db_tracked,
           params  %>%
             dplyr::pull(ref_db) %>%
             unique()%>%
             stringr::str_split(pattern=";", n=Inf) %>% 
             unlist()
  ),
  tar_file(blast_db_tracked,
           params  %>%
             dplyr::pull(blast_db) %>%
             unique() %>%
             stringr::str_split(pattern=";", n=Inf) %>% 
             unlist()
  ),

## IDTAXA -------------------------------------------------------------------
 tar_target(tax_idtaxa,{ 
   process <- temp_samdf3 %>%
     dplyr::select(-one_of("target_gene", "ref_db"))%>%
     dplyr::left_join(params %>% dplyr::select(pcr_primers, target_gene, ref_db)) %>%
     tidyr::separate_rows(ref_db, sep=";") %>%
     dplyr::group_by(target_gene, pcr_primers, ref_db) %>%
     tidyr::nest()  %>% 
     dplyr::mutate(ref_db2 = purrr::map(ref_db, ~{
       ref_db_tracked[stringr::str_detect(ref_db_tracked, .x)]
     }))  %>%
     unnest(ref_db2)%>%
     dplyr::mutate(filtered_seqtab = purrr::map(pcr_primers, ~{
       readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path, .x)])
     }))  %>%
     dplyr::mutate(idtaxa = purrr::pmap(list(target_gene, pcr_primers, filtered_seqtab, ref_db2),
                                 .f = ~step_idtaxa(
                                   seqtab = ..3,
                                   database = ..4,
                                   ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species"),
                                   output = paste0("output/rds/",..2,"_",basename(..4) %>% stringr::str_remove("\\..*$"),"_taxtab.rds"),
                                   qc_dir = "output/logs/",
                                   threshold = 60,
                                   multithread = FALSE, 
                                   quiet = FALSE,
                                   plot=TRUE)
     ))%>%
     tidyr::unnest(data)
   out <- paste0("output/rds/",unique(process$pcr_primers),"_",unique(basename(process$ref_db)  %>% stringr::str_remove("\\..*$")),"_taxtab.rds")
   return(out)
   }, format="file", iteration = "vector"),

## BLAST -------------------------------------------------------------------
 tar_target(tax_blast,
            {
           process <- temp_samdf3 %>%
             dplyr::select(-one_of("target_gene", "blast_db"))%>%
             dplyr::left_join(params %>% dplyr::select(pcr_primers, target_gene, blast_db)) %>%
             tidyr::separate_rows(blast_db, sep=";") %>%
             dplyr::group_by(target_gene, pcr_primers, blast_db) %>%
             tidyr::nest()  %>% 
             dplyr::mutate(blast_db2 = purrr::map(blast_db, ~{
               blast_db_tracked[stringr::str_detect(blast_db_tracked, .x)]
             }))  %>%
             unnest(blast_db2) %>%
             dplyr::mutate(filtered_seqtab = purrr::map(pcr_primers, ~{
               readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path, .x)])
             }))  %>%
             dplyr::mutate(blast = purrr::pmap(list(target_gene, pcr_primers, filtered_seqtab, blast_db2),
                                        .f = ~step_blast_tophit(
                                          seqtab = ..3,
                                          database = ..4,
                                          ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                          output = paste0("output/rds/",..2,"_",basename(..4) %>% stringr::str_remove("\\..*$"),"_blast.rds"),
                                          qc_dir = "output/logs/",
                                          identity = 97,  
                                          coverage=95,
                                          evalue=1e06,
                                          max_target_seqs=5,
                                          max_hsp=5, 
                                          multithread = FALSE, 
                                          quiet = FALSE)
           ))
          out <- paste0("output/rds/",unique(process$pcr_primers), "_",unique(basename(process$blast_db)  %>% stringr::str_remove("\\..*$")),"_blast.rds")
          return(out)
        }, format="file", iteration = "vector"),
 
### Aggregate taxonomic assignment methods-----------------------------------------------
# Need to edit this to detect the PCR_primers as well
 tar_target(joint_tax,
            {
             process <- temp_samdf3 %>%
              dplyr::select(-one_of("target_gene"))%>%
              dplyr::left_join(params %>% dplyr::select(pcr_primers, target_gene, ref_db, blast_db)) %>%
              dplyr::group_by(pcr_primers) %>%
              tidyr::nest() %>% 
              dplyr::mutate(idtaxa = purrr::map2(data, pcr_primers, ~{
                ref_dbs <- .x %>%
                  tidyr::separate_rows(ref_db, sep=";") %>%
                  dplyr::pull(ref_db) %>%
                  unique() %>%
                  basename() %>% 
                  stringr::str_remove("\\..*$")
                taxtabs <- tax_idtaxa[stringr::str_detect(tax_idtaxa, ref_dbs)& stringr::str_detect(tax_idtaxa, .y)] %>%
                  purrr::map(readRDS)
                if(length(taxtabs) == 1){
                  out <- taxtabs[[1]]
                } else if(length(taxtabs) == 2){
                  out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                } else if(length(taxtabs) == 3){
                  temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                  out <- coalesce_tax(temptax, taxtabs[[3]])
                }
                return(out)
              }),
              blast = purrr::map2(data,pcr_primers, ~{
                blast_dbs <- .x %>%
                  tidyr::separate_rows(blast_db, sep=";") %>%
                  dplyr::pull(blast_db) %>%
                  unique() %>%
                  basename() %>% 
                  stringr::str_remove("\\..*$")
                taxtabs <- tax_blast[stringr::str_detect(tax_blast, blast_dbs) & stringr::str_detect(tax_blast, .y)]%>%
                  purrr::map(readRDS)
                if(length(taxtabs) == 1){
                  out <- taxtabs[[1]]
                } else if(length(taxtabs) == 2){
                  out <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                } else if(length(taxtabs) == 3){
                  temptax <- coalesce_tax(taxtabs[[1]], taxtabs[[2]])
                  out <- coalesce_tax(temptax, taxtabs[[3]])
                }
                return(out)
              })) %>%
              dplyr::mutate(joint_tax = purrr::pmap(list(pcr_primers, idtaxa, blast),
                                                     .f = ~step_join_tax_blast(
                                                       tax = ..2,
                                                       blast_spp = ..3,
                                                       ranks = c("Root","Kingdom", "Phylum","Class", "Order", "Family", "Genus","Species") ,
                                                       output = paste0("output/rds/",..1,"_taxblast.rds"),
                                                       propagate_tax = TRUE)
              )) %>%
              tidyr::unnest(data)
            out <- paste0("output/rds/",unique(process$pcr_primers),"_taxblast.rds")
            return(out)
           }, format="file", iteration = "vector"),

### Merge taxonomy tables  -----------------------------------------------------
 tar_target(merged_tax,
            {
            process <- temp_samdf3 %>%
            dplyr::select(-one_of("target_gene"))%>%
            dplyr::left_join(params %>% dplyr::select(pcr_primers, target_gene)) %>%
            tidyr::nest(data=everything())%>%
            dplyr::mutate(merged_tax = purrr::map(data,
                                          .f = ~{
                                            taxtabs <- joint_tax
                                            taxtabs <- taxtabs[taxtabs %>%
                                                                 purrr::map_lgl(function(y){
                                                                   any(stringr::str_detect(y, unique(.x$pcr_primers)))
                                                                 })]
                                            tax <- taxtabs %>%
                                              purrr::map(readRDS) 
                                            tax <- tax[sapply(tax, nrow) >0 ] %>%
                                              dplyr::bind_rows() %>%
                                              as.matrix()
                                            saveRDS(tax, "output/rds/final_tax.rds")
                                            return(tax)
                                          })) %>%
          tidyr::unnest(data)            
          out <- "output/rds/final_tax.rds"
          return(out)
          }, format="file", iteration = "vector"),


### Assignment plot ---------------------------------------------------------

tar_target(assignment_plot, {
  temp_samdf3 %>%
    dplyr::select(-one_of("target_gene", "ref_db"))%>%
    dplyr::left_join(params %>% dplyr::select(pcr_primers, target_gene, ref_db, blast_db)) %>%
    tidyr::separate_rows(blast_db, sep=";") %>%
    dplyr::group_by(target_gene, pcr_primers, ref_db, blast_db) %>%
    tidyr::nest()  %>% 
    dplyr::mutate(filtered_seqtab = purrr::map(pcr_primers, ~{
      readRDS(filtered_seqtab_path[stringr::str_detect(filtered_seqtab_path, .x)])
    }))   %>%
    dplyr::mutate(tax = purrr::map(pcr_primers, ~{
      readRDS(joint_tax[stringr::str_detect(joint_tax, .x)])%>% 
        seqateurs::unclassified_to_na(rownames=FALSE) %>%
        dplyr::mutate(lowest = seqateurs::lowest_classified(.)) 
    })) %>%
    dplyr::mutate(blast = purrr::pmap(list(pcr_primers, filtered_seqtab, blast_db),
                               .f = ~{
                                 seqs <- colnames(..2)
                                 names(seqs) <- colnames(..2)
                                 if(length(seqs) > 0){
                                   blast_top_hit(
                                   query = seqs,
                                   db = ..3,
                                   identity=60,
                                   coverage=80) %>% 
                                   dplyr::select(OTU = qseqid, acc, blastspp = Species, pident, length, evalue, qcovs) 
                                 } else {
                                   tibble::enframe(seqs, name=NULL, value="OTU") %>%
                                     dplyr::mutate(acc = NA_character_,
                                                   blastspp = NA_character_,
                                                   pident = NA_real_, 
                                                   length = NA_integer_,
                                                   evalue = NA_real_,
                                                   qcovs = NA_integer_) 
                                 }
    })) %>%
    dplyr::mutate(joint = purrr::pmap(list(blast, tax),
                               .f = ~{
                                 if(nrow(..1) > 0 & nrow(..2) > 0){
                                 ..1 %>%
                                   dplyr::left_join(..2, by="OTU")
                                 } else {
                                  NULL
                                 }
    })) %>%
    dplyr::mutate(plot = purrr::pmap(list(pcr_primers, joint, ref_db, blast_db),
                              .f= ~{
                                # ADD TITLES HERE!
                                if(!is.null(..2)){
                                  ..2 %>%
                                    dplyr::select(pident, rank = lowest) %>%
                                    dplyr::mutate(rank = factor(rank, levels = c("Root","Kingdom","Phylum","Class","Order","Family","Genus","Species"))) %>%
                                    ggplot(aes(x=pident, fill=rank))+ 
                                    geom_histogram(colour="black", binwidth = 1, position = "stack") + 
                                    labs(title = paste0(..1, " Top hit identity distribution"),
                                         subtitle = paste0("IDTAXA database:", ..3, " BLAST database:", ..4),
                                         x = "BLAST top hit % identity",
                                         y = "OTUs") + 
                                    scale_x_continuous(breaks=seq(60,100,2)) +
                                    scale_fill_brewer(name = "Taxonomic \nAssignment", palette = "Spectral")+
                                    theme_bw()+
                                    theme(
                                      strip.background = element_rect(colour = "black", fill = "lightgray"),
                                      strip.text = element_text(size=9, family = ""),
                                      plot.background = element_blank(),
                                      text = element_text(size=9, family = ""),
                                      axis.text = element_text(size=8, family = ""),
                                      legend.position = "right",
                                      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                                      panel.grid = element_line(size = rel(0.5)),
                                    ) 
                                  } else {
                                    NULL
                                  }
                              }))
}, iteration = "vector"),

# Write out assignment plot
tar_target(write_assignment_plot, {
  pdf(file=paste0("output/logs/taxonomic_assignment_summary.pdf"), width = 11, height = 8 , paper="a4r")
  print(unique(assignment_plot$plot))
  try(dev.off(), silent=TRUE)
  out <- paste0("output/logs/taxonomic_assignment_summary.pdf")
  return(out)
}, format="file", iteration = "vector"),

# Create phyloseq object --------------------------------------------------
  tar_target(ps,{
    process <- step_phyloseq(
       seqtab = merged_seqtab_path,
       taxtab = merged_tax,
       samdf = temp_samdf2,
       seqs_path=NULL,
       phy_path=NULL)
    out <- "output/rds/ps.rds"
    saveRDS(process, out)
    return(out)
  }, format="file", iteration = "vector"),

## Output unfiltered results -----------------------------------------------
  tar_target(ps_summary, {
    ps_obj <- ps %>%
      readRDS()
    
    out <- c(step_output_summary(ps_obj, out_dir="output/results/unfiltered", rank="Species", type="unfiltered"),
             step_output_ps(ps_obj, out_dir="output/results/unfiltered"))
    return(out)
  }, format="file", iteration = "vector"),

  
## Filter phyloseq ---------------------------------------------------------
# Taxonomic and minimum abundance filtering
tar_target(ps_filtered,{
     # if multiple primers were used - split ps into different primers
    process <-  params %>%
       dplyr::mutate(ps_obj = purrr::map(pcr_primers,
              .f = ~{
                    physeq <- ps %>%
                      readRDS()
                    new_samdat <- as(phyloseq::sample_data(physeq), "data.frame") %>%
                      dplyr::filter(pcr_primers == .x)
                    if(nrow(new_samdat) > 0){
                      phyloseq::sample_data(physeq) <- phyloseq::sample_data(new_samdat)
                      return(physeq)
                    } else{
                      return(NULL)
                    }
                      })) %>%
    dplyr::mutate(ps_filt = purrr::pmap(dplyr::select(., ps_obj, target_kingdom, target_phylum, target_class,
                                               target_order, target_family, target_genus, target_species, 
                                               min_sample_reads, min_taxa_reads, min_taxa_ra),
            .f = ~{
              if(!is.null(..1)){
                ..1 %>%
                step_filter_phyloseq(
                  kingdom = ..2,
                  phylum = ..3,
                  class = ..4,
                  order = ..5,
                  family = ..6,
                  genus = ..7,
                  species = ..8,
                  min_sample_reads=..9, 
                  min_taxa_reads = ..10,
                  min_taxa_ra = ..11,
                  quiet=FALSE
                )
              }else {
                return(NULL)
              }
              }))
    ps_merged <- merge_phyloseq_new(process$ps_filt)
    out <- "output/rds/ps_filtered.rds"
    saveRDS(ps_merged, out)
    return(out)
  }, format="file", iteration = "vector"),

## Output filtered results -------------------------------------------------
tar_target(ps_filt_summary, {
  ps_obj <- ps_filtered %>%
    readRDS()
  
  out <- c(step_output_summary(ps_obj, out_dir="output/results/unfiltered", rank="Species", type="unfiltered"),
           step_output_ps(ps_obj, out_dir="output/results/unfiltered"))
  return(out)
}, format="file", iteration = "vector"),

# Read tracking plot ------------------------------------------------------
tar_target(read_tracking, {
  ps_obj <- readRDS(ps)
  phyloseq::tax_table(ps_obj) <- phyloseq::tax_table(ps_obj) %>%
    as("data.frame") %>%
    seqateurs::unclassified_to_na() %>%
    as.matrix()
    
  read_tracker <- temp_samdf1 %>%
    dplyr::select(sample_name, fcid, pcr_primers) %>%
    dplyr::left_join(primer_trim %>%
                tidyr::unnest(primer_trim) %>%
                dplyr::select(sample_name, fcid, trimmed_input, trimmed_output, fwd_out) %>%
                dplyr::mutate(sample_id = basename(fwd_out) %>% stringr::str_remove("_S[0-9]+_R[1-2]_.*$")) %>%
                dplyr::select(sample_name, sample_id, fcid, input_reads=trimmed_input, trimmed=trimmed_output), by = c("sample_name", "fcid")) %>%
    dplyr::left_join(read_filter %>%
                tidyr::unnest(read_filter) %>%
                dplyr::select(sample_name, sample_id, fcid, filtered = filter_output),
              by = c("sample_name", "fcid", "sample_id")) %>%
    dplyr::left_join(dada %>% 
                tidyr::unnest(dada2) %>% 
                dplyr::select(fcid, sample_id, denoised=merged),
              by = c("fcid", "sample_id")
              ) %>%
    dplyr::left_join(filtered_seqtab %>% 
                dplyr::select(sample_id, sample_name, fcid, chimerafilt=reads_chimerafilt,
                              lengthfilt= reads_lengthfilt, phmmfilt=reads_phmmfilt, framefilt = reads_framefilt),
              by = c("sample_name", "fcid", "sample_id")
              ) %>%
    dplyr::left_join(psmelt(ps_obj) %>% # Could replace this with seqateurs::lowest_classified?
                dplyr::filter(Abundance > 0) %>%
                dplyr::select(sample_id, fcid, Abundance, any_of(colnames(phyloseq::tax_table(ps_obj)))) %>%
                pivot_longer(cols=any_of(colnames(phyloseq::tax_table(ps_obj))), 
                             names_to = "rank",
                             values_to="name") %>%
                dplyr::group_by(sample_id, rank) %>%
                summarise(Abundance = sum(Abundance)) %>%
                pivot_wider(names_from="rank",
                            values_from="Abundance")%>%
                rename_with(~stringr::str_to_lower(.), everything()) %>%
                rename_with(~stringr::str_c("classified_", .), -sample_id), by="sample_id")  %>%
    dplyr::select(any_of(c(
      "sample_name","sample_id", "pcr_primers", "fcid", "input_reads", "trimmed", "filtered",
      "denoised", "chimerafilt", "lengthfilt", "phmmfilt", "framefilt", 
      "classified_root", "classified_kingdom", "classified_phylum","classified_class",
      "classified_order", "classified_family", "classified_genus", "classified_species"
    )))
  
  write_csv(read_tracker, "output/logs/read_tracker.csv")

  gg.read_tracker <- read_tracker %>%
    pivot_longer(cols = -c("sample_name", "sample_id", "fcid", "pcr_primers"),
                 names_to = "step",
                 values_to="reads") %>%
    dplyr::mutate(step = factor(step, levels=c(
      "input_reads", "trimmed", "filtered",
      "denoised", "chimerafilt", "lengthfilt", "phmmfilt", "framefilt", 
      "classified_root", "classified_kingdom", "classified_phylum","classified_class",
      "classified_order", "classified_family", "classified_genus", "classified_species"
    ))) %>%
    ggplot(aes(x = step, y = reads, fill=pcr_primers))+
    geom_col() +
    scale_y_continuous(labels = scales::label_number_si())+
    facet_grid(fcid~.)+
    theme_bw()+
    theme(
      strip.background = element_rect(colour = "black", fill = "lightgray"),
      strip.text = element_text(size=9, family = ""),
      axis.text.x =element_text(angle=45, hjust=1, vjust=1),
      plot.background = element_blank(),
      text = element_text(size=9, family = ""),
      axis.text = element_text(size=8, family = ""),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      panel.grid = element_line(size = rel(0.5)),
    ) +
    labs(x = "Pipeline step",
         y = "Reads retained",
         fill = "PCR primers")
  pdf(file=paste0("output/logs/read_tracker.pdf"), width = 11, height = 8 , paper="a4r")
    print(gg.read_tracker)
  try(dev.off(), silent=TRUE)
  
  return(c("output/logs/read_tracker.csv", "output/logs/read_tracker.pdf"))
}, format="file", iteration = "vector")

)

