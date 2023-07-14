library(readr)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggrepel)
library(readxl)
library(tidyr)
library(plyr)
library(cowplot)
library(GenomicRanges)
library(ggrepel)
library(here)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

generate_csv <- FALSE

all_files <- list("IgG_rep1" = "IgG_rep1_closest.bed",
                  "IgG_rep2" = "IgG_rep2_closest.bed",
                  "Rbfox2_rep1" = "Rbfox2_rep1_closest.bed",
                  "Rbfox2_rep2" = "Rbfox2_rep2_closest.bed",
                  "IgG_merged" = "IgG_merged_closest.bed",
                  "RBFOX2_merged" = "RBFOX2_merged_closest.bed")

fold_enrichment_cutoff <- 2

if(generate_csv){
  all_bed_files <- lapply(names(all_files), function(x){
    sample_bed <- all_files[[x]]
    closest <- read_delim(here("revision/SMI_results/genes_bed",
                               sample_bed), 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
    
    colnames(closest) <- c("peak_chr", "peak_start", "peak_end", "log10p",
                           "pvalue", "log2_fold", "peak_strand", "peak_counts_clip", "peak_counts_input", "total_clip_counts", "total_input_counts", 
                           "gene_chr", "gene_start", "gene_end",
                           "gene_ens_id", "gene_score",
                           "gene_strand", "gene_source",
                           "gene_type", "gene_phase", "other_info")
    
    closest <- closest %>%
      dplyr::filter(gene_type == "gene")
    
    full_df <- lapply(1:nrow(closest), function(x){
      full_row <- closest  %>%
        dplyr::slice(x)
      keep_cols <- full_row %>%
        dplyr::select(!"other_info")
      new_info <- strsplit(as.character(full_row$other_info), split = "; ")
      
      named_vect <- lapply(new_info[[1]], function(y){
        name_obj <- strsplit(y, split = " ")
        name_df <- data.frame(name_obj[[1]][2])
        colnames(name_df) <- name_obj[[1]][1]
        return(name_df)
      })
      new_col <- do.call(cbind, named_vect)
      keep_cols <- cbind(keep_cols, new_col)
      return(keep_cols)
    })
    
    eCLIP <- do.call(rbind.fill, full_df)
    
    eCLIP$gene_name <- gsub('"', '', eCLIP$gene_name)
    
    write.csv(eCLIP, file=here("revision/SMI_results/genes_bed",
                               paste0(x, "_eCLIP.csv")))
    
    return(eCLIP)
  })
  
} else {
  all_bed_files <- lapply(names(all_files), function(x){
    read_csv(here("revision/SMI_results/genes_bed",
                  paste0(x, "_eCLIP.csv")))
  })
}




names(all_bed_files) <- names(all_files)

no_filter_bed_files <- all_bed_files

#### Investigate Replicates and Combined Data ----------------------
# 1. Identify peaks that are enriched in the RBFOX2 or IgG eCLIP over the SMI.

all_bed_files <- lapply(all_bed_files, function(x){
  # Select only significant (p < 0.05) peaks that are enriched in the eCLIP
  # over the control (log2FC > 0) and a minimum read length
  x  <- x %>% 
    filter(log2_fold > 0) %>% 
    mutate(Fold_Enrichment = 2^log2_fold)
  
  # Fix the gene_name region
  x$gene_name <- gsub('"', '', x$gene_name)
  
  return(x)
})

# 2. Determine peak width distribution across all samples. 

sig_reads <- lapply(all_bed_files, function(x){
  # filter for only enriched and significant peaks
  x <- x %>% 
    filter(pvalue < (10^-3)) %>% 
    filter(Fold_Enrichment >= fold_enrichment_cutoff) 
  
  return(x)
  
})

names(sig_reads) <- names(all_files)


# write a function to generate a histogram of read length distribution
peak_length <- function(data, plot_name) {
  data_length <- data %>%
    dplyr::mutate(peak_length = peak_end - peak_start) %>%
    dplyr::select(peak_length)
  
  dataset_name <- gsub("_", " ", deparse(substitute(data)))
  mean_value <- mean(data_length$peak_length)
  min_value <- min(data_length$peak_length)
  max_value <- max(data_length$peak_length)
  
  p <- ggplot(data_length, aes(x = peak_length)) +
    geom_histogram() +
    geom_vline(aes(xintercept = mean_value, color = "red"),
               linetype = "dashed",
               show.legend = FALSE) +
    labs(title = dataset_name, x = "Peak Length", y = "Counts") +
    theme_classic() +
    scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 50))
  
  p <- p +
    annotate("text", x = 300, y = 0, vjust = 1, hjust = 1,
             label = paste("Mean:", round(mean_value, 0))) +
    annotate("text", x = 100, y = 0, vjust = 1, hjust = 1,
             label = paste("Min:", min_value)) +
    annotate("text", x = 500, y = 0, vjust = 1, hjust = 1,
             label = paste("Max:", max_value)) +
    ggplot2::ggtitle(plot_name)
  
  return(p)
} 

all_histograms <- lapply(names(sig_reads), function(x){
  return(peak_length(sig_reads[[x]], x))
})


# combine the plots
cowplot::plot_grid(plotlist = all_histograms)


# 3. Identify peaks within genes and further filter based on peak width and fold enrichment. 
# Identify where RBFOX2 eCLIP peaks are within genes using peakAnno.

filterAndAntiJoin <- function(data) {
  NC <- data[grepl("Rik", data$gene_name), ]
  Gm <- data[grepl("Gm", data$gene_name), ]
  peaks_gm <- anti_join(data, Gm)
  peaks_rm <- anti_join(peaks_gm, NC)
  
  return(peaks_rm)
}

sig_reads <- lapply(sig_reads, function(x){
  return(filterAndAntiJoin(x))
})

# Get annotation with chip seeker
all_annotation <- lapply(names(sig_reads), function(x){
  
  all_granges <- makeGRangesFromDataFrame(sig_reads[[x]],
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          seqinfo=NULL,
                                          seqnames.field=c("peak_chr"),
                                          start.field="peak_start",
                                          end.field=c("peak_end"),
                                          strand.field="peak_strand",
                                          starts.in.df.are.0based=FALSE)
  # use peakAnno to annotate where in the gene peaks are
  peakAnno <- annotatePeak(all_granges, TxDb = txdb, annoDb = "org.Mm.eg.db")
  #annotation_types <- gsub("\\s*\\(.*\\)", "", peakAnno$annotation)
  
  return_data <- peakAnno@annoStat
  
  colnames(return_data) <- c("Feature", paste0("Frequency_", x))
  
  full_data <- as.data.frame(peakAnno)
  
  return(list("stats" = return_data, "full" = full_data))
})

stats <- lapply(all_annotation, function(x){
  return(x$stats)
})

stats <- Reduce(merge, stats)

stats <- stats %>%
  tidyr::pivot_longer(cols = starts_with("Frequency"),
                      values_to = "Frequency", names_to = "type") %>%
  dplyr::mutate(type = sub("Frequency_", "", type))

ggplot2::ggplot(stats, ggplot2::aes(x = type, y = Frequency, fill = Feature)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

full_anno <- lapply(all_annotation, function(x){
  anno_peaks <- x$full
  anno_peaks_3UTR_peaks <- anno_peaks %>% filter(annotation == "3' UTR")
  anno_peaks_Intron_peaks <- dplyr::filter(anno_peaks, grepl('Intron', annotation))
  anno_peaks_5UTR_peaks <- anno_peaks %>% filter(annotation == "5' UTR")
  anno_peaks_Exon_peaks <- dplyr::filter(anno_peaks, grepl('Exon', annotation))
  all_peaks <- do.call(rbind, list(anno_peaks_Exon_peaks, anno_peaks_5UTR_peaks,
                      anno_peaks_Intron_peaks, anno_peaks_3UTR_peaks))
  return(all_peaks)
})

names(full_anno) <- names(sig_reads)


peak_exon_overlap <- function(rMATS_data, eCLIP_data,
                              nt_start = 0,
                              nt_end = 300, by = 50,
                              bootstrap = FALSE,
                              nbootstrap = 1000){
  
  # We only want to run this once
  # filter rMATS_data for spliced exons (both included and skipped) and not
  # spliced exons
  
  psis_inclusion <- dplyr::filter(rMATS_data, FDR < 0.05 & 
                                    IncLevelDifference > 0.01)
  psis_exclusion <- dplyr::filter(rMATS_data, FDR < 0.05 & 
                                    IncLevelDifference < -0.01)
  psis_insensitive <- dplyr::filter(rMATS_data, FDR >= 0.05) 
  
  
  # make a GRanages object for each of the new dataframes
  psis_inclusion_gr <- makeGRangesFromDataFrame(psis_inclusion,
                                                keep.extra.columns=TRUE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c("chr"),
                                                start.field="exonStart_0base",
                                                end.field=c("exonEnd"),
                                                strand.field="strand",
                                                starts.in.df.are.0based=FALSE)
  
  psis_exclusion_gr <- makeGRangesFromDataFrame(psis_exclusion,
                                                keep.extra.columns=TRUE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c("chr"),
                                                start.field="exonStart_0base",
                                                end.field=c("exonEnd"),
                                                strand.field="strand",
                                                starts.in.df.are.0based=FALSE)
  
  psis_insensitive_gr <- makeGRangesFromDataFrame(psis_insensitive,
                                                  keep.extra.columns=TRUE,
                                                  ignore.strand=FALSE,
                                                  seqinfo=NULL,
                                                  seqnames.field=c("chr"),
                                                  start.field="exonStart_0base",
                                                  end.field=c("exonEnd"),
                                                  strand.field="strand",
                                                  starts.in.df.are.0based=FALSE)
  
  
  # Find overlaps based on a range of nucelotides
  overlapping_list <- lapply(seq(nt_start, nt_end, by), function(x){
    return_list <- find_overlaps(psis_inclusion_gr, psis_exclusion_gr,
                                 psis_insensitive_gr, eCLIP_data, nt = x,
                                 bin_width = by, bootstrap = bootstrap,
                                 nbootstrap = nbootstrap) 
    return(return_list)
  })
  
  # Combine specific parts --> probably not the most efficient
  all_plotting_df <- lapply(overlapping_list, function(x){
    return(x$plotting_df)
  })
  
  all_plotting_df <- do.call(rbind, all_plotting_df)
  
  all_overlap_dfs <- lapply(overlapping_list, function(x){
    return(x$list_of_df)
  })
  
  if(bootstrap){
    all_bootstrapping_df <- lapply(overlapping_list, function(x){
      return(x$boostraped_df)
    })
    
    all_bootstrapping_df <- do.call(rbind, all_bootstrapping_df)
    return(list("list_of_df" = all_overlap_dfs, "plotting_df" = all_plotting_df,
                "boostraped_df" = all_bootstrapping_df))
    
  } else {
    return(list("list_of_df" = all_overlap_dfs, "plotting_df" = all_plotting_df))
  }  
}


find_overlaps <- function(psis_inclusion_gr, psis_exclusion_gr,
                          psis_insensitive_gr, eCLIP_data, nt,
                          bin_width, bootstrap = FALSE,
                          nbootstrap = 2){
  
  # Move into separate function - nt will be an argument here
  # make a GRanages object for each of the new dataframes  
  psis_inclusion_gr_up <- flank(psis_inclusion_gr, width = nt, start = TRUE)
  psis_inclusion_gr_up <- resize(psis_inclusion_gr_up, width = bin_width,
                                 fix = "start")
  psis_inclusion_gr_down <- flank(psis_inclusion_gr, width = nt,
                                  start = FALSE)
  psis_inclusion_gr_down <- resize(psis_inclusion_gr_down, width = bin_width,
                                   fix = "end")
  
  psis_exclusion_gr_up <- flank(psis_exclusion_gr, width = nt, start = TRUE)
  psis_exclusion_gr_up <- resize(psis_exclusion_gr_up, width = bin_width,
                                 fix = "start")
  psis_exclusion_gr_down <- flank(psis_exclusion_gr, width = nt, start = FALSE)
  psis_exclusion_gr_down <- resize(psis_exclusion_gr_down, width = bin_width,
                                   fix = "end")
  
  psis_insensitive_gr_up <- flank(psis_insensitive_gr, width = nt, start = TRUE)
  psis_insensitive_gr_up <- resize(psis_insensitive_gr_up, width = bin_width,
                                   fix = "start")
  psis_insensitive_gr_down <- flank(psis_insensitive_gr, width = nt,
                                    start = FALSE)
  psis_insensitive_gr_down <- resize(psis_insensitive_gr_down, width = bin_width,
                                     fix = "end")
  
  testing_list <- c("inclusion_up" = psis_inclusion_gr_up,
                    "inclusion_down" = psis_inclusion_gr_down,
                    "exclusion_up" = psis_exclusion_gr_up,
                    "exclusion_down" = psis_exclusion_gr_down,
                    "insensitive_up" = psis_insensitive_gr_up,
                    "insensitive_down" = psis_insensitive_gr_down)
  
  # find overlap
  all_gr_df <- lapply(testing_list, function(x){
    
    overlap_df <- plyranges::join_overlap_intersect(x, eCLIP_data) %>%
      data.frame() %>%
      dplyr::rename(chr = seqnames ,
                    overlap_start = start, overlap_end = end,
                    overlap_width = width, overlap_strand = strand)
    
    return(overlap_df)
  })
  
  names(all_gr_df) <- names(testing_list)
  
  # Run bootstrapping
  if(bootstrap){
    distance_from_exon <- nt - bin_width
    
    # Do it for inclusion exclusion, up down
    return_vals <- lapply(names(testing_list[1:4]), function(test_type){
      
      # Figure out how many exons to test against
      subsample_to <- length(testing_list[[test_type]]) # This might break
      if (grepl('down', test_type)){
        insensitive_gr <- psis_insensitive_gr_down
      } else if (grepl('up', test_type)){
        insensitive_gr <- psis_insensitive_gr_up
      }
      #print(subsample_to)
      
      # Decide if you want to pull out the overlapping list for the test type
      # compare_to_overlaps <- nrow(all_gr_df[test_type])
      
      single_val <- pbapply::pblapply(1:nbootstrap, function(x){
        
        # Downsample randomly to the correct length
        # Look up this section
        # x = all rows
        # size = subsample_to
        # return = what rows to keep
        set.seed(x)
        keep_rows = sample(x = length(insensitive_gr),
                           size = subsample_to,
                           replace = FALSE)
        insensitive_gr <- insensitive_gr[keep_rows, ] # Double check that this worked
        
        # Run the overlap
        results <- plyranges::join_overlap_intersect(insensitive_gr, eCLIP_data)
        # Return results
        #print(results)
        #print(test_type)
        #print(overlaps)
        #print(x)
        #print(distance_from_exon)
        
        return(data.frame("test" = test_type, overlaps = length(results), x = x,
                          "distance_from_exon" = distance_from_exon)) # add extra column for what you are testing
      })
      
      single_val <- do.call(rbind, single_val)
      
    })
    return_vals <- do.call(rbind, return_vals)
  }
  # If boostrap{}
  # Downsample insensitive to # of exclusion x1000
  # Return only the length of the overlaps
  # Downsample insensitive to # of inclusion x1000
  # Return only the lenght of the overlaps
  # Final df from this section
  # col 1 up/down, col 2 inclusion/exclusion, col 3 # overlaps
  # 4000 rows
  # Add this as new item to the return list
  
  # make the plotting df that calculates normalized values
  return_df <- lapply(names(testing_list), function(x){
    number_overlapping_peaks <- nrow(all_gr_df[[x]])
    number_exons <- length(testing_list[[x]])
    normalized_count <- number_overlapping_peaks / number_exons * 100
    return(data.frame("exon_type" = x,
                      "distance_from_exon" = nt - bin_width,
                      "number_overlapping_peaks" = number_overlapping_peaks, # What we care about for insensitive and inclusion
                      "number_of_exons" = number_exons,
                      "normalized_count" = normalized_count))
  })
  
  return_df <- do.call(rbind, return_df)
  
  names(all_gr_df) <- paste0(names(all_gr_df), "_", nt)
  
  if(bootstrap){
    return(list("list_of_df" = all_gr_df, "plotting_df" = return_df,
                "boostraped_df" = return_vals))
  } else {
    return(list("list_of_df" = all_gr_df, "plotting_df" = return_df))
  }
}

calculate_pval_df <- function(SE_bootstrap) {
  bootstrap_df_summary <- SE_bootstrap$boostraped_df %>%
    tidyr::unite(combined, test, distance_from_exon, sep = "_", remove = FALSE) %>%
    dplyr::select(c(combined, overlaps)) %>%
    dplyr::group_by(combined) %>%
    dplyr::summarise(MEAN = mean(overlaps), SD = sd(overlaps))
  
  plotting_df <- SE_bootstrap$plotting_df %>%
    unite(combined, exon_type, distance_from_exon, sep = "_", remove = FALSE)
  
  combined_df <- inner_join(bootstrap_df_summary, plotting_df)
  
  pval_df <- combined_df %>%
    mutate(pval = scales::scientific(pnorm(number_overlapping_peaks, MEAN, 
                                           SD, lower.tail = FALSE, 
                                           log.p = FALSE), digits = 3),
           distance_from_exon = distance_from_exon + 50)
  
  return(pval_df)
}

# Prep input data
SE <- read_delim(here("revision/SE.MATS.JC.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
MXE <- read_delim(here("revision/MXE.MATS.JC.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)

# write function to "clean-up" the data and name samples. 
Clean_1 <- function(data,event)
{data <- data %>% mutate(Event = event) %>% 
  dplyr::rename(ID = ID...1 , NT_included = IJC_SAMPLE_1, NT_skipped = SJC_SAMPLE_1, KD_included = IJC_SAMPLE_2, KD_skipped = SJC_SAMPLE_2) %>% 
  dplyr::select(GeneID,ID,geneSymbol, Event, PValue, FDR, IncLevelDifference, NT_included, NT_skipped, KD_included, KD_skipped, chr, strand, exonStart_0base, exonEnd)
return (data)}


Clean_2 <- function(data,event)
{data <- data %>% mutate(Event = event) %>% 
  dplyr::rename(ID = ID...1 , NT_included = IJC_SAMPLE_1, NT_skipped = SJC_SAMPLE_1, KD_included = IJC_SAMPLE_2, KD_skipped = SJC_SAMPLE_2, exonStart_0base = `1stExonStart_0base`, exonEnd = `1stExonEnd`) %>% 
  dplyr::select(GeneID,ID,geneSymbol, Event, PValue, FDR, IncLevelDifference, NT_included, NT_skipped, KD_included, KD_skipped, chr, strand, exonStart_0base, exonEnd)
return (data)}

# write function to "clean-up" the data and name samples. 
Clean_3 <- function(data,event)
{data <- data %>% mutate(Event = event) %>% 
  dplyr::rename(ID = ID...1 , NT_included = IJC_SAMPLE_1, NT_skipped = SJC_SAMPLE_1, KD_included = IJC_SAMPLE_2, KD_skipped = SJC_SAMPLE_2, exonStart_0base = `2ndExonStart_0base`, exonEnd = `2ndExonEnd`) %>% 
  dplyr::select(GeneID,ID,geneSymbol, Event, PValue, FDR, IncLevelDifference, NT_included, NT_skipped, KD_included, KD_skipped, chr, strand, exonStart_0base, exonEnd)
return (data)}

SE_All  <- Clean_1(SE, "SE")
MXE_1_All <- Clean_2(MXE, "MXE")
MXE_2_All <- Clean_3(MXE, "MXE")

# keep only events with 20+ informative reads for each replicate (this results in removing events with low number of total reads, new total is 38,599)
SE_All <- SE_All %>% 
  separate(., col = NT_included, into = c('NT_I1', 'NT_I2', 'NT_I3', 'NT_I4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = NT_skipped, into = c('NT_S1', 'NT_S2', 'NT_S3', 'NT_S4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = KD_included, into = c('KD_I1', 'KD_I2', 'KD_I3', 'KD_I4'), sep = ',', remove = T, convert = T) %>%
  separate(., col = KD_skipped, into = c('KD_S1', 'KD_S2', 'KD_S3', 'KD_S4'), sep = ',', remove = T, convert = T) %>% 
  mutate(., NT_1_counts = NT_I1 + NT_S1) %>%
  mutate(., NT_2_counts = NT_I2 + NT_S2) %>%
  mutate(., NT_3_counts = NT_I3 + NT_S3) %>%
  mutate(., NT_4_counts = NT_I4 + NT_S4) %>%
  mutate(., KD_1_counts = KD_I1 + KD_S1) %>%
  mutate(., KD_2_counts = KD_I2 + KD_S2) %>%
  mutate(., KD_3_counts = KD_I3 + KD_S3) %>%
  mutate(., KD_4_counts = KD_I4 + KD_S4) %>%
  filter(., NT_1_counts >= 20 & NT_2_counts >= 20 & NT_3_counts >= 20 & NT_4_counts >= 20 & KD_1_counts >= 20 & KD_2_counts >= 20 & KD_3_counts >= 20 & KD_4_counts >= 20) %>%
  dplyr::select(ID, geneSymbol, chr, strand, exonStart_0base, exonEnd, FDR, IncLevelDifference) 



MXE_1_All <- MXE_1_All %>% 
  separate(., col = NT_included, into = c('NT_I1', 'NT_I2', 'NT_I3', 'NT_I4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = NT_skipped, into = c('NT_S1', 'NT_S2', 'NT_S3', 'NT_S4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = KD_included, into = c('KD_I1', 'KD_I2', 'KD_I3', 'KD_I4'), sep = ',', remove = T, convert = T) %>%
  separate(., col = KD_skipped, into = c('KD_S1', 'KD_S2', 'KD_S3', 'KD_S4'), sep = ',', remove = T, convert = T) %>% 
  mutate(., NT_1_counts = NT_I1 + NT_S1) %>%
  mutate(., NT_2_counts = NT_I2 + NT_S2) %>%
  mutate(., NT_3_counts = NT_I3 + NT_S3) %>%
  mutate(., NT_4_counts = NT_I4 + NT_S4) %>%
  mutate(., KD_1_counts = KD_I1 + KD_S1) %>%
  mutate(., KD_2_counts = KD_I2 + KD_S2) %>%
  mutate(., KD_3_counts = KD_I3 + KD_S3) %>%
  mutate(., KD_4_counts = KD_I4 + KD_S4) %>%
  filter(., NT_1_counts >= 20 & NT_2_counts >= 20 & NT_3_counts >= 20 & NT_4_counts >= 20 & KD_1_counts >= 20 & KD_2_counts >= 20 & KD_3_counts >= 20 & KD_4_counts >= 20) %>%
  dplyr::select(ID, geneSymbol, chr, strand, exonStart_0base, exonEnd, FDR, IncLevelDifference) 


MXE_2_All <- MXE_2_All %>% 
  separate(., col = NT_included, into = c('NT_I1', 'NT_I2', 'NT_I3', 'NT_I4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = NT_skipped, into = c('NT_S1', 'NT_S2', 'NT_S3', 'NT_S4'), sep = ',', remove = T, convert = T) %>% 
  separate(., col = KD_included, into = c('KD_I1', 'KD_I2', 'KD_I3', 'KD_I4'), sep = ',', remove = T, convert = T) %>%
  separate(., col = KD_skipped, into = c('KD_S1', 'KD_S2', 'KD_S3', 'KD_S4'), sep = ',', remove = T, convert = T) %>% 
  mutate(., NT_1_counts = NT_I1 + NT_S1) %>%
  mutate(., NT_2_counts = NT_I2 + NT_S2) %>%
  mutate(., NT_3_counts = NT_I3 + NT_S3) %>%
  mutate(., NT_4_counts = NT_I4 + NT_S4) %>%
  mutate(., KD_1_counts = KD_I1 + KD_S1) %>%
  mutate(., KD_2_counts = KD_I2 + KD_S2) %>%
  mutate(., KD_3_counts = KD_I3 + KD_S3) %>%
  mutate(., KD_4_counts = KD_I4 + KD_S4) %>%
  filter(., NT_1_counts >= 20 & NT_2_counts >= 20 & NT_3_counts >= 20 & NT_4_counts >= 20 & KD_1_counts >= 20 & KD_2_counts >= 20 & KD_3_counts >= 20 & KD_4_counts >= 20) %>%
  dplyr::select(ID, geneSymbol, chr, strand, exonStart_0base, exonEnd, FDR, IncLevelDifference) 

igg_names <- names(full_anno)[grepl("IgG", names(sig_reads))]

for(i in igg_names){
  full_anno[[i]] <- NULL
}

all_granges <- lapply(full_anno, function(x){
  
  grange_obj <- makeGRangesFromDataFrame(x,
                                          keep.extra.columns=TRUE,
                                          ignore.strand=FALSE,
                                          seqinfo=NULL,
                                          seqnames.field=c("seqnames"),
                                          start.field="start",
                                          end.field=c("end"),
                                          strand.field="strand",
                                          starts.in.df.are.0based=FALSE)
  

  return(grange_obj)
})

merge_with_by <- function(x, y) {
  merge(x, y, by = c("exon_type", "distance_from_exon"))
}


exon_levels <- seq(0, 300, by = 50)

SE_bootstrap_all <- lapply(names(all_granges), function(x){
  print(x)
  SE_bootstrap <- peak_exon_overlap(SE_All, all_granges[[x]],
                                    nt_start = 0,
                                    nt_end = 300, by = 50,
                                    bootstrap = TRUE,
                                    nbootstrap = 1000)
  
  SE_bootstrap <- calculate_pval_df(SE_bootstrap)
  SE_bootstrap <- SE_bootstrap %>%
    dplyr::mutate(padj = p.adjust(pval, method = "fdr")) %>%
    dplyr::select(exon_type, distance_from_exon, normalized_count,
                  padj)
  colnames(SE_bootstrap) <- c("exon_type", "distance_from_exon",
                              paste0(x, "_normalized_count"),
                              paste0(x, "_padj"))
    return(SE_bootstrap)
})

SE_bootstrap_all <- base::Reduce(merge_with_by,
                                 SE_bootstrap_all)

SE_bootstrap_all$distance_from_exon <- factor(SE_bootstrap_all$distance_from_exon,
                                              levels = exon_levels)

SE_bootstrap_all <- SE_bootstrap_all %>%
  dplyr::arrange(exon_type, distance_from_exon)



MXE1_bootstrap_all <- lapply(names(all_granges), function(x){
  MXE1_bootstrap <- peak_exon_overlap(MXE_1_All, all_granges[[x]],
                                      nt_start = 0,
                                      nt_end = 300, by = 50,
                                      bootstrap = TRUE,
                                      nbootstrap = 1000) 
  
  MXE1_bootstrap <- calculate_pval_df(MXE1_bootstrap)
  MXE1_bootstrap <- MXE1_bootstrap %>%
    dplyr::mutate(padj = p.adjust(pval, method = "fdr")) %>%
    dplyr::select(exon_type, distance_from_exon, normalized_count,
                  pval)
  
  colnames(MXE1_bootstrap) <- c("exon_type", "distance_from_exon",
                              paste0(x, "_normalized_count"),
                              paste0(x, "_padj"))
  return(MXE1_bootstrap)
})

MXE1_bootstrap_all <- base::Reduce(merge_with_by,
                                   MXE1_bootstrap_all)

MXE1_bootstrap_all$distance_from_exon <- factor(MXE1_bootstrap_all$distance_from_exon,
                                              levels = exon_levels)

MXE1_bootstrap_all <- MXE1_bootstrap_all %>%
  dplyr::arrange(exon_type, distance_from_exon)

MXE2_bootstrap_all <- lapply(names(all_granges), function(x){
  MXE2_bootstrap <- peak_exon_overlap(MXE_2_All, all_granges[[x]],
                                      nt_start = 0,
                                      nt_end = 300, by = 50,
                                      bootstrap = TRUE,
                                      nbootstrap = 1000)
  
  MXE2_bootstrap <- calculate_pval_df(MXE2_bootstrap)
  MXE2_bootstrap <- MXE2_bootstrap %>%
    dplyr::mutate(padj = p.adjust(pval, method = "fdr")) %>%
    dplyr::select(exon_type, distance_from_exon, normalized_count,
                  pval)
  colnames(MXE2_bootstrap) <- c("exon_type", "distance_from_exon",
                                paste0(x, "_normalized_count"),
                                paste0(x, "_pval"))
  
  return(MXE2_bootstrap)
  
})

MXE2_bootstrap_all <- base::Reduce(merge_with_by,
                                   MXE2_bootstrap_all)

MXE2_bootstrap_all$distance_from_exon <- factor(MXE2_bootstrap_all$distance_from_exon,
                                              levels = exon_levels)

MXE2_bootstrap_all <- MXE2_bootstrap_all %>%
  dplyr::arrange(exon_type, distance_from_exon)

bootstrap_wb = openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = bootstrap_wb, sheetName = "SE")
openxlsx::writeData(wb = bootstrap_wb, sheet = "SE", x = SE_bootstrap_all)

openxlsx::addWorksheet(wb = bootstrap_wb, sheetName = "MXE1")
openxlsx::writeData(wb = bootstrap_wb, sheet = "MXE1", x = MXE1_bootstrap_all)

openxlsx::addWorksheet(wb = bootstrap_wb, sheetName = "MXE2")
openxlsx::writeData(wb = bootstrap_wb, sheet = "MXE2", x = MXE2_bootstrap_all)

save_name <- paste0("bootstrapped_p_vals_log_", fold_enrichment_cutoff, ".xlsx")

openxlsx::saveWorkbook(wb = bootstrap_wb,
                       file = file.path("revision", save_name),
                       overwrite = TRUE)
