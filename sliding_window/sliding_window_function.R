library(here)
library(GenomicRanges)
library(tidyverse)
library(plyranges)

# Change nt to nt_range

peak_exon_overlap <- function(rMATS_data, eCLIP_data,
                              nt_start = 0,
                              nt_end = 300, by = 50){
  
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
                                 bin_width = by)
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
  
  return(list("list_of_df" = all_overlap_dfs, "plotting_df" = all_plotting_df))
  
}

find_overlaps <- function(psis_inclusion_gr, psis_exclusion_gr,
                          psis_insensitive_gr, eCLIP_data, nt,
                          bin_width){
  
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
                    overlap_width = width, overlap_strand = strand,
                    splicing_ID = ID,   peak_ID =...2)
    
    return(overlap_df)
  })
  
  # make the plotting df that calculates normalized values
  return_df <- lapply(names(testing_list), function(x){
    number_overlapping_peaks <- nrow(all_gr_df[[x]])
    number_exons <- length(testing_list[[x]])
    normalized_count <- number_overlapping_peaks / number_exons * 100
    return(data.frame("exon_type" = x,
                      "distance_from_exon" = nt - bin_width,
                      "number_overlapping_peaks" = number_overlapping_peaks,
                      "number_of_exons" = number_exons,
                      "normalized_count" = normalized_count))
  })
  
  return_df <- do.call(rbind, return_df)
  
  names(all_gr_df) <- paste0(names(all_gr_df), "_", nt)
  
  return(list("list_of_df" = all_gr_df, "plotting_df" = return_df))
}