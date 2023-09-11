```{r, bootstrapping function}
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
      
      single_val <- lapply(1:nbootstrap, function(x){
        
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
```