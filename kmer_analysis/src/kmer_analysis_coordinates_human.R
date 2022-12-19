SE_MATS_JCEC_GSE164416 <- read_delim("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/SE.MATS.JCEC.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

Clean_2 <- function(data,event)
{data <- data %>% mutate(Event = event) %>% 
  dplyr::rename(ND_included = IJC_SAMPLE_1, ND_skipped = SJC_SAMPLE_1, T2D_included = IJC_SAMPLE_2, T2D_skipped = SJC_SAMPLE_2) %>% 
  dplyr::select(GeneID,geneSymbol, Event, PValue, FDR, IncLevelDifference, ND_included, ND_skipped, T2D_included, T2D_skipped, chr, strand, exonStart_0base, exonEnd)
return (data)}

SE_All <- Clean_2(SE_MATS_JCEC_GSE164416, "SE")

# keep only events with 20+ informative reads for each replicate

SE_All <- SE_All %>% 
  separate(., col = ND_included, into = c("ND_I1","ND_I2","ND_I3","ND_I4","ND_I5","ND_I6","ND_I7","ND_I8","ND_I9","ND_I10","ND_I11","ND_I12","ND_I13","ND_I14","ND_I15","ND_I16","ND_I17","ND_I18"), sep = ',', remove = T, convert = T) %>% 
  separate(., col = ND_skipped, into = c("ND_S1","ND_S2","ND_S3","ND_S4","ND_S5","ND_S6","ND_S7","ND_S8","ND_S9","ND_S10","ND_S11","ND_S12","ND_S13","ND_S14","ND_S15","ND_S16","ND_S17","ND_S18"), sep = ',', remove = T, convert = T) %>% 
  separate(., col = T2D_included, into = c("T2D_I1","T2D_I2","T2D_I3","T2D_I4","T2D_I5","T2D_I6","T2D_I7","T2D_I8","T2D_I9","T2D_I10","T2D_I11","T2D_I12","T2D_I13","T2D_I14","T2D_I15","T2D_I16","T2D_I17","T2D_I18","T2D_I19","T2D_I20","T2D_I21","T2D_I22","T2D_I23","T2D_I24","T2D_I25","T2D_I26","T2D_I27","T2D_I28","T2D_I29","T2D_I30","T2D_I31","T2D_I32","T2D_I33","T2D_I34","T2D_I35","T2D_I36","T2D_I37","T2D_I38","T2D_I39"), sep = ',', remove = T, convert = T) %>%
  separate(., col = T2D_skipped, into = c("T2D_S1","T2D_S2","T2D_S3","T2D_S4","T2D_S5","T2D_S6","T2D_S7","T2D_S8","T2D_S9","T2D_S10","T2D_S11","T2D_S12","T2D_S13","T2D_S14","T2D_S15","T2D_S16","T2D_S17","T2D_S18","T2D_S19","T2D_S20","T2D_S21","T2D_S22","T2D_S23","T2D_S24","T2D_S25","T2D_S26","T2D_S27","T2D_S28","T2D_S29","T2D_S30","T2D_S31","T2D_S32","T2D_S33","T2D_S34","T2D_S35","T2D_S36","T2D_S37","T2D_S38","T2D_S39"), sep = ',', remove = T, convert = T) %>% 
  mutate(., ND_1_counts = ND_I1 + ND_S1) %>%
  mutate(., ND_2_counts = ND_I2 + ND_S2) %>%
  mutate(., ND_3_counts = ND_I3 + ND_S3) %>%
  mutate(., ND_4_counts = ND_I4 + ND_S4) %>%
  mutate(., ND_5_counts = ND_I5 + ND_S5) %>%
  mutate(., ND_6_counts = ND_I6 + ND_S6) %>%
  mutate(., ND_7_counts = ND_I7 + ND_S7) %>%
  mutate(., ND_8_counts = ND_I8 + ND_S8) %>%
  mutate(., ND_9_counts = ND_I9 + ND_S9) %>%
  mutate(., ND_10_counts = ND_I10 + ND_S10) %>%
  mutate(., ND_11_counts = ND_I11 + ND_S11) %>%
  mutate(., ND_12_counts = ND_I12 + ND_S12) %>%
  mutate(., ND_13_counts = ND_I13 + ND_S13) %>%
  mutate(., ND_14_counts = ND_I14 + ND_S14) %>%
  mutate(., ND_15_counts = ND_I15 + ND_S15) %>%
  mutate(., ND_16_counts = ND_I16 + ND_S16) %>%
  mutate(., ND_17_counts = ND_I17 + ND_S17) %>%
  mutate(., ND_18_counts = ND_I18 + ND_S18) %>%
  mutate(., T2D_1_counts = T2D_I1 + T2D_S1) %>%
  mutate(., T2D_2_counts = T2D_I2 + T2D_S2) %>%
  mutate(., T2D_3_counts = T2D_I3 + T2D_S3) %>%
  mutate(., T2D_4_counts = T2D_I4 + T2D_S4) %>%
  mutate(., T2D_5_counts = T2D_I5 + T2D_S5) %>%
  mutate(., T2D_6_counts = T2D_I6 + T2D_S6) %>%
  mutate(., T2D_7_counts = T2D_I7 + T2D_S7) %>%
  mutate(., T2D_8_counts = T2D_I8 + T2D_S8) %>%
  mutate(., T2D_9_counts = T2D_I9 + T2D_S9) %>%
  mutate(., T2D_10_counts = T2D_I10 + T2D_S10) %>%
  mutate(., T2D_11_counts = T2D_I11 + T2D_S11) %>%
  mutate(., T2D_12_counts = T2D_I12 + T2D_S12) %>%
  mutate(., T2D_13_counts = T2D_I13 + T2D_S13) %>%
  mutate(., T2D_14_counts = T2D_I14 + T2D_S14) %>%
  mutate(., T2D_15_counts = T2D_I15 + T2D_S15) %>%
  mutate(., T2D_16_counts = T2D_I16 + T2D_S16) %>%
  mutate(., T2D_17_counts = T2D_I17 + T2D_S17) %>%
  mutate(., T2D_18_counts = T2D_I18 + T2D_S18) %>%
  mutate(., T2D_19_counts = T2D_I19 + T2D_S19) %>%
  mutate(., T2D_20_counts = T2D_I20 + T2D_S20) %>%
  mutate(., T2D_21_counts = T2D_I21 + T2D_S21) %>%
  mutate(., T2D_22_counts = T2D_I22 + T2D_S22) %>%
  mutate(., T2D_23_counts = T2D_I23 + T2D_S23) %>%
  mutate(., T2D_24_counts = T2D_I24 + T2D_S24) %>%
  mutate(., T2D_25_counts = T2D_I25 + T2D_S25) %>%
  mutate(., T2D_26_counts = T2D_I26 + T2D_S26) %>%
  mutate(., T2D_27_counts = T2D_I27 + T2D_S27) %>%
  mutate(., T2D_28_counts = T2D_I28 + T2D_S28) %>%
  mutate(., T2D_29_counts = T2D_I29 + T2D_S29) %>%
  mutate(., T2D_30_counts = T2D_I30 + T2D_S30) %>%
  mutate(., T2D_31_counts = T2D_I31 + T2D_S31) %>%
  mutate(., T2D_32_counts = T2D_I32 + T2D_S32) %>%
  mutate(., T2D_33_counts = T2D_I33 + T2D_S33) %>%
  mutate(., T2D_34_counts = T2D_I34 + T2D_S34) %>%
  mutate(., T2D_35_counts = T2D_I35 + T2D_S35) %>%
  mutate(., T2D_36_counts = T2D_I36 + T2D_S36) %>%
  mutate(., T2D_37_counts = T2D_I37 + T2D_S37) %>%
  mutate(., T2D_38_counts = T2D_I38 + T2D_S38) %>%
  mutate(., T2D_39_counts = T2D_I39 + T2D_S39) %>%
  filter(., ND_1_counts >= 20 & ND_2_counts >= 20 & ND_3_counts >= 20 & ND_4_counts >= 20 & ND_5_counts >= 20 & ND_6_counts >= 20 & ND_7_counts >= 20 & ND_8_counts >= 20 & ND_9_counts >= 20 & ND_10_counts >= 20 & ND_11_counts >= 20 & ND_12_counts >= 20 & ND_13_counts >= 20 & ND_14_counts >= 20 & ND_15_counts >= 20 & ND_16_counts >= 20 & ND_17_counts >= 20 & ND_18_counts >= 20 & T2D_1_counts >= 20 & T2D_2_counts >= 20 & T2D_3_counts >= 20 & T2D_4_counts >= 20 & T2D_5_counts >= 20 & T2D_6_counts >= 20 & T2D_7_counts >= 20 & T2D_8_counts >= 20 & T2D_9_counts >= 20 & T2D_10_counts >= 20 & T2D_11_counts >= 20 & T2D_12_counts >= 20 & T2D_13_counts >= 20 & T2D_14_counts >= 20 & T2D_15_counts >= 20 & T2D_16_counts >= 20 & T2D_17_counts >= 20 & T2D_18_counts >= 20 & T2D_19_counts >= 20 & T2D_20_counts >= 20 & T2D_21_counts >= 20 & T2D_22_counts >= 20 & T2D_23_counts >= 20 & T2D_24_counts >= 20 & T2D_25_counts >= 20 & T2D_26_counts >= 20 & T2D_27_counts >= 20 & T2D_28_counts >= 20 & T2D_29_counts >= 20 & T2D_30_counts >= 20 & T2D_31_counts >= 20 & T2D_32_counts >= 20 & T2D_33_counts >= 20 & T2D_34_counts >= 20 & T2D_35_counts >= 20 & T2D_36_counts >= 20 & T2D_37_counts >= 20 & T2D_38_counts >= 20 & T2D_39_counts >= 20) %>%
  select(geneSymbol, chr, strand, exonStart_0base, exonEnd, PValue, FDR, IncLevelDifference)

# separate sensitive (inclusions and exclusion) and insensitive exons

psis_inclusion <- filter(SE_All, PValue < 0.05 & IncLevelDifference > 0.01) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1")
psis_exclusion <- filter(SE_All, PValue < 0.05 & IncLevelDifference < 0.01) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1")
psis_insensitive <- filter(SE_All, PValue >= 0.05) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1") # these crhomosome keys were throwing errors in the python chunk below, these events were removed chrGL456216.1

# get coordinates for sequences 

psis_inclusion %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/Human_T2D_inclusionexon.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)
psis_exclusion %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/Human_T2D_exclusionexon.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)
psis_insensitive %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/Human_T2D_insensitive.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)