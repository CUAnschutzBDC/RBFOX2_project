SE <- read_delim("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/SE.MATS.JC.txt", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)

Clean_2 <- function(data,event)
{data <- data %>% mutate(Event = event) %>% 
  dplyr::rename(ObOb_included = IJC_SAMPLE_1, ObOb_skipped = SJC_SAMPLE_1, NZO_included = IJC_SAMPLE_2, NZO_skipped = SJC_SAMPLE_2) %>% 
  dplyr::select(GeneID,geneSymbol, Event, PValue, FDR, IncLevelDifference, ObOb_included, ObOb_skipped, NZO_included, NZO_skipped, chr, strand, exonStart_0base, exonEnd)
return (data)}

SE_All <- Clean_2(SE, "SE")

# keep only events with 20+ informative reads for each replicate

SE_All <- SE_All %>% 
  separate(., col = ObOb_included, into = c('ObOb_I1', 'ObOb_I2', 'ObOb_I3',"ObOb_I4","ObOb_I5"), sep = ',', remove = T, convert = T) %>% 
  separate(., col = ObOb_skipped, into = c('ObOb_S1', 'ObOb_S2', 'ObOb_S3', "ObOb_S4", "ObOb_S5"), sep = ',', remove = T, convert = T) %>% 
  separate(., col = NZO_included, into = c('NZO_I1', 'NZO_I2', 'NZO_I3',"NZO_I4", "NZO_I5"), sep = ',', remove = T, convert = T) %>%
  separate(., col = NZO_skipped, into = c('NZO_S1', 'NZO_S2', 'NZO_S3', "NZO_S4", "NZO_S5"), sep = ',', remove = T, convert = T) %>% 
  mutate(., ObOb_1_counts = ObOb_I1 + ObOb_S1) %>%
  mutate(., ObOb_2_counts = ObOb_I2 + ObOb_S2) %>%
  mutate(., ObOb_3_counts = ObOb_I3 + ObOb_S3) %>%
  mutate(., ObOb_4_counts = ObOb_I4 + ObOb_S4) %>%
  mutate(., ObOb_5_counts = ObOb_I5 + ObOb_S5) %>%
  mutate(., NZO_1_counts = NZO_I1 + NZO_S1) %>%
  mutate(., NZO_2_counts = NZO_I2 + NZO_S2) %>%
  mutate(., NZO_3_counts = NZO_I3 + NZO_S3) %>%
  mutate(., NZO_4_counts = NZO_I4 + NZO_S4) %>%
  mutate(., NZO_5_counts = NZO_I5 + NZO_S5) %>%
  filter(., ObOb_1_counts >= 20 & ObOb_2_counts >= 20 & ObOb_3_counts >= 20 & ObOb_4_counts >= 20 & ObOb_5_counts >= 20 & NZO_1_counts >= 20 & NZO_2_counts >= 20 & NZO_3_counts >= 20 & NZO_4_counts >= 20 & NZO_5_counts >= 20) %>%
  select(geneSymbol, chr, strand, exonStart_0base, exonEnd, FDR, IncLevelDifference)

# separate sensitive (inclusions and exclusion) and insensitive exons

psis_inclusion <- filter(SE_All, FDR < 0.05 & IncLevelDifference < 0.01) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1")
psis_exclusion <- filter(SE_All, FDR < 0.05 & IncLevelDifference > 0.01) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1")
psis_insensitive <- filter(SE_All, FDR >= 0.05) %>% filter(., chr!= "chrGL456233.1") %>% filter(., chr!= "chrGL456216.1") %>% filter(., chr!= "chrJH584304.1") # these crhomosome keys were throwing errors in the python chunk below, these events were removed chrGL456216.1

# get coordinates for sequences 

psis_inclusion %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/T2D_inclusionexon.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)
psis_exclusion %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/T2D_exclusionexon.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)
psis_insensitive %>% select(., -c(FDR, IncLevelDifference)) %>%
  write.table(., file = 'Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/T2D_insensitive.coordinates.txt', sep = '\t', row.names = F, col.names = F, quote = F)