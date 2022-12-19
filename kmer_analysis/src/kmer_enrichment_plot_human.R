Inclusion_Up <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/inclusionexons.upstream.kmercounts.txt", skip = 1, col.names = c("kmer", "i_up_kmer","i_up_total"))
Inclusion_Down <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/inclusionexons.downstream.kmercounts.txt", skip = 1, col.names = c("kmer", "i_down_kmer","i_down_total"))
Exclusion_Up <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/exclusionexons.upstream.kmercounts.txt", skip = 1, col.names = c("kmer", "e_up_kmer","e_up_total"))
Exclusion_Down <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/GSE183247_rMATS_output_v4.0.2/exclusionexons.downstream.kmercounts.txt", skip = 1, col.names = c("kmer", "e_down_kmer","e_down_total"))
Insense_Up <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/insensitives.upstream.kmercounts.txt", skip = 1, col.names = c("kmer", "insensitive_up_kmer","insensitive_up_total"))
Insense_Down <- read.table("Dropbox/Manuscripts/RBFOX2 Manuscript/Data/rMATS_v4.0.2_output_GSE164416/insensitives.upstream.kmercounts.txt", skip = 1, col.names = c("kmer", "insensitive_down_kmer","insensitive_down_total"))

# Join upstream and downstream kmer counts
# Calculate enrichmnets and p-values

Inc_kmers_up <- left_join(Inclusion_Up, Insense_Up)
Inc_kmers_up_enrichment <- Inc_kmers_up %>% mutate(sensitive_upstream = i_up_kmer/(i_up_kmer+ i_up_total), insensitive_upstream = insensitive_up_kmer/(insensitive_up_kmer+insensitive_up_total)) %>%
  mutate (enrichment = (sensitive_upstream/insensitive_upstream)) %>% 
  mutate(log2_enrichment_inclusion_up = log2(enrichment)) %>%
  rowwise() %>% 
  mutate(., pval_up = binom.test(i_up_kmer, sum(i_up_kmer+ i_up_total), insensitive_upstream)$p.value) %>% 
  mutate(., adj_pval_up = (pval_up * nrow(Inc_kmers_up))) %>%
  mutate(., adj_pval_up = ifelse(adj_pval_up < 1, adj_pval_up, 1)) %>%
  mutate(., pvalue = ifelse(pval_up < 0.05, "<0.05",">=0.05")) %>%
  mutate(., label = ifelse(kmer == "GCAUG", "GCAUG",""))

Inc_kmers_down <- left_join(Inclusion_Down, Insense_Down)
Inc_kmers_down_enrichment <- Inc_kmers_down %>% mutate(sensitive_downstream = i_down_kmer/(i_down_kmer+ i_down_total), insensitive_downstream = insensitive_down_kmer/(insensitive_down_kmer+insensitive_down_total)) %>%
  mutate (enrichment = (sensitive_downstream/insensitive_downstream)) %>% 
  mutate(log2_enrichment_inclusion_down = log2(enrichment)) %>%
  rowwise() %>% 
  mutate(., pval_down = binom.test(i_down_kmer, sum(i_down_kmer+ i_down_total),insensitive_downstream)$p.value) %>% 
  mutate(., adj_pval_down = (pval_down * nrow(Inc_kmers_down))) %>%
  mutate(., adj_pval_down = ifelse(adj_pval_down < 1, adj_pval_down, 1)) %>%
  mutate(., pvalue = ifelse(pval_down < 0.05, "<0.05",">=0.05")) %>%
  mutate(., label = ifelse(kmer == "GCAUG", "GCAUG",""))

Ex_kmers_up <- left_join(Exclusion_Up, Insense_Up)
Ex_kmers_up_enrichment <- Ex_kmers_up %>% mutate(sensitive_upstream = e_up_kmer/(e_up_kmer+ e_up_total), insensitive_upstream = insensitive_up_kmer/(insensitive_up_kmer+insensitive_up_total)) %>%
  mutate (enrichment = (sensitive_upstream/insensitive_upstream)) %>% 
  mutate(log2_enrichment_exclusion_up = log2(enrichment)) %>%
  rowwise() %>% 
  mutate(., pval_up = binom.test(e_up_kmer, sum(e_up_kmer+ e_up_total), insensitive_upstream)$p.value) %>% 
  mutate(., adj_pval_up = (pval_up * nrow(Ex_kmers_up))) %>%
  mutate(., adj_pval_up = ifelse(adj_pval_up < 1, adj_pval_up, 1)) %>%
  mutate(., pvalue = ifelse(pval_up < 0.05, "<0.05",">=0.05")) %>%
  mutate(., label = ifelse(kmer == "GCAUG", "GCAUG",""))

Ex_kmers_down <- left_join(Exclusion_Down, Insense_Down)
Ex_kmers_down_enrichment <- Ex_kmers_down %>% mutate(sensitive_downstream = e_down_kmer/(e_down_kmer+ e_down_total), insensitive_downstream = insensitive_down_kmer/(insensitive_down_kmer+insensitive_down_total)) %>%
  mutate (enrichment = (sensitive_downstream/insensitive_downstream)) %>% 
  mutate(log2_enrichment_exclusion_down = log2(enrichment)) %>%  
  rowwise() %>% 
  mutate(., pval_down = binom.test(e_down_kmer, sum(e_down_kmer+ e_down_total), insensitive_downstream)$p.value) %>% 
  mutate(., adj_pval_down = (pval_down* nrow(Ex_kmers_down))) %>%
  mutate(., adj_pval_down = ifelse(adj_pval_down < 1, adj_pval_down, 1)) %>%
  mutate(., pvalue = ifelse(pval_down < 0.05, "<0.05",">=0.05")) %>%
  mutate(., label = ifelse(kmer == "GCAUG", "RBFOX2", ""))
#mutate(., label = ifelse(kmer == "GCAUG", "GCAUG-RBFOX2", ifelse(kmer == "UUUUU","UUUUU-HNRNPC",ifelse(kmer == "GGGGG","GGGGG-FUS", ifelse(kmer == "CACAC","CACAC-HNRNPL","")))))

# Plot
options(ggrepel.max.overlaps = Inf)

p1 <- ggplot(Inc_kmers_up_enrichment, aes(x = log2_enrichment_inclusion_up, y = -log10(pval_up), col = pvalue, label = label)) +
  geom_point() +
  labs(title = "Upstream Inclusion Sensitive Kmers") +
  geom_text_repel() +
  scale_color_manual(values=c("firebrick4", "black")) +
  theme_minimal()
p2 <- ggplot(Inc_kmers_down_enrichment, aes(x = log2_enrichment_inclusion_down, y = -log10(pval_down), col = pvalue, label = label)) +
  geom_point() +
  labs(title = "Downstream Inclusion Sensitive Kmers") +
  geom_text_repel() +
  scale_color_manual(values=c("firebrick4", "black")) +
  theme_minimal()
p3 <- ggplot(Ex_kmers_up_enrichment, aes(x = log2_enrichment_exclusion_up, y = -log10(pval_up), col = pvalue, label = label)) +
  geom_point() +
  labs(title = "Upstream Exclusion Sensitive Kmers") +
  geom_text_repel() +
  scale_color_manual(values=c("firebrick4", "black")) +
  theme_minimal()

Ex_kmers_down_enrichment <- Ex_kmers_down_enrichment %>% filter(log2_enrichment_exclusion_down > 0)
p4 <- ggplot(Ex_kmers_down_enrichment, aes(x = log2_enrichment_exclusion_down, y = -log10(pval_down), col = pvalue, label = label)) +
  geom_point(size = 1) +
  labs(title = "Downstream Exclusion Sensitive Kmers") +
  scale_color_manual(values=c("steelblue", "black")) +
  theme_classic() +
  #xlim(0, 1) +
  #ylim(0, 6) +
  theme(text = element_text(size = 20)) +
  geom_label_repel(min.segment.length = 0, 
                   nudge_x = .3,nudge_y = 1.8)


# use the miltiplot function to plot multiple charts in a signle figure
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))}
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))}}}

multiplot(p1, p2, p3, p4, cols=2)
