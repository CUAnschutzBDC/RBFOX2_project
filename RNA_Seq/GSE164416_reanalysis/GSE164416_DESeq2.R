#load individual .txt files
GSE16446_count_table <- read_delim("GSE16446_count_table.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

cts <- GSE16446_count_table %>% as.data.frame()

row.names(cts) <- cts$gene

cts <- subset (cts, select = -gene)

cts <- cts %>% select(SRR13380425,SRR13380426,SRR13380427,SRR13380430,SRR13380431,SRR13380443,SRR13380445,SRR13380452,SRR13380457,SRR13380466,SRR13380471,SRR13380477,SRR13380500,SRR13380503,SRR13380511,SRR13380512,SRR13380525,SRR13380531,SRR13380421,SRR13380434,SRR13380435,SRR13380438,SRR13380441,SRR13380446,SRR13380448,SRR13380458,SRR13380470,SRR13380472,SRR13380479,SRR13380480,SRR13380481,SRR13380482,SRR13380483,SRR13380484,SRR13380486,SRR13380488,SRR13380490,SRR13380495,SRR13380497,SRR13380499,SRR13380501,SRR13380502,SRR13380506,SRR13380508,SRR13380513,SRR13380516,SRR13380518,SRR13380520,SRR13380521,SRR13380528,SRR13380529,SRR13380533,SRR13380535,SRR13380537,SRR13380538,SRR13380539,SRR13380541)

# build coldata
col_data <- matrix(c('ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','ND','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D','T2D' ), ncol=1, byrow=TRUE)

col_data_day <- matrix(c("16/01/2014","16/01/2014","16/01/2014","16/01/2014","27/03/2014","27/03/2014","24/07/2015","27/03/2014","30/06/2014","30/06/2014","30/06/2014","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","16/01/2014","16/01/2014","26/02/2013","27/03/2014","27/03/2014","26/02/2013","30/06/2014","30/06/2014","24/07/2015","24/07/2015","27/03/2014","30/06/2014","27/03/2014","30/06/2014","27/03/2014","30/06/2014","27/03/2014","27/03/2014","24/07/2015","27/03/2014","27/03/2014","27/03/2014","27/03/2014","27/03/2014","27/03/2014","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/07/2015","24/11/2016","24/07/2015","24/07/2015","24/07/2015"), ncol=1, byrow=TRUE)
colnames(col_data) <- c('condition')
colnames(col_data_day) <- c('day')
col_data <- cbind(col_data,col_data_day)

rownames(col_data) <- c("SRR13380425","SRR13380426","SRR13380427","SRR13380430","SRR13380431","SRR13380443","SRR13380445","SRR13380452","SRR13380457","SRR13380466","SRR13380471","SRR13380477","SRR13380500","SRR13380503","SRR13380511","SRR13380512","SRR13380525","SRR13380531","SRR13380421","SRR13380434","SRR13380435","SRR13380438","SRR13380441","SRR13380446","SRR13380448","SRR13380458","SRR13380470","SRR13380472","SRR13380479","SRR13380480","SRR13380481","SRR13380482","SRR13380483","SRR13380484","SRR13380486","SRR13380488","SRR13380490","SRR13380495","SRR13380497","SRR13380499","SRR13380501","SRR13380502","SRR13380506","SRR13380508","SRR13380513","SRR13380516","SRR13380518","SRR13380520","SRR13380521","SRR13380528","SRR13380529","SRR13380533","SRR13380535","SRR13380537","SRR13380538","SRR13380539","SRR13380541")


#col_data <- as.table(col_data)

ddsHTSeq <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = col_data,
                                   design = ~ day + condition)

#Run DESeq2 on ddsHTSeq
dds = DESeq(ddsHTSeq)

#Normalize samples
rld = rlogTransformation(dds)

vsd <- vst(dds, blind = F)

mat <- assay(vsd)

# Day is the batch, the model.matrix should include everything but the batch
mat <- limma::removeBatchEffect(mat, vsd$day,
                                design = model.matrix(~vsd$condition))
assay(vsd) <- mat

#make a new PCA plot, note that now the colors identify samples from day 1 vs day 2
plotPCA(vsd, intgroup = c( "condition"),returnData = FALSE) +
  geom_text_repel(aes_string( x = "PC1", y = "PC2", label = "name"), size=3, color = "black") + theme_classic()+ labs(title = "PCA Plot (day batch correction)")

resultsNames(dds) 

res <- results(dds)
resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized =TRUE)), by="row.names", sort = FALSE)

write.csv(resdata, file = "GSE164416_DESeqs.csv")