
# This has turned out to be nightmare fuel after reboot.
# R had problems finding Rcpp, so had to reinstall with dependencies
# Also had problems with limma from Bioconductor, but hopefully that's sorted

library(TCGAbiolinks)
library(DESeq2)

# Download AML RNAseq data
aml.rnaseq.query <- GDCquery(project = "TCGA-LAML", 
                     data.category = "Transcriptome Profiling", 
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "HTSeq - Counts")
GDCdownload(aml.rnaseq.query, method = "client")
aml.rnaseq <- GDCprepare(query = aml.rnaseq.query, 
                         save = TRUE, 
                         save.filename = "tcga-amlRnaseq.rda", 
                         summarizedExperiment = TRUE)

# Download clinical data
aml.clin.query <- GDCquery(project = "TCGA-LAML",
                      data.category = "Clinical")
GDCdownload(aml.clin.query, method="client")
aml.patient <- GDCprepare_clinic(aml.clin.query, clinical.info = "patient")

# Run DESeq2 on AML count data
aml.rnaseq.DE <- DESeqDataSet(aml.rnaseq, design = ~ patient)
aml.rnaseq.DE <- DESeq(aml.rnaseq.DE)