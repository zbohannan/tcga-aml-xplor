
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

aml.patient <- GDCprepare_clinic(aml.clin.query, clinical.info = "patient")

# Convert raw count data into matrix, format patient numbers and convert to row names,
# remove extra clinical entries, and then try doing DESeq2 with that

# CHECK EXTRA PATIENT DATA LATER

aml.rnaseq.matrix <- assay(aml.rnaseq)
colnames(aml.rnaseq.matrix) <- sub("TCGA-AB-","",colnames(aml.rnaseq.matrix))
colnames(aml.rnaseq.matrix) <- sub("-.*","",colnames(aml.rnaseq.matrix))
rownames(aml.patient) <- aml.patient[,"patient_id"]
aml.patient.cut <- aml.patient[colnames(aml.rnaseq.matrix),]

# Double check that data match
all(rownames(aml.patient.cut) %in% colnames(aml.rnaseq.matrix))
all(rownames(aml.patient.cut) == colnames(aml.rnaseq.matrix))

# Create DESeq data set and run DEseq2
aml.rnaseq.DE <- DESeqDataSetFromMatrix(countData = aml.rnaseq.matrix,
                                        colData = aml.patient.cut,
                                        design = ~ acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
aml.rnaseq.DE <- DESeq(aml.rnaseq.DE)

# Run DESeq2 on AML count data
#aml.rnaseq.DE <- DESeqDataSet(aml.rnaseq, design = ~ patient)
#aml.rnaseq.DE <- DESeq(aml.rnaseq.DE)