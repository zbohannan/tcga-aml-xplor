library(TCGAbiolinks)
AMLquery <- GDCquery(project = "TARGET-AML", 
                     data.category = "Transcriptome Profiling", 
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "HTSeq - Counts")
GDCdownload(AMLquery)
aml.rnaseq <- GDCprepare(query = AMLquery, 
                         save = True, 
                         save.filename = "amlRnaseq.rda", 
                         summarizedExperiment = TRUE)
# need to figure out how to open clinical data.
# currently getting errors using both GDCquery_clinic
# and GDCprepare_clinic

# Error in doc_parse_file(con, encoding = encoding, as_html = as_html, options = options) : 
# basic_string::resize