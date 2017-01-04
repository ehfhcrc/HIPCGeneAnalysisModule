# Evan Henrich
# December 8, 2016
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

# PURPOSE: generate DEMO tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

# NOTES: The original code to perform these operations was developed by collaborators of the HIPC Immune
# Signatures project.  This code was developed to provide the same functionality as the original in terms of 
# statistical operation, but use data pulled directly from the ImmuneSpace portal at www.immunespace.org
# instead of data shared among the collaborating labs via a google drive and also to handle all the studies
# used in the meta-analysis in an automated format.

#------TESTING-----------------
 setwd("/home/ehenrich/R/HIPCGeneModuleAnalysis/")


#------Dependencies-------------
# library(ImmuneSpaceR)

#------Helper Functions---------


#------MAIN METHOD---------------
makeDemo <- function(sdy){
  
  # Setup directory vars
  wk_dir <- getwd()
  
  # Get rawdata from ImmuneSpace
  con <- CreateConnection(sdy)
  data <- data.frame(con$getDataset("demographics"))
  
  cols_to_keep <- c("participant_id","age_reported","gender")
  if(sdy %in% c("SDY212","SDY63","SDY404","SDY400")){
    cols_to_keep <- c(cols_to_keep,"phenotype","race","ethnicity")
  }else if(sdy == "SDY67"){
    cols_to_keep <- c(cols_to_keep,"race")
    sdy <- "SDY67-batch2."
  }else if(sdy == "SDY80"){
    cols_to_keep <- c(cols_to_keep,"ethnicity")
    sdy <- "CHI-nih"
  }
  
  data <- data[ , names(data) %in% cols_to_keep]
  
  if(sdy == "SDY67"){
    # data, add exp_samples
  }
  
  write.table(data, file = paste0(wk_dir,"/RDSGen/GE_preproc_data/",sdy,".demographics.EH.txt"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}
  