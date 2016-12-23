# Evan Henrich
# ehenrich@fredhutch.org
# December 16, 2016

# Purpose: Auto the full pipeline of extraction, processing, and analysis of data for the Immune Signatures Project

#------GE Dependencies---------
library(ImmuneSpaceR)
library(httr)
library(R.utils)
library(Rlabkey)
library(tools)
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

#------HAI Dependencies--------



#------Demo Dependencies-------




#------Checks---------------
# check directories are named and located appropriately
full_path <- getwd()
imm_dir <- basename(full_path)
output_dir <- dir.exists(paste0(full_path,"output"))
data_dir <- dir.exists(paste0(full_path,"data"))

if(!output_dir | !data_dir | imm_dir == "ImmSig"){
  stop("File directories not names or located correctly. Please Change before retrying.")
}
  
# check for all necessary files
# file containing all BTM (Blood Transcript Module) gene modules
btm_file <- "BTM_for_GSEA_20131008.GMT"
# A tab-delimited table to map biosample ids to ImmPort subject / participant ids for SDY212
stan_ids <- "SDY212_IDmap.tsv"
yale_one_ids <- "SDY63_IDmap.tsv"
yale_two_ids <- "SDY404_IDmap.tsv"
yale_three_ids <- "SDY400_IDmap.tsv"
id_files <- c(stan_ids, yale_one_ids, yale_two_ids, yale_three_ids)

if(!file.exists(path(full_path,"data/",btm_file))){
  stop("File missing from main directory: BTM_for_GSEA_20131008.GMT")
}

files_present <- list()
for(fl in id_files){
  files_presen[[fl]] <- file.exists(path(full_path,"PreProc_Scripts/",fl))
}

if(FALSE %in% files_present){
  np <- which(files_present == FALSE)
  stop(paste0("ID Mapping Files Missing: ", names(np)))
}

# 
# ***SCRIPTS***
#   --HIPCMetaModuleAnalysis_v2.R
# --getSDY.R
# --makeGEmatrix.R
# --makeHAItable.R
# --makeDemo.R
# --makeRds.R
# check R version

# PROMPT USER FOR USERNAME AND PASSWORD FOR IS

#------Main Methods------------

# Step 1: Extraction and Pre-Processing
studies <- c("SDY212","SDY63","SDY404","SDY400","SDY67","CHI-nih")

source("/PreProc_Scripts/makeGEmatrix.R")
source("/PreProc_Scripts/makeHAItable.R")
source("/PreProc_Scripts/makeDemo.R")

for(sdy in studies){
  print(paste0("Making HAI table for: ",sdy))
  makeHAI(sdy)
  print(paste0("Making GE matrix for: ",sdy))
  makeGE(sdy,user,pwd)
  print(paste0("Making Demo table for: ",sdy))
  makeDemo(sdy)
}

# Step 2: Combination of pre-processed files into Rds (BioConductor eset)
print("Now combining files for each study into rds / eset file")
source("/RDSGen/makeRds.R")

# Step 3: Run meta analysis script
print("Running meta analysis")
source("HIPCMetaModuleAnalysis_v2.R")