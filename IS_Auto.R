# Evan Henrich
# ehenrich@fredhutch.org
# December 16, 2016

# Purpose: Auto the full pipeline of extraction, processing, and analysis of data for the Immune Signatures Project

#------Dependencies---------



#------Checks---------------
# check working directory
# check subdirectories
# check for all necessary files

# ***FILES***
#   --BTM_for_GSEA_20131008.GMT
# >>file containing all BTM (Blood Transcript Module) gene modules
# --bs_to_id_final.tsv
# >>A tab-delimited table to map biosample ids to ImmPort subject / participant ids for SDY212
# 
# ***SCRIPTS***
#   --HIPCMetaModuleAnalysis_v2.R
# --getSDY.R
# --makeGEmatrix.R
# --makeHAItable.R
# --makeDemo.R
# --makeRds.R
# check R version

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
  makeGE(sdy)
  print(paste0("Making Demo table for: ",sdy))
  makeDemo(sdy)
}

# Step 2: Combination of pre-processed files into Rds (BioConductor eset)
print("Now combining files for each study into rds / eset file")
source("/RDSGen/makeRds.R")

# Step 3: Run meta analysis script
print("Running meta analysis")
source("HIPCMetaModuleAnalysis_v2.R")