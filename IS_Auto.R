# Evan Henrich
# ehenrich@fredhutch.org
# December 16, 2016

# PUPROSE: Automate the full pipeline of extraction, processing, and analysis of data for the 
# Immune Signatures Project. 

# NOTE: Short credits for each file are noted here, but please see individual
# files for more information.

#--------TESTING ONLY-------------------------------------
setwd("/home/ehenrich/R/ImmSig/")


#--------HELPERS-----------------------------------------
ftest <- function(file_list){
  if(FALSE %in% file_list){
    np <- which(file_list == FALSE)
    return(np)
  }
}

#------CHECK DIRECTORY STRUCTURE AND FILES ---------------
# need GET() for downloading from ImmuneSpace
library(httr)

# check directories are named and located appropriately
wk_dir <- getwd()
imm_dir <- basename(wk_dir)
output_dir <- file.path(wk_dir,"output")
data_dir <- file.path(wk_dir,"data")
preproc_dir <- file.path(wk_dir,"PreProc")
rds_dir <- file.path(wk_dir,"RDSGen")
hai_pp_data <- file.path(rds_dir,"HAI_PreProc_data")
ge_pp_data <- file.path(rds_dir,"GE_PreProc_data")
dir_list <- list(output_dir, data_dir, preproc_dir,
                  rds_dir, hai_pp_data, ge_pp_data)

dirs_present <- list()
for(dir in dir_list){
  dirs_present <- file.exists(dir)
}
dir_np <- ftest(dirs_present)

if(length(dir_np) == 6){
  dir_chk <- readline(prompt = "Directories absent. Ok to make? [y/n]")
  if(dir_chk == "y"){
    for(dir in dir_list){
      dir.create(dir, showWarnings = F)
    }
  }else{
    stop("User ended process.")
  }
}else if(((length(dir_np) > 0 & length(dir_np) < 6)) | (length(dir_np) > 6)){
  stop("Directories not setup correctly.  Please see readme.txt and fix by hand.")
}else if(length(dir_np) == 0){
  print("All Directories present")
}

# check for all necessary files
files <- list()

# Scripts in main directory, created by Hailong Meng @ Yale
files[["HIPCMetaModuleAnalysis_v2.R"]] <- wk_dir
files[["getSDY.R"]] <- wk_dir

# blood type data file
files[["BTM_for_GSEA_20131008.GMT"]] <- data_dir

# RDS generation file created by Renaud Gaujoux @ Tecnion & Stefan Avey @ Yale
files[["datasets_v2.R"]] <- rds_dir

# A tab-delimited table to map biosample ids to ImmPort subject / participant ids for SDY212
# SDY63, SDY404, SDY400 by Stefan Avey @ Yale, SDY212 and SDY80 maps by Evan Henrich @ Fred Hutch
files[["SDY212_IDmap.tsv"]] <- preproc_dir
files[["SDY63_IDmap.tsv"]] <- preproc_dir
files[["SDY404_IDmap.tsv"]] <- preproc_dir
files[["SDY400_IDmap.tsv"]] <- preproc_dir
files[["SDY80_IDmap.tsv"]] <- preproc_dir

# annotation tables collected / made by Evan Henrich @ Fred Hutch
files[["IlluminaHuman_v4_anno_table.csv"]] <- preproc_dir
files[["IlluminaHuman_450kMethylation_anno_table.txt"]] <- preproc_dir
files[["CHI_nih_gene_map.tsv"]] <- preproc_dir

# pre-processing scripts created by Evan Henrich @ Fred Hutch
files[["EH_makeHAI_PY.R"]] <- preproc_dir
files[["EH_makeGE_PY.R"]] <- preproc_dir
files[["EH_makeDemo.R"]] <- preproc_dir
files[["SDY67_exp_map.tsv"]] <- preproc_dir

files_present <- list()
for(fname in names(files)){
  files_present[[fname]] <- file.exists(file.path(files[[fname]],fname))
}

# downlaod missing files if necessary
files_np <- ftest(files_present)
if(length(files_np) > 0 & length(files_np) <= 16){
  file_chk <- readline(prompt = "Some or all dependency files absent. Ok to download? [y/n]")
  if(file_chk == "y"){
    for(fname in names(files)){
      link_base <- ""
      if(fname %in% files_np){
        GET(paste0(link_base,fname), write_disk(path = file.path(files[[fname]],fname),
                                                overwrite = F))
      }
    }
  }else{
    stop("User ended process.")
  }
}else{
  print("all dependency files found.")
}

#------GE Dependencies---------
print("Loading libraries and other dependencies")
library(ImmuneSpaceR)
# redundant: library(httr)
library(R.utils)
library(Rlabkey)
library(tools)
library(data.table)
library(hash)
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

# The first two tables are provided directly from Illumina for annotation for the corresponding beadchip.
# They are loaded directly because the packages available in bioconductor are either
# defunct (450k Methylation - SDY67) or did not find all probe ids (Humanv4 - Yale Studies).
# Although referenced in makeGE(), they are loaded here once to avoid repetitive loading and 
# slow perfomance.  The last table was generated from the original GEmatrix.txt for the CHI-nih study
# because an annotation table was not available publicly from Affymetrix for the Human Exon 1.0 ST GeneChip
# and no working annotation packages were found in bioconductor.
cat("loading gene annotation tables. May take 5 minutes.")
hmnv4_ann_tbl <- read.csv(paste0(rd_dir,"IlluminaHuman_v4_anno_table.csv"), 
                          stringsAsFactors = F)
hmn450k_ann_tbl <- read.table(paste0(rd_dir,"IlluminaHuman_450kMethylation_anno_table.txt"),
                              stringsAsFactors = F, sep = "\t", header = T)
hmnxn1_ann_tbl <- read.table(paste0(rd_dir,"CHI_nih_gene_map.tsv"),
                             stringsAsFactors = F, sep = "\t", header = T)

#------HAI Dependencies--------
# Redundant: library(ImmuneSpaceR)
# Redundant: library(data.table)
library(dplyr)
library(stringr)


#------Demo Dependencies-------
# Redundant: ImmuneSpaceR


#----GET USER INPUT--------------
go_ahead <- readline(prompt - "Running analysis make take 1 hour. Go ahead? [y/n]")
if(toupper(go_ahead) != "Y"){
  stop("User ended analysis.")
}

user <- readline(prompt = "Username for ImmuneSpace: ")
pwd <- readline(prompt = "Password for ImmuneSpace: ")

#------Main Methods------------
# Step 1: Extraction and Pre-Processing
studies <- c("SDY212","SDY63","SDY404","SDY400","SDY67","SDY80")

source(file.path(preproc_dir, ge_script))
source(file.path(preproc_dir, hai_script))
source(file.path(preproc_dir, demo_script))

for(sdy in studies){
  cat(paste0("Making HAI table for: ", sdy))
  makeHAI(sdy)
  cat(paste0("Making GE matrix for: ", sdy))
  makeGE(sdy, user, pwd)
  cat(paste0("Making Demo table for: ", sdy))
  makeDemo(sdy)
}

# Step 2: Combination of pre-processed files into Rds (BioConductor eset)
cat("Now combining files for each study into rds / eset file")
source(file.path(rds_dir, rds_gen))

# Step 3: Run meta analysis script
cat("Running meta analysis")
source("HIPCMetaModuleAnalysis_v2.R")



