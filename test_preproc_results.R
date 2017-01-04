# Evan Henrich
# ehenrich@fredhutch.org
# December 29, 2016


#-----------PURPOSE--------------------------------
# Testing the results for generating HAI, GE, and Demo txt documents from rawdata
# for the immuneSignatures project ... done from HIPCGeneModuleAnalysis folder

#------GENERAL DEPENDENCIES--------------------------
setwd("/home/ehenrich/R/HIPCGeneModuleAnalysis/")
source("HAIgen/EH_MakeHAI_PY.R")
source("GEmatrixGen/EH_MakeGE_PY.R")
source("DemoGen/EH_MakeDemo.R")

#------GE Dependencies---------
library(ImmuneSpaceR)
library(httr)
library(XML)
library(xml2)
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
# and no annotation packages were found in bioconductor.

wk_dir <- getwd()
mgen_dir <- file.path(wk_dir,"GEmatrixGen")
rd_dir <- file.path(mgen_dir,"rawdata/")

cat("Loading annotation tables .... please be patient")
hmnv4_ann_tbl <- read.csv(paste0(rd_dir,"IlluminaHuman_v4_anno_table.csv"), 
                          stringsAsFactors = F)
hmn450k_ann_tbl <- read.table(paste0(rd_dir,"IlluminaHuman_450kMethylation_anno_table.txt"),
                              stringsAsFactors = F, sep = "\t", header = T)
hmnxn1_ann_tbl <- read.table(paste0(rd_dir,"CHI_nih_gene_map.tsv"),
                             stringsAsFactors = F, sep = "\t", header = T)

user <- "readonly@rglab.org"
pwd <- "JFB4C9hJFB4C9h"

#------HAI Dependencies (Not in GE overlap)--------
# library(ImmuneSpaceR)
# library(data.table)
library(dplyr)
library(stringr)


#------Demo Dependencies-------
# library(ImmuneSpaceR)
# library(hash)

#-----------HELPERS--------------------------------
rd_tbl <- function(fname, order_col){
  res <- read.table(fname,
                    sep = "\t",
                    stringsAsFactors = F,
                    header = T)
  res <- res[order(res[order_col]),]
  return(res)
}

test_ge <- function(sdy){
  print(paste0("Checking GE: ", sdy))
  start_time <- proc.time()
  makeGE(sdy, user, pwd)
  fin_time <- proc.time() - start_time
 
  if(sdy == "SDY67"){
    ofile <- paste0("RDSGen/GE_preproc_data/", sdy ,"-batch2.GEMatrix.txt")
    old_res <- rd_tbl(fname = ofile, order_col = "geneSymbol")
  }else if(sdy == "SDY80"){
    ofile <- "RDSGen/GE_preproc_data/CHI-nih.GEMatrix.txt"
    old_res <- rd_tbl(fname = ofile, order_col = "gs")
  }else{
    ofile <- paste0("RDSGen/GE_preproc_data/", sdy ,".GEMatrix.txt")
    old_res <- rd_tbl(fname = ofile, order_col = "geneSymbol")
  }
  
  if(sdy == "SDY80"){
    nfile <- "RDSGen/GE_preproc_data/CHI-nih.GEMatrix.EH.txt"
  }else{
    nfile <- paste0("RDSGen/GE_preproc_data/", sdy ,".GEMatrix.EH.txt")
  }

  new_res <- rd_tbl(fname = nfile, order_col = "geneSymbol")
  
  res_arr <- list()
  res_arr[["process time"]] <- fin_time
  res_arr[["dims"]] <- c(dim(old_res),dim(new_res))
  res_arr[["dim match"]] <- identical(dim(old_res),dim(new_res))
  
  if(sdy == "SDY80"){
    res_arr[["vals match"]] <- all.equal(old_res$sub_200_d0[2:10],new_res$SUB114450_d0[2:10])
    res_arr[["vals sample"]] <- data.frame(old_res$sub_200_d0[2:10],new_res$SUB114450_d0[2:10])
  }else{
    res_arr[["vals match"]] <- all.equal(old_res[2:10, 3],new_res[2:10, 3])
    res_arr[["vals sample"]] <- data.frame(old_res[2:10, 3],new_res[2:10, 3])
  }
 
  print(res_arr)
  cat("*****************************************")
  cat('\n')
}

test_hai <- function(sdy){
  # HAI
  print(paste0("Checking HAI: ", sdy))
  start_time <- proc.time()
  makeHAI(sdy)
  fin_time <- proc.time() - start_time
  
  if(sdy == "SDY80"){
    ofile <- paste0("RDSGen/HAI_preproc_data/CHI-nih_combined_hai_titer_table.txt")
    nfile <- paste0("RDSGen/HAI_preproc_data/CHI-nih_combined_hai_titer_table_EH.txt")
  }else{
    ofile <- paste0("RDSGen/HAI_preproc_data/", sdy , "_combined_hai_titer_table.txt")
    nfile <- paste0("RDSGen/HAI_preproc_data/", sdy , "_combined_hai_titer_table_EH.txt")
  }
  
  
  old_res <- rd_tbl(fname = ofile, order_col = "subject")
  new_res <- rd_tbl(fname = nfile, order_col =  "subject")
  res_arr <- list()
  res_arr[["process time"]] <- fin_time
  res_arr[["dims"]] <- c(dim(old_res),dim(new_res))
  res_arr[["dim match"]] <- identical(dim(old_res),dim(new_res))
  res_arr[["fc_res_max_d20 match"]] <- all.equal(old_res$fc_res_max_d20, new_res$fc_res_max_d20)
  res_arr[["fc_res_max match"]] <- all.equal(old_res$fc_res_max, new_res$fc_res_max)
  res_arr[["fc_res_max samples"]] <- data.frame(old_res$fc_res_max[1:10], new_res$fc_res_max[1:10])
  print(res_arr)
  cat("*****************************************")
  cat('\n')
}

test_demo <- function(sdy){
  print(paste0("Checking Demo: ", sdy))
  start_time <- proc.time()
  makeDemo(sdy)
  fin_time <- proc.time() - start_time
  
  if(sdy == "SDY80"){
    ofile <- paste0("RDSGen/GE_preproc_data/CHI-nih.demographics.txt")
    nfile <- paste0("RDSGen/GE_preproc_data/CHI-nih.demographics.EH.txt")
  }else if(sdy == "SDY67"){
    ofile <- paste0("RDSGen/GE_preproc_data/", sdy , "-batch2.demographics.txt")
    nfile <- paste0("RDSGen/GE_preproc_data/", sdy , "-batch2.demographics.EH.txt")
  }else{
    ofile <- paste0("RDSGen/GE_preproc_data/", sdy , ".demographics.txt")
    nfile <- paste0("RDSGen/GE_preproc_data/", sdy , ".demographics.EH.txt")
  }
  
  old_res <- rd_tbl(fname = ofile, order_col = 1)
  new_res <- rd_tbl(fname = nfile, order_col = 1)
  res_arr <- list()
  res_arr[["process time"]] <- fin_time
  res_arr[["dims"]] <- c(dim(old_res),dim(new_res))
  res_arr[["dim match"]] <- identical(dim(old_res),dim(new_res))
  res_arr[["data equality"]] <- all.equal(old_res, new_res)
  res_arr[["old sample"]] <-  old_res[1:5, 1:ncol(old_res)]
  res_arr[["new sample"]] <-  new_res[1:5, 1:ncol(new_res)]
  print(res_arr)
  cat("*****************************************")
  cat('\n')
}


test_all <- function(studies){
  for(sdy in studies){
    if(sdy != "SDY400"){
      test_ge(sdy)
    }
    test_hai(sdy)
    test_demo(sdy)
  }
}

#-----------MAIN-----------------------------------

studies <- c("SDY212", "SDY63", "SDY404", "SDY400", "SDY80", "SDY67")
test_all(studies)




