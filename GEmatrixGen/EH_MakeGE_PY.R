#ImmuneSignatures
#Evan Henrich
#Fred Hutch
#ehenrich@fredhutch.org
#dec 7, 2016

#--------GEN NOTES-------------------------------------------------------
# The goal of this program is to preprocess the raw Gene Expression data from either 
# ImmuneSpace.org or ImmPort into the format of GEMatrix.txt file used in the HIPC 
# ImmuneSignatures meta analysis.Much of the code is based on work done by Damian Fermin 
# [dfermin.yale@gmail.com] from Yale in 2014 that was used for the meta analysis paper, 
# particularly the statistical operations. There were some problems encountered due to data  # curation problems.  These include an extra row (gene)  present in SDY212 in ImmPort and 
# ImmuneSpace file (ILMN_2137536) as well as two subjects that are removed (SUB134242 and 
# SUB134267) in the final version.

#-------TESTING NOTES------(uncomment if running alone)--------------------
# Long load time because of loading the two large anno tables
setwd("/home/ehenrich/R/HIPCGeneModuleAnalysis/")
library(ImmuneSpaceR)
library(httr)
library(R.utils)
library(Rlabkey)
library(tools)
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

#-----DIR SETUP / GLOBAL VARS-----------------------------------------------
wk_dir <- getwd()
mgen_dir <- file.path(wk_dir,"GEmatrixGen")
rd_dir <- file.path(mgen_dir,"rawdata/")
dir.create(rd_dir, showWarnings = F)
# These tables are provided directly from Illumina for annotation for the corresponding beadchip.
# I have chosen to load them directly because the packages available in bioconductor are either
# defunct (450k Methylation - SDY67) or did not find all probe ids (Humanv4 - Yale Studies).
hmnv4_ann_tbl <- read.csv(paste0(rd_dir,"IlluminaHuman_v4_anno_table.csv"), stringsAsFactors = F)
hmn450k_ann_tbl <- read.table(paste0(rd_dir,"IlluminaHuman_450kMethylation_anno_table.txt"),
                              stringsAsFactors = F, sep = "\t", header = T)

#------HELPER METHODS ------------------------------------------------------
# Get gene expression file names from ImmuneSpace and return as data frame
get_gef <- function(sdy){
  con <- CreateConnection(sdy)
  gef <- con$getDataset("gene_expression_files")
  gef <- gef[ which(gef$file_info_name != "NA" & gef$study_time_collected == 0),  ]
  return(gef)
}

# download microarray or gene expression files from ImmuneSpace and return list of file paths
get_is_files <- function(gef, sdy, user, pwd){
  inputFiles <- unique(gef$file_info_name)
  
  # download files from IS to local machine
  links <- paste0("https://www.immunespace.org/_webdav/Studies/",
                  sdy,
                  "/%40files/rawdata/gene_expression/", 
                  inputFiles)
  
  # only setting to an object to limit output to console
  dump <- sapply(links, FUN = function(x){
    GET(url = x, 
        write_disk(paste0(rd_dir,basename(x)), 
        overwrite = T), 
        authenticate(user,pwd))
  })
  inputFiles <- file.path(rd_dir, inputFiles)
  return(inputFiles)
}

# download gene expression files from immport for yale studies
get_immport_files <- function(sdy){
  link <- ""
  if(sdy == "SDY63"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59635/suppl/GSE59635_non-normalized.txt.gz"
  }else if(sdy == "SDY404"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59654/suppl/GSE59654_non-normalized.txt.gz"
  }else if(sdy == "SDY400"){
    print("data not available yet")
  }
  inputFile <- paste0(rd_dir,sdy,".txt.gz")
  dump <- GET(url = link, write_disk(inputFile, overwrite = TRUE))
  gunzip(inputFile, overwrite = TRUE)
  inputFile <- paste0(rd_dir,sdy,".txt")
  return(inputFile)
}

# normalize expression values and map heads to correct subject ids
norm_map_mx <- function(sdy, exprs, is_col, smpl_col){
  # normalize expression values
  cnames <- colnames(exprs)
  exprs <- as.matrix(exprs)
  exprs <- preprocessCore::normalize.quantiles(exprs)
  colnames(exprs) <- cnames
  norm_exprs <- log2(pmax(exprs, 1))
  
  # map headers to correct IS subject ids
  exprs_headers <- colnames(norm_exprs)
  id_map <- read.table(paste0("GEmatrixGen/",sdy,"_IDmap.tsv"), 
                       sep = "\t", stringsAsFactors = F)
  
  mapped_headers <- sapply(exprs_headers, FUN = function(x){
    res <- ""  
      if(sdy == "SDY212" | sdy == "SDY67"){
        res <- id_map[[is_col]][which(id_map[[smpl_col]] == x)]
      }else if(sdy %in% c("SDY63","SDY400","SDY404")){
        spl_x <- strsplit(x, split = "_")
        if(sdy == "SDY63"){
          # sample ids were off according to original code, therefore correct
          spl_x[[1]][1] = as.integer(spl_x[[1]][1]) + 91000
        }
        sub_no_date <- id_map[[is_col]][which(id_map[[smpl_col]] == spl_x[[1]][1])]
        res <- paste0(sub_no_date,"_d",spl_x[[1]][3])
      }
    return(res)
    })
  
  # for the case where double of one participant "SUB134307"
  if(sdy == "SDY212"){
    dbl_subid <- which(mapped_headers == "SUB134307_d0")
    mapped_headers[dbl_subid[1]] <- "SUB134307.1_d0"
    mapped_headers[dbl_subid[2]] <- "SUB134307.2_d0"
  }
  
  # change colnames to sub_ids instead of biosample_ids
  colnames(norm_exprs) <- mapped_headers
  return(norm_exprs)
}

# write out final file containing prove ids, gene symbols, and expression values
write_out <- function(probe_ids, gene_syms, norm_map_exprs, sdy){
  out <- data.frame(
    probeID=probe_ids,
    geneSymbol=gene_syms,
    round(norm_map_exprs,6),
    stringsAsFactors=FALSE,
    row.names = NULL)
  
  write.table(out, file = paste0(wk_dir,"/RDSGen/GE_preproc_data/",sdy,".GEMatrix.EH.txt"), 
              sep = "\t", quote=FALSE, row.names=FALSE)
}

#---------MAIN METHOD--------------------------------------------------------
makeGE <- function(sdy, user, pwd){

  # instantiate variables to hold data for output at end
  probe_ids <- ""
  gene_syms <- ""
  final_expr_vals <- ""
  
  # Get raw data and process it uniquely for each study
  if(sdy == "SDY212" | sdy == "SDY67"){
    # Get filenames from Immunespace and then download
    gef <- get_gef(sdy)
    inputFiles <- get_is_files(gef, sdy, user, pwd)
    
    if(sdy == "SDY212"){
      # Clean and Prep
      rawdata <- fread(inputFiles, header = TRUE)
      probe_ids <- rawdata[, PROBE_ID]
      gene_syms <- rawdata[ , SYMBOL]
      sigcols <- grep("Signal", colnames(rawdata), value = TRUE)
      rawdata <- rawdata[, sigcols, with = FALSE]
      setnames(rawdata, gsub(".AVG.*$", "", colnames(rawdata)))
      
      # Norm / Map
      final_expr_vals <- norm_map_mx(sdy, rawdata, "final_id", "BioSampleID")
      write_out(probe_ids, gene_syms, final_expr_vals, sdy)
      
    
    }else if(sdy == "SDY67"){
      # inputFiles here only contain 1 subject worth of data per file, therefore
      # need to iterate through and build a master file with all subjects to return.
      # NOTE: no normalization was done on the original file and that is mimicked her.
      gs_tbl <- fread(inputFiles[[1]])
      gs_tbl$GENE_SYMBOL <- toupper(gs_tbl$GENE_SYMBOL)
      gs_tbl <- gs_tbl[ order(GENE_SYMBOL),]
      gene_syms <- gs_tbl$GENE_SYMBOL
      
      final_expr_vals <- as.data.frame(do.call(cbind,lapply(inputFiles, FUN = function(x){
        fname <- basename(x)
        targ_row <- gef[which(gef$file_info_name == fname),]
        subid <- substr(targ_row$participant_id[[1]],1,9)
        day_coll <- as.integer(targ_row$study_time_collected)
        subid <- paste0(subid,"_d",day_coll)
        
        #read in the table and relabel colnames
        tmp <- fread(x)
        colnames(tmp) <- c("SYMBOL",subid)
        tmp$SYMBOL <- toupper(tmp$SYMBOL)
        tmp <- tmp[ order(SYMBOL),]
        
        # make sure the gene symbol vec for the curr subject matches original
        if(!all.equal(tmp$SYMBOL,gs_tbl$GENE_SYMBOL)){
          stop(paste0("Gene Symbols in ",fname,
                      " do not match those in first file. Please check before re-running."))
        }
        
        # pull out expr values and return as vec
        return(tmp[,2])
      })))
      
      #these subjects were removed from original file, therefore removing here
      sub_rmv <- c("SUB113458_d0","SUB113463_d0","SUB113470_d0","SUB113473_d0","SUB113474_d0",
                    "SUB113476_d0","SUB113483_d0","SUB113487_d0","SUB113490_d0","SUB113494_d0",
                    "SUB113495_d0","SUB113496_d0","SUB113498_d0","SUB113504_d0","SUB113505_d0",
                    "SUB113513_d0","SUB113514_d0","SUB113524_d0","SUB113526_d0","SUB113527_d0",
                    "SUB113532_d0","SUB113535_d0","SUB113537_d0","SUB113545_d0","SUB113548_d0",
                    "SUB113555_d0","SUB113558_d0","SUB113559_d0","SUB113561_d0","SUB113566_d0",
                    "SUB113567_d0","SUB113568_d0","SUB113571_d0","SUB113572_d0","SUB113582_d0",
                    "SUB113583_d0","SUB113588_d0","SUB113595_d0","SUB113610_d0")
      
      final_expr_vals <- final_expr_vals[ , !(names(final_expr_vals) %in% sub_rmv)]
      
      # get ILMN-probes for gene symbols. NOTE: orig file does not actually have probe id vals.
      probe_ids <- sapply(gs_tbl$GENE_SYMBOL, FUN = function(x){
        res <- hmn450k_ann_tbl$ID[which(hmn450k_ann_tbl$Closest_TSS_gene_name == x)]
        if(length(res) > 0){
          return(res[[1]])
        }else{
          return("X")
        }
      })
    }
    
  }else if(sdy %in% c("SDY63","SDY404","SDY400")){
    # download files from ImmPort and read in data
    inputFiles <- get_immport_files(sdy)
    rawdata <- as.data.frame(fread(inputFiles))
    
    # Clean and Prep
    probe_ids <- rawdata$ID_REF
    gene_syms <- sapply(probe_ids, FUN = function(x){
      hmnv4_ann_tbl$ILMN_Gene[which(hmnv4_ann_tbl$ID == x)]
    })
    rawdata <- rawdata[ , grepl("PBMC", names(rawdata))]
    
    # Norm / Map
    final_expr_vals <- norm_map_mx(sdy, rawdata, "Sub.Org.Accession", "User.Defined.ID")
    
  }else{
    # for CHI-nih study
    inputFiles <- XXXX
  }
  write_out(probe_ids, gene_syms, final_expr_vals, sdy)
}
