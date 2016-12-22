#ImmuneSignatures
#Evan Henrich
#Fred Hutch
#ehenrich@fredhutch.org
#dec 7, 2016

#--------NOTES---------------------------------------
# The goal of this program is to preprocess the raw Gene Expression data from either ImmuneSpace.org
# or ImmPort into the format of GEMatrix.txt file used in the HIPC ImmuneSignatures meta analysis.
# Much of the code is based on work done by Damian Fermin [dfermin.yale@gmail.com] from Yale in 2014
# that was used for the meta analysis paper, particularly the statistical operations. 
# There were some problems encountered due to data curation problems.  These include an extra row (gene) 
# present in SDY212 in ImmPort and ImmuneSpace file (ILMN_2137536) as well as two subjects that 
# are removed (SUB134242 and SUB134267) in the final version.

#------helper methods ------------------------------
ge_id_gen <- function(subid, days){
    newstr <- paste0(subid,"_d",days)
    return(newstr)
}

# download microarray or gene expression files from ImmuneSpace and return list of file paths
get_is_files <- function(sdy,user,pwd){
  con <- CreateConnection(sdy)
  gef <- con$getDataset("gene_expression_files")
  gef <- gef[ which(gef$file_info_name != "NA" & gef$study_time_collected == 0),  ]
  inputFiles <- unique(gef$file_info_name)
  
  # download files from IS to local machine
  links <- paste0("https://www.immunespace.org/_webdav/Studies/",
                  sdy,
                  "/%40files/rawdata/gene_expression/", 
                  inputFiles)
  
  # only setting to an object to limit output to console
  dump <- sapply(links, FUN = function(x){
    GET(url = x, write_disk(paste0(targ_dir,basename(x))), authenticate(user,pwd))
  })
  
  inputFiles <- file.path(targ_dir, inputFiles)
  return(inputFiles)
}

# download gene expression files from immport for yale studies
get_immport_files <- function(sdy, user, pwd){
  link <- ""
  if(sdy == "SDY63"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59635/suppl/GSE59635_non-normalized.txt.gz"
  }else if(sdy == "SDY404"){
    link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59654/suppl/GSE59654_non-normalized.txt.gz"
  }
  
}
#---------Main Methods--------------------------------
library(ImmuneSpaceR)
library(httr)
library(Rlabkey)
library(tools)
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")

makeGE <- function(sdy, user, pwd){
  
  # setup target directory
  wk_dir <- getwd()
  targ_dir <- file.path(paste0(wk_dir,"/GEmatrixGen/rawdata/"))
  if(!exists(targ_dir)){dir.create(targ_dir)}
  
  # Get raw data and process it
  if(sdy == "SDY212" | sdy == "SDY67"){
    # download files from ImmuneSpace and return filepaths
    inputFiles <- get_is_files(sdy, user, pwd)
    # process data for each study
    exprs <- data.frame()
    #NTS: each subject had GE for day 0, 3, 28 ... are we just using d0? If so, then
    # need to sort gef to only download those files and build into matrix by mapping
    # to subject id from gef.
    if(sdy == "SDY67"){
      # iterate through inputFiles and bind to one df
      gs_tbl <- fread(inputFiles[[1]])
      exprs <- cbind(gs_tbl$GENE_SYMBOL)
      colnames(exprs) <- c("GENE_SYMBOL")
      lapply(inputFiles, FUN = function(x){
        fname <- basename(x)
        targ_row <- gef[which(gef$file_info_name == fname),]
        subid <- targ_row$participant_id
        tmp <- fread(x)
        colnames(tmp) <- c("GENE_SYMBOL",subid)
        if(tmp$GENE_SYMBOL == exprs[,"GENE_SYMBOL"]){
          exprs <- cbind(exprs,tmp[,2])
        }else{
          stop(paste0("Gene Symbols in ",fname," do not match those in first file. Please check before re-running."))
        }
      })
    }else if(sdy = "SDY212"){
      exprs <- fread(inputFiles, header = TRUE)
      probe_ids <- exprs[, PROBE_ID]
      gene_syms <- exprs[ , SYMBOL]
      
      sigcols <- grep("Signal", colnames(exprs), value = TRUE)
      exprs <- exprs[, sigcols, with = FALSE]
      setnames(exprs, gsub(".AVG.*$", "", colnames(exprs)))
      cnames <- colnames(exprs)
      exprs <- as.matrix(exprs)
    }
    # Back to both SDY212 and SDY67 sharing code
    exprs <- preprocessCore::normalize.quantiles(exprs)
    colnames(exprs) <- cnames
    norm_exprs <- log2(pmax(exprs, 1))
    
    bsmaped <- fread(paste0("GEmatrixGen/",sdy,"_IDmap.tsv"))
    
    # map the bs ids to sub ids
    exprs_headers <- colnames(norm_exprs)
    mapped_headers <- lapply(exprs_headers, FUN = function(x) {bsmaped$final_id[which(bsmaped$BioSampleID == x)]})
    mapped_headers <- unlist(lapply(mapped_headers, as.character))
    
    # for the case where double of one participant "SUB134307"
    if(sdy = "SDY212"){
      dbl_subid <- which(mapped_headers == "SUB134307_d0")
      mapped_headers[dbl_subid[1]] <- "SUB134307.1_d0"
      mapped_headers[dbl_subid[2]] <- "SUB134307.2_d0"
    }
    
    # change colnames to sub_ids instead of biosample_ids
    colnames(norm_exprs) <- mapped_headers
    
    # Add back the gene symbols column and Probe Ids
    out <- data.frame(
      probeID=probe_ids,
      geneSymbol=gene_syms,
      round(norm_exprs,6),
      stringsAsFactors=FALSE)
    
    # write matrix out to file
    write.table(out, file = paste0(sdy,".GEMatrix.EH.txt"), sep = "\t", quote=FALSE, row.names=FALSE)
    
  }else if(sdy %in% c("SDY63","SDY404","SDY400")){
    # download files from ImmPort and return file paths
    inputFiles <- get_immport_files(sdy, user, pwd)
  }else{
    # for CHI-nih study
    inputFiles <- XXXX
  }
}
