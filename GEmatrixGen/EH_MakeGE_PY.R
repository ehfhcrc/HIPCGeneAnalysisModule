#ImmuneSignatures
#Evan Henrich
#Fred Hutch
#ehenrich@fredhutch.org
#dec 7, 2016

#--------NOTES---------------------------------------
# The goal of this program is to preprocess the raw Gene Expression data from either 
# ImmuneSpace.org or ImmPort into the format of GEMatrix.txt file used in the HIPC 
# ImmuneSignatures meta analysis.Much of the code is based on work done by Damian Fermin 
# [dfermin.yale@gmail.com] from Yale in 2014 that was used for the meta analysis paper, 
# particularly the statistical operations. There were some problems encountered due to data  # curation problems.  These include an extra row (gene)  present in SDY212 in ImmPort and 
# ImmuneSpace file (ILMN_2137536) as well as two subjects that are removed (SUB134242 and 
# SUB134267) in the final version.

#----Dependencies, uncomment if running script alone and not part of IS_auto.R
# library(ImmuneSpaceR)
# library(httr)
# library(R.utils)
# library(Rlabkey)
# library(tools)
# library(data.table)
# source("https://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")


#------helper methods ------------------------------
ge_id_gen <- function(subid, days){
    newstr <- paste0(subid,"_d",days)
    return(newstr)
}

# download microarray or gene expression files from ImmuneSpace and return list of file paths
get_is_files <- function(sdy, user, pwd, rd_dir){
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
  
  inputFiles <- file.path(rd_dir, inputFiles)
  return(inputFiles)
}

# download gene expression files from immport for yale studies
get_immport_files <- function(sdy, rd_dir){
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


#---------Main Methods--------------------------------
makeGE <- function(sdy, user, pwd){
  
  # setup target directory
  wk_dir <- getwd()
  mgen_dir <- file.path(wk_dir,"GEmatrixGen")
  rd_dir <- file.path(mgen_dir,"rawdata/")
  dir.create(rd_dir, showWarnings = F)
  norm_map_exprs <- "" # instantiate variable to hold normalized data for output at end
  
  # Get raw data and process it
  if(sdy == "SDY212" | sdy == "SDY67"){
    # download files from ImmuneSpace and return filepaths
    inputFiles <- get_is_files(sdy, user, pwd, rd_dir)
    #NTS: each subject had GE for day 0, 3, 28 for sdy67... are we just using d0? If so, then
    # need to sort gef to only download those files and build into matrix by mapping
    # to subject id from gef.
    if(sdy == "SDY67"){
      # Create a rawdata matrix from individual files
      gs_tbl <- fread(inputFiles[[1]])
      rawdata <- cbind(gs_tbl$GENE_SYMBOL)
      colnames(rawdata) <- c("GENE_SYMBOL")
      lapply(inputFiles, FUN = function(x){
        fname <- basename(x)
        targ_row <- gef[which(gef$file_info_name == fname),]
        subid <- targ_row$participant_id
        tmp <- fread(x)
        colnames(tmp) <- c("GENE_SYMBOL",subid)
        if(tmp$GENE_SYMBOL == exprs[,"GENE_SYMBOL"]){
          exprs <- cbind(exprs,tmp[,2])
        }else{
          stop(paste0("Gene Symbols in ",fname,
                      " do not match those in first file. Please check before re-running."))
        }
      })
    }else if(sdy == "SDY212"){
      # Clean and Prep
      rawdata <- fread(inputFiles, header = TRUE)
      probe_ids <- rawdata[, PROBE_ID]
      gene_syms <- rawdata[ , SYMBOL]
      sigcols <- grep("Signal", colnames(rawdata), value = TRUE)
      rawdata <- rawdata[, sigcols, with = FALSE]
      setnames(rawdata, gsub(".AVG.*$", "", colnames(rawdata)))
      
      # Norm / Map
      norm_map_exprs <- norm_map_mx(sdy, rawdata, "final_id", "BioSampleID")
    }
  }else if(sdy %in% c("SDY63","SDY404","SDY400")){
    # download files from ImmPort and read in data
    inputFiles <- get_immport_files(sdy, rd_dir)
    rawdata <- as.data.frame(fread(inputFiles))
    
    # Clean and Prep
    # NOTE: the full annotation table from Illumina for the expression chip used 
    # in the yale studies is used to map probes to gene symbols because the bioc 
    # library available for IlluminaHumanv4.db was not able to find all genes correctly
    probe_ids <- rawdata$ID_REF
    ilmn_ann_tbl <- read.csv(paste0(rd_dir,"ILMN_ann_tbl_ed.csv"), stringsAsFactors = F)
    gene_syms <- sapply(probe_ids, FUN = function(x){
      ilmn_ann_tbl$ILMN_Gene[which(ilmn_ann_tbl$ID == x)]
    })
    rawdata <- rawdata[ , grepl("PBMC", names(rawdata))]
    
    # Norm / Map
    norm_map_exprs <- norm_map_mx(sdy, rawdata, "Sub.Org.Accession", "User.Defined.ID")
    
  }else{
    # for CHI-nih study
    inputFiles <- XXXX
  }
  
  # Add back the gene symbols and Probe Ids columns
 out <- data.frame(
    probeID=probe_ids,
    geneSymbol=gene_syms,
    round(norm_map_exprs,6),
    stringsAsFactors=FALSE,
    row.names = NULL)
  
  # write matrix out to file
  write.table(out, file = paste0(wk_dir,"/RDSGen/GE_preproc_data/",sdy,".GEMatrix.EH.txt"), 
              sep = "\t", quote=FALSE, row.names=FALSE)
}
