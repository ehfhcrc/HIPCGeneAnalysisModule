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

#-------TESTING------(uncomment if running alone)--------------------
# Setup directory vars.

# setwd("/home/ehenrich/R/HIPCGeneModuleAnalysis/")
# wk_dir <- getwd()
# mgen_dir <- file.path(wk_dir,"GEmatrixGen")
# rd_dir <- file.path(mgen_dir,"rawdata/")
# library(ImmuneSpaceR)
# library(httr)
# library(XML)
# library(xml2)
# library(R.utils)
# library(Rlabkey)
# library(tools)
# library(data.table)
# source("https://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")

# Long load time because of loading the two large anno tables
# Top two tables are provided directly from Illumina for annotation for the corresponding beadchip.
# They are loaded directly because the packages available in bioconductor are either
# defunct (450k Methylation - SDY67) or did not find all probe ids (Humanv4 - Yale Studies). Final table
# was generated from the original tables developed by Yuri Kotliarov at NIH because Affymetrix
# does not provide annotation tables to the public.

# hmnv4_ann_tbl <- read.csv(paste0(rd_dir,"IlluminaHuman_v4_anno_table.csv"),
#                           stringsAsFactors = F)
# hmn450k_ann_tbl <- read.table(paste0(rd_dir,"IlluminaHuman_450kMethylation_anno_table.txt"),
#                               stringsAsFactors = F, sep = "\t", header = T)
# hmnxn1_ann_tbl <- read.table(paste0(rd_dir,"CHI_nih_gene_map.tsv"),
#                              stringsAsFactors = F, sep = "\t", header = T)


#------HELPER METHODS ------------------------------------------------------
# Get gene expression file names from ImmuneSpace and return as data frame
get_gef <- function(sdy){
  con <- CreateConnection(sdy)
  gef <- con$getDataset("gene_expression_files")
  if(sdy != "SDY80"){
    gef <- gef[ which(gef$file_info_name != "NA" & gef$study_time_collected == 0),  ]
  }
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
        authenticate(user,pwd)
        )
    })
  inputFiles <- file.path(rd_dir, inputFiles)
  return(inputFiles)
}

# for SDY80, generate a list of tables by scraping GEO websites
get_geo_files <- function(gef){
  
  # list of accessions to pull
  geo_acc <- gef$geo_accession
  
  # create links with following format and pull data:
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc=GSM1147759&id=13600&db=GeoDb_blob98
  base <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc="
  joint <- "&id="
  tail <- "&db=GeoDb_blob98"
  
  tbl_list <- lapply(geo_acc, FUN = function(x){
    # generate id and link
    digits <- gsub("GSM", "", x)
    id <- as.integer(digits) - 1134159
    link <- paste0(base, x, joint, id, tail)
    
    # scrape link for content and parse to table, return with probe ID
    dump <- GET(link)
    dump_contents <- read_html(dump$content)
    dump_contents_list <- xmlToList(xmlParse(dump_contents))
    raw_table <- dump_contents_list$body$font$pre[[4]]
    tbl_con<- textConnection(raw_table)
    data <- read.table(tbl_con, sep = "\t", stringsAsFactors = F, header = F)
    colnames(data) <- c("probeID", x)
    return(data)
  })
  
  # return list of tables
  return(tbl_list)
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

# Significantly faster than using data.table or which statements for searching paired-values
# NOTE: the key must be presented as a character - it does not auto-caste
hashmap <- function(keys, values, targ_vec){
  tmp_hash <- hash(keys, values)
  mapped_res <- sapply(targ_vec, FUN = function(x){
    val <- tmp_hash[[as.character(x)]]
    val <- ifelse(is.null(val),NA,val)
    names(val) <- x
    return(val)
  })
  return(mapped_res)
}

# Generate subject ID values that can be used downstream and mapped to expr values
map_headers <- function(sdy, exprs_headers, is_col, smpl_col){
  id_map <- read.table(paste0("GEmatrixGen/",sdy,"_IDmap.tsv"), 
                       sep = "\t", stringsAsFactors = F)
  
  if(sdy %in% c("SDY212","SDY67", "SDY80")){
    mapped_headers <- hashmap(id_map[[smpl_col]], id_map[[is_col]], exprs_headers)
  }else if(sdy %in% c("SDY63","SDY400","SDY404")){
    mapped_headers <- sapply(exprs_headers, FUN = function(x){
      spl_x <- strsplit(x, split = "_")
      if(sdy == "SDY63"){
        # sample ids were off according to original code, therefore correct
        spl_x[[1]][1] = as.integer(spl_x[[1]][1]) + 91000
      }
      sub_no_date <- id_map[[is_col]][which(id_map[[smpl_col]] == spl_x[[1]][1])]
      res <- paste0(sub_no_date,"_d",spl_x[[1]][3])
      return(res)
    })
  }
  
  # for the case where two instances / columns of one participant "SUB134307"
  if(sdy == "SDY212"){
    dbl_subid <- which(mapped_headers == "SUB134307_d0")
    mapped_headers[dbl_subid[1]] <- "SUB134307.1_d0"
    mapped_headers[dbl_subid[2]] <- "SUB134307.2_d0"
  }
  return(mapped_headers)
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
  mapped_headers <- map_headers(sdy, exprs_headers, is_col, smpl_col)
  
  # change colnames to sub_ids instead of biosample_ids
  colnames(norm_exprs) <- mapped_headers
  return(norm_exprs)
}


# Remove subjects not found in original files.  Unclear why subjects were removed.
remove_subs <- function(expr_matrix, subs_to_rm){
  expr_matrix <- expr_matrix[ , !(colnames(expr_matrix) %in% subs_to_rm)]
  return(expr_matrix)
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
              sep = "\t", 
              quote=FALSE, 
              row.names=FALSE)
}



#---------MAIN METHOD--------------------------------------------------------
makeGE <- function(sdy, user, pwd){

  # Setup directory vars
  wk_dir <- getwd()
  mgen_dir <- file.path(wk_dir,"GEmatrixGen")
  rd_dir <- file.path(mgen_dir,"rawdata/")
  
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
      # remove gene with NA values. These NAs are found in ImmPort as well.
      # ImmPort staff (Patrick.Dunn@nih.gov) said this was an unresolved issue
      # with ticket "HDB-13" last discussed in October 2015 with study authors.
      rawdata <- rawdata[ -c(which(rawdata$PROBE_ID == "ILMN_2137536")),]
      probe_ids <- rawdata[, PROBE_ID]
      gene_syms <- rawdata[ , SYMBOL]
      sigcols <- grep("Signal", colnames(rawdata), value = TRUE)
      rawdata <- rawdata[, sigcols, with = FALSE]
      setnames(rawdata, gsub(".AVG.*$", "", colnames(rawdata)))
      
      # Norm / Map
      final_expr_vals <- norm_map_mx(sdy, rawdata, "final_id", "BioSampleID")
      
      # remove subjects not found in original file
      subs_rm <- c("SUB134242_d0","SUB134267_d0")
      final_expr_vals <- remove_subs(final_expr_vals, subs_rm)
      
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
      subs_rm <- c("SUB113458_d0","SUB113463_d0","SUB113470_d0","SUB113473_d0","SUB113474_d0",
                    "SUB113476_d0","SUB113483_d0","SUB113487_d0","SUB113490_d0","SUB113494_d0",
                    "SUB113495_d0","SUB113496_d0","SUB113498_d0","SUB113504_d0","SUB113505_d0",
                    "SUB113513_d0","SUB113514_d0","SUB113524_d0","SUB113526_d0","SUB113527_d0",
                    "SUB113532_d0","SUB113535_d0","SUB113537_d0","SUB113545_d0","SUB113548_d0",
                    "SUB113555_d0","SUB113558_d0","SUB113559_d0","SUB113561_d0","SUB113566_d0",
                    "SUB113567_d0","SUB113568_d0","SUB113571_d0","SUB113572_d0","SUB113582_d0",
                    "SUB113583_d0","SUB113588_d0","SUB113595_d0","SUB113610_d0")
      final_expr_vals <- remove_subs(final_expr_vals, subs_rm)
      
      # get ILMN-probes for gene symbols. NOTE: orig file does not actually have probe id vals.
      probe_ids <- hashmap(hmn450k_ann_tbl$Closest_TSS_gene_name, 
                           hmn450k_ann_tbl$ID, 
                           gs_tbl$GENE_SYMBOL)
      probe_ids <- sapply(probe_ids, FUN = function(x){
                            if(is.null(x)){
                              return(NA)
                            }else{
                              return(x)
                            }
                   })
      probe_ids <- unname(unlist(probe_ids))
    }
    
  }else if(sdy %in% c("SDY63","SDY404","SDY400")){
    # download files from ImmPort and read in data
    inputFiles <- get_immport_files(sdy)
    rawdata <- as.data.frame(fread(inputFiles))
    
    # Clean and Prep
    probe_ids <- rawdata$ID_REF
    gene_syms <- hashmap(hmnv4_ann_tbl$ID, hmnv4_ann_tbl$ILMN_Gene, probe_ids)
    rawdata <- rawdata[ , grepl("PBMC", names(rawdata))]
    
    # Norm / Map
    final_expr_vals <- norm_map_mx(sdy, rawdata, "Sub.Org.Accession", "User.Defined.ID")
    
  }else if(sdy == "SDY80"){
    # get data via scraping geo site and parsing into list of tables
    gef <- get_gef(sdy)
    inputFiles <- get_geo_files(gef)
    

    probe_ids <- unname(inputFiles[[1]][,1])
    cnames <- sapply(inputFiles, FUN = function(x){colnames(x)[[2]]})
    expr_vals <- as.data.frame(do.call(cbind,
                                       lapply(inputFiles, FUN = function(x){
                                         return(x[,2])
                                        })
                                ))
    colnames(expr_vals) <- cnames
    rawdata <- cbind(probe_ids, expr_vals)
    
    # map gene symbols to probe ids using a file created by IDs and gs from original
    # CHI-nih.GEmatrix.txt.  To mimic original file, I remove all rows that were not 
    # successfully mapped, which is approximately 45%.  Unsure why this was done.
    gene_syms <- hashmap(hmnxn1_ann_tbl$probeID, hmnxn1_ann_tbl$gs, probe_ids)
    gene_syms <- gene_syms[which(sapply(gene_syms, FUN = function(x){!is.na(x)}))]
    probe_ids <- probe_ids[which(probe_ids %in% names(gene_syms))]
    rawdata <- rawdata[ which(rawdata$probe_ids %in% names(gene_syms)), ]
    rawdata <- rawdata[,-1]
    
    # Setting orig to TRUE ensures that data is not normalized, but maintained in the log2 transform
    # it comes from GEO as.  This was the way the original file was done.  However, by setting orig to
    # false, you can undo the transformation and then send down the same pipeline as the other studies.
    orig <- TRUE
    
    if(orig == T){
      # to mimic original file data, but use ImmuneSpace subjectIDs instead of biosampleIDs
      mapped_headers <- unname(map_headers(sdy, colnames(expr_vals), is_col = "GE_name", smpl_col = "GSM"))
      colnames(rawdata) <- mapped_headers
      final_expr_vals <- rawdata
    }else{
      # log2 untransform rawdata
      rawdata <- 2 ^ rawdata
      final_expr_vals <- norm_map_mx(sdy, rawdata, is_col = "GE_name", smpl_col = "GSM")
    }
    sdy <- "CHI-nih"
  }
  
  write_out(probe_ids, gene_syms, final_expr_vals, sdy)
  
}

#---------------------DEPRECATED---------------------------------------------------------------------
# # for CHI-nih study, the rawdata is taken from the NCBI GEO accession GSE47353 series_matrix file,
# # which as been renamed 'CHI-nih_GE_matrix.txt' for clarity.
# rawdata <- read.table(file.path(rd_dir,"CHI_nih_rawdata.txt"),
#                       sep = "\t", 
#                       stringsAsFactors = F, 
#                       skip = 97, 
#                       fill= T, 
#                       header = F)
# header_row <- read.table(file.path(rd_dir,"CHI_nih_rawdata.txt"),
#                          sep = "\t", 
#                          stringsAsFactors = F, 
#                          skip = 42, 
#                          nrow = 1)
# 
# # Change names to match original file
# header_row[1] <- "ID_REF"
# header_row <- sapply(header_row, FUN = function(x){
#   x <- gsub("PBMC", "sub", x)
#   x <- gsub("day", "d", x)
#   x <- gsub("-","neg",x)
# })
# colnames(rawdata) <- header_row
# 
# # map gene symbols to probe ids // using a file created by IDs and gs from original
# # CHI-nih.GEmatrix.txt.  Remove all rows that were not successfully mapped, which is 
# # approximately 45%.  Unsure why this was done.
# gene_syms <- hashmap(hmnxn1_ann_tbl$probeID, hmnxn1_ann_tbl$gs, rawdata$ID_REF)
# gene_syms <- gene_syms[which(sapply(gene_syms, FUN = function(x){!is.null(x)}))]
# probe_ids <- rawdata$ID_REF[which(rawdata$ID_REF %in% names(gene_syms))]
# rawdata <- rawdata[ which(rawdata$ID_REF %in% names(gene_syms)), ]
# 
# # To match the original file, no normalization is done here.  It would appear that 
# # normalization has already been performed on the data coming from GEO as 'series-matrix.txt'
# final_expr_vals <- rawdata[ , -1]
# gene_syms <- unname(unlist(gene_syms))