#ImmuneSignatures
#Evan Henrich
#Fred Hutch
#ehenrich@fredhutch.org
#dec 7, 2016

#--------NOTES---------------------------------------
# Goal is to reproduce the creation of SDY212 GEMatrix.txt file
# similary to work done by Damian Fermin Damian [dfermin.yale@gmail.com] from Yale in 2014
# that was used for the meta analysis paper.
# Problems encountered include an extra row (gene) present in 
# ImmPort and ImmuneSpace file (ILMN_2137536) as well as
# two subjects that are removed SUB134242 and SUB134267 in the final version.

#------helper method for bs_to_id conv ----------(run separately)
ge_id_gen <- function(subid, days){
    newstr <- paste0(subid,"_d",days)
    return(newstr)
}

#---------Main Methods--------------------------------


# If there is a copy of the ImmPort version of the Microarray data then use this ...
# and skip ahead to 'normalize matrix'
# inputFiles <- "GEmatrixGen/SDY212_WholeBlood_Microarray_update_11242014.389992_ImmPort.txt"

#setup inputs if using files from ImmuneSpace and running on rsT
# NOTES: SDY400 has no GE data?
sdy <- "SDY212"
con <- CreateConnection(sdy)
gef <- con$getDataset("gene_expression_files")
gef <- gef[ which(gef$file_info_name != "NA"),  ]
inputFiles <- unique(gef$file_info_name)
inputFiles <- file.path(paste0("/share/files/Studies/",sdy,"/@files/rawdata/gene_expression"), inputFiles)


# normalize matrix
exprs <- fread(inputFiles, header = TRUE)
sigcols <- grep("Signal", colnames(exprs), value = TRUE)
probe_ids <- exprs[, PROBE_ID]
gene_syms <- exprs[ , SYMBOL]
exprs <- exprs[, sigcols, with = FALSE]
setnames(exprs, gsub(".AVG.*$", "", colnames(exprs)))
cnames <- colnames(exprs)
exprs <- as.matrix(exprs)
#exprs <- normalizeBetweenArrays(exprs, method="quantile")
exprs <- preprocessCore::normalize.quantiles(exprs)
colnames(exprs) <- cnames
norm_exprs <- log2(pmax(exprs, 1))

# load the biosample to subject id mapping file
bsmaped <- read.table("GEmatrixGen/bs_to_id_final.tsv", sep = "\t")

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