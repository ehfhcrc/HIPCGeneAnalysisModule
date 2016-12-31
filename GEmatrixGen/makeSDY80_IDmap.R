gef <- get_gef("SDY80")

subids_out <- sapply(gef$name, FUN = function(x){
  targ_row <- gef[which(gef$name == x),]
  subid <- gsub("\\.80", "", targ_row$participant_id)
  subid <- paste0(subid, "_d", targ_row$study_time_collected)
  subid <- gsub("-","neg", subid)
})

subids <- sapply(gef$name, FUN = function(x){
  targ_row <- gef[which(gef$name == x),]
  subid <- gsub("\\.80", "", targ_row$participant_id)
})

subtypes <- unique(gef$subtype)
subtypes <- subtypes[ !is.na(subtypes)]

subtype_map <- sapply(subtypes, FUN = function(x){
  targ_rows <- gef[which(gef$subtype == x)]
  first_row <- targ_rows[1,]
  subid <- first_row$participant_id
  subid <- gsub("\\.80", "", subid)
})

names_ed <- sapply(subtypes, FUN = function(x){
  return(str_sub(x,1,3))
})

names(subtype_map) <- unname(names_ed)

rev_hash <- hash(unname(subtype_map), names(subtype_map))

final_biosample <- sapply(gef$participant_id, FUN = function(x){
  x <- gsub("\\.80", "", x)
  val <- rev_hash[[x]]
  return(val)
})

sdy80map <- data.frame(participantID = gef$participant_id,
                       GSM = gef$name,
                       GE_name = subids_out,
                       subtype = gef$subtype,
                       bioSampleID = final_biosample,
                       row.names = c())


  
  
