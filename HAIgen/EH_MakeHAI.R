# Evan Henrich
# December 8, 2016
# FHCRC
# ehenrich@fredhutch.org

# PURPOSE: generate HAI tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

#------Dependencies-------------
library(data.table)
library(stringr)
library(ImmuneSpaceR)

#------Helper Functions---------
med_sd_calc <- function(prefix,strains,glob_vals,titer_data){
  suf <- c("med","sd")
  for(virus in strains){
    for(s in suf){
      tar_list <- paste0(prefix,s)
      tar_col <- paste0(prefix,virus)
      if(s == "med"){
        glob_vals[[tar_list]][[virus]] <- median(titer_data[[tar_col]], na.rm = TRUE)
      }else{
        glob_vals[[tar_list]][[virus]] <- sd(titer_data[[tar_col]], na.rm = TRUE)
      }
    }
  }
  return(glob_vals)
}

max_select <- function(subid,trg_col){
  tmp_ls <- list()
  for(virus in gl_strains){
    tmp_ls[virus] <- gl_tdata[gl_tdata$subject == subid, paste0(trg_col,virus)]
  }
  return(max(unlist(tmp_ls),na.rm = TRUE))
}

# calculate the threshold values for discretization of input_col based on percentiles given
thresh_calc <- function(subids, low_perc, input_col, titer_data){
  num_ids <- length(subids)
  low_thresh <- round(num_ids * (low_perc / 100))
  high_thresh <- round(num_ids * ((100-low_perc) / 100))
  ord_tdata <- titer_data[order(titer_data[[input_col]]),]
  low_val <- ord_tdata[[input_col]][[low_thresh]]
  high_val <- ord_tdata[[input_col]][[high_thresh]]

  tmp_vec <- list()
  for(elem in titer_data[[input_col]]){
    tmp_vec <- lapply(titer_data[[input_col]], FUN = function(x){
      if(x <= low_val){
        return(0)
      }else if((low_val < x) && (x <= high_val)){
        return(1)
      }else{
        return(2)
      }
    })
  }
  tmp_vec <- unlist(tmp_vec)
  return(tmp_vec)
}

# This method was created from trial and error to mimic the results of the $fc_res_max column
# in the sdy212 using the $fc_norm_max_ivt and $d0_max columns.  The only guidance from the 
# readme.txt on the google drive was: "decorrelate the fc_norm_max_ivt from the d0_max titers"
# Tried unsuccessfully to contact the author / developer: yuri.kotliarov@nih.gov
decorr <- function(d0_val,ivt_val){
  res <- 0
  slope <- 1.135
  pos <- 0.25
  neg <- -0.20
  if(d0_val <= 40){
    # 8 of the 44 original values <= 40 were tuned up ... can't find patter though
    res <- (ivt_val * slope) + neg
  }else{
    res <- (ivt_val * slope) + pos
  }
  return(res)
}

#-----Main method--------------
# for local version ... set dir / filename
# hai_dir <- "/home/ehenrich/R/HIPCGeneModuleAnalysis/HAIgen/"
# file <- "sdy212_v2.tsv"
# read in raw data
#rawdata <- fread(paste0(hai_dir,file))

# For IS driven version
makeHAI <- function(sdy){
  con <- CreateConnection(sdy)
  rawdata <- con$getDataset("hai")
  
  # Generate vectors from rawdata and assign to global environment for use within helper functions
  subids <- unique(rawdata$participant_id)
  strains <- unique(rawdata$virus)
  days_collected <- c(0,28)
  opts <- c("d0_","fc_")
  
  # parse data into df and make df available in global environment
  titer_data <- as.data.frame(matrix(nrow = length(subids), ncol = 27))
  colnames(titer_data) <- c("subject","fc_B","fc_H1N1","fc_H3N2","d0_B","d0_H1N1","d0_H3N2","d28_B","d28_H1N1","d28_H3N2",
                            "d0_std_norm_B","d0_std_norm_H1N1","d0_std_norm_H3N2","fc_std_norm_B","fc_std_norm_H1N1",
                            "fc_std_norm_H3N2","d0_norm_max","fc_norm_max","fc_norm_max_ivt","d0_max","fc_max","fc_max_4fc",
                            "fc_norm_max_d20","fc_norm_max_d30","fc_res_max","fc_res_max_d20","fc_res_max_d30")
  
  iterator <- 1
  
  for(id in subids){
    titer_data[iterator,1] <- id
    for(virus in strains){
      for(day in days_collected){
        col_to_find <- paste0("d",day,"_",virus)
        rowid <- which(rawdata$participant_id == id & rawdata$study_time_collected == day & rawdata$virus == virus)
        if(length(rowid) != 0){
          target_row <- rawdata[rowid, ]
          titer_data[iterator,col_to_find] <- as.integer(target_row$value_reported)
        }else{
          titer_data[iterator,col_to_find] <- NA
        }
      }
    }
    iterator <- iterator + 1
  }
  
  # setup list of lists to hold median and sd values for later use in global envir accessible var
  glob_names <- c("d0_med", "d0_sd","fc_med", "fc_sd")
  li <- vector("list",length = length(strains))
  names(li) <- strains
  glob_vals <- list()
  for(l in glob_names){
    glob_vals[[l]] <- li
  }
  
  # calculate med and sd for d0 cols
  glob_vals <- med_sd_calc("d0_",strains,glob_vals,titer_data)
  
  # calc fold change
  for(virus in strains){
    d0_col <- paste0("d0_",virus)
    d28_col <- paste0("d28_",virus)
    fc_col <- paste0("fc_",virus)
    titer_data[,fc_col] <- titer_data[,d28_col]/titer_data[,d0_col]
  }
  
  # get fold change med and sd
  glob_vals <- med_sd_calc("fc_",strains,glob_vals,titer_data)
  
  #calc standardized and normalized value for each possibility of (d0,fc) x (B,H1N1,H3N2)
  for(ver in opts){
    for(virus in strains){
      std_norm_col <- paste0(ver,"std_norm_", virus)
      titer_data[std_norm_col] <- lapply(titer_data[paste0(ver,virus)], FUN = function(x){ (x-glob_vals[[paste0(ver,"med")]][[virus]])/glob_vals[[paste0(ver,"sd")]][[virus]]})
    }
  }
  
  # easiest to assign snapshots of variables to global_env for use with mapply,
  # otherwise mapply gets confused if strains / titerdata as arguments and tries to iterate
  # through them too
  assign("gl_tdata", titer_data, envir = .GlobalEnv)
  assign("gl_strains", strains, envir = .GlobalEnv)
  
  # Select maxima for (d0,fc) x ("","std_norm") columns
  for(ver in opts){
    titer_data[[paste0(ver,"max")]] <- mapply(max_select,titer_data$subject,ver)
    titer_data[[paste0(ver,"norm_max")]] <- mapply(max_select,titer_data$subject,paste0(ver,"std_norm_"))
  }
  
  # determine fc_max_4fc
  titer_data$fc_max_4fc <- unlist(lapply(titer_data$fc_max, FUN = function(x){if(x > 4){return(1)}else{return(0)}}))
  
  # fc_norm_max_ivt ... it was unclear what the authors did here, but according to the following paper:
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2921808/
  # the 'P' formula for modified ranking is the most common - see fig 1 in 'distinguishing INTs' section
  ranked <- rank(titer_data$fc_norm_max)
  c <- 3/8
  P <- (ranked - c)/(length(ranked) - (2*c) + 1) 
  normed <- qnorm(P)
  titer_data$fc_norm_max_ivt <- normed
  
  # fc_res_max ... see notes in method call
  titer_data$fc_res_max <- mapply(decorr,titer_data$d0_max,titer_data$fc_norm_max_ivt)
  
  # Categorize maxima for ("fc_norm_max","fc_res_max") x ("d20","d30")
  in_cols <- c("fc_norm_max","fc_res_max")
  in_percs <- c(20,30)
  for(cl in in_cols){
    for(perc in in_percs){
      targ_col <- paste0(cl,"_d",perc)
      titer_data[[targ_col]] <- thresh_calc(subids, perc, cl, titer_data)
    }
  }
  
  # Remove d28 columns b/c not present in original versions
  drops <- c("d28_B","d28_H1N1","d28_H3N2")
  titer_data <- titer_data[ , !(names(titer_data) %in% drops)]
  
  # remove extra sdy info on subids
  titer_data$subject <- unlist(lapply(titer_data$subject, FUN = function(x){
    str_sub(x,1,-5)
  }))
  
  # output the matrix as a tab-delimited file
  fname <- paste0(sdy,"_combined_hai_titer_table.txt")
  write.table(titer_data,fname,sep = "\t", col.names = TRUE)
}

