# Evan Henrich
# December 8, 2016
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

# PURPOSE: generate HAI tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

# NOTES: The original code to perform these operations was developed by Yuri Kotliarov at NIH (yuri.kotliarov@nih.gov).  This version
# was developed to provide the same functionality as the original in terms of statistical operation, but use
# data pulled directly from the ImmuneSpace portal at www.immunespace.org instead of data shared among the 
# collaborating labs via a google drive.

#------Dependencies-------------
library(dtplyr)
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

# Method by Yuri Kotliarov to categorize an observation based on low and high percentiles
discretize <- function(df, input_col, low_perc){
  xq <- quantile(df[[input_col]], c((low_perc/100), 1-(low_perc/100)), na.rm=T)
  xd <- ifelse(is.na(df[[input_col]]),NA,1)
  xd[df[[input_col]] <= xq[1]] <-  0
  xd[df[[input_col]] >= xq[2]] <-  2
  return(xd)
}

# for splitting participant_id string to remove sdy info
sub_split <- function(x){
  tmp <- (strsplit(x, split = "[.]"))
  return(tmp[[1]][1])
}

#-----Main method--------------
# for local version ... set dir / filename
# hai_dir <- "/home/ehenrich/R/HIPCGeneModuleAnalysis/HAIgen/"
# file <- "sdy212_v2.tsv"
# read in raw data
# rawdata <- fread(paste0(hai_dir,file))

# For IS driven version
makeHAI <- function(sdy){
  con <- CreateConnection(sdy)
  rawdata <- con$getDataset("hai")
  
  # Generate vectors from rawdata or instantiate variables based on nomenclature of original column headers
  subids <- unique(rawdata$participant_id)
  strains <- unique(rawdata$virus)
  cohorts <- unique(rawdata$cohort)
  days_collected <- c(0,28)
  opts <- c("d0_","fc_")
  
  # Different nomenclature in sdy212 vs rest of studies in terms of strain naming, therefore need mapping
  strains_adj <- list()
  for(x in strains){
    if(x == "B/Brisbane/60/2008"){
      strains_adj[["B"]] <- x
    }else if(x == "A/California/7/2009"){
      strains_adj[["H1N1"]] <- x
    }else if(x == "A/Perth/16/2009"){
      strains_adj[["H3N2"]] <- x
    }else{
      strains_adj[[x]] <- x
    }
  }
  
  # Setup df with columns, including those that will eventually be removed to mimic original datasets (e.g. d28_x)
  titer_data <- as.data.frame(matrix(nrow = length(subids), ncol = 28))
  colnames(titer_data) <- c("subject","cohort","fc_B","fc_H1N1","fc_H3N2","d0_B","d0_H1N1","d0_H3N2","d28_B","d28_H1N1","d28_H3N2",
                            "d0_std_norm_B","d0_std_norm_H1N1","d0_std_norm_H3N2","fc_std_norm_B","fc_std_norm_H1N1",
                            "fc_std_norm_H3N2","d0_norm_max","fc_norm_max","fc_norm_max_ivt","d0_max","fc_max","fc_max_4fc",
                            "fc_norm_max_d20","fc_norm_max_d30","fc_res_max","fc_res_max_d20","fc_res_max_d30")
  
  # Parse day 0 and day 28 data into df as basis for all future calculations
  iterator <- 1
  for(id in subids){
    titer_data[iterator,1] <- id
    cohort_df <- rawdata[which(rawdata$participant_id == id),]
    cohort_val <- unique(cohort_df$cohort)
    titer_data[iterator,2] <- cohort_val
    for(vir_name in names(strains_adj)){
      for(day in days_collected){
        col_to_find <- paste0("d",day,"_",vir_name)
        # b/c day 28 is not always exactly day 28, it is assumed that any non-zero value represents the 28th day value
        if(day == 0){
          rowid <- which(rawdata$participant_id == id & rawdata$study_time_collected == day & rawdata$virus == strains_adj[[vir_name]])
        }else{
          rowid <- which(rawdata$participant_id == id & rawdata$study_time_collected != 0 & rawdata$virus == strains_adj[[vir_name]])
        }
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
  
  # setup list to hold median and sd values for later use in calculations
  glob_names <- c("d0_med", "d0_sd","fc_med", "fc_sd")
  li <- vector("list",length = length(strains))
  names(li) <- names(strains_adj)
  glob_vals <- list()
  for(l in glob_names){
    glob_vals[[l]] <- li
  }
  
  # calc median and sd for d0 cols
  glob_vals <- med_sd_calc("d0_",names(strains_adj),glob_vals,titer_data)
  
  # calc fold change
  for(vir_name in names(strains_adj)){
    d0_col <- paste0("d0_",vir_name)
    d28_col <- paste0("d28_",vir_name)
    fc_col <- paste0("fc_",vir_name)
    titer_data[,fc_col] <- titer_data[,d28_col]/titer_data[,d0_col]
  }
  
  # calc fold change med and sd
  glob_vals <- med_sd_calc("fc_",names(strains_adj),glob_vals,titer_data)
  
  #calc standardized and normalized value for each possibility of (d0,fc) x (B,H1N1,H3N2)
  for(ver in opts){
    for(vir_name in names(strains_adj)){
      std_norm_col <- paste0(ver,"std_norm_", vir_name)
      titer_data[std_norm_col] <- lapply(titer_data[paste0(ver,vir_name)], FUN = function(x){ 
        (x-glob_vals[[paste0(ver,"med")]][[vir_name]])/glob_vals[[paste0(ver,"sd")]][[vir_name]]})
    }
  }
  
  # Assign snapshots of variables to global_env for use with mapply.
  assign("gl_tdata", titer_data, envir = .GlobalEnv)
  assign("gl_strains", names(strains_adj), envir = .GlobalEnv)
  
  # Select maxima for (d0,fc) x ("","std_norm") columns
  for(ver in opts){
    titer_data[[paste0(ver,"max")]] <- mapply(max_select,titer_data$subject,ver)
    titer_data[[paste0(ver,"norm_max")]] <- mapply(max_select,titer_data$subject,paste0(ver,"std_norm_"))
  }
  
  # determine fc_max_4fc, which is categorization based on fold change > 4
  titer_data$fc_max_4fc <- unlist(lapply(titer_data$fc_max, FUN = function(x){if(x > 4){return(1)}else{return(0)}}))
  
  # Inverse normal transformation of standardized/normalized max fold change column
  # Done by quantile normalization on a modified ranking of observations
  # fc_norm_max_ivt << provided from Yuri Kotliarov
  ranked <- rank(titer_data$fc_norm_max, na.last = "keep")
  P <- ranked / (sum(!is.na(ranked)) + 1) # +1 needed to avoid 0,1 values that generate inf, -inf
  titer_data$fc_norm_max_ivt <- qnorm(P)
  
  # setup meta-list for possible titer tables based on cohorts: young, old, and combined are possible
  submxs <- list()
  submxs[["combined"]] <- titer_data
  
  # Generate subset matrices based on age and perform statistical work on each separately
  # SDY212 and the other studies use different nomenclature for categorization, therefore 
  # need to check against lists
  yng_ls <- c("Cohort_1","Young adults 21-30 years old")
  old_ls <- c("Cohort_2","Older adults >= 65 years old")
  for(coh in cohorts){
    if(coh %in% yng_ls){
      submxs[["young"]]  <- titer_data[which(titer_data$cohort == coh),]
    }else if(coh %in% old_ls){
      submxs[["old"]]  <- titer_data[which(titer_data$cohort == coh),]
    }
  }
  
  # fc_res_max is generated by first binning the subjects' d0_norm_max data 
  # according to a manual selection and then normalizing/standardizing 
  # the inverse normal transformation values within the bins.  This code
  # is based on Yuri Kotliarov's work.
  bins <- c(1,5)
  for(name in names(submxs)){
    df <- submxs[[name]]
    df$bin <- cut(df$d0_norm_max, breaks=c(-Inf, bins, Inf), labels=1:(length(bins)+1))
    df = df %>%
      group_by(bin) %>%
      mutate(fc_res_max = (fc_norm_max_ivt - median(fc_norm_max_ivt, na.rm=T)) / sd(fc_norm_max_ivt, na.rm=T)) %>%
      ungroup()
    
    # discretize for all combinations of ("fc_norm_max","fc_res_max") x ("d20","d30")
    in_cols <- c("fc_norm_max","fc_res_max")
    in_percs <- c(20,30)
    for(cl in in_cols){
      for(perc in in_percs){
        targ_col <- paste0(cl,"_d",perc)
        df[[targ_col]] <- discretize(df, cl, perc)
      }
    }
    
    # Remove d28, bin, and cohort columns b/c not present in original versions
    drops <- c("d28_B","d28_H1N1","d28_H3N2","bin","cohort")
    df <- df[ , !(names(df) %in% drops)]
    
    # remove extra sdy info on subids
    df$subject <- unlist(lapply(df$subject, sub_split))
    
    # remove rows with NA values caused by lack of second titer BUT NOT those that may 
    # have NaN values in other places (e.g. $fc_res_max for SUB120471 in sdy404) due to 
    # other reasons such as extreme d0 titer values
    df <- df[which(!is.na(df$fc_B)),]
    
    # output tibble / df as a tab-delimited file
    base <- "_hai_titer_table_EH.txt"
    fname <- paste0(sdy,"_",name,base)
    write.table(df,fname,sep = "\t", col.names = TRUE)
  }
}
  
  
  

