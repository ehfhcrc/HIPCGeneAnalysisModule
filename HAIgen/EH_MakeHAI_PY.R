# Evan Henrich
# December 8, 2016
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

# PURPOSE: generate HAI tables for downstream HIPC Gene Module Analysis using raw data from Immunespace

# NOTES: The original code to perform these operations was developed by Yuri Kotliarov at NIH (yuri.kotliarov@nih.gov).  This version
# was developed to provide the same functionality as the original in terms of statistical operation, but use
# data pulled directly from the ImmuneSpace portal at www.immunespace.org instead of data shared among the 
# collaborating labs via a google drive and also to handle all the studies used in the meta-analysis in an automated format.

#------Dependencies-------------
# library(ImmuneSpaceR)
# library(data.table)
# library(dplyr)
# library(stringr)

# for local version ... set dir / filename
# hai_dir <- "/home/ehenrich/R/HIPCGeneModuleAnalysis/HAIgen/"
# file <- "sdy212_v2.tsv"
# read in raw data
# rawdata <- fread(paste0(hai_dir,file))


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

# Seems computationally costly, but enough edge cases merit checking for valid data before extraction
sub_check <- function(sub_df, strains){
  result <- list()
  for(vir in strains){
    d_zero <- sub_df[which(sub_df$virus == vir & sub_df$study_time_collected == 0),]
    d_other <- sub_df[which(sub_df$virus == vir & sub_df$study_time_collected != 0),]
    if(nrow(d_zero) > 0 & nrow(d_other) > 0){
      result[[vir]] <- TRUE
    }else{
      result[[vir]] <- FALSE
    }
  }
  final <- !(FALSE %in% result)
  return(final)
}

max_select <- function(subid,trg_col){
  tmp_ls <- list()
  for(virus in gl_strains){
    tmp_ls[virus] <- gl_tdata[gl_tdata$subject == subid, paste0(trg_col,virus)]
  }
  return(max(unlist(tmp_ls),na.rm = TRUE))
}

# Method by Yuri Kotliarov to categorize an observation based on low and high percentiles
# Changed slightly to round fc_res_max values to 7 digits prior to comparison with quantile
# values which are interpolated and therefore can throw off assignment if intention is to 
# use them as if they are nearest order statistic (similar to type = 3 in quantiles args).
discretize <- function(df, input_col, low_perc){
  xq <- quantile(df[[input_col]], 
                 c( (low_perc/100), 1 - (low_perc/100) ), 
                 na.rm=T,
                 type = 7)
  xd <- sapply(df[[input_col]], FUN = function(x){
    if(is.na(x)){
      return(NA)
    }else if(round(x, digits = 7) <= round(xq[[1]], digits = 7)){
      return(0)
    }else if(round(x, digits = 7) >= round(xq[[2]], digits = 7)){
      return(2)
    }else{
      return(1)
    }
  })
}

# for splitting participant_id string to remove sdy info
sub_split <- function(x){
  tmp <- (strsplit(x, split = "[.]"))
  return(tmp[[1]][1])
}

drop_cols <- function(df, cols_to_drop){
  df <- df[ , !(names(df) %in% cols_to_drop)]
  return(df)
}

#-----Main method--------------

# For IS driven version
makeHAI <- function(sdy){
  
  # Setup directory vars
  wk_dir <- getwd()
  hai_dir <- file.path(wk_dir,"RDSGen/HAI_preproc_data")
  
  # Get rawdata from ImmuneSpace
  con <- CreateConnection(sdy)
  if(sdy == "SDY80"){
    rawdata <- con$getDataset("neut_ab_titer")
  }else{
    rawdata <- con$getDataset("hai")
  }
  
  # Generate vectors from rawdata or instantiate variables based on nomenclature of original column headers
  subids <- unique(rawdata$participant_id)
  strains <- unique(rawdata$virus)
  strains <- sapply(strains, FUN = function(x){
    x <- gsub("\\.|\\/| |-", "_", x)
    return(x)
  })
  cohorts <- unique(rawdata$cohort)
  days_collected <- c(0,28)
  opts <- c("d0_","fc_")
  
  # Setup df with columns, including those that will eventually be removed
  str_names <- list()
  str_d28_names <- list()
  for(virus in strains){
    for(opt in opts){
      str_names <- c(str_names, paste0(opt, virus))
      str_names <- c(str_names, paste0(opt, "std_norm_", virus))
    }
    str_d28_names <- c(str_d28_names, paste0("d28_", virus))
  }
  
  first_names <- c("subject","cohort")
  last_names <- c("d0_norm_max","fc_norm_max","fc_norm_max_ivt", "d0_max",
                  "fc_max","fc_max_4fc", "fc_norm_max_d20","fc_norm_max_d30",
                  "fc_res_max","fc_res_max_d20","fc_res_max_d30")
  
  numcol <- length(str_names) + length(str_d28_names) + length(first_names) + length(last_names)
  titer_data <- data.frame(matrix(vector(), 
                                  nrow = 0, 
                                  ncol = numcol), 
                           stringsAsFactors = F)
  colnames(titer_data) <- c(first_names, str_names, str_d28_names, last_names)
  
  # Parse day 0 (initial) and day 28 (follow-up) titer data into df as basis for all future calculations.
  # NOTE: Because the follow-up titer measurement is not always collected on day 28 exactly, 
  # it is assumed that any value greater than 0 represents the 28th day value if there is 
  # not a day 28 value present. E.g. SDY 80 / CHI-nih has negative values possible for days collected
  # as well as possible days 7 and 70. The first of the days over 0 are used because this is what
  # was done in the original file as determined by comparing fold change values for SUB114476.80.
  iterator <- 1
  for(id in subids){
    sub_data <- rawdata[which(rawdata$participant_id == id),]
    valid <- sub_check(sub_data, names(strains))
    # num_strains <- nrow(sub_data)/2
    # second_titer <- max(unique(sub_data$study_time_collected)) > 0
    # # Do not allow those without second titers for each strain to generate record
    if(valid == TRUE){
      titer_data[iterator,1] <- id
      cohort_df <- rawdata[which(rawdata$participant_id == id),]
      cohort_val <- unique(cohort_df$cohort)
      titer_data[iterator,2] <- cohort_val
      for(vir_name in names(strains)){
        for(day in days_collected){
          col_to_find <- paste0("d", day, "_", strains[[vir_name]])
          rowid <- which(rawdata$participant_id == id & 
                           rawdata$study_time_collected == day & 
                           rawdata$virus == vir_name)
          if(length(rowid) == 0){
            rowid <- which(rawdata$participant_id == id & 
                             rawdata$study_time_collected > 0 & 
                             rawdata$virus == vir_name)
          }
          if(length(rowid) > 1){
            rowid <- rowid[[1]]
          }
          target_row <- rawdata[rowid, ]
          
          # for testing:
          # print(paste(id, col_to_find))
          
          if(sdy == "SDY80" & target_row$value_reported == 19){ # To mimic original table
            titer_data[iterator, col_to_find] <- 10
          }else{
            titer_data[iterator, col_to_find] <- as.integer(target_row$value_reported)
          }
          
        }
      }
      iterator <- iterator + 1
    }
  }
  
  # setup list to hold median and sd values for later use in calculations
  glob_names <- c("d0_med", "d0_sd","fc_med", "fc_sd")
  li <- vector("list",length = length(strains))
  names(li) <- strains
  glob_vals <- list()
  for(l in glob_names){
    glob_vals[[l]] <- li
  }
  
  # calc median and sd for d0 cols
  glob_vals <- med_sd_calc("d0_", strains , glob_vals, titer_data)
  
  # calc fold change
  for(virus in strains){
    d0_col <- paste0("d0_", virus)
    d28_col <- paste0("d28_", virus)
    fc_col <- paste0("fc_", virus)
    titer_data[,fc_col] <- titer_data[,d28_col]/titer_data[,d0_col]
  }
  
  # calc fold change med and sd
  glob_vals <- med_sd_calc("fc_", strains , glob_vals, titer_data)
  
  #calc standardized and normalized value for each possibility of (d0,fc) x (B,H1N1,H3N2)
  for(ver in opts){
    for(virus in strains){
      std_norm_col <- paste0(ver, "std_norm_", virus)
      titer_data[std_norm_col] <- lapply(titer_data[paste0(ver, virus)], FUN = function(x){ 
        ( x - glob_vals[[paste0(ver,"med")]][[virus]]) / glob_vals[[paste0(ver,"sd")]][[virus]]} )
    }
  }
  
  # Assign snapshots of variables to global_env for use with mapply.
  assign("gl_tdata", titer_data, envir = .GlobalEnv)
  assign("gl_strains", strains, envir = .GlobalEnv)
  
  # Select maxima for (d0,fc) x ("","std_norm") columns
  for(ver in opts){
    titer_data[[paste0(ver,"max")]] <- mapply(max_select, titer_data$subject, ver)
    titer_data[[paste0(ver,"norm_max")]] <- mapply(max_select, titer_data$subject, paste0(ver,"std_norm_"))
  }
  
  # determine fc_max_4fc, which is categorization based on fold change > 4
  titer_data$fc_max_4fc <- unlist(lapply(titer_data$fc_max, FUN = function(x){
    if(x > 4){return(1)}else{return(0)}
    }))
  
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
  yng_ls <- c("Cohort_1", "Young adults 21-30 years old")
  old_ls <- c("Cohort_2", "Cohort2", "Older adults >= 65 years old", "healthy adults, 50-74 yo")
  for(coh in cohorts){
    if(coh %in% yng_ls){
      submxs[["young"]]  <- titer_data[which(titer_data$cohort == coh),]
    }else if(coh %in% old_ls){
      submxs[["old"]]  <- titer_data[which(titer_data$cohort == coh),]
    }
  }
  
  # fc_res_max is generated by first binning the subjects' d0_norm_max data 
  # according to a manual selection and then normalizing/standardizing 
  # the inverse normal transformation values within bins that have been selected manually.
  # This code is based on Yuri Kotliarov's work.
  SDY404_bins <- list(c(1,5),c(1,5),c(0))
  SDY212_bins <- list(c(0.01),c(0.05),c(0.01))
  SDY400_bins <- list(c(1,5),c(1,5),c(1,5))
  SDY63_bins <- list(c(2,6),c(2),c(2))
  # Following were done by trial and error, b/c not noted in pngs in google drive or Yuri's code
  # Method of finding them was to make data.frame($fc_res_max,$fc_norm_max_ivt,$d0_norm_max),
  # from original file then df <- df[order($fc_norm_max_ivt,$d0_norm_max),] .
  # by looking at where ivt is the same but fc_res_max was different, one can guess that they
  # were therefore in different bins and estimate the cutoff points by looking at the 
  # d0_norm_max column
  SDY67_bins <- list(c(0.1,0.5,1.5,4,6),c(),c(0.1,0.5,1.5,4,6))
  SDY80_bins <- list(c(1,5),c(),c())
  
  bname <- paste0(sdy,"_bins")
  bins <- get(bname)
  names(bins) <- c("combined", "young", "old")
  
  for(name in names(submxs)){
    df <- submxs[[name]]
    df$bin <- cut(df$d0_norm_max, 
                  breaks = c(-Inf, bins[[name]], Inf), 
                  labels = 1:(length(bins[[name]])+1) )
    df = df %>%
      group_by(bin) %>%
      mutate(fc_res_max = 
               (fc_norm_max_ivt - median(fc_norm_max_ivt, na.rm=T)) / sd(fc_norm_max_ivt, na.rm=T)) %>%
      ungroup()
    
    # discretize for all combinations of ("fc_norm_max","fc_res_max") x ("d20","d30")
    in_cols <- c("fc_norm_max", "fc_res_max")
    in_percs <- c(20, 30)
    for(cl in in_cols){
      for(perc in in_percs){
        targ_col <- paste0(cl,"_d", perc)
        df[[targ_col]] <- discretize(df, cl, perc)
      }
    }
    
    # Remove d28, bin, and cohort columns b/c not present in results of original manual versions.
    df <- drop_cols(df, c(str_d28_names, "bin", "cohort"))

    # remove extra sdy info on subids
    df$subject <- unlist(lapply(df$subject, sub_split))
    
    # output tibble / df as a tab-delimited file
    base <- "_hai_titer_table_EH.txt"
    fname <- ""
    if(sdy == "SDY80"){
      fname <- paste0("CHI-nih_", name, base)
    }else{
      fname <- paste0(sdy, "_", name, base)
    }
    
    write.table(df, file = file.path(hai_dir,fname), sep = "\t", col.names = TRUE)
  }
}
  
  
  

