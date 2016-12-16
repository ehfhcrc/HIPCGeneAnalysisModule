Evan Henrich
ehenrich@fredhutch.org
December 16, 2015

Immune Signatures Project - Automated Analysis Script

The script IS_auto.R pulls raw data for gene expression, hai, and demographics from www.immunespace.org for the studies used in the Immune Signatures Meta Analysis project using the ImmuneSpaceR wrapper, then proceeds to process and run the meta-analysis code.

The process of extracting, processing, and analyzing the raw data has the following steps and relies on the indicated code:
1. Extract and preprocess data from ImmuneSpace using 'makeGEmatrix.R','makeHAItable.R', and 'makeDemo.R' scripts
2. Combine data for each study in an rds file that holds a bioconductor eset object (Expression Set) using `makeRds.R`
3. Analyze the rds files with 'HIPCMetaModuleAnalysis_v2.R' to generate the plots and figures seen in the published meta-analysis study

To find information about the scripts, such as the original author and any changes that were made to ensure replication, please see the comments in the script itself.

To run the IS_auto.R script, you first need to do the following:

1. Create a working directory called "ImmSig"

2. Download the following files and scripts:
    ***FILES***
    --BTM_for_GSEA_20131008.GMT
      >>file containing all BTM (Blood Transcript Module) gene modules
    --bs_to_id_final.tsv
      >>A tab-delimited table to map biosample ids to ImmPort subject / participant ids for SDY212
    
    ***SCRIPTS***
    --HIPCMetaModuleAnalysis_v2.R
    --getSDY.R
    --makeGEmatrix.R
    --makeHAItable.R
    --makeDemo.R
    --makeRds.R
    
3. Setup sub-directories and place appropriate files in those directories as seen below
-ImmSig (working directory)
   --data
       --BTM_for_GSEA_20131008.GMT
   --output
   --HIPCMetaModuleAnalysis_v2.R
   --getSDY.R
   --PreProc_Scripts
      --makeGEmatrix.R
      --makeHAItable.R
      --makeDemo.R
      --bs_to_id_final.tsv
      --Yale IDs to immsig
   --RDSGen
      --makeRds.R
      --GE_preproc_data
      --HAI_preproc_data
   --IS_auto.R
      
4. in R console ..
> setwd("/ImmSig")
> source("IS_auto.R")

Please note that package versions and operating systems (32 vs 64 bit) play an important role in producing the same results across different machines.  The work done to pull this code together was done on a Kubuntu v16.04 OS (linux 64bit) and package versions are listed below.  The output files from each script are available to download at:

sessionInfo()

