
The script HIPCMetaModuleAnalysis.R performed gene module meta-analysis on discovery cohorts to detect 
signiture gene modules significantly defferent between high responders and low responders to flu vaccination 
in HIPC signature project. The identifiedy gene modules were then verified by the validation cohort.          

To start using gene module analysis, you need to first download and install R package "qusage" from Bioconductor:
source("https://bioconductor.org/biocLite.R")
biocLite("qusage")

The working directory structure can be set up as follows:

-Working direcory
   --data
       --SDY212.rds
       --SDY63.rds
       --SDY404.rds
       --SDY400.rds
       --CHI-nih.rds
       --SDY67-batch2.rds
       --BTM_for_GSEA_20131008.GMT
       
   --output
   --HIPCMetaModuleAnalysis.R
   --getSDY.R
   
All .rds file can be downloaded from HIPC signature shared folder. 

getSDY.R is the script from Renaud to read R object files.

BTM_for_GSEA_20131008.GMT is the file containing all BTM (Blood Transcript Module) gene modules.

Gene module meta-anlaysis can be performed separately for Young or Older cohorts by set up the parameter 
cohort = 'young' or cohort ='old', at the beginging of HIPCMetaModuleAnalysis.R  

There are three major steps to dicover and validate signiture gene modules:

1. Run gene module anlaysis independently on each data set in discovery cohorts.
2. Run meta-analysis by using the results from step 1.
3. Validate signature gene modules in validation cohort.