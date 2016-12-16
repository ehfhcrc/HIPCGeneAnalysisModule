##############################################################################
##
## Script normalizes SDY212 Illumnia data and creates data.frame
##
##############################################################################


if(file.exists("~/.Rprofile")) source("~/.Rprofile");
library(limma);
library(lumi);
library(annotate);
library(lumiHumanAll.db);

#library(beadarray);




 
dataFile = "SDY212_WholeBlood_Microarray_update_11242014.txt";
targetFile = "subjectID_to_bioSampID.txt"; 

 
if(!exists("x.lumi")) {
  ptm <- proc.time();
  cat("\nReading in raw data\n");
  x.lumi <- lumiR(dataFile, checkDupId=FALSE);
  targets <- read.delim(targetFile, as.is=T);
  deltaTime = proc.time() - ptm;
  cat(paste("Elapsed time:", round(deltaTime[3],2),"seconds\n"));
}

probe_ids <- featureData(x.lumi)$PROBE_ID;
gs <- featureData(x.lumi)$SYMBOL;
map <- data.frame(probe_ids, gs, stringsAsFactors=FALSE);

expr <- log2(exprs(x.lumi));

cat("\nQuantile Normalizing the data\n");
ptm <- proc.time();
expr <- normalizeBetweenArrays(as.matrix(expr), method="quantile");
deltaTime = (proc.time() - ptm);
cat("Elapsed time:", round(deltaTime[3],2)," seconds\n");
rm(ptm,deltaTime);


## There is 1 biological sample (column) that occurs twice.
## The offending sample is "BS694717".
## We'll label the relevant columns twice so others can figure out 
## what to do with it

## data file mapping Stanford ID's to ImmPort IDs
subjects <- read.delim(targetFile, as.is=T);

g <- grep("BS694717", subjects[,2]);
subjects <- rbind(subjects,subjects[g,]);
n <- nrow(subjects);
subjects[g,] <- paste(subjects[g,],".1",sep="")
subjects[n,] <- paste(subjects[n,],".2",sep="")

subjects$idx = NA;

for(i in seq(nrow(subjects))) {
    bs = subjects[i,2];
    sub = subjects[i,1];
    g = grep(bs, colnames(expr));

    if(length(g) == 2) next;    
    if(length(g) == 1) subjects$idx[i] = g;
}


sub2 <- subjects[ !is.na(subjects$idx), ];
expr2 <- expr[ ,sub2$idx ]
colnames(expr2) <- sub2$subject_accession;

out <- data.frame(
	probeID=probe_ids,
	geneSymbol=gs,
	round(expr2,6),
	stringsAsFactors=FALSE);

write.table(out, file="SDY212.GEMatrix.txt", row.names=F, sep="\t", quote=F);





