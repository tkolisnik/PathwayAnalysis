# PathwayAnalysisFromFoldChanges.R

# Created by Tyler Kolisnik and Mark Bieda.
# April 7

# Directions: 
# 1. Install pathview and gage packages from bioconductor if not installed already. 
# 2. Adjust input parameters in the parameter block
# 3. Run the Program

# Important Information:
# This program takes a tab-separated text file as input, with two columns, the first row is a header labeled EntrezIDs /t FoldChange and the subsequent rows being the numerical Entrez Ids \t numerical log2 FoldChange values.
# Output is two lists (One for upregulated pathways, One for downregulated pathways, as well as a folder of image files, KeggNative and Graphviz (when available) Visualizations are included.

########################### Begin Parameter Block ###########################
inputFile <- ("/home/tylerk/corticosteroid562/biopsy_array_data_entrezratio.csv") # Must be a text file with entrezIDs and log2 Fold Change values. 
outputDirectory <- ("/home/tylerk/corticosteroid562/output") # The absolute path files will be outputted to, with no trailing /
runName <- "Corticosteroid562" # The name that will be appended to all output files. 
pValueCutoff <- 0.2 # The pValue threshold for the Fold Change. 
KEGGspeciesCode <- "hsa" # Kegg species code, default = "hsa" (human), more available on bioconductor website. 
pathwayDirectoryName <- "pathway_images_562"
########################### End Parameter Block ###########################

########################### Function Block ###########################
#Max amplitude function
biggestamplitude <- function(x){
  #x is a numeric vector (or single number, but single number case is trivial)
  pos <- max(x)
  neg <- min(x)
  if (abs(pos)>=abs(neg))
  {return(pos)}
  else
  {return(neg)}
}
########################### End Function Block ###########################

setwd(outputDirectory)
#Load Dependencies
library(pathview)
library(gage) 

pathwayDirectory <- paste(outputDirectory, pathwayDirectoryName, sep="/")

# Reads in the file
geneData <- read.delim(inputFile, header=TRUE, sep="\t")
geneDataLog2FoldChange <- log2(geneData$FoldChange)
geneData2 <- geneData
geneData2$FoldChange <- geneDataLog2FoldChange
colnames(geneData2)[2] <- "Log2FoldChange"

geneList <- unique(geneData2$EntrezID)
geneList <- geneList[!is.na(geneList)]
geneListChar <- as.character(geneList)


# Removes blank spaces, ---, entrez ids with /// in them. This program cannot handle these. 
goodGeneListRows <- !grepl("[^0-9]+", geneListChar) & grepl("[0-9]+",geneListChar)
geneListAdj <- geneListChar[goodGeneListRows]
dataMatrix <- matrix(nrow=length(geneListAdj), ncol=1) 
rownames(dataMatrix) <- geneListAdj
colnames(dataMatrix) <- "Log2FoldChange"


# For each unique Entrez ID, finds the highest fold change associated with that gene
for (i in 1:length(geneListAdj)) {
  gene <- as.character(geneListAdj[i])
  biggestampFoldChange <- biggestamplitude(geneData2$Log2FoldChange[which(geneData2$EntrezID == gene)])
  dataMatrix[gene,] <- biggestampFoldChange
}

dataMatrixNotLog2 <- 2^dataMatrix
data(kegg.gs) # load in the KEGG pathway data
gageOutput <- gage(dataMatrixNotLog2, gsets = kegg.gs, ref=NULL, samp=NULL)
#look at what gage expects


# Calculate pathways that were Up-Regulated
pathwayData <- gageOutput$greater[,c("p.val", "q.val")]
pathCode <- substr(rownames(pathwayData), 1, 8)
pathName <- substr(rownames(pathwayData), 9, length(rownames(pathwayData)))
pathwayData <- cbind(pathCode, pathName ,pathwayData)
colnames(pathwayData) <- c("Pathway Code", "Pathway Name", "p Values", "q value")

# Output a text file of upregulated pathways and meaningful data
write.table(pathwayData, paste(runName, "upregulated_pathway_list.txt", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

# Calculate pathways that were Down-Regulated
pathwayDataDownRegulated <- gageOutput$less[,c("p.val", "q.val")] #DO THE LESSER ALSO
pathCodeDown <- substr(rownames(pathwayDataDownRegulated), 1, 8)
pathNameDown <- substr(rownames(pathwayDataDownRegulated), 9, length(rownames(pathwayDataDownRegulated)))
pathwayDataDown <- cbind(pathCodeDown, pathNameDown, pathwayDataDownRegulated)
colnames(pathwayDataDown) <- c("Pathway Code", "Pathway Name", "p Values", "q value")

# Output a text file of downregulated pathways and meaningful data
write.table(pathwayDataDown, paste(runName, "downregulated_pathway_list.txt", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

# Create a folder within the outputDirectory where all of the generated pathway image files will go. Note: Using standard parameters, it is likely several hundred image files will be generated.
system(paste("mkdir ", pathwayDirectory, sep = ""))
setwd(pathwayDirectory)


## Use pathview to create output images

runPathview <- function(pid, upordown){
  # Output Kegg Native Graphs (.png files, works for all Pathways)
  pathview(gene.data=dataMatrix, pathway.id=pid, kegg.native=T, species=KEGGspeciesCode,  low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20),plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="max", out.suffix=paste(runName, "KEGGnative_",upordown,sep=""))
  # Output Graphviz Graphs (.pdf files, doesn't work for all Pathways)
  pathview(gene.data=dataMatrix, pathway.id=pid, kegg.native=F, species=KEGGspeciesCode, low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20), plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="max", out.suffix=paste(runName, "graphviz_",upordown,sep=""))
}

# Select Up Regulated Pathways to be printed based on p value (p.val)
selUp <- gageOutput$greater[,"p.val"] < pValueCutoff & !is.na(gageOutput$greater[,"p.val"])
pathIdsUp <- substr(rownames(gageOutput$greater)[selUp], 1, 8)

# Print Selected Up Regulated Pathways 
# Try Catch allows for error handling if a Graphviz visualization does not exist for the pathway ID entered.
pathListUp <- sapply(pathIdsUp, function(pid) tryCatch(runPathview(pid, "upRegulated"), error=function(e) print("This pathway not found")))

## Select Down Regulated Pathways to be printed based on p value (p.val)
selDown <- gageOutput$less[,"p.val"] < pValueCutoff & !is.na(gageOutput$less[,"p.val"])
pathIdsDown <- substr(rownames(gageOutput$less)[selDown], 1, 8)

# Print Selected Down Regulated Pathways 
# Try Catch allows for error handling if a Graphviz visualization does not exist for the pathway ID entered.
pathListDown <- sapply(pathIdsDown, function(pid) tryCatch(runPathview(pid, "downRegulated"), error=function(e) print("This pathway not found")))

