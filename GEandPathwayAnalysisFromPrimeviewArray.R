#This program annotates an affymetrix primeview gene expression microarray, calculate the gene expression information, and link this
#to visual KEGG pathways

#Created By: Tyler Kolisnik

library(affy)
library(human.db0)
library(AnnotationDbi)
library(AnnotationForge)
library(gage)
library(pathview)
setwd("/home/tylerk/primeview_csv/primeview.db/")
library(primeview.db)
setwd("/home/tylerk/562inputdata")
mytestcond <- "Steroid562"
baseline <- "Placebo"
runName <- "steroid562placebo"
pathwayDirectory <- "/home/tylerk/primeview/pathways"
workingDirectory <- "/home/tylerk/562outputdata"
pValueCutoff <- 0.2
KEGGspeciesCode <- "hsa"
#makeDBPackage("HUMANCHIP_DB", affy = TRUE, prefix="primeview",fileName="PrimeView.na35.annot.csv", baseMapType="gbNRef",outputDir=".", version = "12.1.0", manufacturer="Affymetrix", manufacturerUrl="http://www.affymetrix.com")


##########################################33
dat <- ReadAffy()
eset <- rma(dat)

expressionData <- exprs(eset)


mappedIDs <- select(primeview.db, rownames(expressionData), c("PROBEID","ENTREZID"))
mergedData<- merge(expressionData,mappedIDs,by.x="row.names",by.y="PROBEID")

#xxysmall <- xxy[1:10,]
mergedDataMatrix<-data.matrix(mergedData)
foldchanges<- apply(mergedDataMatrix, 1, function(x) {mean(x[2:12])-mean(x[13:23])})
mergedData$FoldChanges <-foldchanges


mergedData_sorted <- mergedData[order(mergedData$ENTREZID, -abs(mergedData$FoldChanges) ), ]
mergedData_nodupes <- mergedData_sorted[ !duplicated(mergedData_sorted$ENTREZID), ] 
mergedData_noNAs <- mergedData_nodupes[!is.na(mergedData_nodupes$ENTREZID),]
rownames(mergedData_noNAs)<- mergedData_noNAs$ENTREZID
mergedData_forGage <- mergedData_noNAs[,c(-1,-24,-25)]

data_forPathview <- mergedData_noNAs[-1:-24]
dataMatrix <- mergedData_forGage 
dataMatrixNotLog2 <- 2^mergedData_forGage

data(kegg.gs) # load in the KEGG pathway data
gageOutput <- gage(dataMatrix, gsets = kegg.gs, ref=12:22, samp=1:11, compare = "paired")
setwd(workingDirectory)
# Calculate pathways that were Up-Regulated
pathwayData <- gageOutput$greater[,c("p.val", "q.val")]
pathCode <- substr(rownames(pathwayData), 1, 8)
pathName <- substr(rownames(pathwayData), 9, length(rownames(pathwayData)))
pathwayData <- cbind(pathCode, pathName ,pathwayData)
colnames(pathwayData) <- c("Pathway Code", "Pathway Name", "pValues", "qValues")

# Output a text file of upregulated pathways and meaningful data
write.table(pathwayData, paste(runName, "upregulated_pathway_list.txt", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

# Calculate pathways that were Down-Regulated
pathwayDataDownRegulated <- gageOutput$less[,c("p.val", "q.val")] #DO THE LESSER ALSO
pathCodeDown <- substr(rownames(pathwayDataDownRegulated), 1, 8)
pathNameDown <- substr(rownames(pathwayDataDownRegulated), 9, length(rownames(pathwayDataDownRegulated)))
pathwayDataDown <- cbind(pathCodeDown, pathNameDown, pathwayDataDownRegulated)
colnames(pathwayDataDown) <- c("Pathway Code", "Pathway Name", "pValues", "qValues")

# Output a text file of downregulated pathways and meaningful data
write.table(pathwayDataDown, paste(runName, "downregulated_pathway_list.txt", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

# Create a folder within the outputDirectory where all of the generated pathway image files will go. Note: Using standard parameters, it is likely several hundred image files will be generated.
system(paste("mkdir ", pathwayDirectory, sep = ""))
setwd(pathwayDirectory)


## Use pathview to create output images

runPathview <- function(pid, upordown){
  # Output Kegg Native Graphs (.png files, works for all Pathways)
  pathview(gene.data=data_forPathview, pathway.id=pid, kegg.native=T, species=KEGGspeciesCode,  low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20),plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="max", out.suffix=paste(runName, "KEGGnative_",upordown,sep=""))
  # Output Graphviz Graphs (.pdf files, doesn't work for all Pathways)
  #pathview(gene.data=data_forPathview, pathway.id=pid, kegg.native=F, species=KEGGspeciesCode, low = list(gene = "#2947DB", cpd = "blue"), bins=list(gene=20,cpd=20), plot.col.key = TRUE, same.layer=T, split.group=T, node.sum="max", out.suffix=paste(runName, "graphviz_",upordown,sep=""))
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

