library("limma")

# directory where all the data is located
workDir <- "/home/marc/microarray/testData/physco"

# read an annotation table - note that following columns must be defined in there
# Name: a unique sample name
# Group: a name for the sample group
# File: the name of the microarray data file
# Channel: the channel for this sample (Cy3 or Cy5)
annotation <- read.table(file.path(workDir, "annotation.txt"), header = TRUE, sep = '\t', quote = "", stringsAsFactors = FALSE)
annotation$weirdName <- paste(file.path(workDir, annotation$File), annotation$Channel, sep = '_')
if (length(unique(annotation$Name)) != nrow(annotation)) {stop("The names in the annotation file must be unique!")}

# if you like to remove the "swap experiments" (which are in a way technical replicates and maybe not of interest - it is good to check though if sampleA is similar on both channels)
# annotation <- annotation[-grep("swap", annotation$Name),]

# read raw data using read.maimages (from limma) and create a raw data list which mimics a one-channel raw data list
# note that the annotation = c("ProbeName") is important. Make sure that the column ProbeName exists in the microarray data file
twoChannelRawData <- read.maimages(unique(file.path(workDir, annotation$File)), source = "agilent", annotation = c("ProbeName"))
colnames(twoChannelRawData$G) <- paste(colnames(twoChannelRawData$G), ".txt_Cy3", sep = '')
colnames(twoChannelRawData$Gb) <- paste(colnames(twoChannelRawData$Gb), ".txt_Cy3", sep = '')
colnames(twoChannelRawData$R) <- paste(colnames(twoChannelRawData$R), ".txt_Cy5", sep = '')
colnames(twoChannelRawData$Rb) <- paste(colnames(twoChannelRawData$Rb), ".txt_Cy5", sep = '')
greenTargets <- twoChannelRawData$targets; rownames(greenTargets) <- paste(rownames(twoChannelRawData$targets), ".txt_Cy3", sep = '')
redTargets <- twoChannelRawData$targets; rownames(redTargets) <- paste(rownames(twoChannelRawData$targets), ".txt_Cy5", sep = '')
rawDataE <- matrix(0, nrow = nrow(twoChannelRawData$G), ncol = nrow(annotation), dimnames = list(rownames(twoChannelRawData$G), annotation$Name))
rawDataEb  <- matrix(0, nrow = nrow(twoChannelRawData$Gb), ncol = nrow(annotation), dimnames = list(rownames(twoChannelRawData$G), annotation$Name))
rawDataE[,annotation$Name[annotation$Channel == "Cy3"]] <- twoChannelRawData$G[,annotation$weirdName[annotation$Channel == "Cy3"]]
rawDataE[,annotation$Name[annotation$Channel == "Cy5"]] <- twoChannelRawData$R[,annotation$weirdName[annotation$Channel == "Cy5"]]
rawDataEb[,annotation$Name[annotation$Channel == "Cy3"]] <- twoChannelRawData$Gb[,annotation$weirdName[annotation$Channel == "Cy3"]]
rawDataEb[,annotation$Name[annotation$Channel == "Cy5"]] <- twoChannelRawData$Rb[,annotation$weirdName[annotation$Channel == "Cy5"]]
rawData <- list(E = rawDataE, Eb = rawDataEb, targets = rbind(redTargets, greenTargets), source = "agilent", genes = twoChannelRawData$genes)
class(rawData) <- "EListRaw"

# correct, normalize, and extract
bgData <- backgroundCorrect(rawData) # method = "normexp" caused some problems for someone
normData <- normalizeBetweenArrays(bgData, method = "quantile")
normEset <- normData$E; rownames(normEset) <- normData$genes$ProbeName

# load the locus ID to probe name mappings
LOCUStoPNtable <- read.table(file.path(workDir, "agilent_COSMOSStoProbeName.txt"), header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE, row.names = 1); colnames(LOCUStoPNtable) <- c("probeNames")
LOCUStoPN <- sapply(LOCUStoPNtable$probeNames, function(x) strsplit(x, ";", fixed = TRUE)); names(LOCUStoPN) <- rownames(LOCUStoPNtable)

# summarize the probes per locus ID (if there are several, take the mean)
locusData <- list()
for (sample in colnames(normEset)) {
	cat(paste(sample, '\n', sep = ''))
	sub <- normEset[,sample]
	locusData[[sample]] <- sapply(LOCUStoPN, function(x) ifelse(length(x)>1, mean(sub[x]), sub[x])) # switch mean with median if you like
}
locusData <- do.call("cbind", locusData)

# write out
write.table(locusData, file.path(workDir, "COSMOSSphyscoAgilentData_meanWithin.txt"), sep = '\t', quote = FALSE)



