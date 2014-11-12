library("limma")

# directory where all the data is located
workDir <- "/home/marc/microarray/testData/rice"

# read an annotation table - note that following columns must be defined in there
# Name: a unique sample name
# Group: a name for the sample group
# File: the name of the microarray data file
annotation <- read.table(file.path(workDir, "annotation.txt"), header = TRUE, sep = '\t', quote = "", stringsAsFactors = FALSE)
if (length(unique(annotation$Name)) != nrow(annotation)) {stop("The names in the annotation file must be unique!")}

# read raw data using read.maimages (from limma)
# note that the annotation = c("ProbeName") is important. Make sure that the column ProbeName exists in the microarray data file
rawData <- read.maimages(unique(file.path(workDir, annotation$File)), source = "agilent", green.only = TRUE, names = annotation$Name, annotation = c("ProbeName"))

# correct, normalize, and extract
bgData <- backgroundCorrect(rawData, method = "normexp")
normData <- normalizeBetweenArrays(bgData, method = "quantile")
normEset <- normData$E; rownames(normEset) <- normData$genes$ProbeName

# load the locus ID to probe name mappings
LOCUStoPNtable <- read.table(file.path(workDir, "agilent_MSUtoProbeName.txt"), header = FALSE, sep = '\t', quote = "", stringsAsFactors = FALSE, row.names = 1); colnames(LOCUStoPNtable) <- c("probeNames")
LOCUStoPN <- sapply(LOCUStoPNtable$probeNames, function(x) strsplit(x, ";", fixed = TRUE)); names(LOCUStoPN) <- rownames(LOCUStoPNtable)

# summarize the probes per locus ID (if there are several, take the mean)
locusData <- list()
for (sample in colnames(normEset)) {
	cat(paste(sample, '\n', sep = ''))
	sub <- normEset[,sample]
	locusData[[sample]] <- sapply(LOCUStoPN, function(x) ifelse(length(x)>1, mean(sub[x]), sub[x])) # switch mean with median if you like
}
locusData <- do.call("cbind", locusData)

# write out (in case of rice, you may trim the names first)
rownames(locusData) <- substr(rownames(locusData), 5, 14)
write.table(locusData, file.path(workDir, "MSUriceAgilentData_meanWithin.txt"), sep = '\t', quote = FALSE)

