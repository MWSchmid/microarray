import subprocess
import sys

usage = """

python prepareMicroarrayProbes.py MODE ...

there are four modes:

BUILD:
builds a bowtie index and requires two additional arguments:
(i) the cDNA sequence fasta file
(ii) the bowtie index name/path
example:
python prepareMicroarrayProbes.py BUILD myCDNA.fasta myBowtieIndex

TABTOFASTA:
reformats a table of probe names and sequences into a fasta file.
Requires five additional arguments:
(i) file with the tabular probe name and sequence data
(ii) tells if there is a header in the table (1 for yes, 0 for no)
(iii) column number (starting with 1) holding the probe names
(iv) column number (starting with 1) holding the sequences
(v) fasta file where the sequences will be written

ALIGN:
aligns the probes to the cDNA and requires four additional arguments:
(i) the bowtie index name/path (see BUILD above)
(ii) a fasta file with the probe sequences
(iii) a file where the unaligned sequences will be stored
(iv) a file where the aligned sequences will be stored
example:
python prepareMicroarrayProbes.py BUILD myCDNA.fasta myBowtieIndex

EXTRACT:
extracts the mappings from probe names to locus IDs and vice versa.
Requires three additional arguments:
(i) file with the aligned sequences (see ALIGN above)
(ii) file with the probe name to locus ID mappings
(iii) file with the locus ID to probe name mappings

"""

if len(sys.argv) < 2:
	sys.exit(usage)

mode = sys.argv[1]

# build bowtie index
def buildIndex(fastaFile, indexPath):
	command = "bowtie-build -f -q -o 0 %s %s" % (fastaFile, indexPath)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(stdout_data, stderr_data) = process.communicate()
	if stderr_data is not None:
		print >> sys.stderr, "errors raised by bowtie-build: " + stderr_data
	if stdout_data is not None:
		print >> sys.stdout, "bowtie-build message: " + stdout_data
	print >> sys.stderr, "build bowtie index"

# convert the sequence table with probeName, sequence, spot_id to a fasta file
def convertSeqTabToFasta(seqTab, hasHeader, probeColumn, sequenceColumn, fastaFile):
	with open(seqTab, 'r') as infile, open(fastaFile, 'w') as outfile:
		if hasHeader:
			header = infile.readline()
			print >> sys.stderr, "removed this header:", header[:-1]
		for line in infile:
			fields = line[:-1].split('\t')
			pn = fields[probeColumn]
			sequence = fields[sequenceColumn]
			print >> outfile, '>'+pn
			print >> outfile, sequence
	print >> sys.stderr, "converted sequence table to fasta"

# align sequences
def alignSequences(bowtieIndex, sequenceFile, unalignedFile, alignedFile):
	command = "bowtie -v 3 -m 10 -a --best --strata -p 4 -t %s -f %s --un %s %s" % (bowtieIndex, sequenceFile, unalignedFile, alignedFile)
	# ask the user if he wants to run the command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(stdout_data, stderr_data) = process.communicate()
	if stderr_data is not None:
		print >> sys.stderr, "errors raised by bowtie: " + stderr_data
	if stdout_data is not None:
		print >> sys.stdout, "bowtie message: " + stdout_data
	print >> sys.stderr, "aligned sequences"

# get mapping tables
def getMapping(alignedFile, PNtoLOCUSfile, LOCUStoPNfile):
	LOCUStoPN = {}
	PNtoLOCUS = {}
	with open(alignedFile, 'r') as infile:
		for line in infile:
			pn, strand, modelID, pos, sequence, quality, mm, conversions = line[:-1].split('\t')
			geneID = modelID.split('.')[0]
			try:
				if pn not in LOCUStoPN[geneID]:
					LOCUStoPN[geneID].append(pn)
			except KeyError:
				LOCUStoPN[geneID] = [pn]
			try:
				if geneID not in PNtoLOCUS[pn]:
					PNtoLOCUS[pn].append(geneID)
			except KeyError:
				PNtoLOCUS[pn] = [geneID]
	print >> sys.stderr, "extracted mappings"
	with open(PNtoLOCUSfile, 'w') as outfile:
		for pn, locusList in PNtoLOCUS.items():
			print >> outfile, pn + '\t' + ';'.join(locusList)
	with open(LOCUStoPNfile, 'w') as outfile:
		for locus, pnList in LOCUStoPN.items():
			print >> outfile, locus + '\t' + ';'.join(pnList)

# 
mode = sys.argv[1]

if mode == "BUILD":
	fastaFile = sys.argv[2]
	indexPath = sys.argv[3]
	buildIndex(fastaFile, indexPath)
elif mode == "TABTOFASTA":
	seqTab = sys.argv[2]
	hasHeader = bool(int(sys.argv[3]))
	probeColumn = int(sys.argv[4])-1
	sequenceColumn = int(sys.argv[5])-1
	fastaFile = sys.argv[6]
	convertSeqTabToFasta(seqTab, hasHeader, probeColumn, sequenceColumn, fastaFile)
elif mode == "ALIGN":
	bowtieIndex = sys.argv[2]
	sequenceFile = sys.argv[3]
	unalignedFile = sys.argv[4]
	alignedFile = sys.argv[5]
	alignSequences(bowtieIndex, sequenceFile, unalignedFile, alignedFile)
elif mode == "EXTRACT":
	alignedFile = sys.argv[2]
	PNtoLOCUSfile = sys.argv[3]
	LOCUStoPNfile = sys.argv[4]
	getMapping(alignedFile, PNtoLOCUSfile, LOCUStoPNfile)
else:
	print >> sys.stderr, "unkown mode"
	sys.exit(usage)

