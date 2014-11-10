microarray
==========

some sripts for microarray processing

soon available

## some preparation steps (align the probes to cDNA sequences which you are interested in)
The pre-processor script relies on a working [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) installation. The binaries are assumed to be located in your PATH.
To see all options of the python script type:

```shell
python prepareMicroarrayProbes.py
```

### build the bowtie index
To align the probes to the cDNAs of interest, one needs to build a bowtie index first. Note that the locus ID (so the ID that gets the expression value) corresponds to the first field after the arrow (>) in the fasta file (split using space character).

```shell
python prepareMicroarrayProbes.py BUILD cDNA.fasta cDNA_index
```

### reformat the sequence table
Probe sequences are frequently stored in tabular form. This needs to be changed into a fasta file.

```shell
python prepareMicroarrayProbes.py TABTOFASTA seqTable.txt 1 1 2 probes.fasta
```

### align sequences
To align the probes to the cDNAs of interest:

```shell
python prepareMicroarrayProbes.py ALIGN cDNA_index probes.fasta unaligned.txt aligned.txt
```

### extract the mappings
Extract the probe name to locus ID mappings (and vice versa):

```shell
python prepareMicroarrayProbes.py EXTRACT aligned.txt probeNameToID.txt IDtoProbeName.txt
```

## processing in R, 44K agilent single-channel arrays


## processing in R, 44K agilent dual-channel arrays



