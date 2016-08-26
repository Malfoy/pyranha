# pyranha
Clone the project
	
	git clone --recursive https://github.com/Malfoy/pyranha.git

Compile the project: compile the mapping tool

	cd pyranha
	./install.sh

Example commandline for the mapping tool:


/!\ Warning : all fasta files sequences must be written on one line.

* 1- Perform kmer counting with KMC:

	./scripts/kmerCount.py data/1Mrna.fa data/smallref.fa 15 data/smallref15mer

* 2- Map reads on reference:

	./bin/kIWImap -u data/1Mrna.fa -x data/smallref.fa -b data/smallref15mer -k 15 -t 8
 
Example Output

	Filling index of data/smallref.fa
	
	Indexing took : 0 sec
	
	Mapping data/1Mrna.fa
	
	Reads: 1000000
	
	Reads mapped: 403638
	
	Percent Read mapped: 40.3638%
	
	Reads not mapped: 596362
	
	Mapping took : 7 sec
	
	Throughout: 142k read by second or 514M by hour



Example command line for the de novo discovery tool:
	./bin/kIWI ./tools/kIWI_denovo/test/test.fa result.out

Example Output
	ping-pong: AACGTGCCAT
	+: 0, 1, 
	-: 2, 3,

The first line means ping pong has been found in the data with a 10-mer overlap "AACGTGCCAT", the reads participating are reads 0 and 1 of the dataset for the + strand, 2 and 3 for the - strand.
