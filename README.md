# pyranha
Clone the project
	git clone --recursive https://github.com/Malfoy/pyranha.git

Compile the project: compile the mapping tool

	cd pyranha
	./install.sh

Example commandline

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




