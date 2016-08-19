# pyranha
kIWI

Compile the project: compile the mapping tool
	cd pyranha
	cd tools/src/piRANHA_mapping
	make
	cd ../..

Example commandline

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


Example commandline

	./bin/kIWImap -u data/100krna.fa -x data/bigref.fa -b data/bigref15mer -k 15 -t 8 

Example Output

	Filling index of data/bigref.fa
	
	Indexing took : 12 sec
	
	Mapping data/100krna.fa
	
	Reads: 100000
	
	Reads mapped: 90736
	
	Percent Read mapped: 90.736%
	
	Reads not mapped: 9264
	
	Mapping took : 19 sec
	
	Throughout: 5k read by second or 18M by hour



