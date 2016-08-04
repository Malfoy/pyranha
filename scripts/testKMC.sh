#~ rm -rf lol*;
rm ref;
rm allroset;
kmc -k18 -fm -ci1  dmoj-all-chromosome-r1.3.fasta  ref lol1;
kmc -k18 -ci1  smallRNA_hybridA_cutadapt_sup23.fastq allrna lol2;
kmc_tools intersect ref -ci1 allrna -ci1 inter -ci1;
kmc_dump inter kmercount18;                      
