#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil


if len(sys.argv) < 5:
	sys.exit( 'Usage : [Read file] [Reference file] [Kmer size] [Output  file]    NB: fasta file please or modify the script :)' )

scriptDir = os.path.dirname(os.path.realpath(__file__))
binDir = scriptDir + "/../bin/"

debut=time.time()

readFile=sys.argv[1]
refFile=sys.argv[2]
kmerSize=sys.argv[3]
outFile=sys.argv[4]

subprocess.call(["mkdir", "readFolderKMC", "referenceFolderKMC"])

print "Kmercount of the reference file"
subprocess.call([binDir + "kmc", "-fm", "-k"+kmerSize, "-ci1", refFile, "reference", "referenceFolderKMC"])

print "Kmercount of the read file"
subprocess.call([binDir + "kmc", "-fa", "-k"+kmerSize, "-ci1", readFile, "read", "readFolderKMC"])

print "Merging"
subprocess.call([binDir + "kmc_tools", "intersect",   "reference", "read","-ci1", "intersection"])

print "Dumping"
subprocess.call([binDir + "kmc_dump", "intersection", outFile])

print "CLeaning"
subprocess.call(["rm", 'reference.kmc_pre', 'reference.kmc_suf', 'read.kmc_suf', 'read.kmc_pre', 'intersection.kmc_pre', 'intersection.kmc_suf'])

fin=time.time();
print fin-debut
print "secondes"
