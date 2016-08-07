#!/bin/bash


OPTIONS=$(getopt -o q:f:h -n 'parse-options' -- "$@")
if [ $? -ne 0 ]; then
  echo "Option error. Usage: fastq2fasta.sh -q input.fastq -f output.fasta"
  exit 1
fi

echo "$OPTIONS"
eval set -- $OPTIONS

FASTQ=""
FASTA=""
HELP=0

while true; do
  case "$1" in

    -f) FASTA="$2" ; shift ;;
    -q) FASTQ="$2" ; shift ;;
    -h) HELP=1 ;;
    --)        shift ; break ;;
    *)         echo "unknown option: $1. Usage: fastq2fasta.sh -q input.fastq -f output.fasta" ; exit 1 ;;


  esac
  shift
done



if [ $# -ne 0 ]; then
  echo "unknown option(s): $@. Usage: fastq2fasta.sh -q input.fastq -f output.fasta"
  exit 1
fi

if [ $HELP -ne 0 ]; then
  echo "usage: fastq2fasta.sh -q input.fastq -f output.fasta"
  exit 1
fi

echo "from file: $FASTQ"
echo "fasta written in: $FASTA"


awk '{if ( NR%4==1 || NR%4==2){print $0}}' $FASTQ  | sed 's/@/>/g' > $FASTA
