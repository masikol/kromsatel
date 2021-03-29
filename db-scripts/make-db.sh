#!/usr/bin/env bash

if [ $# -lt 1 ]; then
  echo "Please type $0 -h for help"
  exit 1
fi

VERSION='1.0.a'

if [[ $1 == '-h' ]]; then
  echo "make-db.sh version ${VERSION}"
  echo -e "This script creates a BLAST database of amplicons for kromsatel.py.\n"
  echo -e "Depencencies: makeblastdb from BLAST+ toolkit (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).\n"
  echo "Input: fasta file of amplicons,"
  echo "  just like \`amplicons-fasta/nCoV-2019_amplicons.fasta\`.\n"
  echo -n "Output: BLAST-formatted database of amplicons.\n"
  echo "Usage:"
  echo "  bash $0 <input_fasta> <output_dir>"
  echo "Example:"
  echo "  bash $0 amplicons-fasta/nCoV-2019_amplicons.fasta nCoV-2019-database"
  exit 0
fi

infile=$1
outdir=$2

# Check if all args are passed
if [[ (-z "${infile}") || (-z "${outdir}") ]]; then
  echo -e "\aError: not enough arguments."
  echo "Type \`$0 -h\` to see help message."
  exit 1
fi

# Check input file
if [[ ! -f "${infile}" ]]; then
  echo -e "\aError: file \`${infile}\` does not exist!"
  exit 1
fi

# Check if makeblatdb is available
tool='makeblastdb'
if [[ -z `which "${tool}"` ]]; then
  echo -e "\aError: cannot find program \`${tool}\` in your PATH."
  echo 'Please, install it or add to PATH if it is installed.'
  exit 1
fi

# Make outdir
outdir=`realpath "${outdir}"`
if [[ ! -d "${outdir}" ]]; then
  mkdir -p "${outdir}"
  if [[ $? != 0 ]]; then
    echo -e "\aError: cannot create output directory \`${outdir}\`."
    exit 1
  fi
fi

# Configure title of the database: remove `.fasta` extention
title=`basename "${infile}" | sed -E "s|\.fa(sta)?||"`
outprefix="${outdir}/${title}"

echo -e "make-db.sh Version ${VERSION}\n"

echo "Creating BLAST database."
echo "Input file: \`${infile}\`."
echo "Output directory: \`${outdir}\`."

# Proceed
"${tool}" -in "${infile}" -dbtype nucl -parse_seqids -out "${outprefix}" -title "${title}"

if [[ $? != 0 ]]; then
  echo "\aError: cannot create database!"
  exit 1
fi

echo "Path to created database: \`${outprefix}\`"
echo 'Task is completed.'
exit 0
