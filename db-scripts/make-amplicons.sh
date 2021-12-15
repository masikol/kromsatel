#!/usr/bin/env bash

if [ $# -lt 1 ]; then
  echo "Please type $0 -h for help"
  exit 1
fi

VERSION='1.2.a'
HELP=0
primers=''
genome=''
outfile='./my-amplicons.fasta'
seqkit='seqkit'
min_len_minor=25

# Parse command line arguments
while getopts ":p:g:o:s:i:h" opt
do
  case "${opt}" in
    p) primers="${OPTARG}"
    ;;
    g) genome="${OPTARG}"
    ;;
    o) outfile="${OPTARG}"
    ;;
    s) seqkit="${OPTARG}"
    ;;
    i) min_len_minor="${OPTARG}"
    ;;
    h) HELP=1
    break
    ;;
    *) echo -e "\aError: invalid syntax. Here is the help message:" >&2
    HELP=1
    break
    ;;
  esac
done

# Print help message if necessary
if [[ ${HELP} == 1 ]]; then
  echo -e "make-amplicons.sh version ${VERSION}\n"
  echo -e "This script generates fasta file containing amplicons for kromsatel.py.\n"
  echo -e "Depencencies: seqkit (https://github.com/shenwei356/seqkit).\n"
  echo "Input:"
  echo "1. CSV File containing primers (see file \`primers/nCov-2019_primers.csv\` for example)"
  echo "   Important: order of primers is crucial, so keep order exactly as in this example file."
  echo "2. Genome (in fasta format) of damned coronavirus"
  echo "   (I use this one: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512)."
  echo "Output:"
  echo "1. Fasta file containing amplicons and ready for database creation,"
  echo -e "   just like \`amplicons-fasta/nCoV-2019_amplicons.fasta\`.\n"
  echo "Options:"
  echo " -h -- print help message and exit;"
  echo " -p -- CSV file with primers  (mandatory);"
  echo " -g -- fasta file with genome (mandatory);"
  echo " -o -- output fasta file;"
  echo "       default value: \`./my-amplicons.fasta\`"
  echo " -s -- path to \`seqkit\` executable"
  echo -e "       (if \`seqkit\` is in your PATH, just omit this option);"
  echo " -i -- minimum length of a minor amplicon to output;"
  echo -e "       default value: 25.\n"
  echo "Usage:"
  echo -e "  bash $0 <-p primers_csv> <-g genome_fasta> [-o STR] [-s STR] [-i INT]\n"
  echo "Examples:"
  echo " Case 1. \`seqkit\` is in PATH"
  echo "  bash $0 -p my_primers.csv -g Wuhan-Hu-1-compele-genome.fasta -o amplicons-fasta/my_amplicons.fasta"
  echo " Case 2. \`seqkit\` is not in PATH"
  echo "  bash $0 -p my_primers.csv -g Wuhan-Hu-1-compele-genome.fasta -o amplicons-fasta/my_amplicons.fasta -s /home/me/seqkit/bin/seqkit"

  exit 0
fi

# Check mandatory arguments
if [[ -z "${primers}" ]]; then
  echo -e "\aError: option \`-p\` is mandatory!"
  exit 1
fi
if [[ -z "${genome}" ]]; then
    echo -e "\aError: option \`-g\` is mandatory!"
    exit 1
fi

# Check existance of input files
for f in "${primers}" "${genome}"; do
  if [[ ! -f "${f}" ]]; then
    echo -e "\aError: file \`${f}\` does not exits!"
    exit 1
  fi
done

# Create output directory if it does not exist
if [[ ! -d `dirname "${outfile}"` ]]; then
  mkdir -p `dirname "${outfile}"`
fi

# Check if seqkit is available
"${seqkit}" > /dev/null 2> /dev/null
if [[ $? != 0 ]]; then
  echo -e "\aError: cannot find seqkit."
  echo "Possible solutions:"
  echo " 1. Install seqkit and add it to PATH variable."
  echo " 2. Pass correct path to seqkit executable with \`-s\` option."
  exit 1
fi

# Check if mininum length of minor amplicon is positive integer number
if [[ (! ${min_len_minor} =~ [0-9]+) || (${min_len_minor} -le 0) ]]; then
  echo -e "\aError: mininum length of minor amplicon (\`-i\` option)"
  echo '  must be a positive integer number.'
  echo "Your value: \`${min_len_minor}\`."
  exit 1
fi

echo -e "make-amplicons.sh Version ${VERSION}\n"

# Validate and parse primers
primer_seqs=()
echo -n "Valudating primers... "
while read -r line; do
  if [[ ! "${line}" =~ .+,[AGCTagct]+ ]]; then
    echo "Error: invalid line in file \`${primers}\`:"
    echo " \"${line}\""
    exit 2
  fi
  primer_seqs+=( ` echo "${line}" | cut -f2 -d','` )
done < "${primers}"
echo "ok"

n_primer_seqs="${#primer_seqs[@]}"

# Function generates title for a major amplicon by given `i` -- index (in `primer_seqs` array
#   of forward primer sequence forming this amplicon.
get_major_title() {
  i=$1
  let "major_idx = ($i / 2) + 1"
  echo "A${major_idx} ${major_idx}_major"
}

# Function generates title for a minor amplicon by given `i` -- index (in `primer_seqs` array
#   of forward primer sequence forming this amplicon.
get_minor_title() {
  i=$1
  let "minor_idx = ($i - 1) / 2 + 1"
  let "minor_idx_p1 = $minor_idx + 1"
  echo "I${minor_idx} ${minor_idx}-${minor_idx_p1}_minor"
}


# Function that calls seqkit amplicon (inner region), removes fasta header from it's output
#   and finally removed blank lines, leaving only sequence. And returns this sequence.
# $1 -- genome, #2 -- seqkit, $3 -- forward_primer, $4 -- reverse_primer
# This function removes primer sequences from amplicons.
extract_amplicon_seq() {

  # Configure coordinates for `-r` option (see `seqkit amplicon -h`)
  x=`echo ${3} | wc -c`
  y=`echo ${4} | wc -c`

  seq=`cat "$1" | "$2" amplicon --quiet -F "$3" -R "$4" -r ${x}:-${y} | grep -v ">" | tr -d '\n'`

  if [[ $? != 0 ]]; then
    echo "Error occured while extracting amplicon for following primers:"
    echo -en "FORWARD: $3\nREVERSE: $4\n"
    exit 1
  fi

  echo "${seq}"
}


# |== PROCEED ==|

echo -n '' > "${outfile}"

real_outfile=`realpath "${outfile}"`
echo -e "\nExtracting amplicons without primer sequences."
echo -e "Output file:\n  \`${real_outfile}\`"

echo "Extracting major amplicons..."
echo -en '>0'

for (( i=0; i<${n_primer_seqs}; i+=2 )); do
  # Get primers and title for current amplicon
  forward="${primer_seqs[$i]}"
  reverse="${primer_seqs[$i+1]}"
  title=`get_major_title "${i}"`

  # Extract amplicon and write it to outfile
  seq=`extract_amplicon_seq "${genome}" "${seqkit}" "${forward}" "${reverse}"`
  if [[ ! -z ${seq} ]]; then
    echo -en ">${title}\n${seq}\n" >> "${outfile}"
  fi

  let "ndone = ($i/2)+1"
  printf "\b\b\b=>%2d" ${ndone}
done

let "status_value = $n_primer_seqs / 2"
echo "/${status_value}"

echo "Extracting minor amplicons..."
echo -en '>'

for (( i=2; i<${n_primer_seqs}; i+=2 )); do
  # Get primers and title for current amplicon
  forward="${primer_seqs[$i]}"
  reverse="${primer_seqs[$i-1]}"
  title=`get_minor_title "${i}"`

  # Extract amplicon and write it to outfile
  seq=`extract_amplicon_seq "${genome}" "${seqkit}" "${forward}" "${reverse}"`
  if [[ (! -z ${seq}) && (${#seq} -ge ${min_len_minor}) ]]; then
    echo -en ">${title}\n${seq}\n" >> "${outfile}"
  fi

  let "ndone = ($i/2)"
  printf "\b\b\b=>%2d" ${ndone}
done

let "status_value = $n_primer_seqs / 2 - 1"
echo "/${status_value}"

echo -e "\nTask is comleted."
exit 0
