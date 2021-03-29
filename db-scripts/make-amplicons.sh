#!/usr/bin/env bash

if [ $# -lt 1 ]; then
  echo "Please type $0 -h for help"
  exit 1
fi

VERSION='1.1.a'
HELP=0
primers=''
genome=''
outfile_prefix='./my-amplicons'
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
    o) outfile_prefix="${OPTARG}"
    ;;
    s) seqkit="${OPTARG}"
    ;;
    i) min_len_minor="${OPTARG}"
    ;;
    h) HELP=1
    break
    ;;
    *) echo -e "\aError: invalid syntax. Here is help message:" >&2
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
  echo " -o -- prefix of output fasta files;"
  echo "       default value: \`./my-amplicons\`"
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
if [[ ! -d `dirname "${outfile_prefix}"` ]]; then
  mkdir -p `dirname "${outfile_prefix}"`
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

# Define fasta titles for output file
major_titles=( 'A1 1_major' 'A2 2_major' 'A3 3_major' 'A4 4_major'
               'A5 5_major' 'A6 6_major' 'A7 7_major' 'A8 8_major'
               'A9 9_major' 'A10 10_major' 'A11 11_major' 'A12 12_major'
               'A13 13_major' 'A14 14_major' 'A15 15_major' 'A16 16_major'
               'A17 17_major' 'A18 18_major' 'A19 19_major' 'A20 20_major'
               'A21 21_major' 'A22 22_major' 'A23 23_major' 'A24 24_major'
               'A25 25_major' 'A26 26_major' 'A27 27_major' 'A28 28_major'
               'A29 29_major' 'A30 30_major' 'A31 31_major' 'A32 32_major'
               'A33 33_major' 'A34 34_major' 'A35 35_major' 'A36 36_major'
               'A37 37_major' 'A38 38_major' 'A39 39_major' 'A40 40_major'
               'A41 41_major' 'A42 42_major' 'A43 43_major' 'A44 44_major'
               'A45 45_major' 'A46 46_major' 'A47 47_major' 'A48 48_major'
               'A49 49_major' 'A50 50_major' 'A51 51_major' 'A52 52_major'
               'A53 53_major' 'A54 54_major' 'A55 55_major' 'A56 56_major'
               'A57 57_major' 'A58 58_major' 'A59 59_major' 'A60 60_major'
               'A61 61_major' 'A62 62_major' 'A63 63_major' 'A64 64_major'
               'A65 65_major' 'A66 66_major' 'A67 67_major' 'A68 68_major'
               'A69 69_major' 'A70 70_major' 'A71 71_major' 'A72 72_major'
               'A73 73_major' 'A74 74_major' 'A75 75_major' 'A76 76_major'
               'A77 77_major' 'A78 78_major' 'A79 79_major' 'A80 80_major'
               'A81 81_major' 'A82 82_major' 'A83 83_major' 'A84 84_major'
               'A85 85_major' 'A86 86_major' 'A87 87_major' 'A88 88_major'
               'A89 89_major' 'A90 90_major' 'A91 91_major' 'A92 92_major'
               'A93 93_major' 'A94 94_major' 'A95 95_major' 'A96 96_major'
               'A97 97_major' 'A98 98_major' )

minor_titles=( 'I1 1-2_minor' 'I2 2-3_minor' 'I3 3-4_minor' 'I4 4-5_minor'
               'I5 5-6_minor' 'I6 6-7_minor' 'I7 7-8_minor' 'I8 8-9_minor'
               'I9 9-10_minor' 'I10 10-11_minor' 'I11 11-12_minor' 'I12 12-13_minor'
               'I13 13-14_minor' 'I14 14-15_minor' 'I15 15-16_minor' 'I16 16-17_minor'
               'I17 17-18_minor' 'I18 18-19_minor' 'I19 19-20_minor' 'I20 20-21_minor'
               'I21 21-22_minor' 'I22 22-23_minor' 'I23 23-24_minor' 'I24 24-25_minor'
               'I25 25-26_minor' 'I26 26-27_minor' 'I27 27-28_minor' 'I28 28-29_minor'
               'I29 29-30_minor' 'I30 30-31_minor' 'I31 31-32_minor' 'I32 32-33_minor'
               'I33 33-34_minor' 'I34 34-35_minor' 'I35 35-36_minor' 'I36 36-37_minor'
               'I37 37-38_minor' 'I38 38-39_minor' 'I39 39-40_minor' 'I40 40-41_minor'
               'I41 41-42_minor' 'I42 42-43_minor' 'I43 43-44_minor' 'I44 44-45_minor'
               'I45 45-46_minor' 'I46 46-47_minor' 'I47 47-48_minor' 'I48 48-49_minor'
               'I49 49-50_minor' 'I50 50-51_minor' 'I51 51-52_minor' 'I52 52-53_minor'
               'I53 53-54_minor' 'I54 54-55_minor' 'I55 55-56_minor' 'I56 56-57_minor'
               'I57 57-58_minor' 'I58 58-59_minor' 'I59 59-60_minor' 'I60 60-61_minor'
               'I61 61-62_minor' 'I62 62-63_minor' 'I63 63-64_minor' 'I64 64-65_minor'
               'I65 65-66_minor' 'I66 66-67_minor' 'I67 67-68_minor' 'I68 68-69_minor'
               'I69 69-70_minor' 'I70 70-71_minor' 'I71 71-72_minor' 'I72 72-73_minor'
               'I73 73-74_minor' 'I74 74-75_minor' 'I75 75-76_minor' 'I76 76-77_minor'
               'I77 77-78_minor' 'I78 78-79_minor' 'I79 79-80_minor' 'I80 80-81_minor'
               'I81 81-82_minor' 'I82 82-83_minor' 'I83 83-84_minor' 'I84 84-85_minor'
               'I85 85-86_minor' 'I86 86-87_minor' 'I87 87-88_minor' 'I88 88-89_minor'
               'I89 89-90_minor' 'I90 90-91_minor' 'I91 91-92_minor' 'I92 92-93_minor'
               'I93 93-94_minor' 'I94 94-95_minor' 'I95 95-96_minor' 'I96 96-97_minor'
               'I97 97-98_minor' )

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

# Function that calls seqkit amplicon, removes fasta header from it's output
#   and finally removed blank lines, leaving only sequence. And returns this sequence.
# $1 -- genome, #2 -- seqkit, $3 -- forward_primer, $4 -- reverse_primer.
# This function does not remove primer sequences from amplicons.
extract_amplicon_seq_with_primers() {
  seq=`cat "$1" | "$2" amplicon --quiet -F "$3" -R "$4" | grep -v ">" | tr -d '\n'`

  if [[ $? != 0 ]]; then
    echo "Error occured while extracting amplicon for following primers:"
    echo -en "FORWARD: $3\nREVERSE: $4\n"
    exit 1
  fi

  echo "${seq}"
}


# |== PROCEED ==|

outfiles=("${outfile_prefix}.fasta" "${outfile_prefix}_with-primers.fasta")
extract_functions=(extract_amplicon_seq extract_amplicon_seq_with_primers)
strings_to_print=("without" "keeping")

for (( outfile_i=0; outfile_i<${#outfiles[@]}; outfile_i++ )); do

  extract=${extract_functions[$outfile_i]}

  # Init outfile
  outfile=${outfiles[$outfile_i]}
  real_outfile=`realpath "${outfile}"`
  echo -n '' > "${outfile}"

  echo -e "\nExtracting amplicons ${strings_to_print[$outfile_i]} primer sequences."
  echo -e "Output file:\n  \`${real_outfile}\`"

  echo "Extracting major amplicons..."
  echo -en '>0'

  for (( i=0; i<${#primer_seqs[@]}; i+=2 )); do
    # Get primers and title for current amplicon
    forward="${primer_seqs[$i]}"
    reverse="${primer_seqs[$i+1]}"
    title="${major_titles[$i/2]}"

    # Extract amplicon and write it to outfile
    seq=`$extract "${genome}" "${seqkit}" "${forward}" "${reverse}"`
    if [[ ! -z ${seq} ]]; then
      echo -en ">${title}\n${seq}\n" >> "${outfile}"
    fi

    let "ndone = ($i/2)+1"
    printf "\b\b\b=>%2d" ${ndone}
  done

  echo "/${#major_titles[@]}"

  echo "Extracting minor amplicons..."
  echo -en '>'

  for (( i=2; i<${#primer_seqs[@]}; i+=2 )); do
    # Get primers and title for current amplicon
    forward="${primer_seqs[$i]}"
    reverse="${primer_seqs[$i-1]}"
    title="${minor_titles[($i-1)/2]}"

    # Extract amplicon and write it to outfile
    seq=`$extract "${genome}" "${seqkit}" "${forward}" "${reverse}"`
    if [[ (! -z ${seq}) && (${#seq} -ge ${min_len_minor}) ]]; then
      echo -en ">${title}\n${seq}\n" >> "${outfile}"
    fi

    let "ndone = ($i/2)"
    printf "\b\b\b=>%2d" ${ndone}
  done

  echo "/${#minor_titles[@]}"
done

echo -e "\nTask is comleted."
exit 0
