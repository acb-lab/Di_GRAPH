#!/bin/bash

AUTHORS="Lydia Iglesias, Narciso M. Quijada, Andrés Clemente-Blanco"
LASTMODIF="2026-04-24"
VERSION="0.4"
MODIFICATIONS="experiment name dynamically set, COMPLETE, validation 10kb"


DiGRAPH=$0
while [ -h "$DiGRAPH" ]; do # resolve $DiGRAPH until the file is no longer a symlink
  DiGRAPHDIR="$( cd -P "$( dirname "$0" )" && pwd )"
  DiGRAPH="$(readlink "$DiGRAPH")"
  [[ $DiGRAPH != /* ]] && DiGRAPH="$DiGRAPHDIR/$DiGRAPH" # if $DiGRAPH is a symlink, resolve it relative to the path where the symlink file was located
done
DiGRAPHDIR="$( cd -P "$( dirname "$DiGRAPH" )" && pwd )"

# WORKING VARIABLES
CATEGORY_PATH=
MYBLAST=
MYREF=
MYREPORT=
MYWD=
MYSCRIPTS=
STRAIN_REF="no"
NUMBER_EXP=3
N_CPU=$(nproc)
#N_CPU=$(sysctl -n hw.ncpu) #macOS

########## I/O helpers (Di-GRAPH_piped.sh) ##########
# Streaming replacements for "gunzip -k + cat + chains of cutadapt --cut" that
# produce byte-identical FASTQ without writing the intermediate files.

# Decompress to stdout (pigz is multi-threaded; gzip -dc is the fallback)
if command -v pigz > /dev/null 2>&1; then DECOMP=(pigz -dc); else DECOMP=(gzip -dc); fi

# fq_slice SPEC... < FASTQ
# Reads FASTQ on stdin once and writes one output per SPEC, where SPEC is
#   PATH[:TAG]=FILE
# PATH is the sequence of cutadapt --cut steps applied in the original script,
# one digit per step: 1 = "--cut -N" (drop N from the end), 2 = "--cut N" (drop
# N from the start), with N = 75, 38, 19 for steps 1, 2, 3. E.g. "1" = first
# 75 nt, "212" = cut 75 -> cut -38 -> cut 19. TAG renames the read like the
# original sed: TAG is the new suffix plus the mate digit, e.g. "rA1" turns the
# first "/1" in the header into "rA/1" (sed '1~4 s/\/1/rA\/1/').
# FILE may be /dev/stdout. Every FILE is created/truncated even if empty.
fq_slice() {
  awk -v specs="$*" '
    function cn(s, n) { return substr(s, 1, length(s) - n) }   # cutadapt --cut -n
    function cp(s, n) { return substr(s, n + 1) }               # cutadapt --cut n
    function sl(s, p,    i, c, n) {
      for (i = 1; i <= length(p); i++) {
        c = substr(p, i, 1); n = (i == 1 ? 75 : (i == 2 ? 38 : 19))
        s = (c == "1") ? cn(s, n) : cp(s, n)
      }
      return s
    }
    BEGIN {
      n = split(specs, S, " ")
      for (i = 1; i <= n; i++) {
        split(S[i], kv, "="); F[i] = kv[2]
        T[i] = ""; P[i] = kv[1]
        if (index(kv[1], ":")) { split(kv[1], pt, ":"); P[i] = pt[1]; T[i] = pt[2] }
        printf "" > F[i]
      }
    }
    { r = NR % 4 }
    r == 1 { h = $0; next }
    r == 2 { s = $0; next }
    r == 3 { pl = $0; next }
    {
      for (i = 1; i <= n; i++) {
        hh = h
        if (T[i] != "") { m = substr(T[i], length(T[i])); sub("/" m, substr(T[i], 1, length(T[i]) - 1) "/" m, hh) }
        print hh "\n" sl(s, P[i]) "\n" pl "\n" sl($0, P[i]) > F[i]
      }
    }'
}

# Truncate the coverage column (4th) to 5 decimals - same awk as the original script
TRUNC4_AWK='BEGIN{OFS=FS} {split($4,a,"."); if(length(a)>1) $4=a[1]"."substr(a[2],1,5); print}'

# reap_fifo_producer FIFO PRODUCER_PID CONSUMER_STATUS: after the consumer of FIFO has
# finished, reap the producer. If the consumer failed before opening the FIFO,
# the producer may block forever in open(), so kill it.
reap_fifo_producer() {
  local fifo=$1 pid=$2 status=$3
  if [[ $status -ne 0 ]]; then kill "$pid" 2> /dev/null; fi
  wait "$pid" 2> /dev/null
  rm -f "$fifo"
}
# Emit each joined pair line as its two halves with the exact semantics of
# "cut -f1-19" and "cut -f20-38" (a line without a tab is printed whole by cut).
CUT_HALVES_AWK='index($0, "\t") == 0 { print; print; next }
{ r1 = $1; for (i = 2; i <= 19 && i <= NF; i++) r1 = r1 OFS $i
  r2 = ""; for (i = 20; i <= 38 && i <= NF; i++) r2 = (i == 20 ? $i : r2 OFS $i)
  print r1; print r2 }'
#####################################################

current_time=
log_file=
stats=

#MYREF=/home/ibfg/GWSL/RG # Replace with the path to your MYREF directory: path/to/RG
#MYWD=/home/ibfg/GWSL/WD/ # Replace with the path to your MYWD directory: path/to/WD
#CATEGORY_PATH=/home/ibfg/GWSL/Categories # Replace with the path to your CATEGORY_PATH directory: path/to/Categories
#MYBLAST=/home/ibfg/GWSL/BLAST/features_extraction # Replace with the path to your MYBLAST directory: path/to/BLAST
#MYREPORT=/home/ibfg/GWSL/Report_files # Replace with the path to your MYREPORT directory: path/to/Report_files

# Message colors
COL_RESET=$(tput sgr 0)
COL_blue=$(tput setaf 4)
COL_cyan=$(tput setaf 6)
COL_green=$(tput setaf 2)
COL_magenta=$(tput setaf 5)
COL_purple=$(tput setaf 5)
COL_red=$(tput setaf 1)
COL_white=$(tput setaf 7)
COL_yellow=$(tput setaf 3)

usage () {

cat << EOF

=============================================================

██████╗               ██████╗ ██████╗  █████╗ ██████╗ ██╗  ██╗
██╔══██╗  ██╗        ██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██║  ██║
██║  ██║  ╚═╝  ███╗  ██║  ███╗██████╔╝███████║██████╔╝███████║  
██║  ██║  ██╗  ╚══╝  ██║   ██║██╔╚██╗ ██╔══██║██╔═══╝ ██╔══██║  
██████╔╝  ██║        ╚██████╔╝██║ ╚██╗██║  ██║██║     ██║  ██║
╚═════╝   ╚═╝         ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝

      DSB-induced Genome-wide Repair Analysis 
     and Profiling of Homologous recombination
                   ${VERSION}
 ${AUTHORS}
=============================================================

Description:
This tool performs genome-wide analysis of DNA repair following an HO-induced break. 
Di-GRAPH integrates:
  - MAT repair profiling
  - Genome instability evaluation based on coverage 
  - Recombination analysis based on inter-chromosomal discordant reads
Input files are gzipped paired-end fastq files.

OBLIGATORY OPTIONS:
    -b/--blast          Path to BLAST genomic database for cross-validation directory (provided in GitHub)
    -c/--categories     Path to genomic categories and features annotation directory (provided in GitHub)
    -g/--genome         Path to reference genome directory (provided in GitHub)
    -i/--input          Path to input directory (working directory generated by the user)
    -r/--report         Path to report directory (provided in GitHub)
    -S/--scripts        Path to scripts directory (provided in GitHub)


OTHER OPTIONS:
    -h/--help           Shows this help
    -s/--strain_ref     Name of the strain to be used as a reference (example: 1_Wt)
    -n/--number_exp     Number of experiments (default=3)
    -t/--threads        Number of threads to use (default=${N_CPU})
    -v/--version        Show version

EOF
}

if [ $# == 0 ]; then
	usage
	exit 1
fi

# VARIABLE OPTIONS

POSITIONAL=()
while [[ $# -gt 0 ]]
do
ARGS="$1"

case $ARGS in
    -b|--blast)
    if [ "$2" ]; then
        MYBLAST=$2
        shift 2
        if [ ! -d "${MYBLAST}" ]; then
            echo -e "\n${COL_red}ERROR: ${MYBLAST} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${MYBLAST})" ]; then
            echo -e "\n${COL_red}ERROR: ${MYBLAST} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-b/--blast" requires an argument\n'
        exit 1
    fi
    ;;
    -c|--categories)
    if [ "$2" ]; then
        CATEGORY_PATH=$2
        shift 2
        if [ ! -d "${CATEGORY_PATH}" ]; then
            echo -e "\n${COL_red}ERROR: ${CATEGORY_PATH} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${CATEGORY_PATH})" ]; then
            echo -e "\n${COL_red}ERROR: ${CATEGORY_PATH} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-c/--categories" requires an argument\n'
        exit 1
    fi
    ;;
    -i|--input)
    if [ "$2" ]; then
        MYWD=$2
        shift 2
        if [ ! -d "${MYWD}" ]; then
            echo -e "\n${COL_red}ERROR: ${MYWD} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${MYWD})" ]; then
            echo -e "\n${COL_red}ERROR: ${MYWD} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-i/--input" requires an argument\n'
        exit 1
    fi
    ;;
    -g|--genome)
    if [ "$2" ]; then
        MYREF=$2
        shift 2
        if [ ! -d "${MYREF}" ]; then
            echo -e "\n${COL_red}ERROR: ${MYREF} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${MYREF})" ]; then
            echo -e "\n${COL_red}ERROR: ${MYREF} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-g/--genome" requires an argument\n'
        exit 1
    fi
    ;;
    -r|--report)
    if [ "$2" ]; then
        MYREPORT=$2
        shift 2
        if [ ! -d "${MYREPORT}" ]; then
            echo -e "\n${COL_red}ERROR: ${MYREPORT} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${MYREPORT})" ]; then
            echo -e "\n${COL_red}ERROR: ${MYREPORT} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-r/--report" requires an argument\n'
        exit 1
    fi
    ;;
    -S|--scripts)
    if [ "$2" ]; then
        MYSCRIPTS=$2
        shift 2
        if [ ! -d "${MYSCRIPTS}" ]; then
            echo -e "\n${COL_red}ERROR: ${MYSCRIPTS} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [ ! "$(ls -A ${MYSCRIPTS})" ]; then
            echo -e "\n${COL_red}ERROR: ${MYSCRIPTS} is empty!${COL_RESET} Please check \n"
            exit 1
        fi
    else
        echo -e '\nERROR: "-R/--Rscripts" requires an argument\n'
        exit 1
    fi
    ;;
    -s|--strain_ref)
    if [ "$2" ]; then
        STRAIN_REF=$2
        #echo STRAIN_REF=$STRAIN_REF
        shift 2
        if [[ ! -d "${MYWD}/$STRAIN_REF" ]]; then 
            echo -e "\n${COL_red}ERROR: ${STRAIN_REF} doesn't exist!${COL_RESET} Please check \n"
            exit 1
        fi
        if [[ -d "${MYWD}/$STRAIN_REF" ]]; then 
            echo -e "\n${COL_cyan}REFERENCE STRAIN: ${STRAIN_REF} will be used as reference strain!${COL_RESET} \n"
        fi
    else
        echo -e '\nERROR: "-r/--report" requires an argument\n'
        exit 1
    fi
    ;;
    -n|--number_exp)
    if [ "$2" ]; then
		if [ "$2" -eq "$2"  ] 2>/dev/null ; then
            NUMBER_EXP=$2
            shift 2
        else
            echo -e '\nERROR: "-n/--number_exp" requires a numeric argument'
            echo -e "argument parsed: $2 \n"
            exit 1
        fi
	else
		echo -e '\nERROR: "-n/--number_exp" requires a numeric argument\n'
        exit 1
    fi
    ;;
    -t|--threads)
	if [ "$2" ]; then
		if [ "$2" -eq "$2"  ] 2>/dev/null ; then
            NCPUS=$2
            shift 2
        else
            echo -e '\nERROR: "-t/--threads" requires a numeric argument'
            echo -e "argument parsed: $2 \n"
            exit 1
        fi
	else
		echo -e '\nERROR: "-t/--threads" requires a numeric argument\n'
        exit 1
    fi
	;;
    -v|--version)
        echo "${COL_yellow}Di-GRAPH version ${VERSION}${COL_RESET}"
        exit 1
    ;;
    -h|--help)
	    usage
	    exit 1
	;;
    -?*)
	    usage
        echo -e "\n${COL_red}ERROR: unknown option: ${1}${COL_RESET}\n"
	    exit 1
	;;
	*)
	    usage
        echo -e "\n${COL_red}ERROR: unknown option: ${1}${COL_RESET}\n"
        exit 1
    ;;
esac
done
set -- "${POSITIONAL[@]}" #restore positional parameters


REQUIRED_VARS=(MYBLAST CATEGORY_PATH MYWD MYREF MYREPORT MYSCRIPTS)
REQUIRED_OPTS=(
  "-b/--blast"
  "-c/--categories"
  "-i/--input"
  "-g/--genome"
  "-r/--report"
  "-S/--scripts"
)

missing=0

for i in "${!REQUIRED_VARS[@]}"; do
    var="${REQUIRED_VARS[$i]}"
    opt="${REQUIRED_OPTS[$i]}"

    if [ -z "${!var}" ]; then
        echo -e "${COL_red}ERROR:${COL_RESET} Required option $opt is missing"
        missing=1
    fi
done

if [ "$missing" -eq 1 ]; then
    usage
    exit 1
fi

if [[ "$STRAIN_REF" == "no" ]]; then
    echo -e "\n${COL_blue}REFERENCE STRAIN: no reference strain will be used${COL_RESET} \n"
fi

if [[ "$NUMBER_EXP" != 3 ]]; then
    echo -e "\n${COL_blue}NUMBER OF EXPERIMENTS: $NUMBER_EXP experiments will be analyzed${COL_RESET} \n"
fi

read -rp "Please, enter experiment names (comma-separated, e.g. E1,E3,E4): " NAME_EXP

IFS=',' read -ra EXP_LIST <<< "$NAME_EXP"


if [ "${#EXP_LIST[@]}" -ne "$NUMBER_EXP" ]; then
    echo "ERROR: You specified $NUMBER_EXP experiments but provided ${#EXP_LIST[@]} names."
    exit 1
fi
#echo "Working"
#################################

######## Starting the analysis ########
# Create a log and stat file
current_time=$(date "+%d-%m-%Y %H:%M:%S")
log_file="${MYWD}/log_file.txt"
stats="${MYWD}/stats.txt"

# Start timer for the entire script processing
start_time_total_script=$SECONDS

echo "CONCORDANT ANALYSIS - initiated at: \"$current_time\"" > "$log_file"
echo "" >> "$log_file"  # Adds a blank line

echo "=========  Bowtie Alignment ==========" >> "$stats"


# Index RG
bowtie-build "${MYREF}/RG_PMV_v9.fasta" "${MYREF}/S_cerevisiae_indexed" #NMQ --> MOVE TO INSTALL!


###########################################
### COVERAGE AND POLYMORPHISMS ANALYSIS ###
###########################################

echo "Coverage and Polymorphisms Analysis - initiated at: \"$current_time\"" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line

# Start timers for Coverage and Polymorphism Analysis processing

start_time_total_cov=$SECONDS

# Loop through each directory inside MYWD
for subdir in "$MYWD"*/; do
  echo "Processing directory: $subdir" >> "$log_file"
  echo "" >> "$log_file"  # Adds a blank line

  # Loop through sample prefixes (T0, T1, T2, etc.)
    for sample in T0 TSG TLG TLR; do
      for experiment in "${EXP_LIST[@]}"; do
        file1_gz="${subdir}/${sample}_${experiment}_R1.fastq.gz"
        file2_gz="${subdir}/${sample}_${experiment}_R2.fastq.gz"
       
        # Check if both compressed files exist
        if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
          echo "Processing files: $file1_gz, $file2_gz" >> "$log_file"
          echo "" >> "$log_file"  # Adds a blank line

          fq="${subdir}/${sample}_${experiment}"

          # [piped] One decompression pass replaces: gunzip -k R1/R2, cat into one FASTQ,
          # 14 cutadapt --cut runs and their 13 intermediate FASTQs. The output order is
          # unchanged (all R1 then all R2; 75nt = first halves then second halves; 18nt =
          # the 8 slices in the original glob order), so fastp sees identical input files.
          # fastp needs regular files as input (it reads them twice for adapter detection).
          "${DECOMP[@]}" "$file1_gz" "$file2_gz" | fq_slice \
            "1=${fq}_75nt_nonfiltered.fastq" \
            "2=${fq}_75nt_part2.tmp" \
            "111=${fq}_18nt_nonfiltered.fastq" \
            "112=${fq}_18nts_1_1_2.fastq" "121=${fq}_18nts_1_2_1.fastq" "122=${fq}_18nts_1_2_2.fastq" \
            "211=${fq}_18nts_2_1_1.fastq" "212=${fq}_18nts_2_1_2.fastq" "221=${fq}_18nts_2_2_1.fastq" \
            "222=${fq}_18nts_2_2_2.fastq"
          cat "${fq}_75nt_part2.tmp" >> "${fq}_75nt_nonfiltered.fastq"
          cat "${fq}_18nts"*".fastq" >> "${fq}_18nt_nonfiltered.fastq"
          rm "${fq}_75nt_part2.tmp" "${fq}_18nts"*".fastq"

          if [[ -s "${fq}_75nt_nonfiltered.fastq" ]]; then

            # Index files used by bedtools: regenerated once per run instead of once per
            # sample (the original overwrote the repo copy with the same content each time)
            if [[ -z "${CHROM_ORDER_DONE:-}" ]]; then
              samtools faidx "${MYREF}/RG_PMV_v9.fasta"
              cut -f1,2 "${MYREF}/RG_PMV_v9.fasta.fai" > "${MYREF}/chrom_order.txt"
              CHROM_ORDER_DONE=1
            fi

            # [piped] fastp writes its filtered reads into a FIFO read directly by bowtie,
            # and bowtie's SAM goes straight into samtools sort: no filtered FASTQ, no SAM.
            # The fastp command line (and so its html/json reports) is unchanged.

            # Filter 75nt reads Q30 -> Bowtie mapping -> sorted BAM (75nt)
            mkfifo "${fq}_75nt.fastq"
            fastp \
              -i "${fq}_75nt_nonfiltered.fastq" \
              -o "${fq}_75nt.fastq" \
              -q 30 -u 0 -e 30 \
              --thread "$N_CPU" \
              --html "${fq}_75nt.fastq.html" \
              --json "${fq}_75nt.fastq.json" &
            fastp_pid=$!
            echo "### $(basename "$subdir")/${sample}/${experiment} 75nt ###" >> "$stats"
            bowtie -p "$N_CPU" -m 1 -v 0 -S -x "${MYREF}/S_cerevisiae_indexed" "${fq}_75nt.fastq" 2>> "$stats" | \
              samtools sort -@ "$N_CPU" -o "${fq}_75nt.bam" -
            reap_fifo_producer "${fq}_75nt.fastq" "$fastp_pid" "${PIPESTATUS[0]}"
            echo "" >> "$stats"
            rm "${fq}_75nt_nonfiltered.fastq"

            # Filter 18nt reads Q30 -> Bowtie mapping in -k mode -> sorted BAM (18nt)
            mkfifo "${fq}_18nt.fastq"
            fastp \
              -i "${fq}_18nt_nonfiltered.fastq" \
              -o "${fq}_18nt.fastq" \
              -q 30 -u 0 -e 30 \
              --thread "$N_CPU" \
              --html "${fq}_18nt.fastq.html" \
              --json "${fq}_18nt.fastq.json" &
            fastp_pid=$!
            echo "### $(basename "$subdir")/${sample}/${experiment} 18nt ###" >> "$stats"
            bowtie -p "$N_CPU" -k 1 -v 0 -S -x "${MYREF}/S_cerevisiae_indexed" "${fq}_18nt.fastq" 2>> "$stats" | \
              samtools sort -@ "$N_CPU" -o "${fq}_18nt.bam" -
            reap_fifo_producer "${fq}_18nt.fastq" "$fastp_pid" "${PIPESTATUS[0]}"
            echo "" >> "$stats"
            rm "${fq}_18nt_nonfiltered.fastq"

            samtools index "${fq}_75nt.bam" "${fq}_75nt.bai"
            samtools index "${fq}_18nt.bam" "${fq}_18nt.bai"

            # Generate and sort BedGraphs coverage files (75nt)
            bamCoverage -b "${fq}_75nt.bam" -o "${fq}_75nt.bedgraph" -of bedgraph -p "$N_CPU" -bs 1 --normalizeUsing RPGC --effectiveGenomeSize 14272230
            bedtools sort -i "${fq}_75nt.bedgraph" -g "${MYREF}/chrom_order.txt" > "${fq}_75nt_sorted.bedgraph"

            # Generate and sort BedGraphs coverage files (18nt)
            bamCoverage -b "${fq}_18nt.bam" -o "${fq}_CHRIII_18nt.bedgraph" -of bedgraph -p "$N_CPU" -bs 1 --normalizeUsing RPGC --effectiveGenomeSize 14272230 -r CHRIII:199953:201553
            bamCoverage -b "${fq}_18nt.bam" -o "${fq}_CHRV_18nt.bedgraph" -of bedgraph -p "$N_CPU" -bs 1 --normalizeUsing RPGC --effectiveGenomeSize 14272230 -r CHRV:289025:290625
            bedtools sort -i "${fq}_CHRIII_18nt.bedgraph" -g "${MYREF}/chrom_order.txt" > "${fq}_CHRIII_18nt_sorted.bedgraph"
            bedtools sort -i "${fq}_CHRV_18nt.bedgraph" -g "${MYREF}/chrom_order.txt" > "${fq}_CHRV_18nt_sorted.bedgraph"

            # For TSG, TLG and TLR, generate BedGraphs coverage files (18nt) nonRPGC for quantifications
            if [[ "$sample" == "TSG" || "$sample" == "TLG" || "$sample" == "TLR" ]]; then
              bamCoverage -b "${fq}_18nt.bam" -o "${fq}_CHRIII_18nt_nonRPGC.bedgraph" -of bedgraph -p "$N_CPU" -bs 1 --effectiveGenomeSize 14272230 -r CHRIII:199953:201553
              bamCoverage -b "${fq}_18nt.bam" -o "${fq}_CHRV_18nt_nonRPGC.bedgraph" -of bedgraph -p "$N_CPU" -bs 1 --effectiveGenomeSize 14272230 -r CHRV:289025:290625
              bedtools sort -i "${fq}_CHRIII_18nt_nonRPGC.bedgraph" -g "${MYREF}/chrom_order.txt" > "${fq}_CHRIII_18nt_nonRPGC_sorted.bedgraph"
              bedtools sort -i "${fq}_CHRV_18nt_nonRPGC.bedgraph" -g "${MYREF}/chrom_order.txt" > "${fq}_CHRV_18nt_nonRPGC_sorted.bedgraph"
            fi

            # [piped] Per-base coverage (75nt) in ONE pass: bedtools coverage -> decimal
            # truncation -> genome-wide _75nt.tsv, plus the CHRIII/CHRV extracts (chr, cov,
            # position) and their MAT windows. Replaces the .tmp + mv rewrite of the
            # genome-wide TSV and the two extra full reads of it (one per chromosome).
            bedtools coverage -a "${fq}_75nt_sorted.bedgraph" -b "${fq}_75nt.bam" -sorted -g "${MYREF}/chrom_order.txt" -d | \
            awk -F'\t' -v out="${fq}_75nt.tsv" \
                -v c3="${fq}_75nt_CHRIII.tsv" -v c5="${fq}_75nt_CHRV.tsv" \
                -v m3="${fq}_75nt_CHRIII_MATa.tsv" -v m5="${fq}_75nt_CHRV_MATa.tsv" '
              BEGIN { OFS = FS; printf "" > c3; printf "" > c5; printf "" > m3; printf "" > m5 }
              { split($4, a, "."); if (length(a) > 1) $4 = a[1] "." substr(a[2], 1, 5); print > out }
              $1 == "CHRIII" { n3++; l = $1 OFS $4 OFS n3; print l > c3; if (n3 >= 199953 && n3 <= 201553) print l > m3 }
              $1 == "CHRV"   { n5++; l = $1 OFS $4 OFS n5; print l > c5; if (n5 >= 289025 && n5 <= 290625) print l > m5 }'

            # [piped] Per-base coverage (18nt), decimals truncated in the same pipe (no .tmp + mv)
            bedtools coverage -a "${fq}_CHRIII_18nt_sorted.bedgraph" -b "${fq}_18nt.bam" -sorted -g "${MYREF}/chrom_order.txt" -d | awk -F'\t' "$TRUNC4_AWK" > "${fq}_CHRIII_18nt.tsv"
            bedtools coverage -a "${fq}_CHRV_18nt_sorted.bedgraph" -b "${fq}_18nt.bam" -sorted -g "${MYREF}/chrom_order.txt" -d | awk -F'\t' "$TRUNC4_AWK" > "${fq}_CHRV_18nt.tsv"

            # For TSG, TLG and TLR, generate TSV files (18nt) nonRPGC for quantifications
            if [[ "$sample" == "TSG" || "$sample" == "TLG" || "$sample" == "TLR" ]]; then
              bedtools coverage -a "${fq}_CHRIII_18nt_nonRPGC_sorted.bedgraph" -b "${fq}_18nt.bam" -sorted -g "${MYREF}/chrom_order.txt" -d | awk -F'\t' "$TRUNC4_AWK" > "${fq}_CHRIII_18nt_nonRPGC.tsv"
              bedtools coverage -a "${fq}_CHRV_18nt_nonRPGC_sorted.bedgraph" -b "${fq}_18nt.bam" -sorted -g "${MYREF}/chrom_order.txt" -d | awk -F'\t' "$TRUNC4_AWK" > "${fq}_CHRV_18nt_nonRPGC.tsv"
            fi

            # Remove columnes 2, 3, 5, 6, and add coordinates to the last column (18nt)
            awk -F'\t' 'BEGIN {OFS="\t"} {coord=199953+NR-1; print $1, $4, coord}' "${fq}_CHRIII_18nt.tsv" > "${fq}_CHRIII_18nt_ordered.tsv"
            awk -F'\t' 'BEGIN {OFS="\t"} {coord=289025+NR-1; print $1, $4, coord}' "${fq}_CHRV_18nt.tsv" > "${fq}_CHRV_18nt_ordered.tsv"

            # For TSG, TLG and TLR, remove columnes 2, 3, 5, 6, and add coordinates to the last column (nonRPGC for quantifications)
            if [[ "$sample" == "TSG" || "$sample" == "TLG" || "$sample" == "TLR" ]]; then
              awk -F'\t' 'BEGIN {OFS="\t"} {coord=199953+NR-1; print $1, $4, coord}' "${fq}_CHRIII_18nt_nonRPGC.tsv" > "${fq}_CHRIII_18nt_nonRPGC_ordered.tsv"
              awk -F'\t' 'BEGIN {OFS="\t"} {coord=289025+NR-1; print $1, $4, coord}' "${fq}_CHRV_18nt_nonRPGC.tsv" > "${fq}_CHRV_18nt_nonRPGC_ordered.tsv"
            fi

            # Clean intermediate files
            rm -f "${fq}_75nt.bam" "${fq}_75nt.bai" "${fq}_18nt.bam" "${fq}_18nt.bai"

          else
          echo "Warning: Decompressed files missing for sample $sample in $subdir!" >> "$log_file"
          echo "" >> "$log_file"
          fi
        else
        echo "Warning: One or both compressed files for sample $sample are missing in $subdir!" >> "$log_file"
        echo "" >> "$log_file"
        fi
      done
    done
done




for strain in "$MYWD"*/; do
  echo "Processing directory: $strain to calculate average coverage" >> "$log_file"
  echo "" >> "$log_file"  # Adds a blank line

  # Loop through sample prefixes (T0, T1, T2, etc.)
    R_SCRIPT="SR_process_cov_III_V.R"
    strain="${strain}"
    root_dir="${strain}"
    echo "Processing strain: $strain" >> "$log_file"
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
            echo "Skipping ${strain} and continuing..." >> "$log_file"
          continue
        fi  
done




### === DATA ORGANIZATION ===

#Crate folder for the Alignment output files
for subdir in "${MYWD}"*/; do
    mkdir -p "${subdir}/Alignment_data"
    mv "${subdir}"/*.bedgraph "${subdir}/Alignment_data/"
done

# Create a folder for the Coverage TSV files
for subdir in "${MYWD}"*/; do
    mkdir -p "${subdir}/Coverage_data"
    mv "${subdir}"/*_Coverage.tsv "${subdir}/Coverage_data/"
done

for subdir in "${MYWD}"*/; do
    mkdir -p "${subdir}/MATs_quant_data"
    for chr in CHRIII CHRV; do
        mv "${subdir}"/*_"${chr}"_18nt_nonRPGC_ordered.tsv "${subdir}/MATs_quant_data/"
    done
done

# Create a folder for the graphs
for subdir in "$MYWD"*/; do
  mkdir -p "${subdir}/Graphs_Coverage"
  mv "${subdir}/plot_"*".svg" "${subdir}/Graphs_Coverage/"
done

# Create a folder for FASTQP files
for subdir in "$MYWD"*/; do
  mkdir -p "${subdir}/FASTQP"
  mv "${subdir}"/*.fastq.html "${subdir}/FASTQP/"
  mv "${subdir}"/*.fastq.json "${subdir}/FASTQP/"
done

# Remove index files
find "$MYREF" -name "*.ebwt" -exec rm {} +
find "$MYREF" -name "*.fai" -exec rm {} +
find "$MYREF" -name "*.txt" -exec rm {} +

# Delete all TSV files in subdir, excluding those in Alignment_data
for subdir in "${MYWD}"*; do
    find "${subdir}" -maxdepth 1 -type f -name "*.tsv" ! -name "*_75nt.tsv" -exec rm {} \;
done



### === MATs quantification for TSG, TLG and TLR ===
# This section calculates the difference in coverage for polymorphism between CHRIII and CHRV

# Coordinates (you can add more coordinates separated by space)
coordinates_CHRIII=(200119 200167 200212 200272 200326 200386 200449 200509 200542 200575 200635 200689 200753 200817 200882 200947 201012 201077 201148 201207 201272 201337 201402)
coordinates_CHRV=(289191 289239 289284 289344 289398 289458 289521 289581 289614 289647 289707 289761 289825 289889 289954 290019 290084 290149 290220 290278 290343 290408 290473)

offset=19 

for subdir in "$MYWD"*/MATs_quant_data/; do
  for sample in TSG TLG TLR; do
    for exp in "${EXP_LIST[@]}"; do
      input_file_CHRIII="${subdir}/${sample}_${exp}_CHRIII_18nt_nonRPGC_ordered.tsv"
      input_file_CHRV="${subdir}/${sample}_${exp}_CHRV_18nt_nonRPGC_ordered.tsv"

      if [[ -f "$input_file_CHRIII" && -f "$input_file_CHRV" ]]; then
        echo "Processing files in $subdir for $sample _ $exp..."

        # Loop through CHRIII coordinates
        for coordinate_CHRIII in "${coordinates_CHRIII[@]}"; do
          line_number_CHRIII=$(awk -v coord="$coordinate_CHRIII" '$3 == coord {print NR}' "$input_file_CHRIII")

          if [[ -n "$line_number_CHRIII" ]]; then
            top_line_CHRIII=$((line_number_CHRIII - offset))
            bottom_line_CHRIII=$((line_number_CHRIII + offset))
            [[ "$top_line_CHRIII" -lt 1 ]] && top_line_CHRIII=1

            top_line_content=$(sed -n "${top_line_CHRIII}p" "$input_file_CHRIII")
            mid_line_content=$(sed -n "${line_number_CHRIII}p" "$input_file_CHRIII")
            bottom_line_content=$(sed -n "${bottom_line_CHRIII}p" "$input_file_CHRIII")

            top_coverage=$(echo "$top_line_content" | awk '{print $2}')
            bottom_coverage=$(echo "$bottom_line_content" | awk '{print $2}')
            baseline_coverage=$(awk -v top="$top_coverage" -v bottom="$bottom_coverage" 'BEGIN {print top + (bottom - top)/2}')

            # Extract fields from mid_line_content
            mid_chr=$(echo "$mid_line_content" | awk '{print $1}')
            mid_coverage=$(echo "$mid_line_content" | awk '{print $2}')
            mid_coordinate=$(echo "$mid_line_content" | awk '{print $3}')

            # Calculate coverage difference (mid - baseline)
            coverage_diff=$(awk -v mid="$mid_coverage" -v base="$baseline_coverage" 'BEGIN {print mid - base}')

            # Build new line with: original chr, diff coverage, and coordinate (tab-separated)
            diff_line=$(echo -e "$mid_chr\t$coverage_diff\t$mid_coordinate")

            echo "CHRIII Coord: $coordinate_CHRIII — Coverage Diff: $coverage_diff"

            output_file_CHRIII="${subdir}/${sample}_${exp}_CHRIII_diff_${coordinate_CHRIII}.tsv"
            {
              echo "$diff_line"
            } > "$output_file_CHRIII"

          else
          echo "Coordinate $coordinate_CHRIII not found in $input_file_CHRIII for $sample $exp."
          fi
        done
    
        # Loop through CHRV coordinates
        for coordinate_CHRV in "${coordinates_CHRV[@]}"; do
          line_number_CHRV=$(awk -v coord="$coordinate_CHRV" '$3 == coord {print NR}' "$input_file_CHRV")

          if [[ -n "$line_number_CHRV" ]]; then
            top_line_CHRV=$((line_number_CHRV - offset))
            bottom_line_CHRV=$((line_number_CHRV + offset))
            [[ "$top_line_CHRV" -lt 1 ]] && top_line_CHRV=1

            top_line_content=$(sed -n "${top_line_CHRV}p" "$input_file_CHRV")
            mid_line_content=$(sed -n "${line_number_CHRV}p" "$input_file_CHRV")
            bottom_line_content=$(sed -n "${bottom_line_CHRV}p" "$input_file_CHRV")

            top_coverage=$(echo "$top_line_content" | awk '{print $2}')
            bottom_coverage=$(echo "$bottom_line_content" | awk '{print $2}')
            baseline_coverage=$(awk -v top="$top_coverage" -v bottom="$bottom_coverage" 'BEGIN {print top + (bottom - top)/2}')

            # Extract fields from mid_line_content
            mid_chr=$(echo "$mid_line_content" | awk '{print $1}')
            mid_coverage=$(echo "$mid_line_content" | awk '{print $2}')
            mid_coordinate=$(echo "$mid_line_content" | awk '{print $3}')

            # Calculate coverage difference (mid - baseline)
            coverage_diff=$(awk -v mid="$mid_coverage" -v base="$baseline_coverage" 'BEGIN {print mid - base}')

            # Build new line with: original chr, diff coverage, and coordinate (tab-separated)
            diff_line=$(echo -e "$mid_chr\t$coverage_diff\t$mid_coordinate")

            echo "CHRV Coord: $coordinate_CHRV — Coverage Diff: $coverage_diff"

            output_file_CHRV="${subdir}/${sample}_${exp}_CHRV_diff_${coordinate_CHRV}.tsv"
            {
              echo "$diff_line"
            } > "$output_file_CHRV"

          else
          echo "Coordinate $coordinate_CHRV not found in $input_file_CHRV for $sample $exp."
          fi
        done
      else
      echo "Files not found for $sample $exp in $subdir"
      fi
    done
  done
done


for strain in "$MYWD"*/; do
  echo "Processing directory: $strain to calculate MATa coverage (r18)" >> "$log_file"
  echo "" >> "$log_file"  # Adds a blank line
  strain="${strain}"
  root_dir="${strain}"

  # Loop through sample prefixes (T0, T1, T2, etc.)
    R_SCRIPT="SR_process_cov_III_V_18nt.R"
    strain="${strain}"
    root_dir="${strain}"
    echo "Processing strain: $strain" >> "$log_file"
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
            echo "Skipping ${strain} and continuing..." >> "$log_file"
          continue
        fi  
done


# Move the plots to the Graphs_Coverage folder
for strain in "$MYWD"*/; do
  for sample in TSG TLG TLR; do
    mv "${strain}/plot_${sample}_18nt_MATa_Coverage_Quant.svg" "${strain}/Graphs_Coverage/"
  done
done

# Calulate total elapsed time
elapsed_time_total_cov=$((( SECONDS - start_time_total_cov )/60))
echo "Coverage and Polymorphisms Analysis completed in ${elapsed_time_total_cov} minutes" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line            


###########################################
####### GENOMIC CATEGORIES ANALYSIS #######
###########################################

current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Genomic Categories Analysis - initiated at: \"$current_time\"" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

# Start timer for the R category analysis
start_time_total_cat=$SECONDS

# Loop through each directory inside MYWD
# [piped] ONE R session per strain for all 13 categories (was one session per category).
# SR_process_categories_all.R reads the 12 genome-wide *_75nt.tsv coverage files once
# instead of 13 times and runs the unchanged per-category computation. It writes a
# status file, from which the same log lines as before are written, in the same order.
R_SCRIPT="SR_process_categories_fast.R"
CATEGORIES=(ORF LTR TEG Ty tRNA rRNA ncRNA snRNA snoRNA ARS Cen Tel Int)
for strain in "$MYWD"*/; do
    root_dir="${strain}"
    category_path="${CATEGORY_PATH}"
    cat_status="${strain}/.categories_status.tsv"

    Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$category_path" "$cat_status" "${CATEGORIES[@]}"

    for category in "${CATEGORIES[@]}"; do
        echo "Processing category: $category" >> "$log_file"
        if [[ "$(awk -F'\t' -v c="$category" '$1 == c {print $2}' "$cat_status" 2> /dev/null)" != "OK" ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} and ${category}" >> "$log_file"
            echo "Skipping ${strain} and ${category} and continuing..." >> "$log_file"
        fi
    done
    rm -f "$cat_status"
done


### ORDERING CATEGORIES DATASET ###

# Loop through each directory inside MYWD
R_SCRIPT="SR_order_categories.R"
for strain in "${MYWD}"*/; do
  echo "Processing directory: $strain"
  for analysis_suffix in "T0vsTSG" "T0vsTLG" "T0vsTLR"; do
    echo "Processing analysis: $analysis_suffix"
    root_dir="${strain}"
    strain="${strain}"
    category_path="${CATEGORY_PATH}"
    analysis_suffix="${analysis_suffix}"

    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$category_path" "$analysis_suffix"; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}, ${category} and ${analysis_suffix}" >> "$log_file"
    echo "Skipping ${strain}, ${category} and ${analysis_suffix} and continuing..." >> "$log_file"
    continue
    fi  
  done
done

# Move tsv category files into a category data folder
for strain in "${MYWD}"*/; do
  mkdir -p "${strain}/Category_Data_TSV_and_Plots"
  mv "${strain}"/*"fingerprint"* "${strain}/Category_Data_TSV_and_Plots/"
  #mv "${strain}"/*"all_data"* "${strain}/Category_Data_TSV_and_Plots/"
  mv "${strain}"/*".svg" "${strain}/Category_Data_TSV_and_Plots/"
done

### PLOTTING CATEGORIES DATA ###
 R_SCRIPT="SR_plot_categories.R"
for strain in "$MYWD"*/; do
  for analysis_suffix in "T0vsTSG" "T0vsTLG" "T0vsTLR"; do
    echo "Processing analysis: $analysis_suffix"
    root_dir="${strain}"
    strain="${strain}"
    analysis_suffix="${analysis_suffix}"
    echo "Processing time point: $analysis_suffix"
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$analysis_suffix"; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} and ${analysis_suffix}" >> "$log_file"
    echo "Skipping ${strain} and ${analysis_suffix} and continuing..." >> "$log_file"
    continue
    fi  

  done
done

# R plots for gal vs raf coverage
# Start timer for R processing

# Loop inside each subdirectory of MYWD
R_SCRIPT="SR_plot_gal_vs_raf.R"
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      strain="${strain}"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Plotting gal vs raf coverage" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
        echo "Skipping ${strain} and continuing..." >> "$log_file"
        continue
        fi
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total gal vs raf plotting completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done


# Move tsv category files into a category data folder
for subdir in "${MYWD}"*/; do
  mkdir -p "${subdir}/Sorted_datasets"
  mv "${subdir}"/*"sorted"* "${subdir}/Sorted_datasets"
  mv "${subdir}"/*"Template"* "${subdir}/Sorted_datasets"
  mv "${subdir}"/*"UP_DOWN"* "${subdir}/Sorted_datasets"
  mv "${subdir}"/*".svg" "${subdir}/Plots_Analysis/"
done

# Move tsv files into a alignmentfolder
for subdir in "${MYWD}"*/; do
  mkdir -p "${subdir}/Alignment_data/TSV_files"
  mv "${subdir}"/*".tsv" "${subdir}/Alignment_data/TSV_files"
done

# Calulate total elapsed time
elapsed_time_total_cat=$((( SECONDS - start_time_total_cat )/60))
echo "Genomic Categories Analysis completed in ${elapsed_time_total_cat} minutes" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line

##########################################
########### MUTAGENIC RATE AT HO ########## #Updated on 30/09/2025 #R scripts updated on 07/10/2026 #Python scripts updated on 11/04/2026
###########################################

current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Mutagenic Rate at HO Analysis - initiated at: \"$current_time\"" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

#Add a line to separate bowtie and bwa alignments in the stat file
echo "=========  BWA Alignment ==========" >> "$stats"
echo "=========  Mutagenic rate at HO - BWA Alignment ==========" >> "$log_file"

# Set locale
export LC_ALL=C
export LANG=C

# Start timer for the entire script processing
start_time_total_mut=$SECONDS

# Index RG
bwa index "${MYREF}/RG_PMV_v9_CHRIII.fasta"

# Loop through each directory inside MYWD for alignments
for subdir in "$MYWD"*/; do
  echo "Processing directory: $subdir"
    for sample in T0 TLG TLR; do
      for exp in "${EXP_LIST[@]}"; do
        file1_gz="${subdir}/${sample}_${exp}_R1.fastq.gz"
        file2_gz="${subdir}/${sample}_${exp}_R2.fastq.gz"
       
        # Check if both compressed files exist
        if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
          echo "Processing files: $file1_gz, $file2_gz"

          # [piped] Decompress straight into the concatenated FASTQ (no R1/R2 .fastq copies).
          # fastp reads its input twice (adapter detection), so this one stays a file.
          "${DECOMP[@]}" "$file1_gz" "$file2_gz" > "${subdir}/${sample}_${exp}.fastq"

          #Filtering and alinging reads
          if [[ -f "${subdir}/${sample}_${exp}.fastq" ]]; then

            # Detecting number of available CPUs
            N_CPU=$(nproc)

            # [piped] fastp -> FIFO -> bwa mem -> header + CHRIII filter -> samtools sort.
            # Replaces Q30.fastq, Q30.sam, cabecera.sam, Q30_CHRIII.sam and
            # Q30_CHRIII_header.sam. awk 'NR<=2 || $3=="CHRIII"' = the original
            # "sed -n 1,2p" header + "awk '$3 == "CHRIII"'" records, in the same order.
            mkfifo "${subdir}/${sample}_${exp}_Q30.fastq"
            fastp \
            -i "${subdir}/${sample}_${exp}.fastq" \
            -o "${subdir}/${sample}_${exp}_Q30.fastq" \
            -q 30 -u 0 -e 30 \
            --thread "$N_CPU" \
            --html "${subdir}/${sample}_${exp}_Q30.fastq.html" \
            --json "${subdir}/${sample}_${exp}_Q30.fastq.json" &
            fastp_pid=$!
            bwa mem -t "$N_CPU" "${MYREF}/RG_PMV_v9_CHRIII.fasta" "${subdir}/${sample}_${exp}_Q30.fastq" | \
              awk 'NR <= 2 || $3 == "CHRIII"' | \
              samtools sort -@ "$N_CPU" -o "${subdir}/${sample}_${exp}_processed.bam" -
            reap_fifo_producer "${subdir}/${sample}_${exp}_Q30.fastq" "$fastp_pid" "${PIPESTATUS[0]}"
            samtools index "${subdir}/${sample}_${exp}_processed.bam" "${subdir}/${sample}_${exp}_processed.bai"

            # Extact reads mapping to the HO site ([piped] no HOs.sam)
            samtools view "${subdir}/${sample}_${exp}_processed.bam" "CHRIII:200753-200753" -h -O SAM | \
              samtools sort -@ "$N_CPU" -o "${subdir}/${sample}_${exp}_HOs.bam" -
            samtools index "${subdir}/${sample}_${exp}_HOs.bam" "${subdir}/${sample}_${exp}_HOs.bai"
            echo "" >> "$stats"
            echo "### $(basename "$subdir")/${sample}/${exp} 150nt ###" >> "$stats"
            echo "BWA Alignment: Reads mapping the HO site" >> "$stats"
            samtools flagstat "${subdir}/${sample}_${exp}_HOs.bam" | head -n 1 >> "$stats"

            # Remove intermediate files
            rm -f "${subdir}"/*.fastq "${subdir}"/*.sam
            
          else
          echo "Warning: Decompressed files missing for sample $sample in $subdir!" >> "$log_file"
          echo "" >> "$log_file"
          fi
        else
        echo "Warning: One or both compressed files for sample $sample are missing in $subdir!" >> "$log_file"
        echo "" >> "$log_file"
        fi
      done
    done
done

# Remove intermediate files
find "$MYWD" -type f \( -name "*.sam" -o -name "*.fastq" \) -exec rm {} + # intermediate alignmet files

### DATA ANALYSIS ###

# Start timer for data analysis and plotting

### Invoque python script to extract reads with A or G at postion 200689 ###
PYTHON_SCRIPT="MUT_extract_200689.py"

for subdir in "$MYWD"*/; do
  for sample in T0 TLG TLR; do
    for exp in "${EXP_LIST[@]}"; do
      bam_file="${subdir}/${sample}_${exp}_HOs.bam"
      output_bam_G="${subdir}/${sample}_${exp}_200689_G.bam"
      output_bam_A="${subdir}/${sample}_${exp}_200689_A.bam"
      if [[ -f "$bam_file" ]]; then
        echo "Processing BAM file: $bam_file"
      else
        echo "Warning: BAM file $bam_file does not exist!"
      fi

      if ! python3 ${MYSCRIPTS}/python_files/${PYTHON_SCRIPT} "$bam_file" "$output_bam_A" "$output_bam_G" 2>> "$log_file"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR processing ${bam_file}" >> "$log_file"
        echo "Skipping file and continuing..." >> "$log_file"
        continue
      fi
  
        #Generate bai index for the new bam files
        samtools index "${output_bam_A}" "${output_bam_A}.bai"
        samtools index "${output_bam_G}" "${output_bam_G}.bai"

        # Add stats to the log file
        echo "Read count for A at 200689: $(samtools view -c "${output_bam_A}")" 
        #samtools view -c "${output_bam_A}" > "${subdir}/${sample}_${exp}_200689A_row_count.tsv"
        echo "Read count for G at 200689: $(samtools view -c "${output_bam_G}")" 
        samtools view -c "${output_bam_G}" > "${subdir}/${sample}_${exp}_200689G_row_count.tsv"
    done
  done
done

### Invoque python script to extract reads with T or C at postion 200753 from reads with A at position 200689 ###
PYTHON_SCRIPT="MUT_extract_200753.py"

for subdir in "$MYWD"*/; do
  for sample in T0 TLG TLR; do
    for exp in "${EXP_LIST[@]}"; do
      bam_file="${subdir}/${sample}_${exp}_200689_A.bam"
      output_bam_T="${subdir}/${sample}_${exp}_200689_A_200753_T.bam"
      output_bam_C="${subdir}/${sample}_${exp}_200689_A_200753_C.bam"
      if [[ -f "$bam_file" ]]; then
        echo "Processing BAM file: $bam_file"
      else
        echo "Warning: BAM file $bam_file does not exist!"
      fi

      if ! python3 ${MYSCRIPTS}/python_files/${PYTHON_SCRIPT} "$bam_file" "$output_bam_T" "$output_bam_C" 2>> "$log_file"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR processing ${bam_file}" >> "$log_file"
        echo "Skipping file and continuing..." >> "$log_file"
        continue
      fi
  


    #Generate bai index for the new bam files
    samtools index "${output_bam_T}" "${output_bam_T}.bai"
    samtools index "${output_bam_C}" "${output_bam_C}.bai"
    # Add stats to the log file
    echo "Read count for T at 200753 and A at 200689: $(samtools view -c "${output_bam_T}")" 
    samtools view -c "${output_bam_T}" > "${subdir}/${sample}_${exp}_200689AT_row_count.tsv"
    echo "Read count for C at 200753 and A at 200689: $(samtools view -c "${output_bam_C}")"
    samtools view -c "${output_bam_C}" > "${subdir}/${sample}_${exp}_200689AC_row_count.tsv"



    # Merge the two bam files
    samtools merge -f "${subdir}/${sample}_${exp}_200689_A_CT.bam" "${output_bam_C}" "${output_bam_T}"
    samtools view -c "${subdir}/${sample}_${exp}_200689_A_CT.bam" > "${subdir}/${sample}_${exp}_200689ACT_row_count.tsv"
    samtools index "${subdir}/${sample}_${exp}_200689_A_CT.bam" "${subdir}/${sample}_${exp}_200689_A_CT.bai"
    samtools merge -f "${subdir}/${sample}_${exp}_200689_A_CT_G.bam" "${subdir}/${sample}_${exp}_200689_A_CT.bam" "${subdir}/${sample}_${exp}_200689_G.bam"
    samtools view -c "${subdir}/${sample}_${exp}_200689_A_CT_G.bam" > "${subdir}/${sample}_${exp}_200689ACTG_row_count.tsv"
    samtools index "${subdir}/${sample}_${exp}_200689_A_CT_G.bam" "${subdir}/${sample}_${exp}_200689_A_CT_G.bai"
    
    ### Assesing IGVtools for Del and Ins ###
    
    # Extract base scores usign bcftools
    ref_fasta="${MYREF}/RG_PMV_v9_CHRIII.fasta"

    # Define the path to igvtools
    bam_file_igv_A_CT="${subdir}/${sample}_${exp}_200689_A_CT.bam"
    bam_file_igv_A_CT_G="${subdir}/${sample}_${exp}_200689_A_CT_G.bam"
    bam_file_igv_A_200753_T="${subdir}/${sample}_${exp}_200689_A_200753_T.bam"
    bam_file_igv_A_200753_C="${subdir}/${sample}_${exp}_200689_A_200753_C.bam"
    bam_file_igv_G="${subdir}/${sample}_${exp}_200689_G.bam"
    ref_fasta="${MYREF}/RG_PMV_v9_CHRIII.fasta"

    # Create a loop to assess 200689  ACTG, ATG, A_T, A_C or G polymorphisms
      ### Assesisng IGVtools - Del and Ins ###
      for poly in A_CT A_CT_G A_200753_T A_200753_C G ; do
        wig_output="${subdir}/${sample}_${exp}_200689_${poly}.wig"
        var_name="bam_file_igv_${poly}"
        bam_file="${!var_name}" # indirect expansion
      
        # Run igvtools to generate the WIG file
        igvtools count -w 1 --bases --query CHRIII:200743-200762 "$bam_file" "$wig_output" "$ref_fasta"


        ### Assesing nucleotide variants by BCFtools ###

        tsv_output_BCF="${subdir}/${sample}_${exp}_${poly}_BCF.tsv"
        bam_file="${subdir}/${sample}_${exp}_200689_${poly}.bam"      
        vcf_file="${subdir}/${sample}_${exp}_${poly}_variants.vcf.gz"

        # Variant calling
        bcftools mpileup -Ou -f "$ref_fasta" -r CHRIII:200742-200762 "$bam_file" | \
        bcftools call -mA --ploidy 1 -Oz -o "$vcf_file"
        bcftools index "$vcf_file"

        # Fill in missing INFO tags like AF, AN, AC
        vcf_with_af="${vcf_file%.vcf.gz}_with_AF.vcf.gz"
        bcftools +fill-tags "$vcf_file" -Oz -o "$vcf_with_af" -- -t AF,AN,AC
        bcftools index "$vcf_with_af"

        # Extract the region and required fields
        filtered_vcf="${vcf_with_af%.vcf.gz}.filtered.vcf"
        bcftools view -r CHRIII:200743-200762 "$vcf_with_af" > "$filtered_vcf"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\n' "$vcf_file" | \
        awk 'BEGIN { OFS="\t" }
        {
            split($6, dp4, ",");
            ref_count = dp4[1] + dp4[2];
            alt_count = dp4[3] + dp4[4];
            total = ref_count + alt_count;
            af = (total == 0) ? 0 : int((alt_count / total) * 100);
            print $1, $2, $3, $4, $5, af
        }' > "$tsv_output_BCF"

        find "$MYWD" -type f \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.csi" \) -exec rm {} + # BCF intermediate files
         done
    done
  done
done

### Calculate repair pathway frequency
R_SCRIPT="MUT_calc_rep_path_freq.R"

for strain in "${MYWD}"*/; do

    echo "Processing directory: ${strain}" >> "$log_file"

    root_dir="${strain}"
    wd_dir="${MYWD}"

    echo "$(basename "$strain")" >> "$log_file"
    echo "" >> "$log_file"
    echo "Plotting repair_pathway for ${strain}" >> "$log_file"

    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
        echo "Skipping ${strain} and continuing..." >> "$log_file"
        continue
    fi

done


### Plot repair pathway comparison
R_SCRIPT="MUT_plot_rep_path_comparison.R"

echo "Processing directory: ${MYWD}" >> "$log_file"
root_dir="${MYWD}"
echo "Plotting repair_pathway comparison" >> "$log_file"
start_time=$SECONDS
if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir"; then
  echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
  echo "Skipping ${root_dir} and continuing..." >> "$log_file"
fi
      

# Calculate total elapsed time
elapsed_time_total_mut=$((( SECONDS - start_time_total_mut )/60))
echo "Mutagenic Rate at HO Analysis completed in ${elapsed_time_total_mut} minutes" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line



# === ORGANIZING DATA ===

# Remove intermediate files
find "$MYWD" -type f \( -name "*_BFC.tsv" \) -exec rm {} + # BCF.tsv intermediate files

# Move data into subfolders
for subdir in "$MYWD"*/; do
  mkdir -p "${subdir}/HO_mut_rate/FASTQP" "${subdir}/HO_mut_rate/SVG" "${subdir}/HO_mut_rate/TSV"
  mv "${subdir}"/*.fastq.html "${subdir}/HO_mut_rate/FASTQP/" 2> /dev/null
  mv "${subdir}"/*.fastq.json "${subdir}/HO_mut_rate/FASTQP/" 2> /dev/null
  mv "${subdir}"/*.svg "${subdir}/HO_mut_rate/SVG/" 2> /dev/null
  mv "${subdir}"/*.tsv "${subdir}/HO_mut_rate/TSV/" 2> /dev/null
done

# Remove unecesary files
find "$MYREF" \( -name "*.amb" -o -name "*.ann" -o -name "*.bwt" -o -name "*.pac" -o -name "*.sa" -o -name "*.fai" \) -exec rm {} + # index files
find "$MYWD" -type f \( -name "*.bai" -o -name "*.bam" \) -exec rm {} + # intermediate alignmet files
find "$MYWD" -type f \( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.csi" \) -exec rm {} + # BCF intermediate files
find "$MYWD" -type f \( -name "*.wig" \) -exec rm {} + # WIG intermediate files



# Calculate total elapsed time
elapsed_time_total_mut=$((( SECONDS - start_time_total_mut )/60))
echo "HO Mutagenic Rate Analysis completed in ${elapsed_time_total_mut} minutes" >> "$log_file"
echo "" >> "$log_file"


# Create a log and stat file
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "DISCORDANT ANALYSIS - initiated at: \"$current_time\"" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line

# Start timer for the entire script processing
#start_time_total_script=$SECONDS
start_time_total_script_2=$SECONDS

########### MODULE CORE - START ###########
# ############## Alignment and processing of 75bp reads ############

# Start timer for the alignment and processing section
start_time_total_alignment=$SECONDS

# Index RG for bowtie alignment
bowtie-build "${MYREF}/RG_PMV_v9.fasta" "${MYREF}/S_cerevisiae_indexed"

# Index RG for bowtie2 alignment
bowtie2-build "${MYREF}/RG_PMV_v9.fasta" "${MYREF}/S_cerevisiae_indexed"

# Loop through each directory inside MYWD

echo "================================" >> "$log_file"
echo "### 75bp analysis ###" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### Alignments and SAM processing ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file" 
echo "" >> "$log_file"  # Adds a blank line

for strain in "$MYWD"*/; do
  echo "Processing directory: $strain" >> "$log_file"

  # Loop through sample prefixes (T0, TSG, TLG, TLR.)
    for sample in T0 TSG TLG TLR; do
      for experiment in "${EXP_LIST[@]}"; do
        file1_gz="${strain}/${sample}_${experiment}_R1.fastq.gz"
        file2_gz="${strain}/${sample}_${experiment}_R2.fastq.gz"
       
        # Check if both compressed files exist
        if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
          echo "Processing files: $file1_gz, $file2_gz" >> "$log_file"
          echo "" >> "$log_file"  # Adds a blank line

          fq="${strain}/${sample}_${experiment}"

          # [piped] no gunzip -k: the R1/R2 FASTQ are streamed from the .gz files below
          if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
            # Determine the number of cores available
            N_CPU=$(nproc)

            # Start timers for total processing of a strain/sample/experiment
            start_time_total=$SECONDS

            # Trim reads using cutadapt
            echo "### $(basename "$strain")/${sample} ${experiment} ###" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line
            start_time=$SECONDS

            # [piped] r_R1.fastq = first 75 nt of every R1 read renamed .../rA/1, then the
            # last 75 nt renamed .../rB/1 (same for R2). Identical to the original
            # gunzip + 4x cutadapt + 4x sed + cat, without the 10 intermediate FASTQs.
            # (Two decompression passes per mate keep the halves in the original order.)
            { "${DECOMP[@]}" "$file1_gz" | fq_slice "1:rA1=/dev/stdout"
              "${DECOMP[@]}" "$file1_gz" | fq_slice "2:rB1=/dev/stdout"; } > "${fq}_r_R1.fastq" &
            pid_r1=$!
            { "${DECOMP[@]}" "$file2_gz" | fq_slice "1:rA2=/dev/stdout"
              "${DECOMP[@]}" "$file2_gz" | fq_slice "2:rB2=/dev/stdout"; } > "${fq}_r_R2.fastq" &
            pid_r2=$!
            wait "$pid_r1" "$pid_r2"

            # Filter 75nt reads Q30
            fastp \
              -i "${strain}/${sample}_${experiment}_r_R1.fastq" \
              -o "${strain}/${sample}_${experiment}_r_R1_filtered.fastq" \
              -I "${strain}/${sample}_${experiment}_r_R2.fastq" \
              -O "${strain}/${sample}_${experiment}_r_R2_filtered.fastq" \
              -q 30 -u 0 -e 30 \
              --thread "$N_CPU" \
              --html "${strain}/${sample}_${experiment}_75nt.fastq.html" \
              --json "${strain}/${sample}_${experiment}_75nt.fastq.json"

            # Remove intermediate 75nt fastq files
            rm "${strain}/${sample}_${experiment}_r_R1.fastq" "${strain}/${sample}_${experiment}_r_R2.fastq"

            elapsed_time=$(( SECONDS - start_time ))
            echo "Cutting and Filtering Reads completed in ${elapsed_time} seconds" >> "$log_file"

            # [piped] Bowtie SR mapping (-m1 -v0) streamed straight to the sorted list of
            # mapped read names: no SAM, mapped SAM, name list or unsorted list on disk.
            start_time=$SECONDS
            echo "### $(basename "$strain")/${sample} 75nt R1 ###" >> "$stats"
            echo "${experiment}" >> "$stats"
            bowtie -p "$N_CPU" -m 1 -v 0 -S "${MYREF}/S_cerevisiae_indexed" "${fq}_r_R1_filtered.fastq" 2>> "$stats" | \
              samtools view -F 4 - | cut -f1 | sed 's/\/1//' | sort > "${fq}_r_R1_mapped_sorted.txt"
            echo "" >> "$stats"  # Adds a blank line
            elapsed_time=$(( SECONDS - start_time ))
            echo "Bowtie SR alignment 75nt for R1 reads completed in ${elapsed_time} seconds" >> "$log_file"

            start_time=$SECONDS
            echo "### $(basename "$strain")/${sample} 75nt R2 ###" >> "$stats"
            echo "${experiment}" >> "$stats"
            bowtie -p "$N_CPU" -m 1 -v 0 -S "${MYREF}/S_cerevisiae_indexed" "${fq}_r_R2_filtered.fastq" 2>> "$stats" | \
              samtools view -F 4 - | cut -f1 | sed 's/\/2//' | sort > "${fq}_r_R2_mapped_sorted.txt"
            echo "" >> "$stats"  # Adds a blank line
            elapsed_time=$(( SECONDS - start_time ))
            echo "Bowtie SR alignment 75nt for R2 reads completed in ${elapsed_time} seconds" >> "$log_file"

            comm -12 "${fq}_r_R1_mapped_sorted.txt" "${fq}_r_R2_mapped_sorted.txt" > "${fq}_r_R12_mapped.txt"
            rm "${fq}_r_R1_mapped_sorted.txt" "${fq}_r_R2_mapped_sorted.txt"

            # Bowtie2 mapping (in a forced -m1 -v0 mode)
            start_time=$SECONDS
            echo "### $(basename "$strain")/${sample} 75nt paired-end Bowtie2 alignment ###" >> "$stats"
            echo "${experiment}" >> "$stats"
            bowtie2 \
              -p "$N_CPU" \
              --no-1mm-upfront \
              --score-min C,0,0 \
              -N 0 \
              --end-to-end \
              --fr \
              -x "${MYREF}/S_cerevisiae_indexed" \
              -1 "${strain}/${sample}_${experiment}_r_R1_filtered.fastq" \
              -2 "${strain}/${sample}_${experiment}_r_R2_filtered.fastq" \
              -S "${strain}/${sample}_${experiment}.sam" 2>> "$stats"

            echo "" >> "$stats"  # Adds a blank line
            elapsed_time=$(( SECONDS - start_time ))
            echo "Bowtie2 alignment 75nt for paired-end reads completed in ${elapsed_time} seconds" >> "$log_file"
            echo "" >> "$log_file"

            # Remove intermediate files
            rm "${strain}/${sample}_${experiment}_r_R1_filtered.fastq"
            rm "${strain}/${sample}_${experiment}_r_R2_filtered.fastq"

            # Header kept as a small file (the original "sed -n 1,19p")
            head -n 19 "${fq}.sam" > "${fq}_header.sam"

            #########
            # Bowtie2 SAM processing (concordant read-pairs for % recombination)
            start_time=$SECONDS
            echo "### $(basename "$strain")/${sample} Bowtie2 SAM processing - concordant read pairs ###" >> "$log_file"
            echo "${experiment}" >> "$log_file"

            # [piped] Only the row count of concordant_pairs_unique.sam is kept, and a
            # count of unique lines does not depend on line order. So the whole chain
            # (alignment.sam, alignment_processed*.sam x6, concordant_pairs*.sam x2) is one
            # stream; the R1/R2 halves of each joined line are emitted with the exact
            # semantics of "cut -f1-19" and "cut -f20-38".
            tail -n +20 "${fq}.sam" | \
            awk -vOFS='\t' '{$11="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; print}' | \
            cut -f1-19 | \
            sed '$!N;s/\n/ /' | \
            awk '($3 != $7 && $7 == "=")' | \
            awk '($7 != "*")' | \
            awk '$6 == "75M"' | \
            awk '$25 == "75M"' | \
            awk '$13 == "XN:i:0"' | \
            awk '$32 == "XN:i:0"' | \
            sed 's/ /\t/g' | \
            awk '($19 != "YT:Z:UP")' | \
            awk '($38 != "YT:Z:UP")' | \
            awk -F'\t' -v OFS='\t' "$CUT_HALVES_AWK" | \
            cat "${fq}_header.sam" - | \
            samtools view -N "${fq}_r_R12_mapped.txt" - | \
            cat "${fq}_header.sam" - | \
            awk '!seen[$0]++' | \
            tail -n +20 | wc -l | awk '{print $1}' > "${strain}/${sample}_${experiment}_concordant_pairs_unique_row_count.tsv"

            echo "" >> "$log_file"  # Adds a blank line
            elapsed_time=$(( SECONDS - start_time ))
            echo "Concordant read pairs processing completed in ${elapsed_time} seconds" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line


            ########

            # Bowtie2 SAM processing (inter-chromosomal discordant read-pairs)
            start_time=$SECONDS
            echo "### $(basename "$strain")/${sample} Bowtie2 SAM processing - inter-chromosomal discordant read pairs ###" >> "$log_file"
            echo "${experiment}" >> "$log_file"

            # Processing (set pairs in the same row, extract only those aligning both R1 and R2 completely)
            # Consider for further processing only reads pair aligning uniquely
            # [piped] read from the bowtie2 SAM directly (no alignment.sam copy). The
            # filtered result is small and is read twice, so it stays a file.
            tail -n +20 "${fq}.sam" | \
            awk -vOFS='\t' '{$11="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; print}' | \
            cut -f1-19 | \
            sed '$!N;s/\n/ /' | \
            awk '($3 != $7 && $7 != "=")' | \
            awk '($7 != "*")' | \
            awk '$6 == "75M"' | \
            awk '$25 == "75M"' | \
            awk '$13 == "XN:i:0"' | \
            awk '$32 == "XN:i:0"' | \
            sed 's/ /\t/g' | \
            awk '($19 != "YT:Z:UP")' | \
            awk '($38 != "YT:Z:UP")' > "${fq}_alignment_processed.sam"

            # [piped] R1 halves then R2 halves, header added, filtered by the SR-mapped
            # names and de-duplicated in one stream. The original concatenated the
            # samtools output twice before "awk '!seen[$0]++'", which only removed that
            # duplication, so the unique result is identical.
            { cat "${fq}_header.sam"
              { cat "${fq}_header.sam"
                cut -f1-19 "${fq}_alignment_processed.sam"
                cut -f20-38 "${fq}_alignment_processed.sam"
              } | samtools view -N "${fq}_r_R12_mapped.txt" -
            } | awk '!seen[$0]++' > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam"

            # Remove intermediate files
            rm "${fq}_alignment_processed.sam" "${fq}_r_R12_mapped.txt" "${fq}.sam" "${fq}_header.sam"

            # Prepare tsv file for R processing
            tail -n +20 "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam" | \
            # Extract relevant columns (1, 3, 4, 7, 8, 10) and convert to tsv format
            # Note: $1 is read name, $3 is reference sequence name, $4 is position, $7 is mate reference sequence name, $8 is mate position, $10 is sequence
            awk 'BEGIN{OFS="\t"} {print $1, $3, $4, $7, $8, $10}' > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.tsv"
            
            echo "" >> "$log_file"  # Adds a blank line
            elapsed_time=$(( SECONDS - start_time ))
            echo "Inter-chromosomal discordant read pairs processing completed in ${elapsed_time} seconds" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line 

            # Calculate total elapsed time
            elapsed_time_total=$((( SECONDS - start_time_total )/60))
            echo "Total processing completed in ${elapsed_time_total} minutes" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line 

          else
          echo "Warning: Decompressed files missing for sample $sample in $strain!" >> "$log_file"
          echo "" >> "$log_file"
          fi
        else
        echo "Warning: One or both compressed files for sample $sample are missing in $strain!" >> "$log_file"
        echo "" >> "$log_file"
        fi
      done
    done
done

# Calculate alignment and processing total elapsed time
elapsed_time_total_alignment=$((( SECONDS - start_time_total_alignment )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "ALIGNMENT AND SAM PROCESSING COMPLETED IN ${elapsed_time_total_alignment} MINUTES" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line


########### MODULE CORE - END ###########


##### MODULE MATs - START #####
############## MATa-MATa' analysis of 75bp reads ############

# MATa-MATa' analysis

echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### MATa-MATa' analysis, 75bp analysis ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 
start_time_total_MATs=$SECONDS

# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in T0 TSG TLG TLR; do
    for experiment in "${EXP_LIST[@]}"; do
      file_alignment="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam"
      if [[ -f "$file_alignment"  ]]; then
        echo "$(basename "$strain") ${sample} ${experiment}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing MATs for $file_alignment" >> "$log_file" 
        start_time_MATs=$SECONDS
        
        # Determine the number of cores available
        N_CPU=$(nproc)
        
        # Extract inter-chromosomal discordant read pairs for MATa-MATa' analysis
        # Remove header for processing
        sed -n -e '1,19p' "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam"

        # Extract reads and cat them to a new file
        awk ' ($3=="CHRIII" && $7=="CHRV" && $8>289178 && $8<290483) ' "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R1.sam"
        awk ' ($3=="CHRV" && $7=="CHRIII" && $4>289178 && $4<290483) ' "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.sam" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R2.sam"
        cat "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R1.sam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R2.sam" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.sam"
        #cat "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R1.sam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R2.sam" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.tsv"
        #rm  "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam"
        #rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R1.sam"
        #rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R2.sam"

        # convert SAM to sorted BAM (75nt)
        start_time=$SECONDS
        samtools sort -@ "$N_CPU" -o "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.sam"
        samtools index "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bai"
        elapsed_time=$(( SECONDS - start_time ))
        echo "Bam and Bai 75nt MAT reads completed in ${elapsed_time} seconds" >> "$log_file"

        # Generate BedGraphs coverage files (75nt)
        start_time=$SECONDS
        bamCoverage -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" -o "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bedgraph" -of bedgraph -p "$N_CPU"  -bs 1 --smoothLength 1
        elapsed_time=$(( SECONDS - start_time ))
        echo "Bedgraphs 75nt reads completed in ${elapsed_time} seconds" >> "$log_file"

        # Sorting Bedgraphs files (75nt)
        samtools faidx "${MYREF}/RG_PMV_v9.fasta" # Create index file
        cut -f1,2 "${MYREF}/RG_PMV_v9.fasta.fai" > "${MYREF}/chrom_order.txt"
        bedtools sort -i "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bedgraph" -g "${MYREF}/chrom_order.txt" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_sorted.bedgraph"
             
        
        # Generate TSV files (75nt)
        start_time=$SECONDS
        #bedtools coverage -a "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bedgraph" -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" -d > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.tsv"
        bedtools coverage -a "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_sorted.bedgraph" -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" -sorted -g "${MYREF}/chrom_order.txt" -d > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_sorted.tsv"
            
        # Keep only coverage for chromosomes CHRIII and CHRV
        awk '$1=="CHRIII" || $1=="CHRV"' "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_sorted.tsv" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_filtered.tsv"

        # Remove intermediate files   
        rm "${strain}/${sample}_${experiment}_"*.bai
        rm "${strain}/${sample}_${experiment}_"*.bam
        rm "${strain}/${sample}_${experiment}_"*.bedgraph
        rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_sorted.tsv"
        rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R1.sam"
        rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_R2.sam"
        rm "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam"

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time_MATs )))
        echo "Total MATs SAM processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file_alignment" >> "$log_file"
      echo "" >> "$log_file"
      fi
    done
  done
done


# R plots for MATa-MATa' coverage (from inter-chromosomal discordant read pairs)
# Start timer for R processing

R_SCRIPT="DISC_plot_MATa_cov.R"

# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Plotting MATa-MATa' coverage" >> "$log_file"
        start_time=$SECONDS
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
            echo "Skipping ${strain} and continuing..." >> "$log_file"
        continue
        fi

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done

# Calculate total elapsed time
elapsed_time_total=$((( SECONDS - start_time_total_MATs )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "TOTAL MATa-MATa' PROCESSING COMPLETED IN ${elapsed_time_total} MINUTES" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line    


############## MATa-MATa' analysis of 18bp reads ############

# MATa-MATa' analysis

echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### MATa-MATa' analysis, 18bp analysis ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 
start_time_total_MATs_18=$SECONDS

# Start timer for the alignment and processing section
start_time_total_alignment=$SECONDS

# Index RG for bowtie alignment
bowtie-build "${MYREF}/RG_PMV_v9.fasta" "${MYREF}/S_cerevisiae_indexed"

# Index RG for bowtie2 alignment
bowtie2-build "${MYREF}/RG_PMV_v9.fasta" "${MYREF}/S_cerevisiae_indexed"

# Loop through each directory inside MYWD

echo "### Alignments and SAM processing ###" >> "$log_file"
echo "Initiated at ${current_time}" >> "$log_file" 
echo "" >> "$log_file"  # Adds a blank line

# Define the pairs
declare -a R1_ids=("R1_1_1_1" "R1_1_1_2" "R1_1_2_1" "R1_1_2_2" "R1_2_1_1" "R2_1_1_1" "R2_1_1_2" "R2_1_2_1" "R2_1_2_2" "R2_2_1_1")
declare -a R2_ids=("R1_1_2_2" "R1_2_1_1" "R1_2_1_2" "R1_2_2_1" "R1_2_2_2" "R2_1_2_2" "R2_2_1_1" "R2_2_1_2" "R2_2_2_1" "R2_2_2_2")


for strain in "$MYWD"*/; do
  echo "Processing directory: $strain" >> "$log_file"

  # Loop through sample prefixes (TLG)
    for sample in TSG TLG; do
      for experiment in "${EXP_LIST[@]}"; do
        file1_gz="${strain}/${sample}_${experiment}_R1.fastq.gz"
        file2_gz="${strain}/${sample}_${experiment}_R2.fastq.gz"
       
        # Check if both compressed files exist
        if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
          echo "Processing files: $file1_gz, $file2_gz" >> "$log_file"
          echo "" >> "$log_file"  # Adds a blank line

          fq="${strain}/${sample}_${experiment}"

          # [piped] no gunzip -k: the R1/R2 FASTQ are streamed from the .gz files below
          if [[ -f "$file1_gz" && -f "$file2_gz" ]]; then
            # Determine the number of cores available
            N_CPU=$(nproc)

            # Start timers for total processing of a strain/sample/experiment
            start_time_total=$SECONDS

            # Trim reads using cutadapt
            echo "### $(basename "$strain")/${sample} ${experiment} ###" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line
            start_time=$SECONDS

            # [piped] One pass per mate writes the 8 final 18 nt slices directly
            # (R1_a_b_c = cut 75 -> cut 38 -> cut 19, as in the original). Replaces 28
            # cutadapt runs and the 12 intermediate 75/37 nt FASTQs. fastp reads each
            # slice as a regular file, as before.
            for mate in R1 R2; do
              gz="${strain}/${sample}_${experiment}_${mate}.fastq.gz"
              specs=()
              for p in 111 112 121 122 211 212 221 222; do
                specs+=("${p}=${fq}_${mate}_${p:0:1}_${p:1:1}_${p:2:1}.fastq")
              done
              "${DECOMP[@]}" "$gz" | fq_slice "${specs[@]}" &
            done
            wait

            for i in "${!R1_ids[@]}"; do
              R1_file=${R1_ids[$i]}
              R2_file=${R2_ids[$i]}
              pair="${fq}_${R1_file}_${R2_file}"
              echo "Processing files: "${strain}/${sample}_${experiment}_${R1_file}.fastq" "${strain}/${sample}_${experiment}_${R2_file}.fastq"" >> "$log_file"
              echo "" >> "$log_file"

               # Filter 18nt reads Q30
              fastp \
                -i "${strain}/${sample}_${experiment}_${R1_file}.fastq" \
                -o "${strain}/${sample}_${experiment}_${R1_file}_filtered.fastq" \
                -I "${strain}/${sample}_${experiment}_${R2_file}.fastq" \
                -O "${strain}/${sample}_${experiment}_${R2_file}_filtered.fastq" \
                -q 30 -u 0 -e 30 \
                --thread "$N_CPU" \
                --html "${strain}/${sample}_${experiment}_${R1_file}${R2_file}_18nt.fastq.html" \
                --json "${strain}/${sample}_${experiment}_${R1_file}${R2_file}_18nt.fastq.json"


              elapsed_time=$(( SECONDS - start_time ))
              echo "Cutting and Filtering Reads completed in ${elapsed_time} seconds" >> "$log_file"
              echo ""  >> "$log_file"

              # [piped] Bowtie SR mapping (-m1 -v0) streamed to the sorted list of mapped
              # read names (no SAM, mapped SAM, name list or unsorted list on disk)
              start_time=$SECONDS
              echo "### $(basename "$strain")/${sample} 18nt R1 ###" >> "$stats"
              echo "${experiment}" >> "$stats"
              echo "${R1_file}" >> "$stats"
              bowtie -p "$N_CPU" -m 1 -v 0 -S "${MYREF}/S_cerevisiae_indexed" "${fq}_${R1_file}_filtered.fastq" 2>> "$stats" | \
                samtools view -F 4 - | cut -f1 | sed 's/\/[12]//' | sort > "${fq}_${R1_file}_mapped_sorted.txt"
              echo "" >> "$stats"  # Adds a blank line
              elapsed_time=$(( SECONDS - start_time ))
              echo "Bowtie SR alignment 18nt for R1 reads completed in ${elapsed_time} seconds" >> "$log_file"
              echo "" >> "$log_file"

              start_time=$SECONDS
              echo "### $(basename "$strain")/${sample} 18nt R2 ###" >> "$stats"
              echo "${experiment}" >> "$stats"
              echo "${R2_file}" >> "$stats"
              bowtie -p "$N_CPU" -m 1 -v 0 -S "${MYREF}/S_cerevisiae_indexed" "${fq}_${R2_file}_filtered.fastq" 2>> "$stats" | \
                samtools view -F 4 - | cut -f1 | sed 's/\/[12]//' | sort > "${fq}_${R2_file}_mapped_sorted.txt"
              echo "" >> "$stats"  # Adds a blank line
              elapsed_time=$(( SECONDS - start_time ))
              echo "Bowtie SR alignment 18nt for R2 reads completed in ${elapsed_time} seconds" >> "$log_file"
              echo "" >> "$log_file"

              comm -12 "${fq}_${R1_file}_mapped_sorted.txt" "${fq}_${R2_file}_mapped_sorted.txt" > "${pair}_mapped.txt"
              rm "${fq}_${R1_file}_mapped_sorted.txt" "${fq}_${R2_file}_mapped_sorted.txt"

              # Bowtie2 mapping (in a forced -m1 -v0 mode)
              start_time=$SECONDS
              echo "### $(basename "$strain")/${sample} 18nt paired-end Bowtie2 alignment ###" >> "$stats"
              echo "${experiment}" >> "$stats"
              echo "${R1_file} and ${R2_file}" >> "$stats"
              bowtie2 \
                -p "$N_CPU" \
                --no-1mm-upfront \
                --score-min C,0,0 \
                -N 0 \
                --end-to-end \
                --fr \
                -x "${MYREF}/S_cerevisiae_indexed" \
                -1 "${strain}/${sample}_${experiment}_${R1_file}_filtered.fastq" \
                -2 "${strain}/${sample}_${experiment}_${R2_file}_filtered.fastq" \
                -S "${strain}/${sample}_${experiment}_${R1_file}_${R2_file}.sam" 2>> "$stats"

              echo "" >> "$stats"  # Adds a blank line
              elapsed_time=$(( SECONDS - start_time ))
              echo "Bowtie2 alignment 18nt for paired-end reads completed in ${elapsed_time} seconds" >> "$log_file"
              echo "" >> "$log_file"

              # Remove intermediate files
              rm "${strain}/${sample}_${experiment}_${R1_file}_filtered.fastq"
              rm "${strain}/${sample}_${experiment}_${R2_file}_filtered.fastq"

              # Bowtie2 SAM processing (inter-chromosomal discordant read-pairs)
              start_time=$SECONDS
              echo "### $(basename "$strain")/${sample} Bowtie2 SAM processing - inter-chromosomal discordant read pairs ###" >> "$log_file"
              echo "${experiment}" >> "$log_file"
              echo "${R1_file} and ${R2_file}" >> "$log_file"

              # Header kept as a small file (the original "sed -n 1,19p")
              head -n 19 "${pair}.sam" > "${pair}_header.sam"

              # Processing (set pairs in the same row, extract only those aligning both R1 and R2 completely)
              # Consider for further processing only reads pair aligning uniquely
              # [piped] read from the bowtie2 SAM directly (no alignment.sam copy)
              tail -n +20 "${pair}.sam" | \
              awk -vOFS='\t' '{$11="AAAAAAAAAAAAAAAAAA"; print}' | \
              cut -f1-19 | \
              sed '$!N;s/\n/ /' | \
              awk '($3 != $7 && $7 != "=")' | \
              awk '($7 != "*")' | \
              awk '$6 == "18M"' | \
              awk '$25 == "18M"' | \
              awk '$13 == "XN:i:0"' | \
              awk '$32 == "XN:i:0"' | \
              sed 's/ /\t/g' | \
              awk '($19 != "YT:Z:UP")' | \
              awk '($38 != "YT:Z:UP")' > "${pair}_alignment_processed.sam"

              # [piped] R1 halves then R2 halves + header -> filter by SR-mapped names ->
              # de-duplicate -> MATa-MATa' extraction, in one stream. Replaces 9
              # intermediate SAMs per pair. The extraction awk writes the header (first
              # 19 lines, as "sed -n 1,19p") and the MAT reads: CHRIII->CHRV reads first,
              # then CHRV->CHRIII reads (buffered), i.e. the original "cat R1_MAT R2_MAT".
              { cat "${pair}_header.sam"
                { cat "${pair}_header.sam"
                  cut -f1-19 "${pair}_alignment_processed.sam"
                  cut -f20-38 "${pair}_alignment_processed.sam"
                } | samtools view -N "${pair}_mapped.txt" -
              } | awk '!seen[$0]++' | \
              awk -v hdr="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam" \
                  -v out="${pair}_inter_discordant_pairs_unique_MAT_r18.sam" '
                BEGIN { printf "" > hdr; printf "" > out }
                NR <= 19 { print > hdr }
                ($3=="CHRIII" && $7=="CHRV" && $8>289178 && $8<290483) { print > out }
                ($3=="CHRV" && $7=="CHRIII" && $4>289178 && $4<290483) { b[++n] = $0 }
                END { for (i = 1; i <= n; i++) print b[i] > out }'

              # Remove intermediate files
              rm "${pair}_alignment_processed.sam" "${pair}_mapped.txt" "${pair}.sam" "${pair}_header.sam"

              echo "" >> "$log_file"  # Adds a blank line
              elapsed_time=$(( SECONDS - start_time ))
              echo "Inter-chromosomal discordant read pairs processing completed in ${elapsed_time} seconds" >> "$log_file"
              echo "" >> "$log_file"  # Adds a blank line

            done
            # Calculate elapsed time
            elapsed_time=$((( SECONDS - start_time_total )))
            echo "Total MATs processing completed in ${elapsed_time} seconds" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line

            # Remove the 18 nt slice FASTQs
            for mate in R1 R2; do
              for p in 111 112 121 122 211 212 221 222; do
                rm -f "${fq}_${mate}_${p:0:1}_${p:1:1}_${p:2:1}.fastq"
              done
            done

            # [piped] header + every pair's MAT reads (same glob order) in one write,
            # without the intermediate MAT_noheader_r18.sam
            cat "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_header.sam" \
                "${strain}/${sample}_${experiment}_"*"_inter_discordant_pairs_unique_MAT_r18.sam" \
              > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.sam"
            # Original cleanup, kept verbatim: besides the per-pair files (no longer
            # created), "*_pairs_unique.sam" also deletes the 75 nt
            # ${sample}_${experiment}_inter_discordant_pairs_unique.sam, as before.
            rm -f ${strain}/${sample}_${experiment}_*_pairs_unique.sam
            rm -f ${strain}/${sample}_${experiment}_*_pairs_unique_R1.sam
            rm -f ${strain}/${sample}_${experiment}_*_pairs_unique_R2.sam
            rm -f ${strain}/${sample}_${experiment}_*_unique_header.sam
            rm -f ${strain}/${sample}_${experiment}_*_pairs_unique_R2.sam
            rm -f ${strain}/${sample}_${experiment}_*_unique_MAT_noheader_r18.sam
            ###### remove this too
            rm -f ${strain}/${sample}_${experiment}_R*_MAT.sam
            rm -f ${strain}/${sample}_${experiment}_*_header.sam
            rm -f ${strain}/${sample}_${experiment}_*_noheader_r18.sam
            rm -f ${strain}/${sample}_${experiment}_R*_unique.sam

          else
          echo "Warning: Decompressed files missing for sample $sample in $strain!" >> "$log_file"
          echo "" >> "$log_file"
          fi
        else
        echo "Warning: One or both compressed files for sample $sample are missing in $strain!" >> "$log_file"
        echo "" >> "$log_file"
        fi
        
      done
    done
done



#Coverage for MATa-MATa' discordant read pairs

# MATa-MATa' analysis


# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG; do
    for experiment in "${EXP_LIST[@]}"; do
      file_alignment="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.sam"
      if [[ -f "$file_alignment"  ]]; then
        echo "$(basename "$strain") ${sample} ${experiment}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing MATs for $file_alignment" >> "$log_file" 
        start_time_MATs=$SECONDS
        
        # Determine the number of cores available
        N_CPU=$(nproc)
        
        # convert SAM to sorted BAM (75nt)
        start_time=$SECONDS
        samtools sort -@ "$N_CPU" -o "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.sam"
        samtools index "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bam" "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bai"
        elapsed_time=$(( SECONDS - start_time ))
        echo "Bam and Bai 18nt MAT reads completed in ${elapsed_time} seconds" >> "$log_file"

        # Generate BedGraphs coverage files (75nt)
        start_time=$SECONDS
        bamCoverage -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bam" -o "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bedgraph" -of bedgraph -p "$N_CPU"  -bs 1 --smoothLength 1
        elapsed_time=$(( SECONDS - start_time ))
        echo "Bedgraphs 18nt reads completed in ${elapsed_time} seconds" >> "$log_file"

        # Sorting Bedgraphs files (75nt)
        samtools faidx "${MYREF}/RG_PMV_v9.fasta" # Create index file
        cut -f1,2 "${MYREF}/RG_PMV_v9.fasta.fai" > "${MYREF}/chrom_order.txt"
        bedtools sort -i "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bedgraph" -g "${MYREF}/chrom_order.txt" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18_sorted.bedgraph"
             
        
        # Generate TSV files (75nt)
        start_time=$SECONDS
        #bedtools coverage -a "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bedgraph" -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.bam" -d > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT.tsv"
        bedtools coverage -a "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18_sorted.bedgraph" -b "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18.bam" -sorted -g "${MYREF}/chrom_order.txt" -d > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18_sorted.tsv"
            
        # Keep only coverage for chromosomes CHRIII and CHRV
        awk '$1=="CHRIII" || $1=="CHRV"' "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18_sorted.tsv" > "${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_MAT_complete_r18_filtered.tsv"

        # Remove intermediate files   
        rm "${strain}/${sample}_${experiment}_"*.bai
        rm "${strain}/${sample}_${experiment}_"*.bam
        rm "${strain}/${sample}_${experiment}_"*.bedgraph
        rm -f ${strain}/${sample}_${experiment}_*_MAT_complete_r18_sorted.tsv
       

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time_MATs )))
        echo "Total MATs processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file_alignment" >> "$log_file"
      echo "" >> "$log_file"
      fi
    done
  done
done

# elapsed_time_total=$(( SECONDS - start_time ))
# echo "Inter-chromosomal discordant read pairs R processing completed in ${elapsed_time} seconds" >> "$log_file"

# Calculate total elapsed time
elapsed_time_total=$((( SECONDS - start_time_total_MATs_18 )/60))
echo "Total MATs processing completed in ${elapsed_time_total} minutes" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line  



# R plots for MATa-MATa' discordant read pairs

R_SCRIPT="DISC_plot_MATa_cov_r18.R"

# Start timer for R processing
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Plotting MATa-MATa' coverage" >> "$log_file"
        start_time=$SECONDS
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
            echo "Skipping ${strain} and continuing..." >> "$log_file"
        continue
        fi


      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done

# Calculate alignment and processing total elapsed time
elapsed_time_total_r18=$((( SECONDS - start_time_total_MATs_18 )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "ALIGNMENT AND SAM PROCESSING COMPLETED IN ${elapsed_time_total_r18} MINUTES" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line


##### MODULE MATs - DATA ORGANIZATION #####
# Move data into subfolders
for subdir in "$MYWD"*/; do
  mkdir -p "${subdir}/Discordant_alignments" "${subdir}/Discordant_MAT_analysis"
  mkdir -p "${subdir}/Discordant_alignments/FASTQP_75nt" "${subdir}/Discordant_alignments/FASTQP_18nt" "${subdir}/Discordant_alignments/SAM_75nt" "${subdir}/Discordant_alignments/SAM_18nt" 
  mkdir -p "${subdir}/Discordant_MAT_analysis/MAT_Data_75nt" "${subdir}/Discordant_MAT_analysis/MAT_Plots_75nt" "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt" "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt" 
  
  mv "${subdir}"/*_75nt.fastq.html "${subdir}/Discordant_alignments/FASTQP_75nt/" 2> /dev/null
  mv "${subdir}"/*_75nt.fastq.json "${subdir}/Discordant_alignments/FASTQP_75nt/" 2> /dev/null
  mv "${subdir}"/*_18nt.fastq.html "${subdir}/Discordant_alignments/FASTQP_18nt/" 2> /dev/null
  mv "${subdir}"/*_18nt.fastq.json "${subdir}/Discordant_alignments/FASTQP_18nt/" 2> /dev/null

  mv "${subdir}"/*_inter_discordant_pairs_unique_MAT_filtered.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_MAT.sam "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null

  mv "${subdir}"/*_inter_discordant_pairs_unique_MAT_complete_r18_filtered.tsv "${subdir}/Discordant_alignments/SAM_18nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_MAT_complete_r18.sam "${subdir}/Discordant_alignments/SAM_18nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_MAT_r18.sam "${subdir}/Discordant_alignments/SAM_18nt/" 2> /dev/null

  mv "${subdir}"/*_MAT_coverage_r18_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_count_HOinc_reads_summary_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_count_total_reads_summary_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_GC_CO_summary_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_groups_summary_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_r18_processed_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_18nt/" 2> /dev/null

  mv "${subdir}"/MATa_MATa_GC_CO_r18_*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt/" 2> /dev/null
  mv "${subdir}"/MATa_MATa_GC_group_r18_*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt/" 2> /dev/null
  mv "${subdir}"/MATa_MATa_pairs_r18_*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt/" 2> /dev/null
  mv "${subdir}"/MATa_MATa_poly_*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt/" 2> /dev/null
  mv "${subdir}"/Average_MATa_MATa_coverage_r18*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_18nt/" 2> /dev/null

  mv "${subdir}"/*_MAT_coverage_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_MAT_pairs_summary_df.tsv "${subdir}/Discordant_MAT_analysis/MAT_Data_75nt/" 2> /dev/null

  mv "${subdir}"/Average_MATa_MATa_coverage_*.svg "${subdir}/Discordant_MAT_analysis/MAT_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/MATa_MATa_pairs_*.svg  "${subdir}/Discordant_MAT_analysis/MAT_Plots_75nt/" 2> /dev/null

done
  

##### MODULE MATs - END #####


##### MODULE WHOLE GENOME ANALYSIS - START #####

############## R processing of 75bp reads ############

# R processing for inter-chromosomal discordant read pairs
# After processing all samples, all discordant pairs are associated to their genomic features
# and only discordant pairs not present at T0 are kept for further BLAST validation.

start_time_total_R=$SECONDS


echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### R processing for new discordant read pairs, 75bp analysis ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

# Loop inside each subdirectory of MYWD
R_SCRIPT="DISC_process_inter_discordant_pairs_10kb_fast.R"

for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for experiment in "${EXP_LIST[@]}"; do
      echo "### $(basename "$strain")/${sample} ${experiment} ###" >> "$log_file"
      echo "" >> "$log_file"  # Adds a blank line
      strain="${strain}"
      sample="${sample}"
      experiment="${experiment}"
      category_path="${CATEGORY_PATH}"
      file1="${strain}/T0_${experiment}_inter_discordant_pairs_unique.tsv"
      file2="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique.tsv"
      if [[ -f "$file1" && -f "$file2" ]]; then
        # echo "$(basename "$strain") ${sample} ${experiment}" >> "$log_file"
        # echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant tsv files for $file1 and $file2" >> "$log_file"
        start_time=$SECONDS
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$strain" "$sample" "$experiment" "$file1" "$file2" "$category_path"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${experiment}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${experiment} and continuing..." >> "$log_file"
        continue
        fi

      
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file1 or $file2" >> "$log_file"
      echo "" >> "$log_file"
      fi
    done
  done
done


# Calculate total elapsed time
elapsed_time_total_R=$((( SECONDS - start_time_total_R )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "TOTAL INTER-CHROMOSOMAL DISCORDANT READ PAIRS R PROCESSING COMPLETED IN ${elapsed_time_total_R} MINUTES" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
          


############## BLAST processing of 75bp reads ############

# BLAST processing for inter-chromosomal discordant read pairs
# For each inter-chromosomal discordant read pair not present in T0,
# a BLAST cross-validation is performed to discard false positives arising from 
# PCR mutation after library amplification or sequencing errors.

start_time_total_blast=$SECONDS

# [piped] Build each chromosome BLAST database ONCE for the whole run. The original
# rebuilt the database of the read's chromosome for every single read and deleted it
# right after (~0.25 s per read, ~50% of this section). Same names and location as
# before (MYBLAST/<CHR>.n*), because blastn prints the -db name in its "# Database:"
# line, which is kept in the processed results.
(
  cd "$MYBLAST" || exit
  for ref_fasta in *_v9.fasta; do
    [[ -f "$ref_fasta" ]] || continue
    makeblastdb -in "$ref_fasta" -dbtype nucl -out "${ref_fasta%_v9.fasta}"
  done
)

# Loop through each directory inside MYWD

echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### BLAST cross-validation for inter-chromosomal discordant read pairs ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

for strain in "$MYWD"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
  start_time_total=$SECONDS

  # Loop through sample prefixes (TSG, TLG, TLR)
    for sample in TSG TLG TLR; do
      for experiment in "${EXP_LIST[@]}"; do
        for pos in chromosome; do
          input_file="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_blast_${pos}.tsv"
          if [[ -f "$input_file" ]]; then
            echo "### $(basename "$strain")/${sample} ${experiment} ###" >> "$log_file"
            echo "### ${pos} ###" >> "$log_file"
            echo "" >> "$log_file"
            echo "Processing: Input file ${input_file} " >> "$log_file"
            echo "" >> "$log_file"
            cp "$input_file" "$MYBLAST" 
            # Start timer for blast processing
            start_time=$SECONDS
            file_name=$(basename "${input_file%.tsv}")
            output_dir="${MYBLAST}/1.blast_results"
            # Create output directory if it doesn't exist
            mkdir -p "$output_dir"
            (
                cd "$MYBLAST" || exit

                    # [piped] One processed line per read is appended straight to the
                    # combined file. The original wrote each read's BLAST table, three
                    # grep/tr/tr temp files and a per-read result file, then gathered the
                    # 1000s of per-read files with find (in arbitrary directory order).
                    # Rows are identical; their order is now the input order.
                    combined_file="${output_dir}/${file_name}_combined_results.tsv"
                    : > "$combined_file"

                    while IFS=$'\t' read -r nombre_read query_A reference_B || [ -n "$nombre_read" ]; do # -n to read last row also if it doesn't end with a newline
                        if [[ -z "$nombre_read" || -z "$query_A" || -z "$reference_B" ]]; then
                            echo "WARNING: Malformed line in $file_name.tsv - skipping" >> "$log_file"
                            continue
                        fi
                        reference_ext="_v9.fasta"
                        reference_name=$reference_B
                        reference_file_path="${reference_name}${reference_ext}"
                        read_name=$nombre_read
                        query=$query_A

                        if [[ ! -f "$reference_file_path" ]]; then
                            echo "WARNING: Reference file $reference_file_path not found - skipping" >> "$log_file"
                            continue
                        fi
                        # [piped] database built once per run (see blast_db_build above);
                        # query on stdin, BLAST table straight into the processing pipe
                        echo -e ">$read_name\n$query" | \
                        blastn -db "$reference_name" -query - -dust no -outfmt "7 qseqid sseqid mismatch length sstart send qseq sseq" | \
                        grep -E 'Query|Database|hits|E25|V35' | \
                        tr ' ' '\t' | \
                        tr '\n ' '\t' | \
                        cut -f3,6,8,13- >> "$combined_file"

                    done < "$input_file"

            cp "$combined_file" "${strain}"
            rm -r "${output_dir}"
            rm "${MYBLAST}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_blast_${pos}.tsv"
            
            echo "" >> "$log_file"
            echo "Blast processing for $input_file completed and results saved to ${strain}" >> "$log_file"
            echo "" >> "$log_file"

            )  # End subshell to change directory 
            
            
            # End timer for blast processing
            elapsed_time=$((( SECONDS - start_time )/60))
            echo "Blast processing completed in ${elapsed_time} minutes" >> "$log_file"
            echo "" >> "$log_file"  # Adds a blank line
          else
            echo "Warning: Input file $input_file does not exist!" >> "$log_file"
          fi

        done
      done
    done
  # End timer for the current strain processing
  elapsed_time_total=$((( SECONDS - start_time_total )/60))
  echo "" >> "$log_file"
  echo "Processing for strain $(basename "$strain") completed in ${elapsed_time_total} minutes" >> "$log_file"
  echo "" >> "$log_file"  # Adds a blank line
done  

# [piped] Remove the BLAST databases (same file types the original deleted per read)
(
  cd "$MYBLAST" || exit
  rm -f -- *.nhr *.nin *.nsq *.ndb *.njs *.not *.ntf *.nto
)

# Calculate total elapsed time in BLAST processing
elapsed_time_total_blast=$((( SECONDS - start_time_total_blast )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "BLAST CROSS-VALIDATION FOR INTER-CHROMOSOMAL DISCORDANT READ PAIRS COMPLETED IN ${elapsed_time_total_blast} MINUTES" >> "$log_file"
echo "" >> "$log_file"


############## R processing of 75bp VALID reads ############

# R processing for inter-chromosomal discordant read pairs after BLAST cross-validation
# After processing all samples, only valid discordant reads are kept for further analysis.
# With this script, valid reads are identified and saved in a df.
# Start timer for R processing

start_time_total_R_valid=$SECONDS

echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### R processing for new valid discordant read pairs, 75bp analysis ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

# Loop inside each subdirectory of MYWD
R_SCRIPT="DISC_validate_inter_discordant_pairs_10kb.R"

for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for experiment in "${EXP_LIST[@]}"; do

      strain="${strain}"
      sample="${sample}"
      experiment="${experiment}"
      category_path="${CATEGORY_PATH}"

      file1="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv"
      
      fileorig="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed.tsv"
      filecontrol="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_control.tsv"

      if [[ -f "$file1" && -f "$fileorig" && -f "$filecontrol" ]]; then
        echo "$(basename "$strain") ${sample} ${experiment}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant tsv files for $file1 and $fileorig" >> "$log_file"
        start_time=$SECONDS
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$strain" "$sample" "$experiment" "$file1" "$fileorig" "$filecontrol" "$category_path"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${experiment}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${experiment} and continuing..." >> "$log_file"
        continue
        fi

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file1 or $fileorig" >> "$log_file"
      echo "" >> "$log_file"
      fi
    done
  done
done

        

# R processing for inter-chromosomal discordant read pairs global distribution
# Start timer for R processing
R_SCRIPT="DISC_calculate_global_distribution.R"

# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    root_dir="${strain}"
    strain="${strain}"
    sample="${sample}"
    
    echo "$(basename "$strain") ${sample}" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line
    echo "Processing inter-discordant global distribution tsv files for ${strain} , ${sample}" >> "$log_file"
    start_time=$SECONDS

    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample}" >> "$log_file"
            echo "Skipping ${strain} ${sample} and continuing..." >> "$log_file"
        continue
    fi   
    # Calculate elapsed time
      elapsed_time=$((( SECONDS - start_time )/60))
      echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
      echo "" >> "$log_file"  # Adds a blank line 

  done
done



# R processing for inter-chromosomal discordant read pairs category distribution


# Loop inside each subdirectory of MYWD
R_SCRIPT="DISC_calculate_category_distribution.R"

for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    root_dir="${strain}"
    strain="${strain}"
    sample="${sample}"
      
    echo "$(basename "$strain") ${sample}" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line
    echo "Processing inter-discordant category distribution tsv files for ${strain} , ${sample}" >> "$log_file"
    start_time=$SECONDS

    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample}" >> "$log_file"
            echo "Skipping ${strain} ${sample} and continuing..." >> "$log_file"
        continue
    fi   
    
  # Calculate elapsed time
    elapsed_time=$((( SECONDS - start_time )/60))
    echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line 
      

  done
done



# R processing for inter-chromosomal discordant read pairs matrix (merge all experiments)

R_SCRIPT="DISC_calculate_matrix_all_exp.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    root_dir="${strain}"
    strain="${strain}"
    sample="${sample}"
    category_path="${CATEGORY_PATH}"
      
    echo "$(basename "$strain") ${sample}" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line
    echo "Processing inter-discordant valid reads tsv files for ${strain} , ${sample}" >> "$log_file"
    start_time=$SECONDS
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$category_path"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample}" >> "$log_file"
            echo "Skipping ${strain} ${sample} and continuing..." >> "$log_file"
        continue
    fi  
        
    # Calculate elapsed time
    elapsed_time=$((( SECONDS - start_time )))
    echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line 
     

  done
done



##### Calculate % of discordant reads in each Category_A and plot bar plot with SD
R_SCRIPT="DISC_calculate_distrib_categoryA.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      strain="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing discordant reads distribution by Category_A" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
            echo "Skipping ${strain} and continuing..." >> "$log_file"
          continue
        fi
        

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done

# R processing for inter-chromosomal discordant read pairs hotspots (all experiments)
R_SCRIPT="DISC_identify_hotspots.R"
# Loop inside each subdirectory of MYWD

for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
      freqs="${strain}/${sample}_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv"
      root_dir="${strain}"
      strain="${strain}"
      sample="${sample}"

      
      if [[ -f "$freqs" ]]; then
        echo "$(basename "$strain") ${sample}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant valid reads tsv files for ${strain} , ${sample}" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$freqs"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} , ${sample}" >> "$log_file"
            echo "Skipping ${strain} ${sample} and continuing..." >> "$log_file"
          continue
        fi


      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample: $file1 or $file2 or $file3" >> "$log_file"
      echo "" >> "$log_file"
      fi

  done
done


##### R processing for recombination rate

R_SCRIPT="DISC_calculate_recombination_rate.R"
# Get all processed_control.tsv files recursively in root folder

echo "Calculating recombination rate" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "Calculating HO-pGAL read count" >> "$log_file"
start_time=$SECONDS

echo "Processing directory: ${MYWD}" >> "$log_file"
      root_dir="${MYWD}"
      reference_strain="${STRAIN_REF}"
      
        echo "" >> "$log_file"  # Adds a blank line
        echo "Finding all processed_control.tsv files..." >> "$log_file"
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$reference_strain"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${root_dir}" >> "$log_file"
            echo "Skipping ${root_dir} and continuing..." >> "$log_file"
        fi
        
        



# Loop inside each subdirectory of MYWD
R_SCRIPT="DISC_calculate_recombination_rate_conc_disc.R"
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      wd_dir="${MYWD}"
      strain="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Calculating recombination rate" >> "$log_file"
        #start_time=$SECONDS
        # Export variables for R access
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${root_dir} ${strain}" >> "$log_file"
            echo "Skipping ${root_dir} ${strain} and continuing..." >> "$log_file"
          continue
        fi

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done


### R script to calculate Category_A read number and plot all timepoints###
R_SCRIPT="DISC_calculate_distrib_read_number_categoryA.R"
# Loop inside each subdirectory of MYWD 
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      wd_dir="${MYWD}"
      strain="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Calculating Category_A read number" >> "$log_file"
        start_time=$SECONDS
       if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${root_dir} ${strain}" >> "$log_file"
            echo "Skipping ${root_dir} ${strain} and continuing..." >> "$log_file"
          continue
        fi
        
done

# R processing for error rate
R_SCRIPT="DISC_calculate_error_rate.R"
# Get all processed_valid_error_rate.tsv files recursively in root folder

echo "Calculating error rate" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
start_time=$SECONDS

# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  echo "Processing directory: ${strain}" >> "$log_file"
      root_dir="${strain}"
      wd_dir="${MYWD}"
      strain="${strain}"
      #echo "$root_dir"
        echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        
        start_time=$SECONDS
        # Export variables for R access
        
      #echo "$root_dir"
        #echo "$(basename "$strain")" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Finding all processed_valid_error_rate.tsv files..." >> "$log_file"
        
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${root_dir} ${strain}" >> "$log_file"
            echo "Skipping ${root_dir} ${strain} and continuing..." >> "$log_file"
          continue
        fi
    

        
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 

done



### Added BLAST levels of validation
############## R processing of 75bp VALID reads ############

# R processing for inter-chromosomal discordant read pairs after BLAST cross-validation
# After processing all samples, only valid discordant reads are kept for further analysis.
# With this script, valid reads are identified and saved in a df.
# We have defined 5 different levels of validation for each read pair:
# Start timer for R processing

start_time_total_R_valid=$SECONDS
R_SCRIPT="DISC_validate_inter_discordant_pairs_levels_1feature_10kb.R"

echo "================================" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line
echo "### R processing for new valid discordant read pairs, 75bp analysis, different BLAST levels ###" >> "$log_file"
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Initiated at ${current_time}" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 

# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for experiment in "${EXP_LIST[@]}"; do
      file1="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv"
      
      fileorig="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed.tsv"
      filecontrol="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_control.tsv"
      if [[ -f "$file1" && -f "$fileorig" && -f "$filecontrol" ]]; then
        echo "$(basename "$strain") ${sample} ${experiment}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant tsv files for $file1 and $fileorig" >> "$log_file"
        start_time=$SECONDS
        strain="${strain}"
        sample="${sample}"
        experiment="${experiment}"
        category_path="${CATEGORY_PATH}"
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$strain" "$sample" "$experiment" "$category_path" "$file1" "$fileorig" "$filecontrol"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${experiment}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${experiment} and continuing..." >> "$log_file"
          continue
        fi

      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing for different blast options completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file1 or $fileorig" >> "$log_file"
      echo "" >> "$log_file"
      fi
    done
  done
done


# R script to repeat discordant analysis with all blast options
R_SCRIPT="DISC_process_valid_inter_discordant_pairs_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for experiment in "${EXP_LIST[@]}"; do
        for blast_option in option1 option2 option3 option4 option5 ; do
        file1="${strain}/${sample}_${experiment}_${blast_option}_inter_discordant_pairs_unique_processed_valid.tsv"
        file2="${strain}/${sample}_${experiment}_${blast_option}_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"
        filecontrol="${strain}/${sample}_${experiment}_inter_discordant_pairs_unique_processed_control.tsv"
      if [[ -f "$file1" && -f "$file2" ]]; then
        echo "$(basename "$strain") ${sample} ${experiment} ${blast_option}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant tsv files for $file1 , $file2" >> "$log_file"
        start_time=$SECONDS
        strain="${strain}"
        sample="${sample}"
        experiment="${experiment}"
        category_path="${CATEGORY_PATH}"
        blast_option="${blast_option}"
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$strain" "$sample" "$experiment" "$blast_option" "$category_path" "$file1" "$file2" "$filecontrol"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${experiment} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${experiment} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi
        
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      else
      echo "Missing files for $strain/$sample/$experiment: $file1 or $file2 or $filecontrol" >> "$log_file"
      echo "" >> "$log_file"
      fi
      done
    done
  done
done


# R processing for inter-chromosomal discordant read pairs global distribution
# Start timer for R processing
R_SCRIPT="DISC_calculate_global_distribution_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for blast_option in option1 option2 option3 option4 option5 ; do
      root_dir="${strain}"
      strain="${strain}"
      sample="${sample}"
      blast_option="${blast_option}"
        echo "$(basename "$strain") ${sample}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant global distribution tsv files for $strain , $sample and $blast_option" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$blast_option"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
      
    done
  done
done



# R processing for inter-chromosomal discordant read pairs category distribution

R_SCRIPT="DISC_calculate_category_distribution_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for blast_option in option1 option2 option3 option4 option5 ; do

      root_dir="${strain}"
      strain="${strain}"
      sample="${sample}"
      blast_option="${blast_option}"
        echo "$(basename "$strain") ${sample}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant category distribution tsv files for $strain $sample $blast_option" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$blast_option"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )/60))
        echo "Total R processing completed in ${elapsed_time} minutes" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
   
    done
  done
done



# R processing for inter-chromosomal discordant read pairs matrix (merge all experiments)

R_SCRIPT="DISC_calculate_matrix_all_exp_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for blast_option in option1 option2 option3 option4 option5 ; do
        root_dir="${strain}"
        strain="${strain}"
        sample="${sample}"
        blast_option="${blast_option}"
        category_path="${CATEGORY_PATH}"
        echo "$(basename "$strain") ${sample}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant valid reads tsv files for $strain , $sample and $blast_option" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$blast_option" "$category_path"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi    
        # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
        
    done
  done
done



##### Calculate % of discordant reads in each Category_A and plot bar plot with SD
R_SCRIPT="DISC_calculate_distrib_categoryA_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
    for blast_option in option1 option2 option3 option4 option5 ; do
        echo "Processing directory: ${strain}" >> "$log_file"
        root_dir="${strain}"
        strain="${strain}"
        blast_option="${blast_option}"
        #echo "$root_dir"
        echo "$(basename "$strain") ${blast_option}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing discordant reads distribution by Category_A" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$blast_option"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi     
      # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
    done
done

# R processing for inter-chromosomal discordant read pairs hotspots (all experiments)
R_SCRIPT="DISC_identify_hotspots_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    for blast_option in option1 option2 option3 option4 option5 ; do
        root_dir="${strain}"
        strain="${strain}"
        sample="${sample}"
        blast_option="${blast_option}"    
        freqs="${strain}/${sample}_${blast_option}_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv"

        echo "$(basename "$strain") ${sample} ${blast_option}" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line
        echo "Processing inter-discordant valid reads tsv files" >> "$log_file"
        start_time=$SECONDS
        if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample" "$blast_option" "$freqs"; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample} ${blast_option}" >> "$log_file"
            echo "Skipping ${strain} ${sample} ${blast_option} and continuing..." >> "$log_file"
          continue
        fi 
        # Calculate elapsed time
        elapsed_time=$((( SECONDS - start_time )))
        echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
        echo "" >> "$log_file"  # Adds a blank line 
    done
  done
done


##### Hotspots distribution TSG, TLG and TLR
R_SCRIPT="DISC_calculate_hotspots_distribution_levels.R"
root_dir="${MYWD}"
echo "" >> "$log_file"  # Adds a blank line
echo "Plotting hotspots distribution..." >> "$log_file"
start_time=$SECONDS

if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir"; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script" >> "$log_file"
    echo "Skipping and continuing..." >> "$log_file"
fi 
elapsed_time=$((( SECONDS - start_time )))
echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line 


##### Radar plot TSG,TLG,TLR same plot
R_SCRIPT="DISC_calculate_radar_perc_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for blast_option in option1 option2 option3 option4 option5 ; do
    echo "Processing directory: ${strain}" >> "$log_file"
    root_dir="${strain}"
    wd_dir="${MYWD}"
    strain="${strain}"
    blast_option="${blast_option}"

    echo "$(basename "$strain")" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line
    echo "Calculating radar plot for $blast_option" >> "$log_file"
    start_time=$SECONDS
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir" "$blast_option"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${blast_option}" >> "$log_file"
        echo "Skipping ${strain} ${blast_option} and continuing..." >> "$log_file"
        continue
    fi 
    elapsed_time=$((( SECONDS - start_time )))
    echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line   
  done
done


##### Radar plots TSG, TLG and TLR (read_number)
R_SCRIPT="DISC_calculate_radar_number_levels.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
    echo "Processing directory: ${strain}" >> "$log_file"
    root_dir="${strain}"
    wd_dir="${MYWD}"
    strain="${strain}"
    echo "$(basename "$strain")" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line
    echo "Plotting radar plots TSG, TLG and TLR - all blast options..." >> "$log_file"
    start_time=$SECONDS
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$wd_dir"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain}" >> "$log_file"
        echo "Skipping ${strain} and continuing..." >> "$log_file"
        continue
    fi                
    # Calculate elapsed time
    elapsed_time=$((( SECONDS - start_time )))
    echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line 
done

# Discordant networks
R_SCRIPT="DISC_calculate_discordant_network.R"
# Loop inside each subdirectory of MYWD
for strain in "${MYWD}"*/; do
  for sample in TSG TLG TLR; do
    echo "Processing network for sample: ${sample}" >> "$log_file"
    root_dir="${strain}"
    strain="${strain}"
    sample="${sample}"

    echo "" >> "$log_file"  # Adds a blank line
    echo "Calculating discordant reads network" >> "$log_file"
    start_time=$SECONDS
    if ! Rscript "${MYSCRIPTS}/R_files/${R_SCRIPT}" "$root_dir" "$strain" "$sample"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') ERROR running R script for ${strain} ${sample}" >> "$log_file"
        echo "Skipping ${strain} ${sample} and continuing..." >> "$log_file"
        continue
    fi                
    # Calculate elapsed time
    elapsed_time=$((( SECONDS - start_time )))
    echo "Total R processing completed in ${elapsed_time} seconds" >> "$log_file"
    echo "" >> "$log_file"  # Adds a blank line 
          
  done
done

########
# Calculate total elapsed time
elapsed_time_total=$((( SECONDS - start_time_total_R_valid )/60))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "TOTAL VALID INTER-CHROMOSOMAL DISCORDANT READ PAIRS R PROCESSING COMPLETED IN ${elapsed_time_total} MINUTES" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line   

##### MODULE WHOLE GENOME ANALYSIS - DATA ORGANIZATION #####
# Move data into subfolders
for subdir in "$MYWD"*/; do
  mkdir -p "${subdir}/Discordant_global_analysis" 
  mkdir -p "${subdir}/Discordant_global_analysis/Global_Data_75nt" "${subdir}/Discordant_global_analysis/Global_Plots_75nt" 
  
  
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_chromosome.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  # mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_next.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  # mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_prev.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_control.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique.sam "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  # mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_next_combined_results.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null
  # mv "${subdir}"/*_inter_discordant_pairs_unique_processed_blast_prev_combined_results.tsv "${subdir}/Discordant_alignments/SAM_75nt/" 2> /dev/null


  mv "${subdir}"/*_recombination_df.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_recombination_summary.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_valid_counts_summary_ordered.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_concordant_pairs_unique_row_count.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid_category_distribution.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid_counts.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid_discordant_mat*.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid_error_rate.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid_global_distribution.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed_valid.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique_processed.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_inter_discordant_pairs_unique.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_category_A_valid_counts_summary_df.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_hotspots.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  
  mv "${subdir}"/Discordant_reads_distribution_plot_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Inter_chromosomal_discordant_hotspots_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Inter_chromosomal_discordant_matrix_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Inter_chromosomal_discordant_pairs_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Recombination_rate_plot_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Recombination_rate_ratio_plot_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Category_A_read_number_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/Radar_plot_*.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null

  mv "${subdir}"/*_Network_plot.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/*_Network_quantification_plot.svg  "${subdir}/Discordant_global_analysis/Global_Plots_75nt/" 2> /dev/null
  mv "${subdir}"/*_hotspots_matrix.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_hotspots_processed_df.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  mv "${subdir}"/*_Network_summary_positions_strain_sep.tsv  "${subdir}/Discordant_global_analysis/Global_Data_75nt/" 2> /dev/null
  

done


#### MODULE WHOLE GENOME ANALYSIS - END #####


# Calculate full script total elapsed time
elapsed_time_total_2=$((( SECONDS - start_time_total_script_2 )/3600))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "DISCORDANT ANALYSIS - COMPLETED IN ${elapsed_time_total_2} HOURS" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line

elapsed_time_total=$((( SECONDS - start_time_total_script )/3600))
echo "================================" >> "$log_file"
echo "" >> "$log_file"
echo "SCRIPT COMPLETED IN ${elapsed_time_total} HOURS" >> "$log_file"
echo "" >> "$log_file"  # Adds a blank line


##### REPORT GENERATION #####

# Generate report file
current_time=$(date "+%d-%m-%Y %H:%M:%S")
echo "Generating report file at: \"$current_time\"" 

# Paths
RMD_FILE="Di-GRAPH_report.Rmd"
R_SCRIPT="generate_report.R"

# Copy Rmd to folder2
cp "${MYREPORT}/${RMD_FILE}" "${MYWD}/"

# Move into folder2
cd "${MYWD}"

# Run R script
Rscript "${MYREPORT}/${R_SCRIPT}"

# Remove copied Rmd after successful render
rm "${RMD_FILE}"
