#! /bin/bash

if [[ $# != 1 ]]; then
    echo "ERROR: argument needed!"
    exit 1
fi

metadata=$1
expected_header=$'sample\ttime\texperiment\tpath_R1\tpath_R2'
expected_times="T0 TSG TLG TR"

detect_blanks() {
    awk -F'\t' 'NR==1 { next } {
        for (i=1; i<=NF; i++) {
            if ($i ~ /^[[:space:]]*$/) {
                print "ERROR: Empty field detected at line", NR, "column", i
                exit 1
            }
        }
    }' "$1"
}

duplicate_check() {
    declare -A seen
    for sample in "$@"; do
        if [[ -n "${seen[$sample]}" ]]; then
            return 0   # check passed
        fi
        seen[$sample]=1
    done
    return 1   # Check failed: duplicated samples
}

detect_times() {
    awk -F'\t'  -v allowed="$expected_times" ' BEGIN {
        n = split(allowed, a, " ")   # divide la cadena en palabras
        for (i=1; i<=n; i++) {
            permitted[a[i]] = 1
        }
    }
    NR==1 { next } {
        t=$2
        seen[t]=1
    }
    END {
        for (t in permitted) {
            if (!(t in seen)) {
                print "ERROR: there are no samples in the metadata file with time =", t
                print "Please check"
                exit 1
            }
        }
        for (t in seen) {
            if (!(t in permitted)) {
                print "ERROR: Unexpected time value in metadata:", t
                exit 1
            }
        }
    }
    ' $1
}

# Check header
read -r header < "$metadata" # read first line
if [[ "$header" != "$expected_header" ]]; then
    echo -e "ERROR: Metadata file is incorrect!\nThe header of the file must be: ${expected_header}\nArgument parsed: ${header}"
    exit 1
fi

# check duplicate sample names
mapfile -t samples < <(awk -F'\t' 'NR>1 {print $1}' "$metadata")
if duplicate_check "${samples[@]}"; then
    echo -e "ERROR: Some samples have the same name!\nPlease check"
    exit 1
fi

# check empty fields
detect_blanks ${metadata} || exit 1

# check that all times are present
detect_times ${metadata} || exit 1

# chect that are >=3 experiments/replicates per time/condition
for time2search in $expected_times; do
    if [[ $(grep -cw ${time2search} <(cut -f 2 $metadata)) <3 ]]; then
        echo -e "ERROR: at least 3 replicates are needed for ${time2search}\nPlease check"
        exit 1
    fi
done

# check if files exist
for R1 in $(cut -f 4 $metadata | tail -n+2); do
    if [[ ! -s ${R1} ]]; then
        echo -e "ERROR: ${R1} file does not exist!\nPlease check"
        exit 1
    fi
done
for R2 in $(cut -f 5 $metadata | tail -n+2); do
    if [[ ! -s ${R2} ]]; then
        echo -e "ERROR: ${R2} file does not exist!\nPlease check"
        exit 1
    fi
done

# Compute stats of reads and prompt error if not fastq -> Do in one step after moving the files!
# seqkit stat -a -T $file 2>>/dev/null | cut -f 1,2,4,6-8,15-17  # seqkit >2.9; 2.13 compatible with digraph

# END
echo "ALL GOOD!"
# sed "s/metadata check\tPENDING/metadata check\tDONE"
