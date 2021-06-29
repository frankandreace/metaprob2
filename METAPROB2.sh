#!/bin/bash

SCRIPT_NAME=$(basename $0)
DIR_NAME=$(dirname $0)
#TODO: PLEASE SET HERE YOUR MINIMAP2,MINIASM,METAPROB,PYTHON AND ENVIRONMENT FOLDERS
MINIMAP2_DIR=/nfsd/bcb/bcbg/andreace/tools/minimap2 #SET MINIMAP2 DIRECTORY
MINIASM_DIR=/nfsd/bcb/bcbg/andreace/tools/miniasm #SET MINIASM DIRECTORY
PYTHON_DIR=/nfsd/bcb/bcbg/andreace/python #SET PYTHON SCRIPT DIRECTORY
METAPROB_DIR=/nfsd/bcb/bcbg/andreace/tools/MetaProb/Release #SET METAPROB DIRECTORY
CONDA_ACT=/nfsd/bcb/bcbg/andreace/python/env/bin/activate #SET CONDA ENVIRONMENT

#SETTING DEFAULT PARAMETERS
DEFAULT_KSIZE=15
DEFAULT_WSIZE=10
DEFAULT_THREADS=4
DEFAULT_MAXL=20000
DEFAULT_OPT=0.001
DEFAULT_MINM=40

USAGE=$'\nUsage: '"${SCRIPT_NAME}"' [-s NUM SPECIES] [-k KMER-SIZE] [-w WINDOW-SIZE] [-m MAX-CHAINED-UTG-LENGTH] [-o OPT-PARAMETER-MODULARITY] [-l SKIP-READS-LEFT-OUT] [-r MIN-LENGTH MAX-LENGTH] <input_file> <output_folder> <name>

Arguments:
     -h              print this help and exit
     -s              number of species (default:estimated by MetaProb2)
     -k              kmer size for minimap2 (default:15)
     -w              window size for minimap2 (default:10)
     -t              number of cores (default:4)
     -c              discard chains with score < c for minimap2 (default:40)
     -m              max length chained utgs during MetaProb2 (default:20000)
     -o              optimization parameter for modularity optimization (default:0.001)
     -l              skip the reads left out by miniasm (default:False)
     -p              keep the temporary folder with all the output files(default:False)

Positional arguments:
     <input_file>         input in FASTA/Q format
     <output_folder>      output_folder
     <name>               name for output files
'


#PARSING INPUT OPTIONS

species=""
kmerlen=${DEFAULT_KSIZE}
wlen=${DEFAULT_WSIZE}
thr=${DEFAULT_THREADS}
chainv=${DEFAULT_MINM}
maxlen=${DEFAULT_MAXL}
optpar=${DEFAULT_OPT}
leftout=""
keeptemp=false

while getopts "s:k:w:t:c:m:o:lhp" flag; do
    case "${flag}" in
	h) $(>&2 echo "${USAGE}")
	   exit 0
	   ;;
	s) species="-numSp ${OPTARG}"
	   ;;
	k) kmerlen=${OPTARG}
	   ;;
	w) wlen=${OPTARG}
	   ;;
   t) thr=${OPTARG}
	   ;;
   c) chainv=${OPTARG}
	   ;;
	m) maxlen=${OPTARG}
	   ;;
	o) optpar=${OPTARG}
	   ;;
	l) leftout="--lout"
	   ;;
	p) keeptemp=true
	   ;;
    esac
done

#ADDING INPUT FILE AND OUTPUT FOLDER

input_file=""
output_dir=""
name=""

if [[ $# -lt $((${OPTIND} + 2)) ]]
then
    (>&2 echo "ERROR: Wrong number of arguments.")
    (>&2 echo "")
    (>&2 echo "${USAGE}")
    exit 1
fi

input_file=${@:$OPTIND:1}
output_dir=${@:$OPTIND+1:1}
name=${@:$OPTIND+2:1}

#TEMPORARY DIRECTORY
tmp_dir=$output_dir/metaprob2_tmp

#RUNNING METAPROB2
echo "STARTING METAPROB2"

#CREATING OUTPUT AND TEMPORARY FOLDERS
if [ ! -d "$output_dir" ]
then
    mkdir $output_dir
fi
if [ ! -d "$tmp_dir" ]
then
    mkdir $tmp_dir
fi

if [ ! -f $tmp_dir/$name.paf ] && [ ! -f $tmp_dir/$name.gfa ]
then

#RUNNING MINIMAP2
   (>&2 echo "[${SCRIPT_NAME}] Running MINIMAP2")
   $MINIMAP2_DIR/minimap2 -X -sr -B8 -O12,32 -E2,1 -m$chainv -k$kmerlen -w$wlen --end-bonus=100 -t$thr $input_file $input_file > $tmp_dir/$name.paf 

#RUNNING MINIASM
   (>&2 echo "[${SCRIPT_NAME}] Running MINIASM")
   $MINIASM_DIR/miniasm -12 -m1 -g10000 -o1 -d1 -I0.00001 -s1 -i0.00001 -d10000 -F0.01 -c0 -e0 -n0 -r0.011,0.01 -f $input_file $tmp_dir/$name.paf > $tmp_dir/$name.gfa
else
    (>&2 echo  "[${SCRIPT_NAME}] Found MINIMAP2 & MINIASM OUTPUT")
fi


#GETTING NUMBER OF READS IN INPUT FILE
inf=$tmp_dir/info.csv
sh $DIR_NAME/count_reads.sh $input_file p > $inf

#ACTIVATING ENVIRNMENT FOR PYTHON
source $CONDA_ACT

#RUNNING COMMUNITY DETECTION
(>&2 echo "[${SCRIPT_NAME}] Running COMMUNITY DETECTION")
python3 $PYTHON_DIR/community_detection_clear.py --gfa $tmp_dir/$name.gfa --reads $input_file --out $tmp_dir/$name --info $inf --opt $optpar --maxlen $maxlen $leftout

#DEACTIVATING ENVIRNMENT FOR PYTHON
deactivate

#RUNNING METAPROB
(>&2 echo "[${SCRIPT_NAME}] Running METAPROB")
$METAPROB_DIR/MetaProb $species -feature 1 -graphType 1 -si $tmp_dir/$name.fasta -dirOutput $output_dir/


#REMOVE TEMPORARY FOLDER IF NOT SPECIFIED IN THE INPUT
if [ "$keeptemp" = false ]
then
	rm -rf ${tmp_dir}
fi


exit 0