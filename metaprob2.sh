#Since we use an all vs all overlap tool, if you have 2 files with paired end reads, please use "file geenrator" script. If you have an interleaved fasta file, please put ALL first end 
#reads first (.1) and second end reads after (.2). The input file should look like this (if,for example, your sample has 100 paired end reads):
# >read1.1
# >read2.1
# >...
# >read100.1
# >read1.2
# >read2.2
# >...
# >read100.2


python_dir = #SET PYTHON SCRIPT DIRECTORY
minimap2_dir = #SET MINIMAP2 DIRECTORY
miniasm_dir = #SET MINIASM DIRECTORY
metaprob_dir = #SET METAPROB DIRECTORY
input_dir = #SET INPUT FILES DIRECTORY
output_dir = #SET OUTPUT DIRECORY
name = #SET NAME OF OUTPUT FILES (IF YOU ALREADY HAVE A FASTA AS DESCRIBED ABOVE, SKIP THE FIRST PYTHON SCRIPT AND USE THE NAME OF THE FILE AS OUTPUT NAME)
num = # SET NUMBER OF ESTIMATED SPECIES IN THE SAMPLE. IF YOU DON'T HAVE IT OR DON't WANT PROVIDE ONE, PLEASE REMOVE "-numSp $num" FROM METAPROB INPUT PARAMETERS

echo "Starting MetaProb2"

if [ ! -d "$output_dir" ]
then
    mkdir $output_dir
fi
# file generator creates a fasta from the two paired end files. (Note: they must be in fasta format)
# provide the file names as follow: first the sequences that end in .1 or /1, then the ones that end with .2 or /2 and the name of the resulting file.
python3 $python_dir/file_generator.py $input_dir/ $input_dir/ $output_dir/$name.fa

$minimap2_dir/minimap2 -X -sc -t31 $output_dir/$name.fa $output_dir/$name.fa | gzip -1 > $output_dir/$name.paf.gz

$miniasm_dir/miniasm -12 -m1 -o2 -I0.001 -s2 -i0.001 -c1 -e0 -n0 -r0.99,0.01 -f  $output_dir/$name.fa $output_dir/$name.paf.gz > /$output_dir/$name.gfa

python3 $python_dir/community_detection.py $output_dir/$name.gfa $output_dir/$name.fa $output_dir/$name.fasta $output_dir/$name    

$metaprob_dir/MetaProb -numSp $num -feature 1 -si $output_dir/$name.fasta -dirOutput $output_dir/$name
