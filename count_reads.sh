#THIS SCRIPT GETS THE NUMBER OF READS OF A FASTA/Q FILE

#it works like that
#sh count_reads.sh <input-file> <'p' if paired end;'s' if single end>

input_file=$1
form=$2

ext=${input_file#*.}
out_wc=$(wc -l $input_file)
n_lines=${out_wc% *}
if [ "$ext" == "fq" ]
then
	if [ "$form" == "p" ]
	then
		echo $((n_lines/8))
	else
		echo $((n_lines/4))
	fi
	
else
	if [ "$form" == "p" ]
	then
		echo $((n_lines/4))
	else
		echo $((n_lines/2))
	fi
fi


