#!/bin/bash
prodigal -i $1 -o genes -a proteins -f gff -q
if [ $? -eq 0 ]; then
	echo "Prodigal succeeded at predicting genes"
else
	echo "Prodigal failed."
fi
hmmscan --tblout txtt --noali -E 0.01 -o temp All.hmm proteins
if [ $? -eq 0 ]; then
	echo "Hmmscan succeeded at checking proteins for transposases"
else
	echo "Hmmscan failed."
fi
cp proteins $1.faa
tail -n +4 txtt | head -n -10 > hmmscan_clear
grep ">" proteins > headers
python analyzer.py > $1.gdf
grep "Transposase" $1.gdf | awk '{print $1, $2, $3}' > $1.tr.bed
if [ ! -s $1.tr.bed ]; then
	echo "Transposases not found"
else
	z=`cat $1.tr.bed | wc -l`
	echo $z "transposase(s) found in your bacteria"
fi
rm temp
python find_repeat_length_below51-5000.py $1 $1.tr.bed
python transposase_or_not.py $1.gdf repeats_$1.tr.bed
echo "See your transposons in" $1.gdf_Information.csv
echo "You can retrieve inner protein sequences from" $1.faa
