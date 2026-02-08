#!/bin/bash

# << Find the sequence >>

ls

echo "Please choose a FASTA file."
read FASTA
echo "Please choose a list file (csv)."
read LIST

NTH_START="0"
NTH_END=$(($(wc -l $LIST | awk '{print $1}') + 1))

while [ $NTH_START -lt "$NTH_END" ]; do

	NTH_START=$(($NTH_START + 1))

	HEAD=$(head -n $NTH_START $LIST | tail -n 1)

		FIND_HEAD=$(awk "/${HEAD}/{print NR "d"; exit 0}" "$FASTA")

		FIND_TAIL=$(sed "1,${FIND_HEAD}d" "$FASTA" | awk '/fig/{print NR; exit 0}')

		FIND_SEQUENCE=$(sed "1,${FIND_HEAD}d" "$FASTA" | sed "${FIND_TAIL},7202005772004d" | tr "\n" " " | sed "s/ //g")

	firefox --new-tab "https://alphafold.ebi.ac.uk/search/sequence/${FIND_SEQUENCE}"

done
