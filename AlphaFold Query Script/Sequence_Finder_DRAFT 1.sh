#!/bin/bash

# << Reading of the file >>

FASTA=$1
LIST="$2"
NTH_END=$(awk '{n++} {print n}' $2)
NTH_START="0"

while [ "$NTH_START" -lt "85" ]; do

	NTH=$(($NTH + 1))

	RREAD=$(Rscript -e "

		idx <- read.csv(\"$LIST\", header = FALSE)

		print(idx[${NTH},])

"
		)

	echo "$RREAD" | grep "[1]" | awk '{print $2}' | sed "s/\"//g"

done


# << R will find the sequence >>

	RLIST=$(R --vanilla <<EOF

		print("MASDHDADSdaSADJOa")

EOF
		)

	echo "$RLIST" | grep "[1]" | awk '{print $2}' | sed "s/\"//g"
