#!/bin/bash

conda activate vinagpu
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

cat astex_diverse_set_ids.txt | tr "\n" " " > astex_info.txt

n="0"

DIR="/home/li/Downloads/posebusters_paper_data"

mkdir results

while [ $n -lt "85" ]; do

	n=$(($n + 1))
	
	PROT_REF=$(cat astex_info.txt | awk -v x="$n" '{print $x}')
	read LIGAND <<< $PROT_REF

	SPECIFIC_DIR="${DIR}/astex_diverse_set/${PROT_REF}"

		cd $SPECIFIC_DIR
	
		obminimize -ff MMFF94 -n 100000 -h ${LIGAND}_ligand.sdf > ${LIGAND}_ligand.pdb
		obabel -ipdb "${LIGAND}_ligand.pdb" -omol2 -O "${LIGAND}_ligand.mol2"

    export DIR_CARRYOVER="${SPECIFIC_DIR}"
    export LIGAND_CARRYOVER="${LIGAND}_ligand.mol2"

cat <<-EOF > input.pgn
#!home/li/miniconda3/bin/tclsh
cd "$::env(DIR_CARRYOVER)"
mol new "$::env(LIGAND_CARRYOVER)" type mol2 waitfor all
 
puts "Gyration: [ expr [ measure rgyr [ atomselect 0 all ] ] * 2.9 ]"
exit
EOF
    
    GYRATION_OUT=$(prime-run vmd -dispdev none -e input.pgn | grep "Gyration:" | awk {'print $2'})
    read GYRATION <<< $GYRATION_OUT
    
		
	START_LINE=$(awk '/@<TRIPOS>ATOM/{print NR; exit}' ${LIGAND}_ligand.mol2)
	FINISH_LINE=$(awk '/@<TRIPOS>BOND/{print NR; exit}' ${LIGAND}_ligand.mol2)


	awk -v x="$START_LINE" -v y="$FINISH_LINE" 'NR>x&&NR<y' ${LIGAND}_ligand.mol2 > ${LIGAND}.log

	X=$(awk '{x+=$3; n++} END {print x/n}' ${LIGAND}.log)
	Y=$(awk '{x+=$4; n++} END {print x/n}' ${LIGAND}.log)
	Z=$(awk '{x+=$5; n++} END {print x/n}' ${LIGAND}.log)

	#vega ${LIGAND}_protein.pdb -w -f pdb -o ${LIGAND}_vega.pdb
	#grep -vE "HETATM|CONECT" ${LIGAND}_vega.pdb > ${LIGAND}_decat.pdb

	export PROT="${LIGAND}_protein.pdb"
	export PROT_ONLYNAME="${LIGAND}_decat.pdb"
	
cat <<-EOF > prot_isolation.pgn
cd "$::env(DIR_CARRYOVER)"
package require psfgen
mol new "$::env(PROT)" type pdb waitfor all
[ atomselect 0 protein ] writepdb "$::env(PROT_ONLYNAME)"
exit
EOF

    prime-run vmd -dispdev none -e prot_isolation.pgn
	
	pdb_reatom -1 ${LIGAND}_decat.pdb > ${LIGAND}_recat.pdb 
	##pdb_selchain -A ${LIGAND}_recat.pdb > ${LIGAND}_chainonly.pdb
	pdbfixer ${LIGAND}_recat.pdb --add-atoms all --add-residues --replace-nonstandard --output ${LIGAND}_fixed.pdb

	PYTHON_DIR="/home/li/Documents/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
	UTIL_DIR="/home/li/Documents/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"

	$PYTHON_DIR \ $UTIL_DIR/prepare_receptor.py -r ${LIGAND}_fixed.pdb -U -A bonds
	$PYTHON_DIR \ $UTIL_DIR/pdbqs_to_pdbqt.py -s ${LIGAND}_fixed -o ${LIGAND}_pdockready.pdbqt
	$PYTHON_DIR \ $UTIL_DIR/prepare_ligand4.py -l ${LIGAND}_ligand.mol2 -o ${LIGAND}_ldockready.pdbqt
	
	mv ${LIGAND}_ligand.mol2 ${DIR}/results
	cp ${LIGAND}_pdockready.pdbqt ${DIR}/results
	
	VINA_DIR="/home/li/Documents/Vina-GPU-2.1/QuickVina2-GPU-2.1/"

	cd $VINA_DIR

	prime-run ./QuickVina2-GPU-2-1 --receptor "${SPECIFIC_DIR}/${LIGAND}_pdockready.pdbqt" --ligand "${SPECIFIC_DIR}/${LIGAND}_ldockready.pdbqt" --size_x $GYRATION --size_y $GYRATION --size_z $GYRATION --center_x $X --center_y $Y --center_z $Z --thread 8000 --out ${LIGAND}_out.pdbqt

	vina_split --input ${LIGAND}_out.pdbqt

	mv *.pdbqt ${DIR}/results

	cd ${DIR}/results

	#for results in *.pdbqt; do

	#	obabel -ipdbqt $results -omol2 -O ${results%.*}.mol2

	#done

	#rm *.pdbqt

	#for results in ${LIGAND}_out_*; do
			
		#obrms ${LIGAND}_ligand.mol2 $results >> redocking_results.csv
	
	#done

##DOCKER
for results in ${LIGAND}_out_ligand*; do

    podman run --rm -it -v $DIR:/dir:Z --device nvidia.com/gpu=all gnina /bin/bash -c "
        
    cd /dir/results
    
    gnina -r ${LIGAND}_pdockready.pdbqt -l $results --size_x $GYRATION --size_y $GYRATION --size_z $GYRATION --center_x $X --center_y $Y --center_z $Z --minimize -o ${results%.*}_rescored.pdbqt
    "

done

for rescored in *_rescored.pdbqt; do
    
    obabel -ipdbqt $rescored -omol2 -O ${rescored%_rescored.pdbqt}.mol2
    obrms ${LIGAND}_ligand.mol2 ${rescored%_rescored.pdbqt}.mol2 >> redock_rescore.csv
    
done


    OUTPUT=$(cat redock_rescore.csv | awk '{print $3}' | datamash mean 1)
    
    echo ""
    echo "CURRENT AVERAGE RMSD: $OUTPUT" >> RMSD_convergence.csv
    echo ""
    
    rm *.pdbqt    
    
	cd $DIR

done
