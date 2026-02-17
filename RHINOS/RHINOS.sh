#!/bin/env bash

# This script assumes that you already have the binaries installed:
# 1. VMD
# 2. PDBFixer
# 3. PoseBusters Validity Checker
# 4. Conda
# << Setting up of the working directories >>

eval "$(conda shell.bash hook)"

	GNINA=""
	VINA_GPU=""
	NAMD=""
    PDBFIXER_ACTIVATE="conda activate vinagpu"
    MEEKO_ACTIVATE=$(conda activate meeko)
    GNINA_CONDA_ENV=""
    MGL_TOOLS_UTILITIES="/home/li/Documents/pseudobin/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/"
    MGL_TOOLS_PYTHONSH="/home/li/Documents/pseudobin/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"\

# << Exporting conda env var >>

	setLD()
	{
		cat <<-EOF > ld_lib.sh
			export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH"
		EOF

		chmod u+x ld_lib.sh
		source ld_lib.sh
		echo "LD_LIBRARY_PATH is exported."

			rm ld_lib.sh
	}

error_message()
{

echo "ERROR: No input provided. Please follow the format..."
echo "REMEDY: ./Redocking_Script -l <ligand_directory> -p <protein_directory> -f <ligand_file_format> -d <RHINOS_directory>"

}

	isFolder()
	{

	if ls | grep -ox "$1" > /dev/null; then

	echo "Folder found: $1"

	else

	echo "The inputs provided are not folders/directories: $1"
	exit 1

	fi

	}


	while getopts ":l:p:f:d:" options; do
		case ${options} in
			l)
			 LIGAND_FOLDER="${OPTARG}"
			 isFolder "${LIGAND_FOLDER::-1}"
			 LIGAND_REALPATH=$(realpath "${LIGAND_FOLDER}")
			 ;;
			p)
			 PROTEIN_FOLDER="${OPTARG}"
			 isFolder "${PROTEIN_FOLDER::-1}"
			 PROTEIN_REALPATH=$(realpath "${PROTEIN_FOLDER}")
			 ;;
			f)
			 FILE_FORMAT="${OPTARG}"
             ;;
			d)
			 RHINOS_DIR="${OPTARG}"
			 ;;
			:)
			error_message
			 exit 1
			 ;;
			?)
			 echo "Unknown option: -${OPTARG}"
			 exit 1
			 ;;
		esac
	done

	if [ $# -eq 0 ]; then

		error_message

		exit 1

	elif [ -z "$LIGAND_FOLDER" ] || [ -z "$PROTEIN_FOLDER" ] || [ -z "$FILE_FORMAT" ] || [ -z "$RHINOS_DIR" ]; then
		echo ""
		echo "ERROR: Missing an argument."
		error_message
		exit 1

	fi

	adt_util_ligand()
	{
        local LIG_OUT="$2"
		local ADT_UTIL="$MGL_TOOLS_UTILITIES"
		local PYTHONSH="$MGL_TOOLS_PYTHONSH"

		${PYTHONSH} \ ${ADT_UTIL}/${1} -l $LIG_OUT -A '' -U '' -o ${LIG_OUT%.*}.pdbqt

	}

know_gyration()
{

local LIGAND="$1"

BOX_RADIUS=$(python - <<-EOF

from rdkit import Chem
from rdkit.Chem import Descriptors3D
import sys

molecule = Chem.MolFromPDBFile(str("${LIGAND}"))
gyration = Descriptors3D.RadiusOfGyration(molecule)

print(gyration * 2.9)

EOF
)

}

	minconv_ligand_folder()
	{
		local INPUT_FOLDER="$1"

		cd $INPUT_FOLDER
		
		echo "Minimizing Ligands..."
		
		for ligands in *.${FILE_FORMAT}; do

			obminimize -ff UFF -c 1e-6 "${ligands}" > ${ligands%.*}.pdb 2> /dev/null

            know_gyration ${ligands%.*}.pdb

			adt_util_ligand prepare_ligand4.py ${ligands%.*}.pdb > /dev/null
			
			echo "${ligands%.*}.pdbqt ${BOX_RADIUS}" >> input_files.txt
			
			rm ${ligands%.*}.pdb

			done

        echo "Done!"

		cd -
			}

minconv_protein_folder()
{
local INPUT_FOLDER="$1"

cd $INPUT_FOLDER 2> /dev/null

# << Inside the folder >>

    eval $PDBFIXER_ACTIVATE

   for proteins in *.pdb; do
   
   CHAIN=$(eval "echo "$proteins" | cut -d'_' -f2 | cut -d'.' -f1")
   
    echo "Isolating chains..."
      
        if echo "$proteins" | grep -ox "_${CHAIN}.pdb" > /dev/null; then
    
           python3 "${RHINOS_DIR}"/python_scripts/prody_selchain.py "$proteins" "${CHAIN}" "${proteins%.*}_chonly.pdb" 2> /dev/null
        
         else
    
           echo "Assuming the structure is already processed for: ${proteins}"
         
        fi

    done
    
    echo "Preparing proteins..."

    ls *_chonly.pdb | xargs -P 5 -I{} sh -c 'pdbfixer $1 --add-atoms all --keep-heterogens none --replace-nonstandard --add-residues --output "${1%.*}_fixed.pdb" 2> /dev/null' _ {}

    echo "Done!"
    
    ls *_fixed.pdb | xargs -P 5 -I{} sh -c 'python "${0}"/python_scripts/minimize_briefly.py $1 ${1%.*}_min.pdb' "${RHINOS_DIR}" {}
}

	
# << Execution >>

minconv_ligand_folder "${LIGAND_REALPATH}"
minconv_protein_folder "${PROTEIN_REALPATH}"