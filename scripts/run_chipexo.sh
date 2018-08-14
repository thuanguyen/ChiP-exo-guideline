#!/bin/sh
SCRIPTS_CHIPEXO=(\
	"run_chipexo_yieP1"
	"run_chipexo_yieP2" \
	)

SCRIPTS_CHIPPEAK=(\
	"run_chipexo_peak"\
	)

SCRIPTS=("${SCRIPTS_CHIPEXO[@]}" "${SCRIPTS_CHIPPEAK[@]}")

echo "### number of scripts to run: "${#SCRIPTS[*]}
mkdir -p output
mkdir -p ../alignment
mkdir -p ../gff
mkdir -p ../peak

for SCRIPT in "${SCRIPTS[@]}"
do
	SCRIPT_SHELL=$SCRIPT'.sh'
	SCRIPT_OUTPUT='output/'$SCRIPT'.output'
	echo "### running a script: "$SCRIPT_SHELL
	bash $SCRIPT_SHELL 2>&1 | tee $SCRIPT_OUTPUT
done
