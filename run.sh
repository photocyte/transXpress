#!/bin/bash

set -e

if [ $# -lt 1 ]; then
  echo Usage:
  echo "./run.sh <ASSEMBLER> <PROFILE>" 
  echo "where <ASSEMBLER> is either 'trinity' or 'rnaspades', and <PROFILE> are the available profiles in nextflow.config"
  exit 1
fi

if [ $# -gt 1 ]; then
 THEPROFILE=$2
else
 THEPROFILE="notStrandSpecific_local"
fi

ASSEMBLER=$1

export NXF_ANSI_LOG='false'
GIT_DIR=$(dirname $(readlink -f ./transXpress.nf))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
OUTFILE="transXpress-$ASSEMBLER.stdout.log"
ERRFILE="transXpress-$ASSEMBLER.stderr.log"
echo "$(date)" | tee -a $OUTFILE
echo "transXpress-nextflow now running. git hash: "${GIT_HASH} | tee -a $OUTFILE
echo "Nextflow profile set to:$THEPROFILE" | tee -a $OUTFILE
echo "Logs are being written to $OUTFILE and $ERRFILE in the current directory" | tee -a $OUTFILE
echo "Try 'lsof $OUTFILE' if you need to get the process id of the nextflow manager" | tee -a $OUTFILE
echo "'tail -f $OUTFILE' will let you see the output of nextflow in real time" | tee -a $OUTFILE 
echo "transXpress-nextflow dropping to background on host "$HOSTNAME"..." | tee -a $OUTFILE
/lab/solexa_weng/testtube/nextflow/nextflow-19.02.0-edge-all run transXpress.nf -w work-$ASSEMBLER -profile $THEPROFILE --assembler $ASSEMBLER --samples 'samples.txt' --species 'species.txt' -resume 1>>$OUTFILE 2>$ERRFILE &
disown


