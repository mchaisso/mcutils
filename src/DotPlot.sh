#!/usr/bin/env bash
if [ "$#" -eq 0 ]; then
   echo "Usage: DotPlot.sh a.fasta b.fasta k OutName.pdf Title tmp"
   exit 0
fi

blDotPlot $1 $2 $3 --maxCount 10 > $6
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo Rscript $DIR/RenderPlot.R --dots $6 --out $4 --title "$5"
Rscript $DIR/RenderPlot.R --dots $6 --out $4 --title $5 
