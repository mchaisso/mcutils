#!/usr/bin/env bash
blDotPlot $1 $2 $3 --maxCount 10 > tmp.dots
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

Rscript $DIR/RenderPlot.R --dots tmp.dots --out $4 --title $5 
