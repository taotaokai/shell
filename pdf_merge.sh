#!/bin/bash

#======
usage() {
cat<<EOF
NAME
  
  pdf_merge.sh - merge pdf files using gs

SYNOPSIS

  pdf_merge.sh <out file> <1.pdf> [<2.pdf> ...]

DESCRIPTION

  http://stackoverflow.com/questions/2507766/merge-convert-multiple-pdf-files-into-one-pdf
  
  Try the good ghostscript:
  
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=merged.pdf mine1.pdf mine2.pdf
  
  or even this way for an improved version for low resolution PDFs (thanks to Adriano for pointing this out):
  
  gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=merged.pdf mine1.pdf mine2.pdf
  
  In both cases the ouput resolution is much higher and better than this way using convert:
  
  convert -density 300x300 -quality 100 mine1.pdf mine2.pdf merged.pdf
  
  In this way you wouldn't need to install anything else, just work with what you already have installed in your system (at least both come by default in my rhel).
  
  Hope this helps,
  
  UPDATE: first of all thanks for all your nice comments!! just a tip that may work for you guys, after googling, I found a superb trick to shrink the size of PDFs, I reduced with it one PDF of 300 MB to just 15 MB with an acceptable resolution! and all of this with the good ghostscript, here it is:
  
  gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=output.pdf input.pdf

EOF
}

#====== main
if [[ $# -lt 2 ]]
then
  usage
  echo "[ERROR] need at least 2 inputs"
  exit -1
fi

outfile=$1

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress \
  -sOutputFile=$outfile ${@:2}