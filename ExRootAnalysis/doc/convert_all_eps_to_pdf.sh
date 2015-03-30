#! /bin/sh

for file in *.eps
do
  echo ">> Converting file $file" 
  ./doc/epstosmth --gsdev=pdfwrite $file
done
