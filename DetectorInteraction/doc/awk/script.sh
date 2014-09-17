#! /bin/sh

if [ $# -ne 1 ]
then
  echo " Usage: $0 output_file"
  echo " output_file - output file in HTML format."
  exit
fi

awk -f branches_gen.awk \
  ../../test/ExRootLHEFConverter.cpp \
  ../../test/ExRootSTDHEPConverter.cpp \
  ../../test/ExRootLHCOlympicsConverter.cpp > $1
awk -f classes_gen.awk ../../ExRootAnalysis/ExRootClasses.h >> $1

