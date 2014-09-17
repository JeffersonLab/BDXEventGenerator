#! /bin/sh

for file in *.h *.cc
do
  echo $file
  sed 's|CLHEPNAMESPACE HepLorentzVector|CLHEP::HepLorentzVector|g' $file > ${file}.tmp
  mv ${file}.tmp ${file}
done
