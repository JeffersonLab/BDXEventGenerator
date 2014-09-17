#! /bin/sh

for file in *.eps
do
  echo ">> Converting file $file"
  ./doc/epstosmth --gsopt='-r60x60 -dGraphicsAlphaBits=4' --gsdev=jpeg $file
done
