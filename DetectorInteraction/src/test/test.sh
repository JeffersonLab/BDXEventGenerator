#! /bin/sh

root -l -b <<- EOF
  gSystem->Load("../lib/libExRootAnalysis.so");
  .X Example.C("test.list");
  .q
EOF

