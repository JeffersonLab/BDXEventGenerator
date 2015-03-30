#! /bin/sh

root -l -b <<- EOF
  gSystem->Load("../libExRootAnalysis.so");
  .X Example.C("test.list");
  .q
EOF

