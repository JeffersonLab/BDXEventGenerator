##############################################################################
# Input collections (specify ROOT tree name and file listing input ROOT files)
##############################################################################

add InputCollection {LHEF events.list}

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  classification analysis
}

#################
# Analysis module
#################

module MadGraphAnalysis analysis {
  set InputArray "classification/particles"
  set OutputFile "madgraph_plots.root"
}

########################
# Particle class filters
########################

module MadGraphClassFilter classification {

  add ParticleStatus 1 2

  add ClassParticles jet {1 -1 2 -2 3 -3 21}
  add ClassParticles c {4 -4}
  add ClassParticles b {5 -5}

  add ClassParticles t {6}
  add ClassParticles #bar{t} {-6}

#  add ClassParticles lepton {11 -11 13 -13 15 -15}
#  add ClassParticles neutrino {12 -12 14 -14 16 -16}

  add ClassParticles e {11 -11}
  add ClassParticles #mu {13 -13}
  add ClassParticles #tau {15 -15}

  add ClassParticles #nu_{e} {12 -12}
  add ClassParticles #nu_{#mu} {14 -14}
  add ClassParticles #nu_{#tau} {16 -16}

#  add ClassParticles g {21}

  add ClassParticles #gamma {22}
  add ClassParticles Z^{0} {23}
  add ClassParticles W^{+} {24}
  add ClassParticles W^{-} {-24}
  add ClassParticles h^{0} {25}
}

