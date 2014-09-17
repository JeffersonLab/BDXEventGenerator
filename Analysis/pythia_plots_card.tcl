##############################################################################
# Input collections (specify ROOT tree name and file listing input ROOT files)
##############################################################################

add InputCollection {STDHEP pythia_events.list}

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  jets
}

############
# Jet finder
############

module MadGraphKtJetFinder jets {
  set MaxParticleEta 5.0
  set ParticleThreshold 0.5

  set CollisionType 4
  set DistanceScheme 1
  set RecombinationScheme 1
  set ParameterR 1.0
}
