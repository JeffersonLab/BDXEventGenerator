
#
# Makefile for ExRootAnalysis
#
# Author: P. Demin - UCL, Louvain-la-Neuve
#
# multi-platform configuration is taken from ROOT (root/test/Makefile)
#

include doc/Makefile.arch

ifeq ($(ARCH),macosx64)
UNDEFOPT = dynamic_lookup
endif

SrcSuf = cc

CXXFLAGS += $(ROOTCFLAGS) -DDROP_CGAL -I. -Itcl -Imcfio -Istdhep -ICDFCones -ICDFCones/CDFcode
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)
	
###

STATIC = lib/libExRootAnalysisPGS.$(LibSuf)
SHARED = lib/libExRootAnalysis.$(DllSuf)
SHAREDLIB = lib/libExRootAnalysis.lib

all:


ExRootSTDHEPConverter$(ExeSuf): \
	tmp/test/ExRootSTDHEPConverter.$(ObjSuf)

tmp/test/ExRootSTDHEPConverter.$(ObjSuf): \
	test/ExRootSTDHEPConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootUtilities.h \
	ExRootAnalysis/ExRootProgressBar.h
STDHEP_EXECUTABLE =  \
	ExRootSTDHEPConverter$(ExeSuf)

STDHEP_EXECUTABLE_OBJ =  \
	tmp/test/ExRootSTDHEPConverter.$(ObjSuf)

ExRootMain$(ExeSuf): \
	tmp/test/ExRootMain.$(ObjSuf)

tmp/test/ExRootMain.$(ObjSuf): \
	test/ExRootMain.cpp \
	ExRootAnalysis/ExRootAnalysis.h
ExRootLHEFConverter$(ExeSuf): \
	tmp/test/ExRootLHEFConverter.$(ObjSuf)

tmp/test/ExRootLHEFConverter.$(ObjSuf): \
	test/ExRootLHEFConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootUtilities.h \
	ExRootAnalysis/ExRootProgressBar.h
ExRootLHCOlympicsConverter$(ExeSuf): \
	tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf)

tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf): \
	test/ExRootLHCOlympicsConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootUtilities.h \
	ExRootAnalysis/ExRootProgressBar.h
ExRootHEPEVTConverter$(ExeSuf): \
	tmp/test/ExRootHEPEVTConverter.$(ObjSuf)

tmp/test/ExRootHEPEVTConverter.$(ObjSuf): \
	test/ExRootHEPEVTConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootUtilities.h
Example$(ExeSuf): \
	tmp/test/Example.$(ObjSuf)

tmp/test/Example.$(ObjSuf): \
	test/Example.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootUtilities.h
EXECUTABLE =  \
	ExRootMain$(ExeSuf) \
	ExRootLHEFConverter$(ExeSuf) \
	ExRootLHCOlympicsConverter$(ExeSuf) \
	ExRootHEPEVTConverter$(ExeSuf) \
	Example$(ExeSuf)

EXECUTABLE_OBJ =  \
	tmp/test/ExRootMain.$(ObjSuf) \
	tmp/test/ExRootLHEFConverter.$(ObjSuf) \
	tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf) \
	tmp/test/ExRootHEPEVTConverter.$(ObjSuf) \
	tmp/test/Example.$(ObjSuf)

tmp/src/ExRootAnalysisDict.$(SrcSuf): \
	src/ExRootAnalysisLinkDef.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootUtilities.h \
	ExRootAnalysis/ExRootClassifier.h \
	ExRootAnalysis/ExRootFilter.h \
	ExRootAnalysis/ExRootProgressBar.h \
	ExRootAnalysis/ExRootConfReader.h \
	ExRootAnalysis/ExRootAnalysis.h \
	ExRootAnalysis/ExRootFactory.h \
	ExRootAnalysis/ExRootModule.h \
	ExRootAnalysis/ExRootTask.h
tmp/modules/ModulesDict.$(SrcSuf): \
	modules/ModulesLinkDef.h \
	modules/MadGraphClassFilter.h \
	modules/MadGraphAnalysis.h \
	modules/MadGraphKtJetFinder.h \
	modules/MadGraphConeJetFinder.h
DICT =  \
	tmp/src/ExRootAnalysisDict.$(SrcSuf) \
	tmp/modules/ModulesDict.$(SrcSuf)

DICT_OBJ =  \
	tmp/src/ExRootAnalysisDict.$(ObjSuf) \
	tmp/modules/ModulesDict.$(ObjSuf)

tmp/pgs/ExRootAnalysisDict.$(SrcSuf): \
	pgs/ExRootAnalysisLinkDef.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h
PGS_DICT =  \
	tmp/pgs/ExRootAnalysisDict.$(SrcSuf)

PGS_DICT_OBJ =  \
	tmp/pgs/ExRootAnalysisDict.$(ObjSuf)

tmp/src/ExRootFilter.$(ObjSuf): \
	src/ExRootFilter.$(SrcSuf) \
	ExRootAnalysis/ExRootFilter.h \
	ExRootAnalysis/ExRootClassifier.h
tmp/src/ExRootConfReader.$(ObjSuf): \
	src/ExRootConfReader.$(SrcSuf) \
	ExRootAnalysis/ExRootConfReader.h \
	tcl/tcl.h
tmp/src/ExRootFactory.$(ObjSuf): \
	src/ExRootFactory.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootFactory.h
tmp/src/ExRootTreeWriter.$(ObjSuf): \
	src/ExRootTreeWriter.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h
tmp/src/ExRootResult.$(ObjSuf): \
	src/ExRootResult.$(SrcSuf) \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootUtilities.h
tmp/src/ExRootModule.$(ObjSuf): \
	src/ExRootModule.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeReader.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootModule.h \
	ExRootAnalysis/ExRootFactory.h
tmp/src/ExRootClasses.$(ObjSuf): \
	src/ExRootClasses.$(SrcSuf) \
	ExRootAnalysis/ExRootClasses.h
tmp/src/ExRootTreeBranch.$(ObjSuf): \
	src/ExRootTreeBranch.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeBranch.h
tmp/src/ExRootTreeReader.$(ObjSuf): \
	src/ExRootTreeReader.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeReader.h
tmp/src/ExRootTask.$(ObjSuf): \
	src/ExRootTask.$(SrcSuf) \
	ExRootAnalysis/ExRootTask.h \
	ExRootAnalysis/ExRootConfReader.h
tmp/src/ExRootAnalysis.$(ObjSuf): \
	src/ExRootAnalysis.$(SrcSuf) \
	ExRootAnalysis/ExRootAnalysis.h \
	ExRootAnalysis/ExRootFactory.h \
	ExRootAnalysis/ExRootConfReader.h \
	ExRootAnalysis/ExRootTreeReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootUtilities.h \
	ExRootAnalysis/ExRootProgressBar.h
tmp/src/ExRootProgressBar.$(ObjSuf): \
	src/ExRootProgressBar.$(SrcSuf) \
	ExRootAnalysis/ExRootProgressBar.h
tmp/src/ExRootUtilities.$(ObjSuf): \
	src/ExRootUtilities.$(SrcSuf) \
	ExRootAnalysis/ExRootUtilities.h
tmp/modules/MadGraphAnalysis.$(ObjSuf): \
	modules/MadGraphAnalysis.$(SrcSuf) \
	modules/MadGraphAnalysis.h \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootClasses.h
tmp/modules/MadGraphConeJetFinder.$(ObjSuf): \
	modules/MadGraphConeJetFinder.$(SrcSuf) \
	modules/MadGraphConeJetFinder.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootFactory.h \
	CDFCones/JetCluAlgorithm.hh \
	CDFCones/MidPointAlgorithm.hh
tmp/modules/MadGraphKtJetFinder.$(ObjSuf): \
	modules/MadGraphKtJetFinder.$(SrcSuf) \
	modules/MadGraphKtJetFinder.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootFactory.h \
	KtJet/KtEvent.h \
	KtJet/KtLorentzVector.h
tmp/modules/MadGraphClassFilter.$(ObjSuf): \
	modules/MadGraphClassFilter.$(SrcSuf) \
	modules/MadGraphClassFilter.h \
	ExRootAnalysis/ExRootResult.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootFilter.h \
	ExRootAnalysis/ExRootClassifier.h
tmp/CDFCones/ClusterComparisons.$(ObjSuf): \
	CDFCones/ClusterComparisons.$(SrcSuf)
tmp/CDFCones/CalTower.$(ObjSuf): \
	CDFCones/CalTower.$(SrcSuf)
tmp/CDFCones/Cluster.$(ObjSuf): \
	CDFCones/Cluster.$(SrcSuf)
tmp/CDFCones/Centroid.$(ObjSuf): \
	CDFCones/Centroid.$(SrcSuf)
tmp/CDFCones/MidPointAlgorithm.$(ObjSuf): \
	CDFCones/MidPointAlgorithm.$(SrcSuf)
tmp/CDFCones/LorentzVector.$(ObjSuf): \
	CDFCones/LorentzVector.$(SrcSuf)
tmp/CDFCones/PhysicsTower.$(ObjSuf): \
	CDFCones/PhysicsTower.$(SrcSuf)
tmp/CDFCones/JetCluAlgorithm.$(ObjSuf): \
	CDFCones/JetCluAlgorithm.$(SrcSuf)
tmp/KtJet/KtEvent.$(ObjSuf): \
	KtJet/KtEvent.$(SrcSuf) \
	KtJet/KtEvent.h \
	KtJet/KtLorentzVector.h \
	KtJet/KtJetTable.h \
	KtJet/KtDistance.h \
	KtJet/KtDistanceInterface.h \
	KtJet/KtRecom.h \
	KtJet/KtRecomInterface.h
tmp/KtJet/KtRecom.$(ObjSuf): \
	KtJet/KtRecom.$(SrcSuf) \
	KtJet/KtRecom.h \
	KtJet/KtUtil.h \
	KtJet/KtRecomInterface.h
tmp/KtJet/KtUtil.$(ObjSuf): \
	KtJet/KtUtil.$(SrcSuf) \
	KtJet/KtUtil.h \
	KtJet/KtLorentzVector.h
tmp/KtJet/KtDistance.$(ObjSuf): \
	KtJet/KtDistance.$(SrcSuf) \
	KtJet/KtDistance.h \
	KtJet/KtUtil.h \
	KtJet/KtDistanceInterface.h
tmp/KtJet/KtJetTable.$(ObjSuf): \
	KtJet/KtJetTable.$(SrcSuf) \
	KtJet/KtJetTable.h \
	KtJet/KtLorentzVector.h \
	KtJet/KtDistanceInterface.h \
	KtJet/KtRecomInterface.h
tmp/KtJet/KtLorentzVector.$(ObjSuf): \
	KtJet/KtLorentzVector.$(SrcSuf) \
	KtJet/KtLorentzVector.h \
	KtJet/KtUtil.h \
	KtJet/KtRecomInterface.h
tmp/CLHEP/src/RotationInterfaces.$(ObjSuf): \
	CLHEP/src/RotationInterfaces.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h
tmp/CLHEP/src/AxisAngle.$(ObjSuf): \
	CLHEP/src/AxisAngle.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/AxisAngle.h
tmp/CLHEP/src/LorentzVectorR.$(ObjSuf): \
	CLHEP/src/LorentzVectorR.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h
tmp/CLHEP/src/LorentzRotationD.$(ObjSuf): \
	CLHEP/src/LorentzRotationD.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzRotation.h
tmp/CLHEP/src/ThreeVectorR.$(ObjSuf): \
	CLHEP/src/ThreeVectorR.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/Rotation.h
tmp/CLHEP/src/BoostZ.$(ObjSuf): \
	CLHEP/src/BoostZ.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/BoostZ.h \
	CLHEP/Vector/Boost.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/Boost.$(ObjSuf): \
	CLHEP/src/Boost.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Boost.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/RotationXYZ.$(ObjSuf): \
	CLHEP/src/RotationXYZ.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/RotationX.h \
	CLHEP/Vector/RotationY.h \
	CLHEP/Vector/RotationZ.h
tmp/CLHEP/src/RotationP.$(ObjSuf): \
	CLHEP/src/RotationP.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h
tmp/CLHEP/src/ZMxpv.$(ObjSuf): \
	CLHEP/src/ZMxpv.$(SrcSuf) \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/RotationX.$(ObjSuf): \
	CLHEP/src/RotationX.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationX.h \
	CLHEP/Vector/AxisAngle.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/RotationA.$(ObjSuf): \
	CLHEP/src/RotationA.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/BoostY.$(ObjSuf): \
	CLHEP/src/BoostY.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/BoostY.h \
	CLHEP/Vector/Boost.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/Rotation.$(ObjSuf): \
	CLHEP/src/Rotation.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/LorentzRotation.$(ObjSuf): \
	CLHEP/src/LorentzRotation.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/BoostX.$(ObjSuf): \
	CLHEP/src/BoostX.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/BoostX.h \
	CLHEP/Vector/Boost.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/RotationE.$(ObjSuf): \
	CLHEP/src/RotationE.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/RotationZ.$(ObjSuf): \
	CLHEP/src/RotationZ.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationZ.h \
	CLHEP/Vector/AxisAngle.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/RotationL.$(ObjSuf): \
	CLHEP/src/RotationL.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/LorentzRotation.h
tmp/CLHEP/src/EulerAngles.$(ObjSuf): \
	CLHEP/src/EulerAngles.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Vector/ThreeVector.h
tmp/CLHEP/src/LorentzVector.$(ObjSuf): \
	CLHEP/src/LorentzVector.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/RotationY.$(ObjSuf): \
	CLHEP/src/RotationY.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationY.h \
	CLHEP/Vector/AxisAngle.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/ThreeVector.$(ObjSuf): \
	CLHEP/src/ThreeVector.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/ZMxpv.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/RotationIO.$(ObjSuf): \
	CLHEP/src/RotationIO.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h
tmp/CLHEP/src/LorentzVectorC.$(ObjSuf): \
	CLHEP/src/LorentzVectorC.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h
tmp/CLHEP/src/SpaceVectorP.$(ObjSuf): \
	CLHEP/src/SpaceVectorP.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/LorentzVectorK.$(ObjSuf): \
	CLHEP/src/LorentzVectorK.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/LorentzRotationC.$(ObjSuf): \
	CLHEP/src/LorentzRotationC.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzRotation.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/RotationC.$(ObjSuf): \
	CLHEP/src/RotationC.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/SpaceVectorD.$(ObjSuf): \
	CLHEP/src/SpaceVectorD.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/LorentzVectorB.$(ObjSuf): \
	CLHEP/src/LorentzVectorB.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/SpaceVectorR.$(ObjSuf): \
	CLHEP/src/SpaceVectorR.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/AxisAngle.h \
	CLHEP/Vector/EulerAngles.h \
	CLHEP/Vector/ZMxpv.h
tmp/CLHEP/src/SpaceVector.$(ObjSuf): \
	CLHEP/src/SpaceVector.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/ZMxpv.h \
	CLHEP/Units/PhysicalConstants.h
tmp/CLHEP/src/LorentzVectorL.$(ObjSuf): \
	CLHEP/src/LorentzVectorL.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/LorentzRotation.h
tmp/CLHEP/src/TwoVector.$(ObjSuf): \
	CLHEP/src/TwoVector.$(SrcSuf) \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/TwoVector.h \
	CLHEP/Vector/ZMxpv.h \
	CLHEP/Vector/ThreeVector.h
tmp/CLHEP/src/ZMinput.$(ObjSuf): \
	CLHEP/src/ZMinput.$(SrcSuf) \
	CLHEP/Vector/defs.h
SOURCE_OBJ =  \
	tmp/src/ExRootFilter.$(ObjSuf) \
	tmp/src/ExRootConfReader.$(ObjSuf) \
	tmp/src/ExRootFactory.$(ObjSuf) \
	tmp/src/ExRootTreeWriter.$(ObjSuf) \
	tmp/src/ExRootResult.$(ObjSuf) \
	tmp/src/ExRootModule.$(ObjSuf) \
	tmp/src/ExRootClasses.$(ObjSuf) \
	tmp/src/ExRootTreeBranch.$(ObjSuf) \
	tmp/src/ExRootTreeReader.$(ObjSuf) \
	tmp/src/ExRootTask.$(ObjSuf) \
	tmp/src/ExRootAnalysis.$(ObjSuf) \
	tmp/src/ExRootProgressBar.$(ObjSuf) \
	tmp/src/ExRootUtilities.$(ObjSuf) \
	tmp/modules/MadGraphAnalysis.$(ObjSuf) \
	tmp/modules/MadGraphConeJetFinder.$(ObjSuf) \
	tmp/modules/MadGraphKtJetFinder.$(ObjSuf) \
	tmp/modules/MadGraphClassFilter.$(ObjSuf) \
	tmp/CDFCones/ClusterComparisons.$(ObjSuf) \
	tmp/CDFCones/CalTower.$(ObjSuf) \
	tmp/CDFCones/Cluster.$(ObjSuf) \
	tmp/CDFCones/Centroid.$(ObjSuf) \
	tmp/CDFCones/MidPointAlgorithm.$(ObjSuf) \
	tmp/CDFCones/LorentzVector.$(ObjSuf) \
	tmp/CDFCones/PhysicsTower.$(ObjSuf) \
	tmp/CDFCones/JetCluAlgorithm.$(ObjSuf) \
	tmp/KtJet/KtEvent.$(ObjSuf) \
	tmp/KtJet/KtRecom.$(ObjSuf) \
	tmp/KtJet/KtUtil.$(ObjSuf) \
	tmp/KtJet/KtDistance.$(ObjSuf) \
	tmp/KtJet/KtJetTable.$(ObjSuf) \
	tmp/KtJet/KtLorentzVector.$(ObjSuf) \
	tmp/CLHEP/src/RotationInterfaces.$(ObjSuf) \
	tmp/CLHEP/src/AxisAngle.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVectorR.$(ObjSuf) \
	tmp/CLHEP/src/LorentzRotationD.$(ObjSuf) \
	tmp/CLHEP/src/ThreeVectorR.$(ObjSuf) \
	tmp/CLHEP/src/BoostZ.$(ObjSuf) \
	tmp/CLHEP/src/Boost.$(ObjSuf) \
	tmp/CLHEP/src/RotationXYZ.$(ObjSuf) \
	tmp/CLHEP/src/RotationP.$(ObjSuf) \
	tmp/CLHEP/src/ZMxpv.$(ObjSuf) \
	tmp/CLHEP/src/RotationX.$(ObjSuf) \
	tmp/CLHEP/src/RotationA.$(ObjSuf) \
	tmp/CLHEP/src/BoostY.$(ObjSuf) \
	tmp/CLHEP/src/Rotation.$(ObjSuf) \
	tmp/CLHEP/src/LorentzRotation.$(ObjSuf) \
	tmp/CLHEP/src/BoostX.$(ObjSuf) \
	tmp/CLHEP/src/RotationE.$(ObjSuf) \
	tmp/CLHEP/src/RotationZ.$(ObjSuf) \
	tmp/CLHEP/src/RotationL.$(ObjSuf) \
	tmp/CLHEP/src/EulerAngles.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVector.$(ObjSuf) \
	tmp/CLHEP/src/RotationY.$(ObjSuf) \
	tmp/CLHEP/src/ThreeVector.$(ObjSuf) \
	tmp/CLHEP/src/RotationIO.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVectorC.$(ObjSuf) \
	tmp/CLHEP/src/SpaceVectorP.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVectorK.$(ObjSuf) \
	tmp/CLHEP/src/LorentzRotationC.$(ObjSuf) \
	tmp/CLHEP/src/RotationC.$(ObjSuf) \
	tmp/CLHEP/src/SpaceVectorD.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVectorB.$(ObjSuf) \
	tmp/CLHEP/src/SpaceVectorR.$(ObjSuf) \
	tmp/CLHEP/src/SpaceVector.$(ObjSuf) \
	tmp/CLHEP/src/LorentzVectorL.$(ObjSuf) \
	tmp/CLHEP/src/TwoVector.$(ObjSuf) \
	tmp/CLHEP/src/ZMinput.$(ObjSuf)

PGS_SOURCE_OBJ =  \
	tmp/src/ExRootClasses.$(ObjSuf) \
	tmp/src/ExRootTreeBranch.$(ObjSuf) \
	tmp/src/ExRootTreeWriter.$(ObjSuf)

tmp/pgs/ExRootAnalysis.$(ObjSuf): \
	pgs/ExRootAnalysis.$(SrcSuf) \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootTreeWriter.h
PGS_OBJ =  \
	tmp/pgs/ExRootAnalysis.$(ObjSuf)

tmp/mcfio/mcf_NTuIOUtils.$(ObjSuf): \
	mcfio/mcf_NTuIOUtils.c
tmp/mcfio/mcf_NTuIOFiles.$(ObjSuf): \
	mcfio/mcf_NTuIOFiles.c
tmp/mcfio/mcfio_SeqDummy.$(ObjSuf): \
	mcfio/mcfio_SeqDummy.c
tmp/mcfio/mcf_evt_xdr.$(ObjSuf): \
	mcfio/mcf_evt_xdr.c
tmp/mcfio/mcfio_UserDictionary.$(ObjSuf): \
	mcfio/mcfio_UserDictionary.c
tmp/mcfio/mcfio_Direct.$(ObjSuf): \
	mcfio/mcfio_Direct.c
tmp/mcfio/mcfio_Util1.$(ObjSuf): \
	mcfio/mcfio_Util1.c
tmp/mcfio/mcfio_Block.$(ObjSuf): \
	mcfio/mcfio_Block.c
tmp/stdhep/mcf_hepup_xdr.$(ObjSuf): \
	stdhep/mcf_hepup_xdr.c
tmp/stdhep/stdhep_internal_utils.$(ObjSuf): \
	stdhep/stdhep_internal_utils.c
tmp/stdhep/stdhep_mcfio.$(ObjSuf): \
	stdhep/stdhep_mcfio.c
tmp/stdhep/mcf_stdcm1_xdr.$(ObjSuf): \
	stdhep/mcf_stdcm1_xdr.c
tmp/stdhep/mcf_Stdhep_xdr.$(ObjSuf): \
	stdhep/mcf_Stdhep_xdr.c
tmp/stdhep/mcf_hepev4_xdr.$(ObjSuf): \
	stdhep/mcf_hepev4_xdr.c
STDHEP_OBJ =  \
	tmp/mcfio/mcf_NTuIOUtils.$(ObjSuf) \
	tmp/mcfio/mcf_NTuIOFiles.$(ObjSuf) \
	tmp/mcfio/mcfio_SeqDummy.$(ObjSuf) \
	tmp/mcfio/mcf_evt_xdr.$(ObjSuf) \
	tmp/mcfio/mcfio_UserDictionary.$(ObjSuf) \
	tmp/mcfio/mcfio_Direct.$(ObjSuf) \
	tmp/mcfio/mcfio_Util1.$(ObjSuf) \
	tmp/mcfio/mcfio_Block.$(ObjSuf) \
	tmp/stdhep/mcf_hepup_xdr.$(ObjSuf) \
	tmp/stdhep/stdhep_internal_utils.$(ObjSuf) \
	tmp/stdhep/stdhep_mcfio.$(ObjSuf) \
	tmp/stdhep/mcf_stdcm1_xdr.$(ObjSuf) \
	tmp/stdhep/mcf_Stdhep_xdr.$(ObjSuf) \
	tmp/stdhep/mcf_hepev4_xdr.$(ObjSuf)

tmp/tcl/tclObj.$(ObjSuf): \
	tcl/tclObj.c
tmp/tcl/tclUtil.$(ObjSuf): \
	tcl/tclUtil.c
tmp/tcl/tclAsync.$(ObjSuf): \
	tcl/tclAsync.c
tmp/tcl/tclPosixStr.$(ObjSuf): \
	tcl/tclPosixStr.c
tmp/tcl/tclCompile.$(ObjSuf): \
	tcl/tclCompile.c
tmp/tcl/tclHistory.$(ObjSuf): \
	tcl/tclHistory.c
tmp/tcl/tclResolve.$(ObjSuf): \
	tcl/tclResolve.c
tmp/tcl/tclParse.$(ObjSuf): \
	tcl/tclParse.c
tmp/tcl/tclVar.$(ObjSuf): \
	tcl/tclVar.c
tmp/tcl/tclIndexObj.$(ObjSuf): \
	tcl/tclIndexObj.c
tmp/tcl/tclCkalloc.$(ObjSuf): \
	tcl/tclCkalloc.c
tmp/tcl/tclListObj.$(ObjSuf): \
	tcl/tclListObj.c
tmp/tcl/tclHash.$(ObjSuf): \
	tcl/tclHash.c
tmp/tcl/tclCmdIL.$(ObjSuf): \
	tcl/tclCmdIL.c
tmp/tcl/tclStringObj.$(ObjSuf): \
	tcl/tclStringObj.c
tmp/tcl/tclAlloc.$(ObjSuf): \
	tcl/tclAlloc.c
tmp/tcl/tclCompExpr.$(ObjSuf): \
	tcl/tclCompExpr.c
tmp/tcl/tclLink.$(ObjSuf): \
	tcl/tclLink.c
tmp/tcl/tclCmdAH.$(ObjSuf): \
	tcl/tclCmdAH.c
tmp/tcl/panic.$(ObjSuf): \
	tcl/panic.c
tmp/tcl/tclBasic.$(ObjSuf): \
	tcl/tclBasic.c
tmp/tcl/tclPreserve.$(ObjSuf): \
	tcl/tclPreserve.c
tmp/tcl/tclGet.$(ObjSuf): \
	tcl/tclGet.c
tmp/tcl/tclNamesp.$(ObjSuf): \
	tcl/tclNamesp.c
tmp/tcl/tclProc.$(ObjSuf): \
	tcl/tclProc.c
tmp/tcl/tclExecute.$(ObjSuf): \
	tcl/tclExecute.c
tmp/tcl/tclCmdMZ.$(ObjSuf): \
	tcl/tclCmdMZ.c
TCL_OBJ =  \
	tmp/tcl/tclObj.$(ObjSuf) \
	tmp/tcl/tclUtil.$(ObjSuf) \
	tmp/tcl/tclAsync.$(ObjSuf) \
	tmp/tcl/tclPosixStr.$(ObjSuf) \
	tmp/tcl/tclCompile.$(ObjSuf) \
	tmp/tcl/tclHistory.$(ObjSuf) \
	tmp/tcl/tclResolve.$(ObjSuf) \
	tmp/tcl/tclParse.$(ObjSuf) \
	tmp/tcl/tclVar.$(ObjSuf) \
	tmp/tcl/tclIndexObj.$(ObjSuf) \
	tmp/tcl/tclCkalloc.$(ObjSuf) \
	tmp/tcl/tclListObj.$(ObjSuf) \
	tmp/tcl/tclHash.$(ObjSuf) \
	tmp/tcl/tclCmdIL.$(ObjSuf) \
	tmp/tcl/tclStringObj.$(ObjSuf) \
	tmp/tcl/tclAlloc.$(ObjSuf) \
	tmp/tcl/tclCompExpr.$(ObjSuf) \
	tmp/tcl/tclLink.$(ObjSuf) \
	tmp/tcl/tclCmdAH.$(ObjSuf) \
	tmp/tcl/panic.$(ObjSuf) \
	tmp/tcl/tclBasic.$(ObjSuf) \
	tmp/tcl/tclPreserve.$(ObjSuf) \
	tmp/tcl/tclGet.$(ObjSuf) \
	tmp/tcl/tclNamesp.$(ObjSuf) \
	tmp/tcl/tclProc.$(ObjSuf) \
	tmp/tcl/tclExecute.$(ObjSuf) \
	tmp/tcl/tclCmdMZ.$(ObjSuf)

ExRootAnalysis/ExRootModule.h: \
	ExRootAnalysis/ExRootTask.h
	@touch $@

CLHEP/Vector/BoostX.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

CLHEP/Vector/EulerAngles.h: \
	CLHEP/Vector/defs.h
	@touch $@

KtJet/KtDistanceInterface.h: \
	KtJet/KtUtil.h
	@touch $@

CLHEP/Vector/BoostY.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

CLHEP/Vector/RotationX.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h
	@touch $@

CLHEP/Vector/Rotation.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/RotationX.h \
	CLHEP/Vector/RotationY.h \
	CLHEP/Vector/RotationZ.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

CLHEP/Vector/RotationY.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h
	@touch $@

CLHEP/Vector/BoostZ.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

CLHEP/Vector/LorentzVector.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h
	@touch $@

KtJet/KtUtil.h: \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

modules/MadGraphKtJetFinder.h: \
	ExRootAnalysis/ExRootModule.h \
	KtJet/KtLorentzVector.h
	@touch $@

CLHEP/Vector/RotationZ.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h
	@touch $@

CLHEP/Vector/ThreeVector.h: \
	CLHEP/Vector/defs.h
	@touch $@

modules/MadGraphClassFilter.h: \
	ExRootAnalysis/ExRootModule.h
	@touch $@

CLHEP/Vector/AxisAngle.h: \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/defs.h
	@touch $@

KtJet/KtLorentzVector.h: \
	KtJet/KtUtil.h
	@touch $@

ExRootAnalysis/ExRootAnalysis.h: \
	ExRootAnalysis/ExRootTask.h
	@touch $@

CLHEP/Vector/TwoVector.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h
	@touch $@

KtJet/KtEvent.h: \
	KtJet/KtUtil.h \
	KtJet/KtDistanceInterface.h \
	KtJet/KtJetTable.h \
	KtJet/KtRecomInterface.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

CLHEP/Vector/LorentzRotation.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/Rotation.h \
	CLHEP/Vector/Boost.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

KtJet/KtDistance.h: \
	KtJet/KtDistanceInterface.h \
	KtJet/KtUtil.h \
	KtJet/KtLorentzVector.h
	@touch $@

KtJet/KtJetTable.h: \
	KtJet/KtLorentzVector.h \
	KtJet/KtDistanceInterface.h \
	KtJet/KtRecomInterface.h
	@touch $@

CLHEP/Vector/Boost.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/RotationInterfaces.h \
	CLHEP/Vector/BoostX.h \
	CLHEP/Vector/BoostY.h \
	CLHEP/Vector/BoostZ.h \
	CLHEP/Vector/LorentzVector.h
	@touch $@

modules/MadGraphConeJetFinder.h: \
	ExRootAnalysis/ExRootModule.h \
	CDFCones/PhysicsTower.hh \
	CDFCones/Cluster.hh
	@touch $@

KtJet/KtRecom.h: \
	KtJet/KtUtil.h \
	KtJet/KtLorentzVector.h \
	KtJet/KtRecomInterface.h
	@touch $@

CLHEP/Units/PhysicalConstants.h: \
	CLHEP/Units/defs.h \
	CLHEP/Units/SystemOfUnits.h
	@touch $@

KtJet/KtRecomInterface.h: \
	KtJet/KtUtil.h
	@touch $@

modules/MadGraphAnalysis.h: \
	ExRootAnalysis/ExRootModule.h
	@touch $@

ExRootAnalysis/ExRootTask.h: \
	ExRootAnalysis/ExRootConfReader.h
	@touch $@

CLHEP/Vector/RotationInterfaces.h: \
	CLHEP/Vector/defs.h \
	CLHEP/Vector/ThreeVector.h \
	CLHEP/Vector/LorentzVector.h \
	CLHEP/Vector/AxisAngle.h
	@touch $@



###

all: $(SHARED) $(EXECUTABLE) $(STDHEP_EXECUTABLE)

static: $(STATIC)

$(STATIC): $(PGS_DICT_OBJ) $(PGS_SOURCE_OBJ) $(PGS_OBJ)
	@mkdir -p $(@D)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

$(SHARED): $(DICT_OBJ) $(SOURCE_OBJ) $(TCL_OBJ)
	@mkdir -p $(@D)
	@echo ">> Building $@"
ifeq ($(ARCH),aix5)
	@$(MAKESHARED) $(OutPutOpt) $@ $(LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	@$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	@ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
	@bindexplib $* $^ > $*.def
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(SHAREDLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
	@$(MT_DLL)
endif
endif
endif

clean:
	@rm -f $(PGS_DICT_OBJ) $(PGS_SOURCE_OBJ) $(PGS_OBJ) $(DICT_OBJ) $(SOURCE_OBJ) $(TCL_OBJ) $(STDHEP_OBJ) core

distclean: clean
	@rm -f $(SHARED) $(STATIC) $(EXECUTABLE) $(STDHEP_EXECUTABLE)

###

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

%Dict.$(SrcSuf):
	@mkdir -p $(@D)
	@echo ">> Generating $@"
	@rootcint -f $@ -c $<
	@echo "#define private public" > $@.arch
	@echo "#define protected public" >> $@.arch
	@mv $@ $@.base
	@cat $@.arch $< $@.base > $@
	@rm $@.arch $@.base

$(SOURCE_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(PGS_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(PGS_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(TCL_OBJ): tmp/%.$(ObjSuf): %.c
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@gcc $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(STDHEP_OBJ): tmp/%.$(ObjSuf): %.c
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@gcc $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(STDHEP_EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(STDHEP_EXECUTABLE): %$(ExeSuf): $(DICT_OBJ) $(SOURCE_OBJ) $(TCL_OBJ) $(STDHEP_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@

$(EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE): %$(ExeSuf): $(DICT_OBJ) $(SOURCE_OBJ) $(TCL_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@

###


