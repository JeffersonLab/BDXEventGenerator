
#
# Makefile for ExRootAnalysis
#
# Author: P. Demin - UCL, Louvain-la-Neuve
#
# multi-platform configuration is taken from ROOT (root/test/Makefile)
#

include doc/Makefile.arch

ROOT_MAJOR := $(shell $(RC) --version | cut -d'.' -f1)

SrcSuf = cc
PcmSuf = _rdict.pcm

CXXFLAGS += $(ROOTCFLAGS) -DDROP_CGAL -I.
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)
	
###

SHARED = libExRootAnalysis.$(DllSuf)
SHAREDLIB = libExRootAnalysis.lib

all:


ExRootLHEFConverter$(ExeSuf): \
	tmp/test/ExRootLHEFConverter.$(ObjSuf)

tmp/test/ExRootLHEFConverter.$(ObjSuf): \
	test/ExRootLHEFConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootLHEFReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootProgressBar.h
ExRootLHCOlympicsConverter$(ExeSuf): \
	tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf)

tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf): \
	test/ExRootLHCOlympicsConverter.cpp \
	ExRootAnalysis/ExRootStream.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
	ExRootAnalysis/ExRootProgressBar.h
ExRootSTDHEPConverter$(ExeSuf): \
	tmp/test/ExRootSTDHEPConverter.$(ObjSuf)

tmp/test/ExRootSTDHEPConverter.$(ObjSuf): \
	test/ExRootSTDHEPConverter.cpp \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootSTDHEPReader.h \
	ExRootAnalysis/ExRootTreeWriter.h \
	ExRootAnalysis/ExRootTreeBranch.h \
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
	ExRootLHEFConverter$(ExeSuf) \
	ExRootLHCOlympicsConverter$(ExeSuf) \
	ExRootSTDHEPConverter$(ExeSuf) \
	ExRootHEPEVTConverter$(ExeSuf) \
	Example$(ExeSuf)

EXECUTABLE_OBJ =  \
	tmp/test/ExRootLHEFConverter.$(ObjSuf) \
	tmp/test/ExRootLHCOlympicsConverter.$(ObjSuf) \
	tmp/test/ExRootSTDHEPConverter.$(ObjSuf) \
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
	ExRootAnalysis/ExRootFactory.h
ExRootAnalysisDict$(PcmSuf): \
	tmp/src/ExRootAnalysisDict$(PcmSuf) \
	tmp/src/ExRootAnalysisDict.$(SrcSuf)
DICT +=  \
	tmp/src/ExRootAnalysisDict.$(SrcSuf)

DICT_OBJ +=  \
	tmp/src/ExRootAnalysisDict.$(ObjSuf)

DICT_PCM +=  \
	ExRootAnalysisDict$(PcmSuf)

tmp/src/ExRootFilter.$(ObjSuf): \
	src/ExRootFilter.$(SrcSuf) \
	ExRootAnalysis/ExRootFilter.h \
	ExRootAnalysis/ExRootClassifier.h
tmp/src/ExRootSTDHEPReader.$(ObjSuf): \
	src/ExRootSTDHEPReader.$(SrcSuf) \
	ExRootAnalysis/ExRootSTDHEPReader.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootFactory.h \
	ExRootAnalysis/ExRootTreeBranch.h
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
tmp/src/ExRootStream.$(ObjSuf): \
	src/ExRootStream.$(SrcSuf) \
	ExRootAnalysis/ExRootStream.h
tmp/src/ExRootClasses.$(ObjSuf): \
	src/ExRootClasses.$(SrcSuf) \
	ExRootAnalysis/ExRootClasses.h
tmp/src/ExRootTreeBranch.$(ObjSuf): \
	src/ExRootTreeBranch.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeBranch.h
tmp/src/ExRootTreeReader.$(ObjSuf): \
	src/ExRootTreeReader.$(SrcSuf) \
	ExRootAnalysis/ExRootTreeReader.h
tmp/src/ExRootProgressBar.$(ObjSuf): \
	src/ExRootProgressBar.$(SrcSuf) \
	ExRootAnalysis/ExRootProgressBar.h
tmp/src/ExRootLHEFReader.$(ObjSuf): \
	src/ExRootLHEFReader.$(SrcSuf) \
	ExRootAnalysis/ExRootLHEFReader.h \
	ExRootAnalysis/ExRootClasses.h \
	ExRootAnalysis/ExRootFactory.h \
	ExRootAnalysis/ExRootStream.h \
	ExRootAnalysis/ExRootTreeBranch.h
tmp/src/ExRootUtilities.$(ObjSuf): \
	src/ExRootUtilities.$(SrcSuf) \
	ExRootAnalysis/ExRootUtilities.h
SOURCE_OBJ =  \
	tmp/src/ExRootFilter.$(ObjSuf) \
	tmp/src/ExRootSTDHEPReader.$(ObjSuf) \
	tmp/src/ExRootFactory.$(ObjSuf) \
	tmp/src/ExRootTreeWriter.$(ObjSuf) \
	tmp/src/ExRootResult.$(ObjSuf) \
	tmp/src/ExRootStream.$(ObjSuf) \
	tmp/src/ExRootClasses.$(ObjSuf) \
	tmp/src/ExRootTreeBranch.$(ObjSuf) \
	tmp/src/ExRootTreeReader.$(ObjSuf) \
	tmp/src/ExRootProgressBar.$(ObjSuf) \
	tmp/src/ExRootLHEFReader.$(ObjSuf) \
	tmp/src/ExRootUtilities.$(ObjSuf)



###

ifeq ($(ROOT_MAJOR),6)
all: $(SHARED) $(DICT_PCM) $(EXECUTABLE)
else
all: $(SHARED) $(EXECUTABLE)
endif

$(SHARED): $(DICT_OBJ) $(SOURCE_OBJ)
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
	@rm -f $(DICT_OBJ) $(SOURCE_OBJ) core
	@rm -rf tmp

distclean: clean
	@rm -f $(SHARED) $(SHAREDLIB) $(DICT_PCM) $(EXECUTABLE)

###

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf) $(PcmSuf)

%Dict.$(SrcSuf):
	@mkdir -p $(@D)
	@echo ">> Generating $@"
	@rootcint -f $@ -c $<
	@echo "#define private public" > $@.arch
	@echo "#define protected public" >> $@.arch
	@mv $@ $@.base
	@cat $@.arch $< $@.base > $@
	@rm $@.arch $@.base

%Dict$(PcmSuf):
	@echo ">> Copying $@"
	@cp $< $@

$(SOURCE_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE): %$(ExeSuf): $(DICT_OBJ) $(SOURCE_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@

###


