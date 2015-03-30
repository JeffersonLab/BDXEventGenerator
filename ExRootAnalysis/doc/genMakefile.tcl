#!/usr/bin/env tclsh

set prefix "tmp/"
set suffix " \\\n\t"

set srcSuf {.$(SrcSuf)}
set objSuf {.$(ObjSuf)}
set pcmSuf {$(PcmSuf)}
set exeSuf {$(ExeSuf)}

proc dependencies {fileName firstLine {force 1} {command {}}} {
  global suffix headerFiles sourceFiles

  if {[info exists sourceFiles($fileName)]} return

  set sourceFiles($fileName) 1

  set list {}
  set fid [open $fileName]
  while {! [eof $fid]} {
    set line [gets $fid]
    if [regexp -- {^\s*#include\s*"((\w+/)+\w+\.(h|hh))"} $line] {
      set elements [split $line {"}]
      set file [lindex $elements 1]
      if [file exists $file] {
        lappend list $file
        set headerFiles($file) 1
      }
    }
  }
  
  if {[llength $list] > 0} {
    puts -nonewline $firstLine
    foreach file $list {puts -nonewline $suffix$file}
    if {$command != {}} {
      puts {}
      puts $command
    }
    puts {}
  } elseif {$force} {
    puts -nonewline $firstLine
    if {$command != {}} {
      puts {}
      puts $command
    }
    puts {}
  }

  close $fid  
}

proc dictDeps {dictVar args} {

  global prefix suffix srcSuf objSuf pcmSuf

  set dict [eval glob -nocomplain $args]
  
  set dictSrcFiles {}
  set dictObjFiles {}

  foreach fileName $dict {
    regsub {LinkDef\.h} $fileName {Dict} dictName
    set dictName $prefix$dictName
  
    lappend dictSrcFiles $dictName$srcSuf
    lappend dictObjFiles $dictName$objSuf
    lappend dictPcmFiles [file tail $dictName$pcmSuf]
  
    dependencies $fileName "$dictName$srcSuf:$suffix$fileName"

    puts -nonewline [file tail $dictName$pcmSuf]:$suffix
    puts -nonewline $dictName$pcmSuf$suffix
    puts -nonewline $dictName$srcSuf
    puts {}
  }
  
  puts -nonewline "${dictVar} += $suffix"
  puts [join $dictSrcFiles $suffix]
  puts {}

  puts -nonewline "${dictVar}_OBJ += $suffix"
  puts [join $dictObjFiles $suffix]
  puts {}

  puts -nonewline "${dictVar}_PCM += $suffix"
  puts [join $dictPcmFiles $suffix]
  puts {}
}

proc sourceDeps {srcPrefix args} {

  global prefix suffix srcSuf objSuf
  
  set source [eval glob -nocomplain $args]
    
  set srcObjFiles {}
  
  foreach fileName $source {
    regsub {\.cc} $fileName {} srcName
    set srcObjName $prefix$srcName
  
    lappend srcObjFiles $srcObjName$objSuf
  
    dependencies $fileName "$srcObjName$objSuf:$suffix$srcName$srcSuf"
  }

  puts -nonewline "${srcPrefix}_OBJ = $suffix"
  puts [join $srcObjFiles $suffix]
  puts {}
}

proc executableDeps {} {

  global prefix suffix objSuf exeSuf
   
  set executable [glob -nocomplain {test/*.cpp}]
  
  set exeFiles {}
  
  foreach fileName $executable {
    regsub {\.cpp} $fileName {} exeObjName
    set exeObjName $prefix$exeObjName
    set exeName [file tail $exeObjName]

    lappend exeFiles $exeName$exeSuf
    lappend exeObjFiles $exeObjName$objSuf
    
    puts "$exeName$exeSuf:$suffix$exeObjName$objSuf"
    puts {}
  
    dependencies $fileName "$exeObjName$objSuf:$suffix$fileName"
  }
  
  if [info exists exeFiles] {
    puts -nonewline "EXECUTABLE = $suffix"
    puts [join $exeFiles $suffix]
    puts {}
  }
  if [info exists exeObjFiles] {
    puts -nonewline "EXECUTABLE_OBJ = $suffix"
    puts [join $exeObjFiles $suffix]
    puts {}
  }

}

proc headerDeps {} {
  global suffix headerFiles
    
  foreach fileName [array names headerFiles] {  
    dependencies $fileName "$fileName:" 0 "\t@touch \$@"
  }
}

puts {
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

}

executableDeps

dictDeps {DICT} {src/*LinkDef.h} {modules/*LinkDef.h}

sourceDeps {SOURCE} {src/*.cc} {modules/*.cc}

headerDeps

puts {

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

}
