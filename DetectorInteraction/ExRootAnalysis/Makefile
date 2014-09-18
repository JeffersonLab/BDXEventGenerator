CC = g++

ROOT_INC = `root-config --cflags`

INC_PATH = -I./include -I../tcl/include

SRC=$(wildcard src/*.cc)
OBJ_DIR=obj
OBJS=$(patsubst src/%,$(OBJ_DIR)/%,$(SRC:.cc=.o))

all: src/ExRootAnalysisDict.cc  $(OBJS)

$(OBJ_DIR)/%.o: src/%.cc
	$(CC) -c -o $@ $^ $(ROOT_INC) $(INC_PATH)

src/ExRootAnalysisDict.cc:
	rootcint -f $@ -c include/ExRootAnalysisLinkDef.h
	mv src/ExRootAnalysisDict.h ./include
	mv src/ExRootAnalysisDict.cc src/ExRootAnalysisDict.cc.base	
	echo "#define private public" > src/ExRootAnalysisDict.cc.arch
	echo "#define protected public" >> src/ExRootAnalysisDict.cc.arch
	cat src/ExRootAnalysisDict.cc.arch include/ExRootAnalysisLinkDef.h src/ExRootAnalysisDict.cc.base > src/ExRootAnalysisDict.cc
	rm src/ExRootAnalysisDict.cc.arch src/ExRootAnalysisDict.cc.base

clean:
	rm -f $(OBJ_DIR)/*.o
	rm src/ExRootAnalysisDict.cc