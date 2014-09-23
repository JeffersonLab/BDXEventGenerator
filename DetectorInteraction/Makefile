CC = g++
MAIN = DetectorInteraction.exe

SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,obj/%.o,$(SRC))

ROOT_INC = `root-config --cflags`
ROOT_LIB = `root-config --glibs`
INC_PATH = -I./ -I./tcl/include -I./ExRootAnalysis/include -I./include
OBJ_TCL = $(patsubst tcl/src/%.c,tcl/obj/%.o,$(wildcard tcl/src/*.c))
OBJ_EXROOTANALYSIS = $(patsubst ExRootAnalysis/src/%.cc,ExRootAnalysis/obj/%.o,$(wildcard ExRootAnalysis/src/*.cc))
OBJ_EXROOTANALYSIS += ExRootAnalysis/obj/ExRootAnalysisDict.o

all: $(OBJ)
	$(MAKE) -C tcl    
	$(MAKE) -C ExRootAnalysis
	$(CC) -o $(MAIN) $(OBJ_TCL) $(OBJ_EXROOTANALYSIS) $(OBJ) $(OBJ_DICT) $(ROOT_LIB) 

obj/%.o: %.cpp
	$(CC) $(ROOT_INC) $(INC_PATH) -c -o $@ $^ 

clean:
	$(MAKE) clean -C tcl
	$(MAKE) clean -C ExRootAnalysis
	rm -f $(MAIN)
	rm -f obj/*.o