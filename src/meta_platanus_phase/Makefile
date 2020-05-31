CXX = g++
CXXFLAGS = -std=c++0x -O3 -funroll-loops -Wall -Wno-sign-compare -fopenmp -finline-limit-50000 -lm -Dnullptr=0


PRG = meta_platanus_phase
OBJ = main.o assemble.o scaffold.o scaffoldGraph.o gapClose.o common.o baseCommand.o seqlib.o mapper.o gapCloseOLC.o merge.o divide.o iterate.o polish.o solveDBG.o pairedDBG.o phase.o consensus.o mgc.o cluster.o clusterScaffold.o


all: $(PRG)

$(PRG): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^
.cpp.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(PRG) $(OBJ)

