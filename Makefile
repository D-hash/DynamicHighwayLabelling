CXX = g++ -std=c++17
DBG = -g
OPT = -Ofast -DNDEBUG -march=native
VALGRIND = -g -DNDEBUG


OPTIONS = -lnetworkit -lboost_serialization -lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lboost_timer

INCLUDEPATH = /usr/include/valgrind
PATHLIB = $(HOME)/networkit/build/

TARGETS = testing_dynamic
OTHERS = highway_cover_labelling.h libraries.h


all:
	$(foreach var,$(TARGETS),$(CXX) $(DBG) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS);)
valgrind:
	$(foreach var,$(TARGETS),$(CXX) $(VALGRIND) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS);)
debug:
	$(foreach var,$(TARGETS),$(CXX) $(DBG) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS);)
release:
	$(foreach var,$(TARGETS),$(CXX) $(OPT) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS);)
clean:
	$(foreach var,$(TARGETS),rm -rf $(var);)

