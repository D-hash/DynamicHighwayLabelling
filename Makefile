CXX = g++ -std=c++17
DBG = -g
OPT = -Ofast -DNDEBUG -march=native
VALGRIND = -g -DNDEBUG


OPTIONS = -lnetworkit -lroutingkit -lboost_serialization -lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lboost_timer

INCLUDEPATH = /usr/include/routingkit
PATHLIB = $(HOME)/Ricerca/RoutingKit/lib

TARGETS = shifting_tests
OTHERS = highway_cover_labelling.h libraries.h


all:
	$(foreach var,$(TARGETS),$(CXX) $(DBG) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS);)
debug:
	$(foreach var,$(TARGETS),$(CXX) $(DBG) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS) -I$(INCLUDEPATH) -L$(PATHLIB);)
release:
	$(foreach var,$(TARGETS),$(CXX) $(OPT) -o $(var) $(var).cpp $(OTHERS) $(OPTIONS) -I$(INCLUDEPATH) -L$(PATHLIB);)
clean:
	$(foreach var,$(TARGETS),rm -rf $(var);)

