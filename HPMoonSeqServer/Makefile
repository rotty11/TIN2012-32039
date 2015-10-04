## @file Makefile
## @author Juan José Escobar Pérez
## @date 15/06/2015
## @brief File with the necessary rules for compile the project. Also, documentation is generated

SRC = src
INC = include
OBJ = obj
BIN = bin
LIB = lib
DOC = doc
GNUPLOT = gnuplot

COMP = g++
CXXFLAGS = -c -Iinclude -Wno-unused-result
OPT = -O3 -funroll-loops -ftree-vectorize -march=native
#OPT =

POPSIZE = -D POPULATION_SIZE=$(POPULATION_SIZE)
NFEATURES = -D N_FEATURES=$(N_FEATURES)
NOBJECTIVES = -D N_OBJECTIVES=$(N_OBJECTIVES)
NINSTANCES = -D N_INSTANCES=$(N_INSTANCES)
VARS = $(POPSIZE) $(NFEATURES) $(NOBJECTIVES) $(NINSTANCES)

INCLUDES = $(INC)/tinyxml2.h $(INC)/xml.h $(INC)/bd.h $(INC)/initialization.h $(INC)/evaluation.h $(INC)/sort.h $(INC)/tournament.h $(INC)/crossover.h
OBJECTS = $(OBJ)/tinyxml2.o $(OBJ)/xml.o $(OBJ)/bd.o $(OBJ)/individual.o $(OBJ)/initialization.o $(OBJ)/evaluation.o $(OBJ)/sort.o $(OBJ)/tournament.o $(OBJ)/crossover.o $(OBJ)/main.o

# ************ Targets ************

all: $(LIB)/libhv.a $(BIN)/hpmoonSeq

# ************ Documentation ************

documentation:
	doxygen $(DOC)/Doxyfile

# ************ Compilation of modules ************

$(OBJ)/tinyxml2.o: $(SRC)/tinyxml2.cpp $(INC)/tinyxml2.h
	$(COMP) $(CXXFLAGS) $(OPT) $(SRC)/tinyxml2.cpp -o $(OBJ)/tinyxml2.o

$(OBJ)/xml.o: $(SRC)/xml.cpp $(INC)/xml.h
	$(COMP) $(CXXFLAGS) $(OPT) $(SRC)/xml.cpp -o $(OBJ)/xml.o
$(OBJ)/bd.o: $(SRC)/bd.cpp $(INC)/bd.h
	$(COMP) $(CXXFLAGS) $(OPT) $(SRC)/bd.cpp -o $(OBJ)/bd.o
$(OBJ)/individual.o: $(SRC)/individual.cpp $(INC)/individual.h
	$(COMP) $(CXXFLAGS) $(NFEATURES) $(NOBJECTIVES) $(OPT) $(SRC)/individual.cpp -o $(OBJ)/individual.o
$(OBJ)/initialization.o: $(SRC)/initialization.cpp $(INC)/initialization.h $(INC)/individual.h
	$(COMP) $(CXXFLAGS) $(NFEATURES) $(NOBJECTIVES) $(OPT) $(SRC)/initialization.cpp -o $(OBJ)/initialization.o
$(OBJ)/evaluation.o: $(SRC)/evaluation.cpp $(INC)/evaluation.h $(INC)/individual.h $(INC)/hv.h
	$(COMP) $(CXXFLAGS) $(NFEATURES) $(NOBJECTIVES) $(NINSTANCES) $(OPT) $(SRC)/evaluation.cpp -o $(OBJ)/evaluation.o
$(OBJ)/sort.o: $(SRC)/sort.cpp $(INC)/sort.h $(INC)/individual.h
	$(COMP) $(CXXFLAGS) $(NFEATURES) $(NOBJECTIVES) $(OPT) $(SRC)/sort.cpp -o $(OBJ)/sort.o
$(OBJ)/tournament.o: $(SRC)/tournament.cpp $(INC)/tournament.h
	$(COMP) $(CXXFLAGS) $(OPT) $(SRC)/tournament.cpp -o $(OBJ)/tournament.o
$(OBJ)/crossover.o: $(SRC)/crossover.cpp $(INC)/crossover.h $(INC)/individual.h
	$(COMP) $(CXXFLAGS) $(NFEATURES) $(NOBJECTIVES) $(OPT) $(SRC)/crossover.cpp -o $(OBJ)/crossover.o

$(OBJ)/main.o: $(SRC)/main.cpp $(INCLUDES)
	$(COMP) $(CXXFLAGS) $(VARS) $(OPT) $(SRC)/main.cpp -o $(OBJ)/main.o

# ************ Libraries ************

$(LIB)/libhv.a:
	@cd hv-1.3-Fonseca; $(MAKE) -s
	@\cp hv-1.3-Fonseca/libhv.a $(LIB)

# ************ Linking and creating executable ************

$(BIN)/hpmoonSeq: $(OBJECTS) $(LIB)/libhv.a
	$(COMP) $(OBJECTS) -o $(BIN)/hpmoonSeq -Llib -lhv

# ************ Cleaning ***************

clean:
	@echo "Partial cleaning..."
	@printf "\t- Binary files\n"
	@printf "\t- .o files\n"
	@printf "\t- ~ files\n"
	@printf "\t- .a libraries\n"
	@printf "\t- Hypervolume project of Fonseca\n"
	@\rm -rf $(OBJ)/* $(BIN)/* $(LIB)/* *~
	@cd hv-1.3-Fonseca; $(MAKE) -s clean

eraseAll: clean
	@echo "Additionally..."
	@printf "\t- gnuplot files\n"
	@printf "\t- Documentation files\n"
	@\rm -rf $(GNUPLOT)/*
	@\rm -rf $(DOC)/html/* $(DOC)/*.db