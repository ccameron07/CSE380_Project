EXEC :=   main.out 
SRC  :=   main.cpp elements.cpp solvers.cpp masa_helper.cpp
OBJ  :=   $(patsubst %.cpp,%.o,$(SRC))
# Options
INCL     := -I$(TACC_GRVY_INC) -I$(TACC_MASA_INC) -I$(HOME)/CSE380/CSE380_Project/FEM/include
CXX      := g++
CXXFLAGS := -O0 -std=c++11 $(INCL)
LINK.o   := $(LINK.cc)
LDFLAGS  := 
LDLIBS   := -L$(TACC_GRVY_LIB) -lgrvy -L$(TACC_MASA_LIB) -lmasa 
# Rules
$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) -o $@ $^

main.o: elements.o solvers.o masa_helper.o

# Useful phony targets
.PHONY: clobber clean neat echo
clobber: clean
	$(RM) $(EXEC)
clean: neat
	$(RM) $(OBJ)
neat:
	$(RM) *~ .*~
echo:
	@echo $(OBJ)
