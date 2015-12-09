SRCDIR := src/
OBJDIR := build/
BINDIR := bin/
EXEC   := main
SRCS :=   main.cpp elements.cpp solvers.cpp masa_helper.cpp

SRC  :=   $(addprefix $(SRCDIR), $(SRCS))
OBJ  :=   $(addprefix $(OBJDIR),$(patsubst %.cpp,%.o,$(SRCS))) 
# Options
INCL     := -I$(TACC_GRVY_INC) -I$(TACC_MASA_INC) -I$(HOME)/CSE380/CSE380_Project/FEM/include
CXX      := g++
CXXFLAGS := -O0 -std=c++11 $(INCL)
LINK.o   := $(LINK.cc)
LDFLAGS  := 
LDLIBS   := -L$(TACC_GRV_LIB) -lgrvy -L$(TACC_MASA_LIB) -lmasa

# Rules

$(BINDIR)$(EXEC): $(OBJ)
	$(LINK.o) $(LDLIBS) $(OBJ) -o $@

#$(OBJ) : $(SRC)
$(OBJ) : $(OBJDIR)%.o : $(SRCDIR)%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Useful phony targets
.PHONY: clobber clean neat echo
clobber: clean
	$(RM) $(EXEC)
clean: neat
	$(RM) $(OBJ)
neat:
	$(RM) *~ .*~
echo:
	@echo $(SRC)
	@echo $(OBJ)
	@echo $(OBJ) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
