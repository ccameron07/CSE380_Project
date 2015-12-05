CC := icpc # This is the main compiler
 # CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
TARGET := bin/a.out

SRCEXT := cpp
SOURCES := $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g # -Wall
LIB := 
INC := -I include/

all:
$(TARGET): $(OBJECTS)
	@echo $(SOURCES)
	  @echo " Linking..."
	    @echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	  @mkdir -p $(BUILDDIR)
	    @echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Tests
tests: 	
	$(CC) $(CFLAGS) test/eigentest.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	  $(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean, list
clean:
	  @echo " Cleaning..."; 
	    @echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

list:
	@echo " $(SOURCES) ";

