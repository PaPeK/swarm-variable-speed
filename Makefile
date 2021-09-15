#HOME	= /home/prom
TARGET	= swarmdyn

ARCH	= $(shell uname -m)
OS		= $(shell uname -s)
LIB_LOC	:= lib/lib64/
LIB		= /lib64/

# Flag explanation:
# -O3 = highest code optimization level
# -Wall = shows much more warnings than normally
# -Wextra = some extra waringns (also -ansi - pedantic)
# -g = debug mode -> retain symbol information in executable
#	g0: no debug info, g1:minimal debug info, g:default debug info, g3:max
C++	= h5c++
CXXFLAGS	= -O3 -Wall -std=c++14

LINKER	= h5c++
LFLAGS	= -lgsl -lgslcblas -lm -lgmp -lboost_system
ifeq ($(OS), Linux)
	LFLAGS	+= -lboost_thread # -lCGAL 
endif
ifeq ($(OS), Darwin)
	LFLAGS	+= -lboost_thread-mt
endif

# defines directories for specific file-type
SRCDIR	= src
OBJDIR	= obj
BINDIR	= .

SRC = $(wildcard $(SRCDIR)/*.cpp)
INC = $(wildcard $(SRCDIR)/*.h)
OBJ = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

$(BINDIR)/$(TARGET): $(OBJ)
	@$(LINKER) $(OBJ) $(LFLAGS) -o $@
	@echo "Linking complete"

$(OBJ): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(C++) $(CXXFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully"

cl:
	rm -f out_*.dat *.bin *.h5
	@echo "removed data-files"

clean: cl
	@rm -f $(BINDIR)/$(TARGET)
	@rm -f $(OBJ)
	@echo "Exec, and Objects removed"
