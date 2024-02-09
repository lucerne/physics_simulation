#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

DEPSDIR := $(shell mkdir -p .deps; echo .deps)

# Get the shell name to determine the OS
UNAME := $(shell uname)

# define the C compiler to use
CC := gcc
ifeq ($(UNAME), Linux)
CXX := g++ -std=gnu++0x
endif
ifeq ($(UNAME), Darwin)
CXX := $(shell for i in 4.8 4.7; do if g++-mp-$$i -v >/dev/null 2>&1; then echo g++-mp-$$i; exit; fi; done; echo false) -std=gnu++0x
SDLOBJS := CS207/SDLMain.o
endif

# define any compile-time flags
CFLAGS := -g -O2 -W -Wall -Wextra #-pedantic
DEPCFLAGS = -MD -MF $(DEPSDIR)/$*.d -MP


# define any directories containing header files other than /usr/include
#   include directories like -Ipath/to/files
INCLUDES += -I.

# define any directories containing libraries other than /usr/lib
#   include directories like -Lpath/to/libraries
ifeq ($(UNAME), Linux)
LIBS := -lSDL -lGL -lGLU
endif
ifeq ($(UNAME), Darwin)
LIBS := -framework SDL -framework OpenGL -framework Cocoa
endif

# define any libraries to link into executable:
#   if I want to link in libraries (libXXX.so or libXXX.a) I use the -lXXX
#   option, something like (this will link in libm.so):
LIBS += -lm

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

# 'make' - default rule
all: shallow_water

# generic rule - make sure object files are up-to-date, then compile MAIN
viewer: viewer.o $(SDLOBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

mass_spring: mass_spring.o $(SDLOBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

subgraph: subgraph.o $(SDLOBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

shortest_path: shortest_path.o $(SDLOBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

shallow_water: shallow_water.o $(SDLOBJS)
	$(CXX) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

# suffix replacement rule for building .o's from .cpp's
#   $<: the name of the prereq of the rule (a .cpp file)
#   $@: the name of the target of the rule (a .o file)
.cpp.o:
	$(CXX) $(CFLAGS) $(DEPCFLAGS) $(DEFS) $(INCLUDES) -c -o $@ $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	$(RM) *.o primes viewer test_nodes test_edges mass_spring subgraph shortest_path shallow_water
	$(RM) -r $(DEPSDIR)

DEPFILES := $(wildcard $(DEPSDIR)/*.d)
ifneq ($(DEPFILES),)
include $(DEPFILES)
endif

# define rules that do not actually generate the corresponding file
.PHONY: clean all always
