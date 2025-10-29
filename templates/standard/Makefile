# ROOT Analysis Framework Makefile

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -fPIC
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
LIBDIR = lib

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cxx)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.h)

# Library name
LIBNAME = libWaveformAnalysis
SHAREDLIB = $(LIBDIR)/$(LIBNAME).so

# Default target
all: $(SHAREDLIB)

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# Compile object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cxx $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -I$(INCDIR) -c $< -o $@

# Create shared library
$(SHAREDLIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS)

initial: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing initial processing macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x InitialProcessing.C");'

calibrate: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing calibration macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Calibration.C");'

background: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing background subtraction macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Background.C");'

plots: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing main plotting macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x Plotting.C");'

PSD: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing pulse shape discrimination macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x PSD.C");'

optimize: $(SHAREDLIB)
	@export ROOT_INCLUDE_PATH="$(PWD)/$(INCDIR):$$ROOT_INCLUDE_PATH" && \
	 export LD_LIBRARY_PATH="$(PWD)/$(LIBDIR):$$LD_LIBRARY_PATH" && \
	 cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing PSD gate optimization macro..." && \
	 root -l -b -q -e 'gSystem->Load("../$(LIBDIR)/$(LIBNAME).so"); gROOT->ProcessLine(".x OptimizeGates.C");'

clean:
	rm -rf $(OBJDIR) $(LIBDIR)


