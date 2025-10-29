# ROOT Analysis Framework Makefile

CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -std=c++17 -Wall -Wextra -O2 -fPIC
ROOTLIBS = $(shell root-config --libs)

# Check if dsp-toolkit is available via pkg-config
ifneq ($(shell pkg-config --exists dsp-toolkit && echo yes),yes)
    $(warning dsp-toolkit not found via pkg-config, using fallback paths)
    TOOLKIT_CFLAGS = 
    TOOLKIT_LIBS = 
    TOOLKIT_LIBDIR = 
else
    TOOLKIT_CFLAGS = $(shell pkg-config --cflags dsp-toolkit)
    TOOLKIT_LIBS = $(shell pkg-config --libs dsp-toolkit)
    TOOLKIT_LIBDIR = $(shell pkg-config --variable=libdir dsp-toolkit)
endif

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
LIBNAME = lib-DSP-Toolkit
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
	$(CXX) $(CXXFLAGS) $(TOOLKIT_CFLAGS) -I$(INCDIR) -c $< -o $@

# Create shared library
$(SHAREDLIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS) $(TOOLKIT_LIBS)

# Helper function to run ROOT macros with proper library loading
define run_root_macro
	@cd macros && \
	 echo "Loading library: ../$(LIBDIR)/$(LIBNAME).so" && \
	 echo "Executing $(1) macro..." && \
	 root -l -b -q -e 'gSystem->AddIncludePath("-I$(PWD)/$(INCDIR)"); \
	                    gSystem->AddLinkedLibs(" $(PWD)/$(LIBDIR)/$(LIBNAME).so"); \
	                    $(if $(TOOLKIT_LIBDIR),gSystem->AddLinkedLibs(" $(TOOLKIT_LIBS)");,) \
	                    gSystem->Load("$(PWD)/$(LIBDIR)/$(LIBNAME).so"); \
	                    $(if $(TOOLKIT_LIBDIR),gSystem->Load("$(TOOLKIT_LIBDIR)/lib-DSP-Toolkit.so");,) \
	                    .L $(2).C+; $(2)();'
endef

initial: $(SHAREDLIB)
	$(call run_root_macro,initial processing,InitialProcessing)

calibrate: $(SHAREDLIB)
	$(call run_root_macro,calibration,Calibration)

background: $(SHAREDLIB)
	$(call run_root_macro,background subtraction,Background)

plots: $(SHAREDLIB)
	$(call run_root_macro,main plotting,Plotting)

PSD: $(SHAREDLIB)
	$(call run_root_macro,pulse shape discrimination,PSD)

optimize: $(SHAREDLIB)
	$(call run_root_macro,PSD gate optimization,OptimizeGates)

clean:
	rm -rf $(OBJDIR) $(LIBDIR)
	rm -f macros/*_C.so macros/*_C.d macros/*_C_ACLiC_dict_rdict.pcm

.PHONY: all initial calibrate background plots PSD optimize clean
