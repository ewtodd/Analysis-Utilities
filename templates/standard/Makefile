CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -std=c++17 -Wall -Wextra -O2 -fPIC
ROOTLIBS = $(shell root-config --libs)

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
LIBDIR = lib

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)

# Library name
LIBNAME = lib-nm-toolkit
SHAREDLIB = $(LIBDIR)/$(LIBNAME).so
STATICLIB = $(LIBDIR)/$(LIBNAME).a

# Default target - build both shared and static libraries
all: $(SHAREDLIB) $(STATICLIB)

# Create directories
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

# Compile object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Create shared library
$(SHAREDLIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS)
	@echo "Built shared library: $(SHAREDLIB)"

# Create static library
$(STATICLIB): $(OBJECTS) | $(LIBDIR)
	ar rcs $@ $^
	@echo "Built static library: $(STATICLIB)"

# Install target (for manual installation, Nix handles this differently)
install: $(SHAREDLIB) $(STATICLIB)
	@echo "Use 'nix build' to install via Nix"

clean:
	rm -rf $(OBJDIR) $(LIBDIR)

.PHONY: all install clean
