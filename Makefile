TARGET = GWA
# compilers
CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

#compiling flags
CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR) -I$(EIGEN)
CPPFLAGS += -O3
CPPFLAGS += -I/opt/homebrew/opt/libomp/include
#CPPFLAGS += -g
#CPPFLAGS += -pg
#CPPFLAGS += -march=native

#linking flags
LDFLAGS := -Xlinker -rpath . $(shell root-config --glibs ) $(STDLIBDIR) -L/opt/homebrew/opt/libomp/lib

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

rm = rm -f

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(LD) $(CPPFLAGS) -o $(BINDIR)/$(TARGET) $(OBJECTS) $(LDFLAGS)
	@echo "Linking complete"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
	$(CXX) $(CPPFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully"
clean:
	$(rm) $(OBJECTS)
	$(rm) libGWA.so  GWADict.cxx GWADict_rdict.pcm
	@echo "Cleanup complete"

