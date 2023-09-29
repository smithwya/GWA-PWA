TARGET = GWA
# compilers
CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

#compiling flags
CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR) -I/usr/include/eigen3
CPPFLAGS += -O3

#linking flags
LDFLAGS := -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
HEADERS := $(filter-out src/Linkdef.h,$(INCLUDES))
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm = rm -f

all:  libGWA.so $(BINDIR)/$(TARGET)

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
	@echo "Cleanup complete"

GWADict.cxx: $(HEADERS) src/Linkdef.h
	rootcling -f $@ -c -p $^

libGWA.so: GWADict.cxx $(SOURCES)
	g++ -shared -o$@ $(CPPFLAGS) -fPIC -I$(ROOTSYS)/include $^ `root-config --ldflags --libs`

