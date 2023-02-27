CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)


CPPFLAGS := $(shell root-config --cflags) $(STDINCDIR) -I/usr/include/eigen3
LDFLAGS := -Xlinker -rpath . $(shell root-config --glibs) $(STDLIBDIR)
READERFLAGS := -lgsl -lgslcblas
#G_OPTS := -malign-double
#G_OPTSPOLE := -ffinite-math-only -fforce-addr -funroll-loops -ffast-math -malign-double -fno-trapping-math -fcaller-saves -fstrength-reduce -funsafe-math-optimizations -ffinite-math-only

#CPPFLAGS += -g
CPPFLAGS += -O3



OBJ = $(SRC:.C=.o)

all : GWA

GWA : GWA.o amplitude.o observable.o channel.o
	$(LD) $(CPPFLAGS) -o GWA  GWA.o amplitude.o observable.o channel.o $(LDFLAGS)


%.o : %.cpp %.h
	$(CXX) $(CPPFLAGS) -o $@ -c $<

clean :
	rm -f *.o 
