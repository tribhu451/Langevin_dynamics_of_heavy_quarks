
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

GSLLIBS      := $(shell gsl-config --libs)

#EXTRA_FLAGS = -D SIMPLE # EoS p=e/3
EXTRA_FLAGS  = -D TABLE # Laine EoS, tabulated

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3 -march=native
LD            = g++
LDFLAGS       = -O3 -march=native

CXXFLAGS     += $(ROOTCFLAGS) $(EXTRA_FLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) $(GSLLIBS)

vpath %.cpp src
objdir     = obj

SRC        = main.cpp bulk.cpp langevin.cpp random.cpp initialisation.cpp hadronise.cpp input_params.cpp

             
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = langevin
#-------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  -fopenmp  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

cfiles : 
	rm -rf outputs/*

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)
	
obj/%.o : %.cpp
	$(CXX) -fopenmp $(CXXFLAGS) -c $< -o $@


