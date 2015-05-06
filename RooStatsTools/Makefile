INCDIR  = include
SRCDIR  = src
OBJDIR  = obj
LIBDIR  = lib
BINDIR  = bin
TESTDIR = test
LIBRARY = CGStatistics

LIBFILE  = $(LIBDIR)/lib$(LIBRARY).so
DICTHEAD = $(SRCDIR)/$(LIBRARY)_Dict.h
DICTCPP  = $(SRCDIR)/$(LIBRARY)_Dict.cxx
DICTOBJ  = $(OBJDIR)/$(LIBRARY)_Dict.o

HLIST = $(wildcard $(INCDIR)/*.h)
CPPLIST = $(filter-out $(DICTCPP),$(wildcard $(SRCDIR)/*.cxx))
OLIST = $(patsubst %.cxx,%.o,$(notdir $(CPPLIST)))

#CXXFLAGS += -I./ -Iinclude/ -Wall -Wno-overloaded-virtual -Wno-unused -pthread -Wno-deprecated-declarations -m64 -I$(shell root-config --incdir) -fPIC -DCG_EXPERIMENTAL
#CXX11    += --std=c++11
CXXFLAGS += -I./ -Iinclude/ -fPIC -DCG_EXPERIMENTAL $(shell root-config --cflags)

LIBS     += $(shell root-config --libs) -lRooFitCore -lRooFit -lRooStats 
LDFLAGS  += $(shell root-config --ldflags)

VPATH += $(SRCDIR) $(OBJDIR) $(TESTDIR)

all: $(LIBFILE) tests

%.o: %.cxx
	@echo "compiling $<"
	@mkdir -p $(OBJDIR)
	@g++ $(CXXFLAGS) -c -o $(OBJDIR)/$(notdir $@) $<

$(DICTCPP): $(HLIST)
	@echo "generating RootCint dictionary"
	rootcint -v -f $@ $(CXXFLAGS) $^ 

$(DICTOBJ): $(DICTCPP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@g++ $(CXXFLAGS) -c -o $@ $<

$(LIBFILE): $(OLIST) $(DICTOBJ)
	@echo "Making shared library: $(LIBFILE)"
	@rm -f $(LIBFILE)
	@mkdir -p $(LIBDIR)
	@g++ $(LDFLAGS) $(LIBS) -shared -fPIC $(addprefix $(OBJDIR)/,$(OLIST)) $(DICTOBJ) -o $@

.PHONY: tests
tests: GaussLimitPlot

.PHONY: GaussLimitPlot
GaussLimitPlot: GaussLimitPlot.o $(LIBFILE)
	@echo "creating test for limits on Gaussian problem"
	@mkdir -p $(BINDIR)
	@g++ $(addprefix $(OBJDIR)/,$(notdir $<)) $(LIBS) -L$(shell pwd)/$(LIBDIR) -l$(LIBRARY) -o $(BINDIR)/GaussLimitPlot

.PHONY: clean
clean:
	rm -f $(LIBFILE)
	rm -f $(DICTHEAD) $(DICTCPP)
	rm -rf $(OBJDIR)
	rm -f $(BINDIR)/*
