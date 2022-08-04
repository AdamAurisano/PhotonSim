SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cxx .o Dict.cxx .so .hh .h

.PHONY: all lib clean bin

ROOTCONFIG   := root-config

OPT2         = -g -O3
COMPAT_OPTS  = -m64 -march=native
#OPT2         = -g -O2
#COMPAT_OPTS  = -m32 -march=i686
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTINC      := $(shell $(ROOTCONFIG) --incdir)
ROOTCINT     := rootcint
ROOTVER      := $(shell $(ROOTCONFIG) --version | root-config --version | sed 's,\([0-9]*\).\([0-9]*\)/\([0-9]*\),\1\2\3,')

CXX          = g++
CXXFLAGS     = $(OPT2) $(COMPAT_OPTS) -Wall -fPIC -D_ROOTVER_=$(ROOTVER)
LD           = g++
LDFLAGS      = $(OPT2) $(COMPAT_OPTS) -flto
SOFLAGS      = -shared
OutPutOpt     = -o # keep whitespace after "-o"

CXXFLAGS     += $(ROOTCFLAGS) 
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lMathMore -lMinuit2 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#$(shell gsl-config --cflags)
#$(shell gsl-config --libs)

WORKDIR = $(shell pwd)
SRCDIR = src
INCDIR = include
LIBDIR = lib
OBJDIR = tmp
DICTDIR = dict
EXEDIR = drivers

CXXFLAGS  += -I$(WORKDIR)/$(INCDIR) -I$(WORKDIR) -I$(ROOTINC)

vpath %.cxx $(SRCDIR)
vpath %.hh $(INCDIR)

vpath %Dict.cxx $(DICTDIR)
vpath %_LinkDef.h $(DICTDIR)
vpath %Dict.h  $(DICTDIR)

vpath %.so $(LIBDIR)
vpath %.o  $(OBJDIR)
vpath %.d  $(OBJDIR)

STANDALONE_DEP= -MM $< > $(OBJDIR)/$(basename $(notdir $<)).d

define postprocess_d
test -f $(OBJDIR)/$(basename $(notdir $<)).d && \
cat $(OBJDIR)/$(basename $(notdir $<)).d | \
sed 's?$*\.o?$(OBJDIR)/$*.o ?g' > \
$(OBJDIR)/dep_tmp.$$$$ ; \
mv $(OBJDIR)/dep_tmp.$$$$ $(OBJDIR)/$(basename $(notdir $<)).d
endef

define cxx_generate_depends
$(CXX) $(CXXFLAGS) $(STANDALONE_DEP) || /bin/rm -f $(OBJDIR)/$(basename $(notdir $<)).d
endef

$(DICTDIR)/%Dict.cxx: $(INCDIR)/%.hh $(DICTDIR)/%_Linkdef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

$(OBJDIR)/%Dict.o: $(DICTDIR)/%Dict.cxx
	@echo Prereqs: $^
	$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	@echo Prereqs: $^
	$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(OBJDIR)/%.d: $(SRCDIR)/%.cxx 
	$(cxx_generate_depends)
	$(postprocess_d)

#creates list of dictionaries
DICTS := $(patsubst $(DICTDIR)/%_Linkdef.h, $(DICTDIR)/%Dict.cxx, $(wildcard $(DICTDIR)/*_Linkdef.h))
#creates list of all source files (normal cxx files)
SRCS  := $(wildcard $(SRCDIR)/*.cxx)
#creates dependencies list ofr normal cxx files
DEPS := $(patsubst $(SRCDIR)/%.cxx, $(OBJDIR)/%.d, $(SRCS))
#list of file which are part of executables
EXESRCS := $(SRCDIR)/GetFiberTime.cxx $(SRCDIR)/GetFiberTime_mp.cxx $(SRCDIR)/PhotonSim.cxx $(SRCDIR)/PhotonSim_mp.cxx
LIBSRCS := $(filter-out $(EXESRCS), $(SRCS))
#creates list of all object files
LIBOBJS  := $(patsubst $(SRCDIR)/%.cxx, $(OBJDIR)/%.o, $(LIBSRCS))
EXEOBJS  := $(patsubst $(SRCDIR)/%.cxx, $(OBJDIR)/%.o, $(EXESRCS))
DICT_OBJS := $(patsubst $(DICTDIR)/%.cxx, $(OBJDIR)/%.o, $(DICTS))

-include $(DEPS)

.SECONDARY: $(DICTS)

all:
	$(MAKE) lib
	$(MAKE) GetFiberTime.x
	$(MAKE) GetFiberTimeMP.x
	$(MAKE) PhotonSim.x
	$(MAKE) PhotonSim_mp.x

lib: $(DICT_OBJS) $(LIBOBJS)
	@echo $^
	$(LD) $(SOFLAGS)  $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$(WORKDIR)/$(LIBDIR)/libPhotonSim.so

GetFiberTime.x: $(OBJDIR)/GetFiberTime.o $(LIBDIR)/libPhotonSim.so
	@echo $^
	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(filter %.so, $^) $(OutPutOpt)$@

GetFiberTimeMP.x: $(OBJDIR)/GetFiberTime_mp.o $(LIBDIR)/libPhotonSim.so
	@echo $^
	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(filter %.so, $^) $(OutPutOpt)$@

PhotonSim.x: $(OBJDIR)/PhotonSim.o $(LIBDIR)/libPhotonSim.so
	@echo $^
	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(filter %.so, $^) $(OutPutOpt)$@

PhotonSim_mp.x: $(OBJDIR)/PhotonSim_mp.o $(LIBDIR)/libPhotonSim.so
	@echo $^
	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(filter %.so, $^) $(OutPutOpt)$@



#PhotonSim.x: $(OBJDIR)/PhotonSim.o $(DICT_OBJS) $(LIBOBJS)
#	@echo $^
#	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(OutPutOpt)$@ $(filter %.so, $^)	

#profileToyMCFits.x: $(OBJDIR)/profileToyMCFits.o $(LIBDIR)/libFitFramework.so
#	@echo $^
#	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(OutPutOpt)$@ $(filter %.so, $^)

#fitHisto.x: $(OBJDIR)/fitHisto.o $(LIBDIR)/libFitFramework.so
#	@echo $^
#	$(LD) $(LDFLAGS) $(filter %.o, $^) $(LIBS) $(OutPutOpt)$@ $(filter %.so, $^)	

clean:
	@rm -rf $(WORKDIR)/$(OBJDIR)/*.o $(WORKDIR)/$(OBJDIR)/*.d $(WORKDIR)/$(LIBDIR)/*.so $(WORKDIR)/$(DICTDIR)/*Dict.h $(WORKDIR)/$(DICTDIR)/*Dict.cxx *.x

