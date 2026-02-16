include Makefile.inc

DIRS = python
SRCDIR = src
INCLUDEDIR = include
ANADIR = analyzers
BINDIR = bin
FASTJET = fastjet-install/bin/fastjet-config

ANALYZERS = $(wildcard $(ANADIR)/*.cc)
ANALYZERSH = $(ANALYZERS:cc=h)
ANALYZERSOBJ = $(ANALYZERS:cc=o)
RUNNERS = $(addprefix $(BINDIR)/Run,$(notdir $(basename $(ANALYZERS))))
RUNNERSCC = $(addsuffix .cc,$(addprefix $(ANADIR)/,$(notdir $(RUNNERS))))
MDSNANO_RUNNER = $(BINDIR)/Runllp_MuonSystem_CA_mdsnano
RUNNERS_GENERIC = $(filter-out $(MDSNANO_RUNNER),$(RUNNERS))
UTILS =$(SRCDIR)/RazorHelper.cc $(SRCDIR)/JetCorrectorParameters.cc \
        $(SRCDIR)/SimpleJetCorrectionUncertainty.cc \
		$(SRCDIR)/JetCorrectionUncertainty.cc \
	       	$(SRCDIR)/CACluster.cc ${SRCDIR}/TreeMuonSystemCombination.cc
UTILSOBJ = $(UTILS:cc=o)
MDSNANO_UTILSOBJ = $(SRCDIR)/RazorHelper.o $(SRCDIR)/JetCorrectorParameters.o \
		$(SRCDIR)/SimpleJetCorrectionUncertainty.o $(SRCDIR)/JetCorrectionUncertainty.o \
		$(SRCDIR)/CACluster.o $(SRCDIR)/TreeMuonSystemCombination.o
EXECUTABLES = $(RUNNERS)
HELPERSCRIPT = python/MakeAnalyzerCode.py


.PHONY: clean all lxplus copy_runners

all: $(FASTJET) copy_runners $(EXECUTABLES)

lxplus: all

clean:
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o
copy_runners:
		@for d in $(subst Run,,$(notdir $(basename $(RUNNERSCC)))); do ( if [ ! -f "src/Run"$$d".cc" ]; then echo $$d" file does not exists, copying"; $(HELPERSCRIPT) $$d; fi ) ; done

$(BINDIR):
	@mkdir -p $(BINDIR)

$(FASTJET):
	$(ANADIR)/fastjet_install.sh

$(SRCDIR)/llp_event.o: $(SRCDIR)/llp_event.C $(INCLUDEDIR)/llp_event.h
	$(CXX) $(SRCDIR)/llp_event.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/merged_event.o: $(SRCDIR)/merged_event.C $(INCLUDEDIR)/merged_event.h
	$(CXX) $(SRCDIR)/merged_event.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
	
$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/merged_event.o $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(RUNNERS): | $(BINDIR)
$(RUNNERS_GENERIC): $(BINDIR)/Run%: $(SRCDIR)/merged_event.o $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/%.o $(SRCDIR)/Run%.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(MDSNANO_RUNNER): $(SRCDIR)/merged_event.o $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.o $(MDSNANO_UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_mdsnano.o $(SRCDIR)/Runllp_MuonSystem_CA_mdsnano.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
