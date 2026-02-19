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
ANALYZER_NAMES = $(notdir $(basename $(ANALYZERS)))
RUNNERS = $(addprefix $(BINDIR)/Run,$(notdir $(basename $(ANALYZERS))))
RUNNER_SRCS = $(addprefix $(SRCDIR)/Run,$(addsuffix .cc,$(ANALYZER_NAMES)))
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
RUNNER_TEMPLATES = $(INCLUDEDIR)/AnalyzerTemplate.txt \
		   $(INCLUDEDIR)/AnalyzerTemplate_mdsnano.txt \
		   $(SRCDIR)/RunAnalyzerTemplate.txt


.PHONY: clean all lxplus copy_runners

all: $(FASTJET) $(RUNNER_SRCS) $(EXECUTABLES)

lxplus: all

clean:
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o $(RUNNERS) $(SRCDIR)/Run*.cc
copy_runners: $(RUNNER_SRCS)

$(RUNNER_SRCS): $(SRCDIR)/Run%.cc: $(ANADIR)/%.cc $(HELPERSCRIPT) $(RUNNER_TEMPLATES)
	@if [ ! -f "$@" ]; then \
		echo "$* file does not exists, copying"; \
		$(HELPERSCRIPT) $*; \
	fi

$(BINDIR):
	@mkdir -p $(BINDIR)

$(FASTJET):
	$(ANADIR)/fastjet_install.sh

# Ensure FastJet is fully installed before any compilation/linking that
# expands fastjet-config flags.
$(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANALYZERSOBJ): | $(FASTJET)

$(SRCDIR)/mdsnano_event.o: $(SRCDIR)/mdsnano_event.cc $(INCLUDEDIR)/mdsnano_event.h
	$(CXX) $(SRCDIR)/mdsnano_event.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
	
$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(RUNNERS): | $(BINDIR) $(FASTJET)
$(RUNNERS_GENERIC): $(BINDIR)/Run%: $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/%.o $(SRCDIR)/Run%.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(MDSNANO_RUNNER): $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(MDSNANO_UTILSOBJ) $(ANADIR)/llp_MuonSystem_CA_mdsnano.o $(SRCDIR)/Runllp_MuonSystem_CA_mdsnano.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
