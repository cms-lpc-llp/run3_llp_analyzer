include Makefile.inc

DIRS = python
SRCDIR = src
INCLUDEDIR = include
ANADIR = analyzers
BINDIR = bin
FASTJET = fastjet-install/bin/fastjet-config

# Add CMSSW package roots for includes like DataFormats/... when building
# this external analyzer with plain g++ (outside scram's build system).
ifneq ($(strip $(CMSSW_BASE)),)
	CXXFLAGS += -I$(CMSSW_BASE)/src
endif
ifneq ($(strip $(CMSSW_RELEASE_BASE)),)
	CXXFLAGS += -I$(CMSSW_RELEASE_BASE)/src
endif

ANALYZERS = $(wildcard $(ANADIR)/*.cc)
ANALYZERSH = $(ANALYZERS:cc=h)
ANALYZERSOBJ = $(ANALYZERS:cc=o)
ANALYZER_NAMES = $(notdir $(basename $(ANALYZERS)))
RUNNERS = $(addprefix $(BINDIR)/Run,$(notdir $(basename $(ANALYZERS))))
RUNNER_SRCS = $(addprefix $(SRCDIR)/Run,$(addsuffix .cc,$(ANALYZER_NAMES)))

MDSNANO_ANALYZER_NAME = llp_MuonSystem_CA_mdsnano
MDSNANO_ANALYZER_DIR = $(ANADIR)/$(MDSNANO_ANALYZER_NAME)
MDSNANO_ANALYZER_INCLUDEDIR = $(MDSNANO_ANALYZER_DIR)/include
MDSNANO_ANALYZER_SRCDIR = $(MDSNANO_ANALYZER_DIR)/src
MDSNANO_ANALYZER_SRC = $(MDSNANO_ANALYZER_DIR)/main.cc
MDSNANO_ANALYZER_HDR = $(MDSNANO_ANALYZER_DIR)/main.h
MDSNANO_ANALYZER_OBJ = $(MDSNANO_ANALYZER_DIR)/main.o
MDSNANO_PHASE_SRCS = $(MDSNANO_ANALYZER_SRCDIR)/MuonSystemSynthesisPhase.cc \
					 $(MDSNANO_ANALYZER_SRCDIR)/MuonSystemCuttingPhase.cc \
					 $(MDSNANO_ANALYZER_SRCDIR)/MuonSystemFillingPhase.cc \
					 $(MDSNANO_ANALYZER_SRCDIR)/MuonSystemSignalScanManager.cc
MDSNANO_PHASE_OBJS = $(MDSNANO_PHASE_SRCS:.cc=.o)
MDSNANO_RUNNER = $(BINDIR)/Run$(MDSNANO_ANALYZER_NAME)
MDSNANO_RUNNER_SRC = $(SRCDIR)/Run$(MDSNANO_ANALYZER_NAME).cc

RUNNERS += $(MDSNANO_RUNNER)
RUNNERS_GENERIC = $(filter-out $(MDSNANO_RUNNER),$(RUNNERS))
UTILS =$(SRCDIR)/RazorHelper.cc $(SRCDIR)/JetCorrectorParameters.cc \
        $(SRCDIR)/SimpleJetCorrectionUncertainty.cc \
					$(SRCDIR)/JetCorrectionUncertainty.cc \
		       	$(SRCDIR)/CACluster.cc ${SRCDIR}/TreeMuonSystemCombination.cc
UTILSOBJ = $(UTILS:cc=o)
MDSNANO_UTILSOBJ = $(SRCDIR)/RazorHelper.o $(SRCDIR)/JetCorrectorParameters.o \
					$(SRCDIR)/SimpleJetCorrectionUncertainty.o $(SRCDIR)/JetCorrectionUncertainty.o \
					$(SRCDIR)/CACluster.o $(SRCDIR)/TreeMuonSystemCombination.o $(MDSNANO_PHASE_OBJS)
EXECUTABLES = $(RUNNERS)
HELPERSCRIPT = python/MakeAnalyzerCode.py
RUNNER_TEMPLATES = $(INCLUDEDIR)/AnalyzerTemplate.txt \
		   $(INCLUDEDIR)/AnalyzerTemplate_mdsnano.txt \
		   $(SRCDIR)/RunAnalyzerTemplate.txt


.PHONY: clean all lxplus copy_runners

all: $(FASTJET) $(RUNNER_SRCS) $(MDSNANO_RUNNER_SRC) $(EXECUTABLES)

lxplus: all

clean:
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o $(ANADIR)/*/*.o $(ANADIR)/*/*/*.o $(RUNNERS) $(filter-out $(MDSNANO_RUNNER_SRC),$(wildcard $(SRCDIR)/Run*.cc))
copy_runners: $(RUNNER_SRCS)

$(RUNNER_SRCS): $(SRCDIR)/Run%.cc: $(ANADIR)/%.cc $(HELPERSCRIPT) $(RUNNER_TEMPLATES)
	@if [ ! -f "$@" ]; then \
		echo "$* file does not exists, copying"; \
		$(HELPERSCRIPT) $*; \
	fi

$(MDSNANO_RUNNER_SRC):
	@if [ ! -f "$@" ]; then \
		echo "[ERROR]: missing $@; restore it from git" ; \
		exit 1; \
	fi

$(BINDIR):
	@mkdir -p $(BINDIR)

$(FASTJET):
	$(ANADIR)/fastjet_install.sh

# Ensure FastJet is fully installed before any compilation/linking that
# expands fastjet-config flags.
$(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANALYZERSOBJ) $(MDSNANO_ANALYZER_OBJ) $(MDSNANO_PHASE_OBJS): | $(FASTJET)

$(SRCDIR)/mdsnano_event.o: $(SRCDIR)/mdsnano_event.cc $(INCLUDEDIR)/mdsnano_event.h
	$(CXX) $(SRCDIR)/mdsnano_event.cc $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
		
$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(MDSNANO_ANALYZER_OBJ): $(MDSNANO_ANALYZER_SRC) $(MDSNANO_ANALYZER_HDR)
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) -I$(MDSNANO_ANALYZER_INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(MDSNANO_PHASE_OBJS): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) -I$(MDSNANO_ANALYZER_INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(RUNNERS): | $(BINDIR) $(FASTJET)
$(RUNNERS_GENERIC): $(BINDIR)/Run%: $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/%.o $(SRCDIR)/Run%.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(MDSNANO_RUNNER): $(SRCDIR)/mdsnano_event.o $(SRCDIR)/RazorAnalyzer.o $(MDSNANO_UTILSOBJ) $(MDSNANO_ANALYZER_OBJ) $(MDSNANO_RUNNER_SRC)
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
