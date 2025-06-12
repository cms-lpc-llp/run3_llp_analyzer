# run3_llp_analyzer



Class for analyzing the 2015 razor ntuples

Setup (for el9 machines)
-------------

    cmsrel CMSSW_14_1_0_pre4
    cd CMSSW_14_1_0_pre4/src
    cmsenv
    source /cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.30.07-024df6516c17fd2edef848a927a788f1/bin/thisroot.sh
    git clone git@github.com:cms-lpc-llp/run3_llp_analyzer.git run3_llp_analyzer
    cd run3_llp_analyzer
    source set_SSL_vars.sh
    make
  
Defining a new analysis
-------------
1) Copy analyzers/DummyAnalyzer.cc and replace each instance of "DummyAnalyzer" with the name of your desired analyzer.
   Modify the body of the Analyze function to define your analyzer's behavior.
   DO NOT need to write a header file for the analyzer class; the Makefile will generate one for you automatically.  

2) Do `make`.  This will create an executable `bin/Run<name of your analyzer>`. You can execute your analysis using this program directly or by calling it via the `RazorRun` script. 

Running
------------
After compiling, 

    ./RazorRun <list of input files> <name of your analyzer> <options>
  

The "options" are the following:
    
    -d   --isData
    -f=  --outputFile=<output filename> (optional)
    -n=  --optionNumber=<option number> (optional)
    -l=  --optionLabel=<option Label> (optional)
    -h   --help


## Run the llp_analyzer
    ./RazorRun <list of input files> llp_MuonSystem -d=${isData} -n=${option} -f=${outputfile} -l=${tag}
* ```isData``` is ```yes``` or ```no```
* ```option``` is currently not used
* ```tag``` is currently not used
* list of input files are stored in ```lists```


### Submit condor jobs on tier2
Before submitting jobs, make sure proxy and CMSSW environment is setup.

* run the ```llp_MuonSystem``` analyzer for data:
	* ```python scripts_condor/submit_condor_caltech.py```
	* increment `analyzer_version` everytime the analyzer is being rerun

Normalizing the processed ntuples
------------
The NormalizeNtuple macro opens a specified set of files and adds a 'weight' branch to each TTree in each file.  The value of 'weight' is the same for all events in a tree and is equal to lumi * CrossSection/NEvents, where NEvents is the total number of events processed for the given dataset, and lumi is the luminosity normalized to.  The cross sections can be found in the file ```data/xSections.dat```.  To run NormalizeNtuple:

    ./NormalizeNtuple <input file list> [lumi]

* Make sure the dataset being processed have xSections in ```data/xSections.dat```
  
* Normalize the ntuples with condor job
  * ```python submit_normalize_caltech.py```
  * the script hadd the condor jobs and normalize the ntuple to cross section for signal, only hadd is done for data
  


### Filter good lumi events for data

https://github.com/RazorCMS/RazorCommon
