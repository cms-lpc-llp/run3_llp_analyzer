# Post Analyzer Selection + Plotting for Tau HNL Search (AlexAlbert41)



Notebooks for applying cuts to the output of HNL Analyzers

Jupyter Notebook Kernel Setup (for cmslpc el9 on VSCode)
--------------
I tried to get the Jupyter Notebook to run with ROOT on VSCode. I found that after adding ROOT to the PYTHONPATH in the kernel env, it seem to lose access to the local site packages (numpy, pandas, etc), so they have to be re-added to the PYTHONPATH specified in the kernel config file (*kernel.json*) to have access to these packages. I recommend using the supplied kernel.json file here if initializing your own kernel (please change the few local paths with */uscms/amalbert/*), though the script find_pkg_path.py can be used to append package paths to the PYTHONPATH of an existing kernel config file. This kernel config will get you almost all of the required packages, but pip install may be occasionally required within the notebook.

Generating a Kernel (after cmsenv)
-------------

    /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_14_1_0_pre4/external/el8_amd64_gcc12/bin/python3 -m ipykernel install --user --name Pyroot_Local --display-name "PyROOT + Local Site Packages"
    #You should see *Installed kernelspec Pyroot_Local in <directory>*
    cp kernel.json <directory>/kernel.json #remember to modify local paths in PYTHONPATH in kernel.json

Then, restart VSCode, and select *PyROOT + Local Site Packages* from the kernel manager when opening a notebook.
