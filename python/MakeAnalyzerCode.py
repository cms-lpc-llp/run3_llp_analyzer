#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    print("Please specify an analyzer name!")
    sys.exit()

analyzer = sys.argv[1]

if "mdsnano" in analyzer or "merge" in analyzer: #added so that analyzer header files calls correct RazorAnalyzer
    inNames = ['include/AnalyzerTemplate_mdsnano.txt','src/RunAnalyzerTemplate.txt']
else:
     inNames = ['include/AnalyzerTemplate.txt','src/RunAnalyzerTemplate.txt']
outNames = ['analyzers/'+analyzer+'.h','src/Run'+analyzer+'.cc']

print(inNames)
print(outNames)

if analyzer.find("Run1") > 0:
	inNames = ['include/AnalyzerTemplateRun1.txt','src/RunAnalyzerTemplateRun1.txt']

if analyzer.find("UpgradeTiming") > 0:
	inNames = ['include/AnalyzerTemplateUpgradeTiming.txt','src/RunAnalyzerTemplateUpgradeTiming.txt']

for i in range(len(inNames)):
    with open(inNames[i]) as inF:
        with open(outNames[i],'w') as outF:
            for line in inF:
                outF.write( line.replace('%ANALYZER%',analyzer) )       
