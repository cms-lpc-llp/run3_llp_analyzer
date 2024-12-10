import ROOT as rt
import correctionlib
files = ['Summer22', 'Summer22EE', 'Summer23','Summer23BPix']
for name in files:
    ceval = correctionlib.CorrectionSet.from_file(name+".json")
    outFile = rt.TFile("PileupReweight_" + name + '.root', 'RECREATE')


    for corr in ceval.values():
        key = corr.name
        print(f"Correction {corr.name} has {len(corr.inputs)} inputs")
        for ix in corr.inputs:
            print(f"   Input {ix.name} ({ix.type}): {ix.description}")
            print(len(ix.description))
    print(key)
    for weight in ['nominal','up','down']:
        hist =  rt.TH1F("npu_{}".format(weight),"npu_{}".format(weight),100,-0.5,99.5)
        for i in range(hist.GetNbinsX()):hist.SetBinContent(i+1, ceval[key].evaluate(hist.GetBinCenter(i+1), weight))
    
        outFile.WriteTObject(hist, "npu_{}".format(weight), "WriteDelete");
    
    outFile.Close();
