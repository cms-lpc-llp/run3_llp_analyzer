//C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>

#include <llp_MuonSystem_SUEP.h>

using namespace std;

std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
        {
          if ( tmp.find( "=" )  != std::string::npos ) return tmp.substr( tmp.find_last_of("=") + 1 );
	  if ( tmp.find( "--" ) != std::string::npos ) return "yes";
	}
    }

  return "";
};

void usage()
{
  std::cerr << "Usage: Runllp_MuonSystem_SUEP  <input list>  [options]\n[options]:\n"
	    << "-d  --isData\n"
	    << "-w  --isWeighted\n"
	    << "-f=  --outputFile=<output filename> (optional)\n"
	    << "-n=  --optionNumber=<option number> (optional)\n"
	    << "-l=  --optionLabel=<option Label> (optional)\n"
	    << "-h  --help"
	    << std::endl;
};

int main(int argc, char* argv[]){

  //get input files and analysis type from command line
  if ( ParseCommandLine( argc, argv, "--help" ) != ""  || ParseCommandLine( argc, argv, "-h" ) != ""  || argc < 2 )
    {
      usage();
      return -1;
    }

  //----------------------------------------
  //Getting <input list> and <analysis type>
  //----------------------------------------
  string inputFileName(argv[1]);

  //--------------------------------
  //G e t t i n g   d a t a  f l a g
  //--------------------------------
  std::string _isData = ParseCommandLine( argc, argv, "--isData=" );
  std::string _isWeighted = ParseCommandLine( argc, argv, "--isWeighted=" );
  std::string _d = ParseCommandLine( argc, argv, "-d=" );
  std::string _w = ParseCommandLine( argc, argv, "-w=" );
  bool isData = false;
  if ( _isData == "yes" || _d == "yes" ) isData = true;
  bool isWeighted = false;
  cout << _isWeighted << endl;
  if ( _isWeighted == "yes" || _w == "yes" ) isWeighted = true;

  //---------------------------------------------
  //G e t t i n g   o u t p u t F i l e   N a m e
  //---------------------------------------------
  std::string _outFile = ParseCommandLine( argc, argv, "--outputFile=" );
  std::string _f = ParseCommandLine( argc, argv, "-f=" );
  string outputFileName = "";
  if ( _outFile != "" )
    {
      outputFileName = _outFile;
    }
  else if ( _f != "" )
    {
      outputFileName = _f;
    }
  else
    {
      std::cerr << "[WARNING]: output ROOT file not provided, using default output" << std::endl;
    }

  //-----------------------------------------
  //G e t t i n g   o p t i o n   n u m b e r
  //-----------------------------------------
  int option = -1;
  std::string _optionNumber = ParseCommandLine( argc, argv, "--optionNumber=" );
  std::string _n = ParseCommandLine( argc, argv, "-n=" );
  if ( _optionNumber != "" )
    {
      option = atoi( _optionNumber.c_str() );
    }
  else if ( _n != "" )
    {
      option = atoi( _n.c_str() );
    }
  else
    {
      std::cerr << "[WARNING]: option number not provided, using default option number" << std::endl;
    }

  string label = "";
  std::string _optionLabel = ParseCommandLine( argc, argv, "--optionLabel=" );
  std::string _l = ParseCommandLine( argc, argv, "-l=" );
  if ( _optionLabel != "" )
    {
      label = _optionLabel;
    }
  else if ( _l != "" )
    {
      label = _l;
    }
  else
    {
      std::cerr << "[WARNING]: optional label not provided, using default optional label" << std::endl;
    }

  std::cout << "[INFO]: <input list> --> " << inputFileName << std::endl;
  std::cout << "[INFO]: isData --> " << isData << std::endl;
  std::cout << "[INFO]: isWeighted --> " << isWeighted << std::endl;
  std::cout << "[INFO]: outputFileName --> " << outputFileName << std::endl;
  std::cout << "[INFO]: option --> " << option << std::endl;
  std::cout << "[INFO]: optionalLabel --> " << label << std::endl;

  int mh = 0;
  int mx = 0;
  int ctau = 0;
  int T = 0;

  //build the TChain
  //tree name is set give the structure in the first root file, see while loop below
  TChain* theChain = new TChain();
  string curFileName;
  ifstream inputFile(inputFileName.c_str());
  int NFilesLoaded = 0;
  if ( !inputFile ){
    cerr << "Error: input file not found!" << endl;
    return -1;
  }


  while ( getline(inputFile, curFileName) )
    {
      cout << "Here check if weighted or not! " << endl;
      string outFileString;
      if(isWeighted)
	{
	  string outDir = (curFileName.substr(0, curFileName.find_last_of("/"))) + "/weighted/";
	  outFileString = outDir+(curFileName.substr(curFileName.find_last_of("/")+1,-1));
	}
      else
	{
	  outFileString = curFileName;
	}
      std::cout << "Opening " << outFileString << endl;

      if ( NFilesLoaded == 0 )
	{
	  /*
	    checks root file structure and add first file
	  */
	  string mh_substring = outFileString.substr(outFileString.find("mMed-")+5);
	  mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
	  string mx_substring = outFileString.substr(outFileString.find("mDark-")+6);
	  mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
	  string T_substring = outFileString.substr(outFileString.find("temp-")+5);
	  T = stoi(T_substring.substr(0,T_substring.find('_')));
	  string ctau_substring = outFileString.substr(outFileString.find("ctau-")+5);
	  ctau = stoi(ctau_substring.substr(0,ctau_substring.find('p')));

	  TFile* f_0 = TFile::Open( outFileString.c_str() );
	  //TH1F *h_NEventsGen = (TH1F*) f_0->Get("ntuples/NEventsGen");
	  //h_NEventsGen->SetDirectory(0);
	  //double NEventsGen = h_NEventsGen->GetBinContent(1);
	  //cout << "NEventsGen: " << NEventsGen << endl;//LB
	  //TH1F *h_NEventsGenPass = (TH1F*) f_0->Get("ntuples/NEventsGenPass");
	  //h_NEventsGenPass->SetDirectory(0);
	  //double NEventsGenPass = h_NEventsGenPass->GetBinContent(1);
	  //cout << "NEventsGenPass: " << NEventsGenPass << endl;//LB

	  if(isWeighted)
	    {
	      theChain->SetName("llp");
	      std::cout << "[INFO]: weighted tree configuration for tchain" << std::endl;
	    }
	  else
	    {
	      if( f_0->GetDirectory("ntuples") )
		{
		  theChain->SetName("ntuples/llp");
		  std::cout << "[INFO]: default configuration for tchain" << std::endl;
		}
	      else
		{
		  theChain->SetName("Events");
		  std::cout << "[INFO]: alternative configuration for tchain" << std::endl;
		}
	    }
	  theChain->Add( outFileString.c_str() );
	  delete f_0;
	}
      else
	{
	  theChain->Add( outFileString.c_str() );
	}
      NFilesLoaded++;
    }
  std::cout << "Loaded Total of " << NFilesLoaded << " files\n";
  if ( theChain == NULL ) return -1;

  llp_MuonSystem_SUEP analyzer(theChain);

  //------ EXECUTE ------//
  cout << "Executing llp_MuonSystem_SUEP..." << endl;
  analyzer.EnableAll();
  analyzer.Analyze(isData, option, outputFileName, label, mh, mx, ctau, T);
  cout << "Process completed!" << endl;
  cerr << "------------------------------" << endl; //present so that an empty .err file corresponds to a failed job

  return 0;
}
