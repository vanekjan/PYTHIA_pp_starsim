#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include<fstream>
#include<sstream>
#include "TBenchmark.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TString.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TList.h"
#include "TFileInfo.h"

//#include "genevents_new.h"

//class genevents_new;

void run_new()   //0 = D0, 1 = D+-
{
	//TFile *_file0 = TFile::Open(in);


  TString out;

  ifstream fileList;
   
  out = "./output/output_Lambda_pp_510_MB_1M_events_daughter_open_cuts.root";
  fileList.open("/star/u/vanekjan/pwg/vanekjan/myPYTHIA_8_pp/production/2022-10-14_13-34__pp_510_MB_Lambda_filt_open_cuts/fileList.list");
  //fileList.open("/star/u/vanekjan/pwg/vanekjan/myPYTHIA_8_pp/production/2020-08-12_06-39/simdata/pythia8/sngdmeson/outFileList.list");


  TChain *myChain = new TChain("genevents");

  string fileFromList;

  while(getline(fileList, fileFromList))
  {
    myChain->Add(fileFromList.c_str());
  }

	//TTree *m = (TTree*)_file0->Get("genevents");
	gROOT->LoadMacro("genevents_new.C++");

  gSystem->Load("genevents_new_C.so");
	//genevents t(m);
  genevents_new t(myChain);

  //genevents *t = new genevents(myChain);

	t.setoutputFileName(out);
	t.Loop();

  fileList.close();

  cout<<"end"<<endl;  

  return;

}
