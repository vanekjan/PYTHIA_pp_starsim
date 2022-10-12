#define genevents_new_cxx
#include "genevents_new.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>

void genevents_new::Loop() 
{
	//   In a ROOT session, you can do:
	//      Root > .L genevents.C
	//      Root > genevents t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	TString oFile(outputFileName);
	TFile *fout = new TFile(oFile,"recreate");

  TH1D* D0_pt;
	TH1D* D0_mass;

  const int D0PDGid = 421;

  TH1D* Dpm_pt;
	TH1D* Dpm_mass;

  const int DpmPDGid = 411;


  D0_pt = new TH1D("D0_pt","D0_pt",100,0,10);
  D0_mass = new TH1D("D0_mass","D0_mass",100,1.8,2.0);


  Dpm_pt = new TH1D("Dpm_pt","Dpm_pt",100,0,10);
  Dpm_mass = new TH1D("Dpm_pt_mass","Dpm_mass",100,1.8,2.0);


	Long64_t nentries = fChain->GetEntriesFast();

  float lastEventNo = -1;
  Long64_t nEvents = 0;

	Long64_t nbytes = 0, nb = 0;

  cout<<"Nuber of entries: "<<nentries<<endl;

	for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

    if(jentry % 10000 == 0)
    {
      cout<<"Working on event: "<<jentry<<endl;
    }

    nEvents++;

		for (int i = 0; i < mParticles_; i++)  
    {
			if (abs(mParticles_mId[i]) == D0PDGid && mParticles_mStatus[i] == 2) 
      {
				D0_mass->Fill(mParticles_mMass[i]);
				double pt = TMath::Sqrt(mParticles_mPx[i]*mParticles_mPx[i] + mParticles_mPy[i]*mParticles_mPy[i]);
				D0_pt->Fill(pt);
			}

      if (abs(mParticles_mId[i]) == DpmPDGid && mParticles_mStatus[i] == 2) 
      {
				Dpm_mass->Fill(mParticles_mMass[i]);
				double pt = TMath::Sqrt(mParticles_mPx[i]*mParticles_mPx[i] + mParticles_mPy[i]*mParticles_mPy[i]);
				Dpm_pt->Fill(pt);
			}

		}

	}

  D0_pt->Sumw2();
  D0_pt->Scale(1./nEvents); //scale pT spectrum by number of events

  Dpm_pt->Sumw2();
  Dpm_pt->Scale(1./nEvents); //scale pT spectrum by number of events



  cout<<"Nuber of events: "<<nEvents<<endl;

  fout->cd();
  fout->Write();
  fout->Close();

  cout<<"end loop"<<endl;

}
