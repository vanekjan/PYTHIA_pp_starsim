#define genevents_new_cxx
#include "genevents_new.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>


void genevents_new::Loop(int daughterCuts, float energy) //0 - no additional cuts on dacay daughters, 1 - apply additional cuts
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

  
  const int nPtBins = 8;
  float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

  const int nEtaBins = 3;
  float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

  const int L0PDGid = 3122;
  const int L0barPDGid = -3122;


  //variables for event and particle loop
  TLorentzVector L_fourmom;
  TLorentzVector L_fourmom_reverse; //for proton boost

  TLorentzVector p_fourmom;
  TLorentzVector pi_fourmom;


  //histograms
  TH1D* L0_pt_hist = new TH1D("L0_pt_hist","L0_pt_hist",100,0,10);
	TH1D* L0_mass_hist = new TH1D("L0_mass_hist","L0_mass_hist",100,1.,1.2);
	
	TH2D* L0_pT_vs_L0_y = new TH2D("L0_pT_vs_L0_y", "L0_pT_vs_L0_y", 100, 0, 10, 100, -1, 1);
	
  TH1D *L0_thetaProdPlane[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_y_vs_p_eta[nPtBins+1];
  TH2D *L0_y_vs_pi_eta[nPtBins+1];
  
  TH1D *L0_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0_xF[nPtBins+1][nEtaBins+1];
  
  TH2D *L0_pT_vs_L0_pz[nEtaBins+1];
  
  
  TH2D *L0_p_eta_vs_pi_eta[nPtBins+1];


  TH1D* L0bar_pt_hist = new TH1D("L0bar_pt_hist","L0bar_pt_hist",100,0,10);
	TH1D* L0bar_mass_hist = new TH1D("L0bar_mass_hist","L0bar_mass_hist",100,1.,1.2);
	
	TH2D* L0bar_pT_vs_L0bar_y = new TH2D("L0bar_pT_vs_L0bar_y", "L0bar_pT_vs_L0bar_y", 100, 0, 10, 100, -1, 1);

  TH1D *L0bar_thetaProdPlane[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_y_vs_p_eta[nPtBins+1];  
  TH2D *L0bar_y_vs_pi_eta[nPtBins+1];
  
  TH1D *L0bar_pz[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_xF[nPtBins+1][nEtaBins+1];
  
  TH2D *L0bar_pT_vs_L0bar_pz[nEtaBins+1];
  
  TH2D *L0bar_p_eta_vs_pi_eta[nPtBins+1];
  
  

  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
    L0_y_vs_p_eta[pTbin] = new TH2D(Form("L0_y_vs_p_eta_pT_%i", pTbin), Form("L0_y_vs_p_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0_y_vs_pi_eta[pTbin] = new TH2D(Form("L0_y_vs_pi_eta_pT_%i", pTbin), Form("L0_y_vs_pi_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    L0_p_eta_vs_pi_eta[pTbin] = new TH2D(Form("L0_p_eta_vs_pi_eta_L0_pT_%i", pTbin), Form("L0_p_eta_vs_pi_eta_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
    
    L0bar_y_vs_p_eta[pTbin] = new TH2D(Form("L0bar_y_vs_p_eta_pT_%i", pTbin), Form("L0bar_y_vs_p_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    L0bar_y_vs_pi_eta[pTbin] = new TH2D(Form("L0bar_y_vs_pi_eta_pT_%i", pTbin), Form("L0bar_y_vs_pi_eta_pT_%i", pTbin), 100, -1, 1, 500, -5, 5);
    
    L0bar_p_eta_vs_pi_eta[pTbin] = new TH2D(Form("L0bar_p_eta_vs_pi_eta_L0_pT_%i", pTbin), Form("L0bar_p_eta_vs_pi_eta_L0_pT_%i", pTbin), 500, -5, 5, 500, -5, 5);
    
    
  
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {
      if(pTbin == 0)
      {
        L0_pT_vs_L0_pz[etaBin] = new TH2D(Form("L0_pT_vs_L0_pz_eta_%i", etaBin), Form("L0_pT_vs_L0_pz_eta_%i", etaBin), 100, 0, 10, 200, -10, 10);
        
        L0bar_pT_vs_L0bar_pz[etaBin] = new TH2D(Form("L0bar_pT_vs_L0bar_pz_eta_%i", etaBin), Form("L0bar_pT_vs_L0bar_pz_eta_%i", etaBin), 100, 0, 10, 200, -10, 10);    
      }
    

      L0_thetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0_cosThetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0_pz[pTbin][etaBin] = new TH1D(Form("L0_pz_pT_%i_eta__%i", pTbin, etaBin), Form("L0_pz_pT_%i_eta_%i", pTbin, etaBin), 200,-10, 10);
      L0_xF[pTbin][etaBin] = new TH1D(Form("L0_xF_pT_%i_eta__%i", pTbin, etaBin), Form("L0_xF_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);

      L0bar_thetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0bar_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
      L0bar_cosThetaProdPlane[pTbin][etaBin] = new TH1D(Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
      
      L0bar_pz[pTbin][etaBin] = new TH1D(Form("L0bar_pz_pT_%i_eta__%i", pTbin, etaBin), Form("L0bar_pz_pT_%i_eta_%i", pTbin, etaBin), 200, -10, 10);
      L0bar_xF[pTbin][etaBin] = new TH1D(Form("L0bar_xF_pT_%i_eta__%i", pTbin, etaBin), Form("L0bar_xF_pT_%i_eta_%i", pTbin, etaBin), 100, 0, 0.01);
       
    }

  }



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

    
    //loop over particles in event
		for (int i = 0; i < mParticles_; i++)  
    {     

      //select only Lambda and Lambda-bar
      //first daughter is proton, second is pion
			if( fabs(mParticles_mId[i]) == L0PDGid && mParticles_mStatus[i] == 2 )
      {        

        //get Ids of daughters in mParticles array
        int daughter1_Id = mParticles_mDaughter[i][0];
        int daughter2_Id = mParticles_mDaughter[i][1];	

        //Lambda fourmomentum
        L_fourmom.SetPxPyPzE(mParticles_mPx[i], mParticles_mPy[i], mParticles_mPz[i], mParticles_mEnergy[i]);
        
        double L_xF = fabs(L_fourmom.Pz())/energy/2.; //proton in beam has momentum mCmsEnergy/2.
                
        //if( fabs(L_fourmom.Eta()) >= 0.5 ) continue;
        
        L_fourmom_reverse.SetPxPyPzE(-mParticles_mPx[i], -mParticles_mPy[i], -mParticles_mPz[i], mParticles_mEnergy[i]);
                

        p_fourmom.SetPxPyPzE(mParticles_mPx[daughter1_Id], mParticles_mPy[daughter1_Id], mParticles_mPz[daughter1_Id], mParticles_mEnergy[daughter1_Id]);
        pi_fourmom.SetPxPyPzE(mParticles_mPx[daughter2_Id], mParticles_mPy[daughter2_Id], mParticles_mPz[daughter2_Id], mParticles_mEnergy[daughter2_Id]); 
        
        
        
        if( daughterCuts == 1 )
        {
          if( fabs(p_fourmom.Eta()) >= 1. || fabs(pi_fourmom.Eta()) >= 1. ) continue;
          
          if( p_fourmom.Pt() < 0.15 || p_fourmom.Pt() > 20. ) continue;
          if( pi_fourmom.Pt() < 0.15 || pi_fourmom.Pt() > 20. ) continue;
        }
        
        

        TLorentzVector p_fourmom_star = p_fourmom;
        p_fourmom_star.Boost(L_fourmom_reverse.BoostVector());  

        TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
        TVector3 mProdPlane = beamVector.Cross(L_fourmom.Vect());
        mProdPlane = ( mProdPlane )*(1./mProdPlane.Mag() );

        float mThetaProdPlane = mProdPlane.Angle(p_fourmom_star.Vect());


        
        //fill all histograms for all pT and centrality bins
        int pT_bin = -1;

        //find pT bin of Lambda
        for(int j = 0; j < nPtBins; j++) //loop over pT bins
        {
          if(L_fourmom.Pt() > pT_bins[j] && L_fourmom.Pt() <= pT_bins[j+1])
          {
            pT_bin = j;
            break; //stop after pT bin is found
          }
        }

        if( pT_bin == -1 ) continue;


        //fill all histograms for all eta and centrality bins
        int eta_bin = -1;

        //find eta bin of Lambda
        for(int j = 0; j < nEtaBins; j++) //loop over eta bins
        {
          if(L_fourmom.Eta() > eta_bins[j] && L_fourmom.Eta() <= eta_bins[j+1])
          {
            eta_bin = j;
            break; //stop after eta bin is found
          }
        }

        if( eta_bin == -1 ) continue;


        //Lambda
        if( mParticles_mId[i] == L0PDGid )
        {
          L0_mass_hist->Fill(mParticles_mMass[i]);
          L0_pt_hist->Fill(L_fourmom.Pt());
          
          L0_pT_vs_L0_y->Fill(L_fourmom.Pt(), L_fourmom.Rapidity());

          L0_thetaProdPlane[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0_thetaProdPlane[nPtBins][eta_bin]->Fill(mThetaProdPlane);  //pT integrated, eta bins
          L0_thetaProdPlane[pT_bin][nEtaBins]->Fill(mThetaProdPlane);  //pT bins, -1 < eta < 1
          L0_thetaProdPlane[nPtBins][nEtaBins]->Fill(mThetaProdPlane); //pT integrated and -1 < eta < 1


          L0_cosThetaProdPlane[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0_cosThetaProdPlane[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          L0_y_vs_p_eta[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_y_vs_p_eta[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0_y_vs_pi_eta[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0_p_eta_vs_pi_eta[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0_p_eta_vs_pi_eta[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          
          L0_pT_vs_L0_pz[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0_pT_vs_L0_pz[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());          
          
          L0_pz[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0_pz[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0_pz[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());          
          
          L0_xF[pT_bin][eta_bin]->Fill(L_xF);
          L0_xF[nPtBins][eta_bin]->Fill(L_xF);
          L0_xF[pT_bin][nEtaBins]->Fill(L_xF);
          L0_xF[nPtBins][nEtaBins]->Fill(L_xF);

        }

        //Lambda-bar
        if( mParticles_mId[i] == L0barPDGid )
        {
          L0bar_mass_hist->Fill(mParticles_mMass[i]);
          L0bar_pt_hist->Fill(L_fourmom.Pt());
          
          L0bar_pT_vs_L0bar_y->Fill(L_fourmom.Pt(), L_fourmom.Rapidity());

          L0bar_thetaProdPlane[pT_bin][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[nPtBins][eta_bin]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[pT_bin][nEtaBins]->Fill(mThetaProdPlane);
          L0bar_thetaProdPlane[nPtBins][nEtaBins]->Fill(mThetaProdPlane);

          L0bar_cosThetaProdPlane[pT_bin][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[nPtBins][eta_bin]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[pT_bin][nEtaBins]->Fill(cos(mThetaProdPlane));
          L0bar_cosThetaProdPlane[nPtBins][nEtaBins]->Fill(cos(mThetaProdPlane));
          
          L0bar_y_vs_p_eta[pT_bin]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta[pT_bin]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          L0bar_y_vs_p_eta[nPtBins]->Fill(L_fourmom.Rapidity(), p_fourmom.Eta());
          L0bar_y_vs_pi_eta[nPtBins]->Fill(L_fourmom.Rapidity(), pi_fourmom.Eta());
          
          
          L0bar_p_eta_vs_pi_eta[pT_bin]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          L0bar_p_eta_vs_pi_eta[nPtBins]->Fill(p_fourmom.Eta(), pi_fourmom.Eta());
          
          L0bar_pT_vs_L0bar_pz[eta_bin]->Fill(L_fourmom.Pt(), L_fourmom.Pz());
          L0bar_pT_vs_L0bar_pz[nEtaBins]->Fill(L_fourmom.Pt(), L_fourmom.Pz());          
          
          L0bar_pz[pT_bin][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz[nPtBins][eta_bin]->Fill(L_fourmom.Pz());
          L0bar_pz[pT_bin][nEtaBins]->Fill(L_fourmom.Pz());
          L0bar_pz[nPtBins][nEtaBins]->Fill(L_fourmom.Pz());          
          
          L0bar_xF[pT_bin][eta_bin]->Fill(L_xF);
          L0bar_xF[nPtBins][eta_bin]->Fill(L_xF);
          L0bar_xF[pT_bin][nEtaBins]->Fill(L_xF);
          L0bar_xF[nPtBins][nEtaBins]->Fill(L_xF);
          
        }
        
			}      

		}

	}

  L0_pt_hist->Sumw2();
  L0_pt_hist->Scale(1./nEvents); //scale pT spectrum by number of events

  L0bar_pt_hist->Sumw2();
  L0bar_pt_hist->Scale(1./nEvents); //scale pT spectrum by number of events



  cout<<"Nuber of events: "<<nEvents<<endl;

  fout->cd();
  fout->Write();
  fout->Close();

  cout<<"end loop"<<endl;

}
