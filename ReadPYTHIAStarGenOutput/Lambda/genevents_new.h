//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  5 13:09:34 2020 by ROOT version 5.34/38
// from TTree genevents/TTree containing event generator information
// found on file: pythia8tune_sngd0_12099044_4580001.starsim.root
//////////////////////////////////////////////////////////

#ifndef genevents_new_h
#define genevents_new_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <vector>
#include "TBenchmark.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TMath.h"
#include <vector>

using namespace std;


// Header file for the classes stored in the TTree if any.
//#include "/home/vanekjan/root5/root/core/base/inc/TObject.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxmParticles = 1000;

class genevents_new {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Declaration of leaf types
		//StarGenEvent    *primaryEvent;
		UInt_t          fUniqueID;
		UInt_t          fBits;
		TString         mName;
		TString         mTitle;
		Int_t           mParticles_;
		UInt_t          mParticles_fUniqueID[kMaxmParticles];   //[mParticles_]
		UInt_t          mParticles_fBits[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mStatus[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mId[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mMother[kMaxmParticles][2];   //[mParticles_]
		Int_t           mParticles_mDaughter[kMaxmParticles][2];   //[mParticles_]
		Float_t         mParticles_mPx[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mPy[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mPz[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mEnergy[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mMass[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mVx[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mVy[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mVz[kMaxmParticles];   //[mParticles_]
		Float_t         mParticles_mTof[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mIndex[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mStack[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mPrimaryKey[kMaxmParticles];   //[mParticles_]
		Int_t           mParticles_mGeneratorId[kMaxmParticles];   //[mParticles_]
		Int_t           mGeneratorId;
		Int_t           mProcessId;
		Int_t           mOffset;
		Int_t           mEventNumber;
		Int_t           mRunNumber;
		Int_t           mDaqRunNumber;
		Int_t           mDaqFileNumber;
		Int_t           mBlueId;
		Int_t           mYellId;
		Double_t        mCmsEnergy;
		Int_t           mNumRejected[3];
		UInt_t          mFilterResult;
		vector<Double_t> mWeights;
		Int_t           mNumParticles;
		//StarGenPPEvent  *Pythia8;
		//StarGenPPEvent  *StarGenEvent;
		Int_t           idBlue;
		Int_t           idYell;
		Int_t           process;
		Int_t           subprocess;
		Int_t           idParton1;
		Int_t           idParton2;
		Double_t        xParton1;
		Double_t        xParton2;
		Double_t        xPdf1;
		Double_t        xPdf2;
		Double_t        Q2fac;
		Double_t        Q2ren;
		Bool_t          valence1;
		Bool_t          valence2;
		Double_t        sHat;
		Double_t        tHat;
		Double_t        uHat;
		Double_t        ptHat;
		Double_t        thetaHat;
		Double_t        phiHat;
		Double_t        weight;

		// List of branches
		TBranch        *b_primaryEvent_fUniqueID;   //!
		TBranch        *b_primaryEvent_fBits;   //!
		TBranch        *b_primaryEvent_mName;   //!
		TBranch        *b_primaryEvent_mTitle;   //!
		TBranch        *b_primaryEvent_mParticles_;   //!
		TBranch        *b_mParticles_fUniqueID;   //!
		TBranch        *b_mParticles_fBits;   //!
		TBranch        *b_mParticles_mStatus;   //!
		TBranch        *b_mParticles_mId;   //!
		TBranch        *b_mParticles_mMother;   //!
		TBranch        *b_mParticles_mDaughter;   //!
		TBranch        *b_mParticles_mPx;   //!
		TBranch        *b_mParticles_mPy;   //!
		TBranch        *b_mParticles_mPz;   //!
		TBranch        *b_mParticles_mEnergy;   //!
		TBranch        *b_mParticles_mMass;   //!
		TBranch        *b_mParticles_mVx;   //!
		TBranch        *b_mParticles_mVy;   //!
		TBranch        *b_mParticles_mVz;   //!
		TBranch        *b_mParticles_mTof;   //!
		TBranch        *b_mParticles_mIndex;   //!
		TBranch        *b_mParticles_mStack;   //!
		TBranch        *b_mParticles_mPrimaryKey;   //!
		TBranch        *b_mParticles_mGeneratorId;   //!
		TBranch        *b_primaryEvent_mGeneratorId;   //!
		TBranch        *b_primaryEvent_mProcessId;   //!
		TBranch        *b_primaryEvent_mOffset;   //!
		TBranch        *b_primaryEvent_mEventNumber;   //!
		TBranch        *b_primaryEvent_mRunNumber;   //!
		TBranch        *b_primaryEvent_mDaqRunNumber;   //!
		TBranch        *b_primaryEvent_mDaqFileNumber;   //!
		TBranch        *b_primaryEvent_mBlueId;   //!
		TBranch        *b_primaryEvent_mYellId;   //!
		TBranch        *b_primaryEvent_mCmsEnergy;   //!
		TBranch        *b_primaryEvent_mNumRejected;   //!
		TBranch        *b_primaryEvent_mFilterResult;   //!
		TBranch        *b_primaryEvent_mWeights;   //!
		TBranch        *b_primaryEvent_mNumParticles;   //!
		TBranch        *b_Pythia8_idBlue;   //!
		TBranch        *b_Pythia8_idYell;   //!
		TBranch        *b_Pythia8_process;   //!
		TBranch        *b_Pythia8_subprocess;   //!
		TBranch        *b_Pythia8_idParton1;   //!
		TBranch        *b_Pythia8_idParton2;   //!
		TBranch        *b_Pythia8_xParton1;   //!
		TBranch        *b_Pythia8_xParton2;   //!
		TBranch        *b_Pythia8_xPdf1;   //!
		TBranch        *b_Pythia8_xPdf2;   //!
		TBranch        *b_Pythia8_Q2fac;   //!
		TBranch        *b_Pythia8_Q2ren;   //!
		TBranch        *b_Pythia8_valence1;   //!
		TBranch        *b_Pythia8_valence2;   //!
		TBranch        *b_Pythia8_sHat;   //!
		TBranch        *b_Pythia8_tHat;   //!
		TBranch        *b_Pythia8_uHat;   //!
		TBranch        *b_Pythia8_ptHat;   //!
		TBranch        *b_Pythia8_thetaHat;   //!
		TBranch        *b_Pythia8_phiHat;   //!
		TBranch        *b_Pythia8_weight;   //!

		genevents_new(TTree *tree=0);
		virtual ~genevents_new();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop(int daughterCuts, float energy);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		TString outputFileName;
		void setoutputFileName(TString o);

};

void genevents_new::setoutputFileName(TString o) {
        outputFileName = o;
}


#endif

#ifdef genevents_new_cxx
genevents_new::genevents_new(TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) 
  {
		TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pythia8tune_sngd0_12099044_4580001.starsim.root");
		if (!f || !f->IsOpen())
    {
			f = new TFile("pythia8tune_sngd0_12099044_4580001.starsim.root");
		}
		f->GetObject("genevents",tree);
	}

  outputFileName = "out.root";
	Init(tree);
}

genevents_new::~genevents_new()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t genevents_new::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t genevents_new::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void genevents_new::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_primaryEvent_fUniqueID);
	fChain->SetBranchAddress("fBits", &fBits, &b_primaryEvent_fBits);
	fChain->SetBranchAddress("mName", &mName, &b_primaryEvent_mName);
	fChain->SetBranchAddress("mTitle", &mTitle, &b_primaryEvent_mTitle);
	fChain->SetBranchAddress("mParticles", &mParticles_, &b_primaryEvent_mParticles_);
	fChain->SetBranchAddress("mParticles.fUniqueID", mParticles_fUniqueID, &b_mParticles_fUniqueID);
	fChain->SetBranchAddress("mParticles.fBits", mParticles_fBits, &b_mParticles_fBits);
	fChain->SetBranchAddress("mParticles.mStatus", mParticles_mStatus, &b_mParticles_mStatus);
	fChain->SetBranchAddress("mParticles.mId", mParticles_mId, &b_mParticles_mId);
	fChain->SetBranchAddress("mParticles.mMother[2]", mParticles_mMother, &b_mParticles_mMother);
	fChain->SetBranchAddress("mParticles.mDaughter[2]", mParticles_mDaughter, &b_mParticles_mDaughter);
	fChain->SetBranchAddress("mParticles.mPx", mParticles_mPx, &b_mParticles_mPx);
	fChain->SetBranchAddress("mParticles.mPy", mParticles_mPy, &b_mParticles_mPy);
	fChain->SetBranchAddress("mParticles.mPz", mParticles_mPz, &b_mParticles_mPz);
	fChain->SetBranchAddress("mParticles.mEnergy", mParticles_mEnergy, &b_mParticles_mEnergy);
	fChain->SetBranchAddress("mParticles.mMass", mParticles_mMass, &b_mParticles_mMass);
	fChain->SetBranchAddress("mParticles.mVx", mParticles_mVx, &b_mParticles_mVx);
	fChain->SetBranchAddress("mParticles.mVy", mParticles_mVy, &b_mParticles_mVy);
	fChain->SetBranchAddress("mParticles.mVz", mParticles_mVz, &b_mParticles_mVz);
	fChain->SetBranchAddress("mParticles.mTof", mParticles_mTof, &b_mParticles_mTof);
	fChain->SetBranchAddress("mParticles.mIndex", mParticles_mIndex, &b_mParticles_mIndex);
	fChain->SetBranchAddress("mParticles.mStack", mParticles_mStack, &b_mParticles_mStack);
	fChain->SetBranchAddress("mParticles.mPrimaryKey", mParticles_mPrimaryKey, &b_mParticles_mPrimaryKey);
	fChain->SetBranchAddress("mParticles.mGeneratorId", mParticles_mGeneratorId, &b_mParticles_mGeneratorId);
	fChain->SetBranchAddress("mGeneratorId", &mGeneratorId, &b_primaryEvent_mGeneratorId);
	fChain->SetBranchAddress("mProcessId", &mProcessId, &b_primaryEvent_mProcessId);
	fChain->SetBranchAddress("mOffset", &mOffset, &b_primaryEvent_mOffset);
	fChain->SetBranchAddress("mEventNumber", &mEventNumber, &b_primaryEvent_mEventNumber);
	fChain->SetBranchAddress("mRunNumber", &mRunNumber, &b_primaryEvent_mRunNumber);
	fChain->SetBranchAddress("mDaqRunNumber", &mDaqRunNumber, &b_primaryEvent_mDaqRunNumber);
	fChain->SetBranchAddress("mDaqFileNumber", &mDaqFileNumber, &b_primaryEvent_mDaqFileNumber);
	fChain->SetBranchAddress("mBlueId", &mBlueId, &b_primaryEvent_mBlueId);
	fChain->SetBranchAddress("mYellId", &mYellId, &b_primaryEvent_mYellId);
	fChain->SetBranchAddress("mCmsEnergy", &mCmsEnergy, &b_primaryEvent_mCmsEnergy);
	fChain->SetBranchAddress("mNumRejected[3]", mNumRejected, &b_primaryEvent_mNumRejected);
	fChain->SetBranchAddress("mFilterResult", &mFilterResult, &b_primaryEvent_mFilterResult);
	fChain->SetBranchAddress("mWeights", &mWeights, &b_primaryEvent_mWeights);
	fChain->SetBranchAddress("mNumParticles", &mNumParticles, &b_primaryEvent_mNumParticles);
	fChain->SetBranchAddress("idBlue", &idBlue, &b_Pythia8_idBlue);
	fChain->SetBranchAddress("idYell", &idYell, &b_Pythia8_idYell);
	fChain->SetBranchAddress("process", &process, &b_Pythia8_process);
	fChain->SetBranchAddress("subprocess", &subprocess, &b_Pythia8_subprocess);
	fChain->SetBranchAddress("idParton1", &idParton1, &b_Pythia8_idParton1);
	fChain->SetBranchAddress("idParton2", &idParton2, &b_Pythia8_idParton2);
	fChain->SetBranchAddress("xParton1", &xParton1, &b_Pythia8_xParton1);
	fChain->SetBranchAddress("xParton2", &xParton2, &b_Pythia8_xParton2);
	fChain->SetBranchAddress("xPdf1", &xPdf1, &b_Pythia8_xPdf1);
	fChain->SetBranchAddress("xPdf2", &xPdf2, &b_Pythia8_xPdf2);
	fChain->SetBranchAddress("Q2fac", &Q2fac, &b_Pythia8_Q2fac);
	fChain->SetBranchAddress("Q2ren", &Q2ren, &b_Pythia8_Q2ren);
	fChain->SetBranchAddress("valence1", &valence1, &b_Pythia8_valence1);
	fChain->SetBranchAddress("valence2", &valence2, &b_Pythia8_valence2);
	fChain->SetBranchAddress("sHat", &sHat, &b_Pythia8_sHat);
	fChain->SetBranchAddress("tHat", &tHat, &b_Pythia8_tHat);
	fChain->SetBranchAddress("uHat", &uHat, &b_Pythia8_uHat);
	fChain->SetBranchAddress("ptHat", &ptHat, &b_Pythia8_ptHat);
	fChain->SetBranchAddress("thetaHat", &thetaHat, &b_Pythia8_thetaHat);
	fChain->SetBranchAddress("phiHat", &phiHat, &b_Pythia8_phiHat);
	fChain->SetBranchAddress("weight", &weight, &b_Pythia8_weight);
	Notify();
}

Bool_t genevents_new::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void genevents_new::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t genevents_new::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
#endif // #ifdef genevents_cxx
