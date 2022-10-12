//This filter is to find D+-->k- pi+ pi+ process. 

#include "StDiDpmPiPiKFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include "StarGenerator/EVENT/StarGenEvent.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#define muMass 0.106
#define piMass 0.139570
#define kMass 0.493677
#define DpmMass 1.869

using namespace std;
//_______________________________________________________________
StDiDpmPiPiKFilter::StDiDpmPiPiKFilter():StarFilterMaker("DiDpmPiPiK")
{

	nEvents = new int;
	(*nEvents) = 0;

	mMotherKineFlag = false;
	mDaughterKineFlag = false;
	mParent1YMin = -1.;
	mParent1YMax = 1.;
	mParent2YMin = -1.;
	mParent2YMax = 1.;

	mMode = 1;

	mPtMin1  = 0.3;
	mPtMax1  = 1e5;
	mEtaMin1 = -1.;
	mEtaMax1 = 1.;
	mPhiMin1 = -TMath::Pi();
	mPhiMax1 = TMath::Pi();
	mPtMin2  = 0.3;
	mPtMax2  = 1e5;
	mEtaMin2 = -1.;
	mEtaMax2 = 1.;
	mPhiMin2 = -TMath::Pi();
	mPhiMax2 = TMath::Pi();
}
//_______________________________________________________________
StDiDpmPiPiKFilter::~StDiDpmPiPiKFilter()
{
	delete nEvents;
}
//_______________________________________________________________
Int_t StDiDpmPiPiKFilter::Filter( StarGenEvent *mEvent )
{
	//if(mEvent->GetNumberOfParticles <= 0) {return kError;}
	if(int(mEvent->GetNumberOfParticles()) <= 0)return -1;

	TIter Iterator = mEvent->IterAll();
	StarGenParticle *p = 0;

	// Find D+- meson ( id = 411 and -411), find their daughters
	int indexDpm(0);
	int DpmDaughter1(0), DpmDaughter2(0);
	int DpmMother1(0), DpmMother2(0);

	int indexDpmbar(0);
	int DpmbarDaughter1(0), DpmbarDaughter2(0);
	int DpmbarMother1(0), DpmbarMother2(0);

	// Flags for mother Dpm mesons
	int mPassDpmMotherKineFlag = 0;
	int mPassDpmbarMotherKineFlag = 0;
	TLorentzVector mot1(0,0,0,0);	
	TLorentzVector mot2(0,0,0,0);	

	while( ( p = (StarGenParticle*)Iterator.Next() ) ){
		if(p->GetId() == 411){ // Dpm
			indexDpm = p->GetIndex();
			cout<<"D+ found!"<<endl;
			//Find hadronization products of c
			DpmDaughter1 = p->GetFirstDaughter();
			DpmDaughter2 = p->GetLastDaughter();
			DpmMother1 = p->GetFirstMother();
			DpmMother2 = p->GetLastMother();

			mot1.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),DpmMass);
			cout<<"Dpm y="<<mot1.Rapidity()<<endl;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<DpmDaughter1<<"\t"<<DpmDaughter2<<"\t"<<DpmMother1<<"\t"<<DpmMother2<<endl;
		}
		if(p->GetId() == -411){ // Dpm bar
			indexDpmbar = p->GetIndex();
			cout<<"D- found!"<<endl;
			//Find hadronization products of c
			DpmbarDaughter1 = p->GetFirstDaughter();
			DpmbarDaughter2 = p->GetLastDaughter();
			DpmbarMother1 = p->GetFirstMother();
			DpmbarMother2 = p->GetLastMother();

			mot2.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),DpmMass);
			cout<<"Dpm_bar y="<<mot2.Rapidity()<<endl;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<DpmbarDaughter1<<"\t"<<DpmbarDaughter2<<"\t"<<DpmbarMother1<<"\t"<<DpmbarMother2<<endl;
		}
	}// particle loop
	if(!indexDpm || !indexDpmbar) return  StarGenEvent::kReject;

	//loop for daughters
	TIter Iterator2 = mEvent->IterAll();
	StarGenParticle *p2 = 0;
	int yesDpm(0);
	int DpmdaughterK(0);
	int DpmdaughterMu1(0);
  int DpmdaughterMu2(0);

	int yesDpmbar(0);
	int DpmbardaughterK(0);
	int DpmbardaughterMu1(0);
  int DpmbardaughterMu2(0);

	// Flags for daughter muons
	int mPassDpmDaughterKineFlag = 0;
	int mPassDpmbarDaughterKineFlag = 0;
	TLorentzVector dau1pi1(0,0,0,0);
  TLorentzVector dau1pi2(0,0,0,0);
	TLorentzVector dau1k (0,0,0,0);	
	TLorentzVector dau2pi1(0,0,0,0);
  TLorentzVector dau2pi2(0,0,0,0);	
	TLorentzVector dau2k (0,0,0,0);


	while( ( p2 = (StarGenParticle*)Iterator2.Next() ) )
	{
		if((p2->GetIndex() >= DpmDaughter1) && (p2->GetIndex() <= DpmDaughter2)) //go over all daughters
		{
			cout<<"Dpm daughter:"<<p2->GetId()<<endl;
			if(p2->GetId() == +211) //pions
			{
        if(DpmdaughterMu1 == 0)
        {
          DpmdaughterMu1=1;
				  cout<<"first pi+ found!"<<endl;
				  dau1pi1.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);        
        }
        else
        {
          DpmdaughterMu2=1;
				  cout<<"second pi+ found!"<<endl;
				  dau1pi2.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);        
        }
				
			}
			if(p2->GetId() == -321) //kaons
			{
				DpmdaughterK=1;
				cout<<"K- found!"<<endl;
				dau1k.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
		if((p2->GetIndex() >= DpmbarDaughter1) && (p2->GetIndex() <= DpmbarDaughter2))
		{
			cout<<"Dpm_bar daughter:"<<p2->GetId()<<endl;
			if(p2->GetId() == -211)
			{
        if(DpmbardaughterMu1 == 0)
        {
          DpmbardaughterMu1=1;
				  cout<<"first pi- found!"<<endl;
				  dau2pi1.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);        
        }
        else
        {
          DpmbardaughterMu2=1;
				  cout<<"second pi- found!"<<endl;
				  dau2pi2.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);        
        }

			}
			if(p2->GetId() == +321)
			{
				DpmbardaughterK=1;
				cout<<"K+ found!"<<endl;
				dau2k.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
	}// particle loop 2

	//loop for daughters
	TIter Iterator3 = mEvent->IterAll();
	StarGenParticle *p3 = 0;
	if(DpmdaughterK && DpmdaughterMu1 && DpmdaughterMu2)// &&mothercb)
		yesDpm=1;
	if(DpmbardaughterK && DpmbardaughterMu1 && DpmbardaughterMu2)// &&mothercb)
		yesDpmbar=1;

	int mothercbDpm(0);
	int mothercbDpmbar(0);
	while( ( p3 = (StarGenParticle*)Iterator3.Next() ) )
	{
		if(yesDpm && (p3->GetIndex() == DpmMother1) )//&& (p3->GetIndex() <= DpmMother2))
			//if(yesDpmbar && (p3->GetIndex() >= DpmMother1) && (p3->GetIndex() <= DpmMother2))
		{
			mothercbDpm = 1;
			cout<<"Dpm mother found = "<<p3->GetId()<<endl;
		}

		//if(yesDpmbar && (p3->GetIndex() >= DpmbarDaughter1) && (p3->GetIndex() <= DpmbarDaughter2))
		//{
		//	cout<<"********************************************"<<endl;
		//	cout<<"daughter found = "<<p3->GetId()<<endl;
		//	cout<<"********************************************"<<endl;
		//}

		if(yesDpmbar && (p3->GetIndex() == DpmbarMother1))// && (p3->GetIndex() <= DpmbarMother2))
			//if(yesDpmbar && (p3->GetIndex() >= DpmbarMother1) && (p3->GetIndex() <= DpmbarMother2))
		{
			mothercbDpmbar = 1;
			cout<<"Dpm_bar mother found = "<<p3->GetId()<<endl;
		}
	}// particle loop 2

	// Check Dpm kinematics
	if(mMotherKineFlag &&
			mot1.Rapidity()>mParent1YMin && mot1.Rapidity()<mParent1YMax
		) mPassDpmMotherKineFlag = 1;

	if(mMotherKineFlag &&
			mot2.Rapidity()>mParent2YMin && mot2.Rapidity()<mParent2YMax 
		) mPassDpmbarMotherKineFlag = 1;

	// Check Dpm daughter kinematics
	if(mDaughterKineFlag){
		if(dau1pi1.Pt()>mPtMin1 && dau1pi1.Pt()<mPtMax1 && dau1pi1.Eta()>mEtaMin1 && dau1pi1.Eta()<mEtaMax1 &&	dau1pi1.Phi()>mPhiMin1 && dau1pi1.Phi()<mPhiMax1 &&
       dau1pi2.Pt()>mPtMin1 && dau1pi2.Pt()<mPtMax1 && dau1pi2.Eta()>mEtaMin1 && dau1pi2.Eta()<mEtaMax1 &&	dau1pi2.Phi()>mPhiMin1 && dau1pi2.Phi()<mPhiMax1 &&
				dau1k.Pt() >mPtMin2 && dau1k.Pt()<mPtMax2 && dau1k.Eta() >mEtaMin2 && dau1k.Eta() <mEtaMax2 && dau1k.Phi() >mPhiMin2 && dau1k.Phi() <mPhiMax2) {
			mPassDpmDaughterKineFlag = 1;
		}
		if(dau2pi1.Pt()>mPtMin1 && dau2pi1.Pt()<mPtMax1 && dau2pi1.Eta()>mEtaMin1 && dau2pi1.Eta()<mEtaMax1 &&	dau2pi1.Phi()>mPhiMin1 && dau2pi1.Phi()<mPhiMax1 &&
       dau2pi2.Pt()>mPtMin1 && dau2pi2.Pt()<mPtMax1 && dau2pi2.Eta()>mEtaMin1 && dau2pi2.Eta()<mEtaMax1 &&	dau2pi2.Phi()>mPhiMin1 && dau2pi2.Phi()<mPhiMax1 &&
				dau2k.Pt() >mPtMin2 && dau2k.Pt()<mPtMax2 && dau2k.Eta() >mEtaMin2 && dau2k.Eta() <mEtaMax2 && dau2k.Phi() >mPhiMin2 && dau2k.Phi() <mPhiMax2)
		{
			mPassDpmbarDaughterKineFlag = 1;
		}
	}

	// No Di Dpm meson found
	if (
			(
			 (mMotherKineFlag && mPassDpmMotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassDpmDaughterKineFlag == 0)  || 
			 mothercbDpm == 0 
			) ||
			(
			 (mMotherKineFlag && mPassDpmbarMotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassDpmbarDaughterKineFlag == 0)  || 
			 mothercbDpmbar == 0
			)
		 ) {
		return StarGenEvent::kReject;
	}
			 //cout << mPassDpmMotherKineFlag <<" "
			 //<<  mPassDpmDaughterKineFlag <<" " 
			 //<< mothercbDpm <<endl;
			 //cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
			 //cout << mPassDpmbarMotherKineFlag <<" "
			 //<<  mPassDpmbarDaughterKineFlag <<" " 
			 //<< mothercbDpmbar <<endl;

		return StarGenEvent::kAccept;
}

void  StDiDpmPiPiKFilter::SetDauKine(double ptMin1, double ptMax1, double etaMin1, double etaMax1, double phiMin1, double phiMax1, double ptMin2, double ptMax2, double etaMin2, double etaMax2, double phiMin2, double phiMax2){
  //pion cuts (same for both pions)
	mPtMin1 = ptMin1; mPtMax1 = ptMax1;
	mEtaMin1 = etaMin1; mEtaMax1 = etaMax1;
	mPhiMin1 = phiMin1; mPhiMax1 = phiMax1;

  //kaon cuts
	mPtMin2 = ptMin2; mPtMax2 = ptMax2;
	mEtaMin2 = etaMin2; mEtaMax2 = etaMax2;
	mPhiMin2 = phiMin2; mPhiMax2 = phiMax2;

	mDaughterKineFlag = true;
}

void StDiDpmPiPiKFilter::SetParentRapidities(double yMin1, double yMax1, double yMin2, double yMax2){
	mParent1YMin = yMin1;
	mParent1YMax = yMax1;
	mParent2YMin = yMin2;
	mParent2YMax = yMax2;
	mMotherKineFlag = true;
}
void StDiDpmPiPiKFilter::SetMode(int val){
	mMode = val;
	cout<<"StDiDpmPiPiKFilter::SetMode to mMode = "<<mMode<<endl;
}
