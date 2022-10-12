//This filter is to find D0->k- + pi + process. 

#include "StDiD0PiKFilter.h"
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
#define d0Mass 1.864

using namespace std;
//_______________________________________________________________
StDiD0PiKFilter::StDiD0PiKFilter():StarFilterMaker("DiD0PiK")
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

	mPtMin1  = 0.2;
	mPtMax1  = 1e5;
	mEtaMin1 = -1.;
	mEtaMax1 = 1.;
	mPhiMin1 = -TMath::Pi();
	mPhiMax1 = TMath::Pi();
	mPtMin2  = 0.2;
	mPtMax2  = 1e5;
	mEtaMin2 = -1.;
	mEtaMax2 = 1.;
	mPhiMin2 = -TMath::Pi();
	mPhiMax2 = TMath::Pi();
}
//_______________________________________________________________
StDiD0PiKFilter::~StDiD0PiKFilter()
{
	delete nEvents;
}
//_______________________________________________________________
Int_t StDiD0PiKFilter::Filter( StarGenEvent *mEvent )
{
	//if(mEvent->GetNumberOfParticles <= 0) {return kError;}
	if(int(mEvent->GetNumberOfParticles()) <= 0)return -1;

	TIter Iterator = mEvent->IterAll();
	StarGenParticle *p = 0;

	// Find D meson ( id = 421 and -421), find their daughters
	int indexD0(0);
	int D0Daughter1(0), D0Daughter2(0);
	int D0Mother1(0), D0Mother2(0);

	int indexD0bar(0);
	int D0barDaughter1(0), D0barDaughter2(0);
	int D0barMother1(0), D0barMother2(0);

	// Flags for mother D0 mesons
	int mPassD0MotherKineFlag = 0;
	int mPassD0barMotherKineFlag = 0;
	TLorentzVector mot1(0,0,0,0);	
	TLorentzVector mot2(0,0,0,0);	

	while( ( p = (StarGenParticle*)Iterator.Next() ) ){
		if(p->GetId() == 421){ // D0
			indexD0 = p->GetIndex();
			cout<<"D0 found!"<<endl;
			//Find hadronization products of c
			D0Daughter1 = p->GetFirstDaughter();
			D0Daughter2 = p->GetLastDaughter();
			D0Mother1 = p->GetFirstMother();
			D0Mother2 = p->GetLastMother();

			mot1.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),d0Mass);
			cout<<"D0 y="<<mot1.Rapidity()<<endl;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<D0Daughter1<<"\t"<<D0Daughter2<<"\t"<<D0Mother1<<"\t"<<D0Mother2<<endl;
		}
		if(p->GetId() == -421){ // D0 bar
			indexD0bar = p->GetIndex();
			cout<<"D0bar found!"<<endl;
			//Find hadronization products of c
			D0barDaughter1 = p->GetFirstDaughter();
			D0barDaughter2 = p->GetLastDaughter();
			D0barMother1 = p->GetFirstMother();
			D0barMother2 = p->GetLastMother();

			mot2.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),d0Mass);
			cout<<"D0_bar y="<<mot2.Rapidity()<<endl;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<D0barDaughter1<<"\t"<<D0barDaughter2<<"\t"<<D0barMother1<<"\t"<<D0barMother2<<endl;
		}
	}// particle loop
	if(!indexD0 || !indexD0bar) return  StarGenEvent::kReject;

	//loop for daughters
	TIter Iterator2 = mEvent->IterAll();
	StarGenParticle *p2 = 0;
	int yesD0(0);
	int D0daughterK(0);
	int D0daughterMu(0);

	int yesD0bar(0);
	int D0bardaughterK(0);
	int D0bardaughterMu(0);

	// Flags for daughter muons
	int mPassD0DaughterKineFlag = 0;
	int mPassD0barDaughterKineFlag = 0;
	TLorentzVector dau1pi(0,0,0,0);	
	TLorentzVector dau1k (0,0,0,0);	
	TLorentzVector dau2pi(0,0,0,0);	
	TLorentzVector dau2k (0,0,0,0);	

	while( ( p2 = (StarGenParticle*)Iterator2.Next() ) )
	{
		if((p2->GetIndex() >= D0Daughter1) && (p2->GetIndex() <= D0Daughter2))
		{
			cout<<"D0 daughter:"<<p2->GetId()<<endl;
			if(p2->GetId() == +211)
			{
				D0daughterMu=1;
				cout<<"pi- found!"<<endl;
				dau1pi.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);
			}
			if(p2->GetId() == -321)
			{
				D0daughterK=1;
				cout<<"K- found!"<<endl;
				dau1k.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
		if((p2->GetIndex() >= D0barDaughter1) && (p2->GetIndex() <= D0barDaughter2))
		{
			cout<<"D0_bar daughter:"<<p2->GetId()<<endl;
			if(p2->GetId() == -211)
			{
				D0bardaughterMu=1;
				cout<<"pi- found!"<<endl;
				dau2pi.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);
			}
			if(p2->GetId() == +321)
			{
				D0bardaughterK=1;
				cout<<"K+ found!"<<endl;
				dau2k.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
	}// particle loop 2

	//loop for daughters
	TIter Iterator3 = mEvent->IterAll();
	StarGenParticle *p3 = 0;
	if(D0daughterK && D0daughterMu)// &&mothercb)
		yesD0=1;
	if(D0bardaughterK && D0bardaughterMu)// &&mothercb)
		yesD0bar=1;

	int mothercbd0(0);
	int mothercbd0bar(0);
	while( ( p3 = (StarGenParticle*)Iterator3.Next() ) )
	{
		if(yesD0 && (p3->GetIndex() == D0Mother1) )//&& (p3->GetIndex() <= D0Mother2))
			//if(yesD0bar && (p3->GetIndex() >= D0Mother1) && (p3->GetIndex() <= D0Mother2))
		{
			mothercbd0 = 1;
			cout<<"D0 mother found = "<<p3->GetId()<<endl;
		}

		//if(yesD0bar && (p3->GetIndex() >= D0barDaughter1) && (p3->GetIndex() <= D0barDaughter2))
		//{
		//	cout<<"********************************************"<<endl;
		//	cout<<"daughter found = "<<p3->GetId()<<endl;
		//	cout<<"********************************************"<<endl;
		//}

		if(yesD0bar && (p3->GetIndex() == D0barMother1))// && (p3->GetIndex() <= D0barMother2))
			//if(yesD0bar && (p3->GetIndex() >= D0barMother1) && (p3->GetIndex() <= D0barMother2))
		{
			mothercbd0bar = 1;
			cout<<"D0_bar mother found = "<<p3->GetId()<<endl;
		}
	}// particle loop 2

	// Check D0 kinematics
	if(mMotherKineFlag &&
			mot1.Rapidity()>mParent1YMin && mot1.Rapidity()<mParent1YMax
		) mPassD0MotherKineFlag = 1;

	if(mMotherKineFlag &&
			mot2.Rapidity()>mParent2YMin && mot2.Rapidity()<mParent2YMax 
		) mPassD0barMotherKineFlag = 1;

	// Check D0 daughter kinematics
	if(mDaughterKineFlag){
		if(dau1pi.Pt()>mPtMin1&&dau1pi.Pt()<mPtMax1 && dau1pi.Eta()>mEtaMin1&&dau1pi.Eta()<mEtaMax1 &&	dau1pi.Phi()>mPhiMin1&&dau1pi.Phi()<mPhiMax1 &&
				dau1k.Pt() >mPtMin2&&dau1k.Pt()<mPtMax2 && dau1k.Eta() >mEtaMin2&&dau1k.Eta() <mEtaMax2 && dau1k.Phi() >mPhiMin2&&dau1k.Phi() <mPhiMax2) {
			mPassD0DaughterKineFlag = 1;
		}
		if(dau2pi.Pt()>mPtMin1&&dau2pi.Pt()<mPtMax1 && dau2pi.Eta()>mEtaMin1&&dau2pi.Eta()<mEtaMax1 &&	dau2pi.Phi()>mPhiMin1&&dau2pi.Phi()<mPhiMax1 &&
				dau2k.Pt() >mPtMin2&&dau2k.Pt()<mPtMax2 && dau2k.Eta() >mEtaMin2&&dau2k.Eta() <mEtaMax2 && dau2k.Phi() >mPhiMin2&&dau2k.Phi() <mPhiMax2)
		{
			mPassD0barDaughterKineFlag = 1;
		}
	}

	// No Di D0 meson found
	if (
			(
			 (mMotherKineFlag && mPassD0MotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassD0DaughterKineFlag == 0)  || 
			 mothercbd0 == 0 
			) 
        && // changed from || -> Finds single D0 or D0-bar
			(
			 (mMotherKineFlag && mPassD0barMotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassD0barDaughterKineFlag == 0)  || 
			 mothercbd0bar == 0
			)
		 ) {
		return StarGenEvent::kReject;
	}
			 //cout << mPassD0MotherKineFlag <<" "
			 //<<  mPassD0DaughterKineFlag <<" " 
			 //<< mothercbd0 <<endl;
			 //cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
			 //cout << mPassD0barMotherKineFlag <<" "
			 //<<  mPassD0barDaughterKineFlag <<" " 
			 //<< mothercbd0bar <<endl;

		return StarGenEvent::kAccept;
}

void  StDiD0PiKFilter::SetDauKine(double ptMin1, double ptMax1, double etaMin1, double etaMax1, double phiMin1, double phiMax1, double ptMin2, double ptMax2, double etaMin2, double etaMax2, double phiMin2, double phiMax2){
	mPtMin1 = ptMin1; mPtMax1 = ptMax1;
	mEtaMin1 = etaMin1; mEtaMax1 = etaMax1;
	mPhiMin1 = phiMin1; mPhiMax1 = phiMax1;
	mPtMin2 = ptMin2; mPtMax2 = ptMax2;
	mEtaMin2 = etaMin2; mEtaMax2 = etaMax2;
	mPhiMin2 = phiMin2; mPhiMax2 = phiMax2;
	mDaughterKineFlag = true;
}

void StDiD0PiKFilter::SetParentRapidities(double yMin1, double yMax1, double yMin2, double yMax2){
	mParent1YMin = yMin1;
	mParent1YMax = yMax1;
	mParent2YMin = yMin2;
	mParent2YMax = yMax2;
	mMotherKineFlag = true;
}
void StDiD0PiKFilter::SetMode(int val){
	mMode = val;
	cout<<"StDiD0PiKFilter::SetMode to mMode = "<<mMode<<endl;
}
