//This filter is to find Lambda->k- + pi + process. 

#include "StLambdaPiPFilter.h"
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
#define LambdaMass 1.115683 //mass in GeV/c^2 from latest PDG

using namespace std;
//_______________________________________________________________
StLambdaPiPFilter::StLambdaPiPFilter():StarFilterMaker("LambdaPiP")
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
StLambdaPiPFilter::~StLambdaPiPFilter()
{
	delete nEvents;
}
//_______________________________________________________________
Int_t StLambdaPiPFilter::Filter( StarGenEvent *mEvent )
{
	//if(mEvent->GetNumberOfParticles <= 0) {return kError;}
	if(int(mEvent->GetNumberOfParticles()) <= 0)return -1;

	TIter Iterator = mEvent->IterAll();
	StarGenParticle *p = 0;

	// Find Lambda baryon ( id = 3122 and -3122), find their daughters
	int indexLambda(0);
	int LambdaDaughter1(0), LambdaDaughter2(0);
	int LambdaMother1(0), LambdaMother2(0);

	int indexLambdabar(0);
	int LambdabarDaughter1(0), LambdabarDaughter2(0);
	int LambdabarMother1(0), LambdabarMother2(0);

	// Flags for mother Lambda mesons
	int mPassLambdaMotherKineFlag = 0;
	int mPassLambdabarMotherKineFlag = 0;
	TLorentzVector mot1(0,0,0,0);	
	TLorentzVector mot2(0,0,0,0);	

	while( ( p = (StarGenParticle*)Iterator.Next() ) ){
		if(p->GetId() == 3122)
    { // Lambda
			indexLambda = p->GetIndex();
			//cout<<"Lambda found!"<<endl;
			//Find hadronization products of c
			LambdaDaughter1 = p->GetFirstDaughter();
			LambdaDaughter2 = p->GetLastDaughter();
			LambdaMother1 = p->GetFirstMother();
			LambdaMother2 = p->GetLastMother();

			mot1.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),LambdaMass);
			//cout<<"Lambda y="<<mot1.Rapidity()<<endl;

			//cout<<"d1\td2\tm1\tm2"<<endl;
			//cout<<LambdaDaughter1<<"\t"<<LambdaDaughter2<<"\t"<<LambdaMother1<<"\t"<<LambdaMother2<<endl;
		}
		if(p->GetId() == -3122)
    { // Lambda bar
			indexLambdabar = p->GetIndex();
			//cout<<"Lambdabar found!"<<endl;
			//Find hadronization products of c
			LambdabarDaughter1 = p->GetFirstDaughter();
			LambdabarDaughter2 = p->GetLastDaughter();
			LambdabarMother1 = p->GetFirstMother();
			LambdabarMother2 = p->GetLastMother();

			mot2.SetXYZM(p->GetPx(),p->GetPy(),p->GetPz(),LambdaMass);
			//cout<<"Lambda_bar y="<<mot2.Rapidity()<<endl;

			//cout<<"d1\td2\tm1\tm2"<<endl;
			//cout<<LambdabarDaughter1<<"\t"<<LambdabarDaughter2<<"\t"<<LambdabarMother1<<"\t"<<LambdabarMother2<<endl;
		}
	}// particle loop
	if(!indexLambda || !indexLambdabar) return  StarGenEvent::kReject;

	//loop for daughters
	TIter Iterator2 = mEvent->IterAll();
	StarGenParticle *p2 = 0;
	int yesLambda(0);
	int LambdadaughterP(0);
	int LambdadaughterPi(0);

	int yesLambdabar(0);
	int LambdabardaughterP(0);
	int LambdabardaughterPi(0);

	// Flags for daughter muons
	int mPassLambdadaughterKineFlag = 0;
	int mPassLambdabardaughterKineFlag = 0;
	TLorentzVector dau1pi(0,0,0,0);	
	TLorentzVector dau1p (0,0,0,0);	
	TLorentzVector dau2pi(0,0,0,0);	
	TLorentzVector dau2p (0,0,0,0);	

	while( ( p2 = (StarGenParticle*)Iterator2.Next() ) )
	{
		if((p2->GetIndex() >= LambdaDaughter1) && (p2->GetIndex() <= LambdaDaughter2))
		{
			//cout<<"Lambda daughter:"<<p2->GetId()<<endl;
			if( fabs(p2->GetId()) == 211)
			{
				LambdadaughterPi=1;
				//cout<<"pi- found!"<<endl;
				dau1pi.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);
			}
			if( fabs(p2->GetId()) == 2212)
			{
				LambdadaughterP=1;
				//cout<<"p found!"<<endl;
				dau1p.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
		if((p2->GetIndex() >= LambdabarDaughter1) && (p2->GetIndex() <= LambdabarDaughter2))
		{
			//cout<<"Lambda_bar daughter:"<<p2->GetId()<<endl;
			if( fabs(p2->GetId()) == 211)
			{
				LambdabardaughterPi=1;
				//cout<<"pi- found!"<<endl;
				dau2pi.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),piMass);
			}
			if( fabs(p2->GetId()) == 2212)
			{
				LambdabardaughterP=1;
				//cout<<"p-bar found!"<<endl;
				dau2p.SetXYZM(p2->GetPx(),p2->GetPy(),p2->GetPz(),kMass);
			}
		}
	}// particle loop 2

	//loop for daughters
	TIter Iterator3 = mEvent->IterAll();
	StarGenParticle *p3 = 0;
	if(LambdadaughterP && LambdadaughterPi)// &&mothercb)
		yesLambda=1;
	if(LambdabardaughterP && LambdabardaughterPi)// &&mothercb)
		yesLambdabar=1;

	int mothercbLambda(0);
	int mothercbLambdabar(0);
	while( ( p3 = (StarGenParticle*)Iterator3.Next() ) )
	{
		if(yesLambda && (p3->GetIndex() == LambdaMother1) )//&& (p3->GetIndex() <= LambdaMother2))
			//if(yesLambdabar && (p3->GetIndex() >= LambdaMother1) && (p3->GetIndex() <= LambdaMother2))
		{
			mothercbLambda = 1;
			//cout<<"Lambda mother found = "<<p3->GetId()<<endl;
		}

		//if(yesLambdabar && (p3->GetIndex() >= LambdabarDaughter1) && (p3->GetIndex() <= LambdabarDaughter2))
		//{
		//	//cout<<"********************************************"<<endl;
		//	//cout<<"daughter found = "<<p3->GetId()<<endl;
		//	//cout<<"********************************************"<<endl;
		//}

		if(yesLambdabar && (p3->GetIndex() == LambdabarMother1))// && (p3->GetIndex() <= LambdabarMother2))
			//if(yesLambdabar && (p3->GetIndex() >= LambdabarMother1) && (p3->GetIndex() <= LambdabarMother2))
		{
			mothercbLambdabar = 1;
			//cout<<"Lambda_bar mother found = "<<p3->GetId()<<endl;
		}
	}// particle loop 2

	// Check Lambda kinematics
	if(mMotherKineFlag &&	mot1.Rapidity()>mParent1YMin && mot1.Rapidity()<mParent1YMax) mPassLambdaMotherKineFlag = 1;

	if(mMotherKineFlag && mot2.Rapidity()>mParent2YMin && mot2.Rapidity()<mParent2YMax) mPassLambdabarMotherKineFlag = 1;

	// Check Lambda daughter kinematics
	if(mDaughterKineFlag)
  {
		if(dau1pi.Pt()>mPtMin1&&dau1pi.Pt()<mPtMax1 && dau1pi.Eta()>mEtaMin1&&dau1pi.Eta()<mEtaMax1 &&	dau1pi.Phi()>mPhiMin1&&dau1pi.Phi()<mPhiMax1 &&
			 dau1p.Pt() >mPtMin2&&dau1p.Pt()<mPtMax2 && dau1p.Eta() >mEtaMin2&&dau1p.Eta() <mEtaMax2 && dau1p.Phi() >mPhiMin2&&dau1p.Phi() <mPhiMax2) 
    {
			mPassLambdadaughterKineFlag = 1;
		}
		if(dau2pi.Pt()>mPtMin1&&dau2pi.Pt()<mPtMax1 && dau2pi.Eta()>mEtaMin1&&dau2pi.Eta()<mEtaMax1 &&	dau2pi.Phi()>mPhiMin1&&dau2pi.Phi()<mPhiMax1 &&
			 dau2p.Pt() >mPtMin2&&dau2p.Pt()<mPtMax2 && dau2p.Eta() >mEtaMin2&&dau2p.Eta() <mEtaMax2 && dau2p.Phi() >mPhiMin2&&dau2p.Phi() <mPhiMax2)
		{
			mPassLambdabardaughterKineFlag = 1;
		}
	}

	// Lambda or Lambda-bar not found
	if (
			(
			 (mMotherKineFlag && mPassLambdaMotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassLambdadaughterKineFlag == 0)  || 
			 mothercbLambda == 0 
			) 
      &&
			(
			 (mMotherKineFlag && mPassLambdabarMotherKineFlag == 0)  || 
			 (mDaughterKineFlag && mPassLambdabardaughterKineFlag == 0)  || 
			 mothercbLambdabar == 0
			)
		 ) 
  {
		return StarGenEvent::kReject;
	}
			 ////cout << mPassLambdaMotherKineFlag <<" "
			 //<<  mPassLambdadaughterKineFlag <<" " 
			 //<< mothercbLambda <<endl;
			 ////cout<<"++++++++++++++++++++++++++++++++++++++++"<<endl;
			 ////cout << mPassLambdabarMotherKineFlag <<" "
			 //<<  mPassLambdabardaughterKineFlag <<" " 
			 //<< mothercbLambdabar <<endl;
  if((mPassLambdaMotherKineFlag == 1) &&
		 (mPassLambdadaughterKineFlag == 1) &&
		 mothercbLambda == 1 )
  {
    cout<<"Good Lambda found!"<<endl;
  }

  if((mPassLambdabarMotherKineFlag == 1) &&
		 (mPassLambdabardaughterKineFlag == 1) &&
		 mothercbLambdabar == 1 )
  {
    cout<<"Good Lambda-bar found!"<<endl;
  }
    

  return StarGenEvent::kAccept;
}

void  StLambdaPiPFilter::SetDauKine(double ptMin1, double ptMax1, double etaMin1, double etaMax1, double phiMin1, double phiMax1, double ptMin2, double ptMax2, double etaMin2, double etaMax2, double phiMin2, double phiMax2)
{
	mPtMin1 = ptMin1; mPtMax1 = ptMax1;
	mEtaMin1 = etaMin1; mEtaMax1 = etaMax1;
	mPhiMin1 = phiMin1; mPhiMax1 = phiMax1;
	mPtMin2 = ptMin2; mPtMax2 = ptMax2;
	mEtaMin2 = etaMin2; mEtaMax2 = etaMax2;
	mPhiMin2 = phiMin2; mPhiMax2 = phiMax2;
	mDaughterKineFlag = true;
}

void StLambdaPiPFilter::SetParentRapidities(double yMin1, double yMax1, double yMin2, double yMax2)
{
	mParent1YMin = yMin1;
	mParent1YMax = yMax1;
	mParent2YMin = yMin2;
	mParent2YMax = yMax2;
	mMotherKineFlag = true;
}

void StLambdaPiPFilter::SetMode(int val)
{
	mMode = val;
	//cout<<"StLambdaPiPFilter::SetMode to mMode = "<<mMode<<endl;
}
