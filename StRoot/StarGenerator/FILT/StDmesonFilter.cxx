//This filter is to find D0 and D+- mesons 

#include "StDmesonFilter.h"
#include "StarGenerator/EVENT/StarGenParticle.h"
#include "StarGenerator/EVENT/StarGenEvent.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#define DpmMass 1.869
#define d0Mass 1.864

using namespace std;
//_______________________________________________________________
StDmesonFilter::StDmesonFilter():StarFilterMaker("Dmeson")
{

  nEvents = new int;
  (*nEvents) = 0;

}
//_______________________________________________________________
StDmesonFilter::~StDmesonFilter()
{
  delete nEvents;
}
//_______________________________________________________________
Int_t StDmesonFilter::Filter( StarGenEvent *mEvent )
{
  //if(mEvent->GetNumberOfParticles <= 0) {return kError;}
  if(int(mEvent->GetNumberOfParticles()) <= 0)return -1;

  TIter Iterator = mEvent->IterAll();
  StarGenParticle *p = 0;

  // Find D meson ( id = 421), find its daughters
  int indexD0(0);
	int D0Daughter1(0), D0Daughter2(0);
	int D0Mother1(0), D0Mother2(0);
  int yesD0(0);

	int indexD0bar(0);
	int D0barDaughter1(0), D0barDaughter2(0);
	int D0barMother1(0), D0barMother2(0);
  int yesD0bar(0);

  int indexDpm(0);
	int DpmDaughter1(0), DpmDaughter2(0);
	int DpmMother1(0), DpmMother2(0);
  int yesDplus(0);

	int indexDpmbar(0);
	int DpmbarDaughter1(0), DpmbarDaughter2(0);
	int DpmbarMother1(0), DpmbarMother2(0);
  int yesDminus(0);

  while( ( p = (StarGenParticle*)Iterator.Next() ) )
  {
    // D0
    if(p->GetId() == 421)
    { 
			indexD0 = p->GetIndex();
			cout<<"D0 found!"<<endl;
			//Find hadronization products of c
			D0Daughter1 = p->GetFirstDaughter();
			D0Daughter2 = p->GetLastDaughter();
			D0Mother1 = p->GetFirstMother();
			D0Mother2 = p->GetLastMother();

      yesD0 = 1;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<D0Daughter1<<"\t"<<D0Daughter2<<"\t"<<D0Mother1<<"\t"<<D0Mother2<<endl;
		}
    // D0 bar
		if(p->GetId() == -421)
    { 
			indexD0bar = p->GetIndex();
			cout<<"D0bar found!"<<endl;
			//Find hadronization products of c
			D0barDaughter1 = p->GetFirstDaughter();
			D0barDaughter2 = p->GetLastDaughter();
			D0barMother1 = p->GetFirstMother();
			D0barMother2 = p->GetLastMother();

      yesD0bar = 1;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<D0barDaughter1<<"\t"<<D0barDaughter2<<"\t"<<D0barMother1<<"\t"<<D0barMother2<<endl;
		}
    
    
    if(p->GetId() == 411)
    { 
			indexDpm = p->GetIndex();
			cout<<"D+ found!"<<endl;
			//Find hadronization products of c
			DpmDaughter1 = p->GetFirstDaughter();
			DpmDaughter2 = p->GetLastDaughter();
			DpmMother1 = p->GetFirstMother();
			DpmMother2 = p->GetLastMother();

      yesDplus = 1;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<DpmDaughter1<<"\t"<<DpmDaughter2<<"\t"<<DpmMother1<<"\t"<<DpmMother2<<endl;
		}
    // D-
		if(p->GetId() == -411)
    { 
			indexDpmbar = p->GetIndex();
			cout<<"D- found!"<<endl;
			//Find hadronization products of c
			DpmbarDaughter1 = p->GetFirstDaughter();
			DpmbarDaughter2 = p->GetLastDaughter();
			DpmbarMother1 = p->GetFirstMother();
			DpmbarMother2 = p->GetLastMother();

      yesDminus = 1;

			cout<<"d1\td2\tm1\tm2"<<endl;
			cout<<DpmbarDaughter1<<"\t"<<DpmbarDaughter2<<"\t"<<DpmbarMother1<<"\t"<<DpmbarMother2<<endl;
		}
  }// particle loop

//do not need cuts on daugter kinematics - want all D0 and D+- mesons

//loop for daughters
  TIter Iterator3 = mEvent->IterAll();
  StarGenParticle *p3 = 0;  
	
  int mothercbd0(0);
  int mothercbd0bar(0);

  int mothercbDplus(0);
  int mothercbDminus(0);

  while( ( p3 = (StarGenParticle*)Iterator3.Next() ) )
  {
	  if(yesD0 && (p3->GetIndex() == D0Mother1) )//&& (p3->GetIndex() <= DmesonMother2))
	  {
		  mothercbd0 = 1;
		  cout<<"D0 meson mother found = "<<p3->GetId()<<endl;	
	  }
    if(yesD0bar && (p3->GetIndex() == D0barMother1) )//&& (p3->GetIndex() <= DmesonMother2))
	  {
		  mothercbd0bar = 1;
		  cout<<"D0bar meson mother found = "<<p3->GetId()<<endl;	
	  }

    if(yesDplus && (p3->GetIndex() == DpmMother1) )//&& (p3->GetIndex() <= DmesonMother2))
	  {
		  mothercbDplus = 1;
		  cout<<"D+ meson mother found = "<<p3->GetId()<<endl;	
	  }
    if(yesDminus && (p3->GetIndex() == DpmbarMother1) )//&& (p3->GetIndex() <= DmesonMother2))
	  {
		  mothercbDminus = 1;
		  cout<<"D- meson mother found = "<<p3->GetId()<<endl;	
	  }
  }// particle loop 2


  // No D meson found
  
  if ( mothercbd0 == 0 && mothercbd0bar == 0 && mothercbDplus == 0 && mothercbDminus == 0)
  {
    return StarGenEvent::kReject;
  }

  return StarGenEvent::kAccept;
}
