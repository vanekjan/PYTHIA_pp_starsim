#ifndef StLambdaPiPFilter_h
#define StLambdaPiPFilter_h

#include <vector>

/*!
  \class StarParticleFilter
  \brief Filter which requires one or more particles in the final state of the event record
 */

#include "StarFilterMaker.h"

class StLambdaPiPFilter : public StarFilterMaker
{
public:
  StLambdaPiPFilter(  );
  virtual ~StLambdaPiPFilter();

  int Filter( StarGenEvent *event = 0 );

  void SetDauKine(double ptMin1, double ptMax1, double etaMin1, double etaMax1, double phiMin1, double phiMax1, double ptMin2, double ptMax2, double etaMin2, double etaMax2, double phiMin2, double phiMax2);

  void SetParentRapidities(double yMin1, double yMax1, double yMin2, double yMax2);

  void SetMode(int val);


private:
protected:

  int *nEvents;

	bool mMotherKineFlag = false;
	bool mDaughterKineFlag = false;
	double mParent1YMin = -1.;
	double mParent1YMax = 1.;
	double mParent2YMin = -1.;
	double mParent2YMax = 1.;

	int mMode = 1;

	double mPtMin1  = 0.2;
	double mPtMax1  = 1e5;
	double mEtaMin1 = -1.;
	double mEtaMax1 = 1.;
	double mPhiMin1 = -TMath::Pi();
	double mPhiMax1 = TMath::Pi();
	double mPtMin2  = 0.2;
	double mPtMax2  = 1e5;
	double mEtaMin2 = -1.;
	double mEtaMax2 = 1.;
	double mPhiMin2 = -TMath::Pi();
	double mPhiMax2 = TMath::Pi();

//#if 1 // we dont really need, but this triggers cons to create dictionary
  ClassDef(StLambdaPiPFilter,0);
//#endif

};

#endif
