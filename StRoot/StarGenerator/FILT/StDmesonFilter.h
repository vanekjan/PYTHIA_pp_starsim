#ifndef StDmesonFilter_h
#define StDmesonFilter_h

#include <vector>

/*!
  \class StarParticleFilter
  \brief Filter which requires one or more particles in the final state of the event record
 */

#include "StarFilterMaker.h"

class StDmesonFilter : public StarFilterMaker
{
public:
  StDmesonFilter(  );
  virtual ~StDmesonFilter();

  int Filter( StarGenEvent *event = 0 );



private:
protected:

  int *nEvents;


//#if 1 // we dont really need, but this triggers cons to create dictionary
  ClassDef(StDmesonFilter,0);
//#endif

};

#endif
