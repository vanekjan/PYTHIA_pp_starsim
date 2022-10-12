/**!
 * Example macro for running an event generator in standalone mode.
 *
 * Use cases:
 *   1) Generate input root files for geant simulation
 *   2) Generate root files for quick analysis.  
 *
 * NOTE:  By default, a number of particles are made stable (e.g. pi0, K0, etc...),
 *        to allow GEANT to handle the decays.  If this is not desired, change the
 *        flag below.
 *
 * Usage:
 *
 * root4star
 * .L standalone.pythia6.C
 * int nevents=100;
 * standalone( nevents )
 */
#include"TMath.h"

const bool GEANT_HANDLES_DECAYS = true;

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarPythia8;
StarPythia8* _pythia8;

class StDmesonFilter;
StDmesonFilter *sngDmesonFilter = 0;

class StDiD0PiKFilter;
StDiD0PiKFilter *DiD0PiKFilter = 0;

class StLambdaPiPFilter;
StLambdaPiPFilter *sngLambdaPiPFilter = 0;

Float_t ptHatMin = 0;
Float_t ptHatMax = 128;

// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; i++ ) {
    chain->Clear();
    chain->Make();
    //_primary -> event() -> Print(); //commented to reduce log file size
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia6( TString mode="pp:W", Int_t tune=320 )
{
  
  //  gSystem->Load( "libStarGeneratorPoolPythia6_4_23.so" );
  gSystem->Load( "libPythia6_4_28.so");

  StarPythia6 *pythia6 = new StarPythia6("pythia6");
  if ( mode=="pp:W" )
  {
    pythia6->SetFrame("CMS", 510.0 );
    pythia6->SetBlue("proton");
    pythia6->SetYell("proton");
    if ( tune ) pythia6->PyTune( tune );

    // Setup pythia process
    PySubs_t &pysubs = pythia6->pysubs();
    int& msel = pysubs.msel;
    msel = 12;
    
    // Setup other common block variables / array elements
    float& ckin3 = pysubs.ckin(3); 
    ckin3 = 4.0;

    //
    // Set particles to be stable so that the decay manager
    // can handle them in the starsim phase
    //
    pythia6 -> SetDecayFlag( +24, 0 ); // W+
    pythia6 -> SetDecayFlag( -24, 0 ); // W-
    pythia6 -> SetDecayFlag( +23, 0 ); // Z0
    pythia6 -> SetDecayFlag( -23, 0 ); // Z0
    pythia6 -> SetDecayFlag( +15, 0 ); // tau+
    pythia6 -> SetDecayFlag( -15, 0 ); // tau-
    

  }
  if ( mode == "pp:minbias" )
  {
    pythia6->SetFrame("CMS", 510.0 );
    pythia6->SetBlue("proton");
    pythia6->SetYell("proton");
    if ( tune ) pythia6->PyTune( tune );
  }
  if ( mode == "ep" )
  {
    Double_t pblue[]={0.,0.,30.0};
    Double_t pyell[]={0.,0.,-320.0};
    pythia6->SetFrame("3MOM", pblue, pyell );
    pythia6->SetBlue("e-");
    pythia6->SetYell("proton");
    if ( tune ) pythia6->PyTune( tune );
  }
    
  _primary->AddGenerator(pythia6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Pythia8( TString config="pp:W", Int_t collEnergy = 200,  const char* _library="libPythia8_1_62.so" )
{
  gSystem->Load( "/star/u/tdrk/software/PDF/LHAPDF-6.1.4/lib/libLHAPDF.so"  ); // LHAPDF needs to be called before PYTHIA8
  gSystem->Load( _library );

  //
  // Create the pythia 8 event generator and add it to 
  // the primary generator
  //
  StarPythia8 *pythia8 = new StarPythia8();    
  _pythia8=pythia8;

  pythia8->SetFrame("CMS", collEnergy);
  pythia8->SetBlue("proton");
  pythia8->SetYell("proton");

  if ( config=="pp:W" )
  {     
    pythia8->Set("WeakSingleBoson:all=off");
    pythia8->Set("WeakSingleBoson:ffbar2W=on");
    pythia8->Set("24:onMode=0");              // switch off all W+/- decaus
    pythia8->Set("24:onIfAny 11 -11");        // switch on for decays to e+/-      
  }

  if ( config=="pp:minbias" )
  {
    pythia8->Set("SoftQCD:minBias = on");
  }

  if ( config=="pp:minbiasLambda")
  {
    pythia8->Set("SoftQCD:minBias = on");

    //pythia8->Set("3122:mayDecay= on");
    //pythia8->Set("3122:onMode=1");
    //pythia8->Set("3122:onMode=0");
    //pythia8->Set("3122:OnIfMatch=2212 -211");
    //pythia8->Set("-3122:onMode=0");
    //pythia8->Set("-3122:OnIfMatch=-2212 211");
  }

 //D+- and D0 settings for PYTHIA8
  if ( config=="pp:dmeson" || config=="pp:sngdmeson"  )
  {
    //ccbar process
    pythia8->Set("SigmaProcess:renormScale2 = 3");
    pythia8->Set("SigmaProcess:factorScale2 = 3");
    pythia8->Set("SigmaProcess:renormMultFac = 2"); //2mT
    pythia8->Set("SigmaProcess:factorMultFac = 2");
    pythia8->Set(Form("PhaseSpace:pTHatMin = %f",ptHatMin));
    pythia8->Set(Form("PhaseSpace:pTHatMax = %f",ptHatMax));

    pythia8->Set("PDF:useLHAPDF = on");
    pythia8->Set("PDF:LHAPDFset = MRSTMCal.LHgrid");
    pythia8->Set("PDF:extrapolateLHAPDF = on");
    pythia8->Set("PartonLevel:MI = on");
    pythia8->Set("PartonLevel:ISR = on");
    pythia8->Set("BeamRemnants:primordialKT = on");
    pythia8->Set("PartonLevel:FSR = on");
    pythia8->Set("StringFlav:mesonCvector = 1.5");
    pythia8->Set("StringFlav:mesonBvector = 3");
    pythia8->Set("4:m0 = 1.43");
    pythia8->Set("5:m0 = 4.30");
    //pythia8->Set("HardQCD:all = on");

    //pythia8->Set("HardQCD:hardccbar = on"); // some how does not switch on qq2ccbar and qqbar2ccbar
    pythia8->Set("HardQCD:gg2ccbar = on");
    pythia8->Set("HardQCD:qqbar2ccbar = on");
    //pythia8->Set("HardQCD:gg2bbbar = on");
    //pythia8->Set("HardQCD:qqbar2bbbar = on");
/*
    pythia8->Set("411:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("411:OnIfMatch=-321 211 211"); // switch on any decay mode into K- pi+ pi+
    pythia8->Set("-411:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("-411:OnIfMatch=321 -211 -211"); // switch on any decay mode into  K+ pi- pi-
*/
  }
  //D+- and D0 settings for PYTHIA8 tune
  if ( config=="pp:dmesontune" || config=="pp:sngdmesontune"   )
  {

    //ccbar process
    pythia8->Set("HardQCD:hardccbar = on");

    pythia8->Set("Tune:pp = 6");

    //http://home.thep.lu.se/~torbjorn/pythia81html/Tunes.html
    //option 6 : "Tune 4Cx", based on tune 4C, but using the x-dependent matter profile,
    //MultipartonInteractions:bProfile = 4 and an increased MultipartonInteractions:pT0Ref [Cor11]. 

    pythia8->Set("SigmaProcess:renormScale2 = 3");
    pythia8->Set("SigmaProcess:factorScale2 = 3");
    pythia8->Set("SigmaProcess:renormMultFac = 2"); //2mT
    pythia8->Set("SigmaProcess:factorMultFac = 2");
    pythia8->Set(Form("PhaseSpace:pTHatMin = %f",ptHatMin));
    pythia8->Set(Form("PhaseSpace:pTHatMax = %f",ptHatMax));

    pythia8->Set("PDF:useLHAPDF = on");
    pythia8->Set("PDF:LHAPDFset = MRSTMCal.LHgrid");
    pythia8->Set("PDF:extrapolateLHAPDF = on");
    pythia8->Set("PartonLevel:MI = on");
    pythia8->Set("PartonLevel:ISR = on");
    pythia8->Set("BeamRemnants:primordialKT = on");
    pythia8->Set("PartonLevel:FSR = on");
    pythia8->Set("StringFlav:mesonCvector = 1.5");
    pythia8->Set("StringFlav:mesonBvector = 3");
    pythia8->Set("4:m0 = 1.43");
    pythia8->Set("5:m0 = 4.30");
    //pythia8->Set("HardQCD:all = on");

    //pythia8->Set("HardQCD:hardccbar = on"); // some how does not switch on qq2ccbar and qqbar2ccbar
    pythia8->Set("HardQCD:gg2ccbar = on");
    pythia8->Set("HardQCD:qqbar2ccbar = on");
    //pythia8->Set("HardQCD:gg2bbbar = on");
    //pythia8->Set("HardQCD:qqbar2bbbar = on");
/*
    pythia8->Set("411:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("411:OnIfMatch=-321 211 211"); // switch on any decay mode into K- pi+ pi+
    pythia8->Set("-411:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("-411:OnIfMatch=321 -211 -211"); // switch on any decay mode into  K+ pi- pi-

*/

  }

  //D+- and D0 settings for PYTHIA8 new tune
  if ( config=="pp:dmesontune_new" || config=="pp:sngdmesontune_new"   )
  {
    //new PYTHIA 8 tune   
    pythia8->Set("PDF:pSet = 17");
	  pythia8->Set("MultipartonInteractions:ecmRef = 200");
	  pythia8->Set("MultipartonInteractions:bprofile = 2");
	  pythia8->Set("MultipartonInteractions:pT0Ref = 0.140");
	  pythia8->Set("MultipartonInteractions:ecmPow  = 0.135");
	  pythia8->Set("MultipartonInteractions:coreRadius = 0.56");
	  pythia8->Set("MultipartonInteractions:coreFraction = 0.78");
	  pythia8->Set("ColourReconnection:range = 5.4");

	  pythia8->Set(Form("PhaseSpace:pTHatMin = %f",ptHatMin));
	  pythia8->Set(Form("PhaseSpace:pTHatMax = %f",ptHatMax));

	  pythia8->Set("HardQCD:hardccbar = on");
	  pythia8->Set("HardQCD:gg2ccbar = on");
	  pythia8->Set("HardQCD:qqbar2ccbar = on");

    pythia8->Set("421:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("421:OnIfMatch=-321 211"); // switch on any decay mode into K- pi+ pi+
    pythia8->Set("-421:onMode=0"); // switch off all decay modes for pdgid
    pythia8->Set("-421:OnIfMatch=321 -211"); // switch on any decay mode into  K+ pi- pi-

  }

  _primary -> AddGenerator( pythia8 );
  
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=1000, Int_t collEnergy = 200, UInt_t rngSeed = 12345, TString config = "pp:dmeson" )
{ 

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "tables nodefault";
    bfc(0, simple );
  }

  cout<<"Load libraries"<<endl;

  gSystem->Load( "libVMC.so");
  gSystem->Load( "St_g2t.so" );
  gSystem->Load( "St_geant_Maker.so" );
 
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "StarGeneratorFilt.so" );

  gSystem->Load( "libMathMore.so"   );  

  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    if(config.Contains("dmesontune_new"))
    {
      _primary -> SetFileName( "pythia8.dmesontune_new.starsim.root");
    }
    else if(config.Contains("dmesontune"))
    {
      _primary -> SetFileName( "pythia8.dmesontune.starsim.root");
    }
    else if( config.Contains("dmeson"))
    {
      _primary -> SetFileName( "pythia8.dmeson.starsim.root");
    }
    else if(config.Contains("minbiasLambda"))
    {
      _primary->SetFileName("pythia8.minbiasLambda.starsim.root");
    }
    else if(config.Contains("minbias"))
    {
      _primary->SetFileName("pythia8.minbias.starsim.root");
    }
    
    //  chain -> AddBefore( "geant", primary );
  }

  //
  // Setup an event generator
  //
  //Pythia8( "pp:minbias" );
  //Pythia8( "pp:dmesontune" );
  Pythia8( config, collEnergy );

  //
  // Initialize random number generator
  //
  StarRandom &random = StarRandom::Instance();
  random.capture(); // maps all ROOT TRandoms to StarRandom
  random.seed( rngSeed );


  //
  // Setup cuts on which particles get passed to geant for
  //   simulation.  (To run generator in standalone mode,
  //   set ptmin=1.0E9.)
  //                    ptmin  ptmax
  _primary->SetPtRange  (1.0E9,  -1.0);         // GeV
  //                    etamin etamax
  _primary->SetEtaRange ( -10.0, +10.0 );
  //                    phimin phimax
  _primary->SetPhiRange ( 0., TMath::TwoPi() );
  
  
  // 
  // Fixed x, y, z vertex
  // 
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0., 0., 0. );

  //set particle filter
  if( config.Contains("dmeson"))
  {
    //sngDmesonFilter = new StDmesonFilter(); //this filter does not have any daughter or mother cuts - selects all D0 and D+- mesons

    DiD0PiKFilter = new StDiD0PiKFilter();
    DiD0PiKFilter->SetDauKine(0.15, 20., -1, 1, 0, TMath::TwoPi(), 0.15, 20., -1, 1, 0, TMath::TwoPi());
    DiD0PiKFilter->SetParentRapidities(-1, 1, -1, 1); 

    //_primary->AddFilter( sngDmesonFilter );
    _primary->AddFilter( DiD0PiKFilter );
    _primary->SetAttr("FilterKeepAll",    int(0));
    _primary->SetAttr("FilterKeepHeader", int(0));
  }
  if(config.Contains("minbiasLambda")) 
  {
    sngLambdaPiPFilter = new StLambdaPiPFilter();
    //StLambdaPiPFilter::SetDauKine(double ptMin1, double ptMax1, double etaMin1, double etaMax1, double phiMin1, double phiMax1, double ptMin2, double ptMax2, double etaMin2, double etaMax2, double phiMin2, double phiMax2)
    //sngLambdaPiPFilter->SetDauKine(0.15, 20., -1, 1, 0, 2*TMath::Pi(), 0.15, 20., -1, 1, 0, 2*TMath::Pi());
    sngLambdaPiPFilter->SetDauKine(0.15, 20., -1, 1, -TMath::Pi(), TMath::Pi(), 0.15, 20., -1, 1, -TMath::Pi(), TMath::Pi());
    sngLambdaPiPFilter->SetParentRapidities(-1, 1, -1, 1); //setter for mother rapidities (Lambda and Lambda-bar)

    _primary->AddFilter( sngLambdaPiPFilter );
    _primary->SetAttr("FilterKeepAll",    int(0));
    _primary->SetAttr("FilterKeepHeader", int(0));
    //_primary->SetAttr("Debug", int(1));
    //_primary->SetAttr("FilterSkipRejects",    int(1));
  }
 
  

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();


  //
  // By default we have configured pythia8 such that the following
  // particles are stable, so that geant is responsible for decaying
  // them.  If you need pythia to decay these for analysis, enable
  // the following code block...
  //
  if ( GEANT_HANDLES_DECAYS ) {
    _pythia8->Set("111:onMode=1"); // pi0 
    _pythia8->Set("211:onMode=1"); // pi+/-                         
    _pythia8->Set("221:onMode=1"); // eta                              
    _pythia8->Set("321:onMode=1"); // K+/-                             
    _pythia8->Set("310:onMode=1"); // K short                                               
    _pythia8->Set("130:onMode=1"); // K long                                               
    //_pythia8->Set("3122:onMode=1"); // Lambda 0
    _pythia8->Set("3122:onMode=0");
    _pythia8->Set("3122:OnIfMatch=2212 -211");
    _pythia8->Set("-3122:onMode=0");
    _pythia8->Set("-3122:OnIfMatch=-2212 211");                                          
    _pythia8->Set("3112:onMode=1"); // Sigma -                                              
    _pythia8->Set("3222:onMode=1"); // Sigma +                                              
    _pythia8->Set("3212:onMode=1"); // Sigma 0                                              
    _pythia8->Set("3312:onMode=1"); // Xi -                                                 
    _pythia8->Set("3322:onMode=1"); // Xi 0                                                 
    _pythia8->Set("3334:onMode=1"); // Omega -              
  }


  //
  // Trigger on nevents
  //
  trig( nevents );

}
// ----------------------------------------------------------------------------

