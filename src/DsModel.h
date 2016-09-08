#ifndef DsModel_h
#define DsModel_h

#include <string>

#include "TRandom3.h"

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooSuperCategory.h"
#include "RooMultiVarGaussian.h"
#include "RooUnblindPrecision.h"
#include "RooUnblindUniform.h"
#include "RooUnblindOffset.h"

#include "Base.h"
#include "Parameters.h"
//#include "MisIDShapes.h"
//#include "PIDCalibration.h"
//#include "PartRecoShapes.h"

class RooGaussian;

class DsModel : public Base {
public :
  DsModel(Parameters*,RooRealVar*,bool);
  ~DsModel(){}
	void DefineRooCategories();
	void DefineModel();
  RooCategory* Cat(){return cat;}
	RooSimultaneous* Pdf(){return sim;}
	void PrintResult();
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > GetResult();
  RooArgSet* GetParameters();
  const RooArgSet& Constraints();

private:
  TRandom3 *rand;
  Parameters* par;
	RooCategory* cat;
	RooSimultaneous* sim;


  RooRealVar* mDs;
  RooCategory type;
  RooCategory mode;
  RooCategory Bmode;
  RooCategory charge;
  RooCategory magnet;
  RooCategory helBin;
  RooCategory DsBDTBin;
  RooCategory PhiBDTBin;
  std::string blindString;
  bool needsBlinding;
  bool isPhi;
  std::string particle_name; 
  std::vector<std::string> typeList;
  std::vector<std::string> modeList; 
  std::vector<std::string> BmodeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> magnetList;
  std::vector<std::string> HelBinList;
  std::vector<std::string> DsBDTBinList;
  std::vector<std::string> PhiBDTBinList;

  public : 

  //Helicity Bin        DsBDTBin            PhiBDTBin            Bmode                Dmode                charge               magnet
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > yield_Ds_peak;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > yield_Ds_comb; 
  
  //Signal Shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map< std::string, std::map< std::string, std::map< std::string, RooRealVar* > > > > > > > mean_Ds;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > sigma_Ds;

  
  //Combinatoric
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > comb_a1;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > comb_a2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > comb_slope;

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > Ds_yield;

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > fit_results_Ds;

};


#endif
