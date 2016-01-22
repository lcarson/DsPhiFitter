#ifndef DsPhiModel_h
#define DsPhiModel_h

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

#include "Base.h"
#include "Parameters.h"
//#include "MisIDShapes.h"
//#include "PIDCalibration.h"
//#include "PartRecoShapes.h"

class RooGaussian;

class DsPhiModel : public Base {
public :
  DsPhiModel(Parameters*,RooRealVar*,bool);
  ~DsPhiModel(){}
	void DefineRooCategories();
	void DefineModel();
  RooCategory* Cat(){return cat;}
	RooSimultaneous* Pdf(){return sim;}
	void PrintResult();
  RooArgSet* GetParameters();
  const RooArgSet& Constraints();

private:
  TRandom3 *rand;
  Parameters* par;
	RooCategory* cat;
	RooSimultaneous* sim;


  RooRealVar* mB;
  RooCategory type;
  RooCategory mode;
  RooCategory charge;
  RooCategory magnet;
  	std::string blindString;
  	bool needsBlinding;
	
  std::vector<std::string> typeList;
	std::vector<std::string> modeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> magnetList;

  public :  
	std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsd0;
       	std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_comb;
	//std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsd0st;
	//std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsstd0;
};


#endif
