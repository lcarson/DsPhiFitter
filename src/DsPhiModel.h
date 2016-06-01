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
#include "RooUnblindPrecision.h"
#include "RooUnblindUniform.h"
#include "RooUnblindOffset.h"

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
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > GetResult();
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
  RooCategory Bmode;
  RooCategory charge;
  RooCategory magnet;
  RooCategory helBin;
  RooCategory DsBDTBin;
  RooCategory PhiBDTBin;
  std::string blindString;
  bool needsBlinding;
	
  std::vector<std::string> typeList;
  std::vector<std::string> modeList; 
  std::vector<std::string> BmodeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> magnetList;
  std::vector<std::string> HelBinList;
  std::vector<std::string> DsBDTBinList;
  std::vector<std::string> PhiBDTBinList;

  public : 
  /* 
  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsd0;
  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_comb;
  //std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsd0st;
  //std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > yield_dsstd0;  
  std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > yield_dsd0st;
  std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > yield_dsstd0;

  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > PR_total_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > frac;
  
  //Signal Shape
  std::map< std::string, std::map< std::string, std::map< std::string, RooRealVar* > > > mean_B;
  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > sigma_dsd0;
  std::map<std::string,RooRealVar*> ncb_dsd0;
  std::map<std::string,RooRealVar*> alpha_dsd0;
  
  //Combinatoric
  std::map<std::string,RooRealVar*> comb_slope_dsd0;

  //PartReco BG RooHORNSDini
  std::map<std::string,RooRealVar*> HORNS_a;       //first horn peak mass
  std::map<std::string,RooRealVar*> HORNS_b;       //second horn peak mass
  std::map<std::string,RooRealVar*> HORNS_csi;     //slope parameter, <1 if peaks of unequal height
  std::map<std::string,RooRealVar*> HORNS_shift;   //shift along B mass - rigid shift of full shape
  std::map<std::string,RooRealVar*> HORNS_sigma;   //width of core resolution Gaussian
  std::map<std::string,RooRealVar*> HORNS_R;       //ratio of widths of two Gaussians
  std::map<std::string,RooRealVar*> HORNS_f;       //between 0 and 1, fraction of total area in core Gauss
  
  //PartReco BG RooHILLSDini
  std::map<std::string,RooRealVar*> HILL_a;       
  std::map<std::string,RooRealVar*> HILL_b;      
  std::map<std::string,RooRealVar*> HILL_csi;    
  std::map<std::string,RooRealVar*> HILL_shift;   
  std::map<std::string,RooRealVar*> HILL_sigma; 
  std::map<std::string,RooRealVar*> HILL_R;       
  std::map<std::string,RooRealVar*> HILL_f; 
  */
  //Helicity Bin        DsBDTBin            PhiBDTBin            Bmode                Dmode                charge               magnet
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > yield_peak;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > yield_comb; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > yield_XXst;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > yield_XstX;

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > PR_total_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > frac;
  
  //Signal Shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map< std::string, std::map< std::string, std::map< std::string, RooRealVar* > > > > > > > mean_B;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > sigma;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > ncb;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > alpha;
  
  //Combinatoric
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > comb_slope;

  //PartReco BG RooHORNSDini
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_a;       //first horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_b;       //second horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_csi;     //slope parameter, <1 if peaks of unequal height
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_shift;   //shift along B mass - rigid shift of full shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_sigma;   //width of core resolution Gaussian
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_R;       //ratio of widths of two Gaussians
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HORNS_f;       //between 0 and 1, fraction of total area in core Gauss
  
  //PartReco BG RooHILLSDini
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_a;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_b;      
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_csi;    
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_shift;   
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_sigma; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_R;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > HILL_f; 

  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindPrecision*> > > > B_yield;
  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindOffset*> > > > B_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > B_yield;
  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindUniform*> > > > B_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > yield_double;


  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > fit_results;

};


#endif
