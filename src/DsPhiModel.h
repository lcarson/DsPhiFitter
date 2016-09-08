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
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > > GetResult();
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
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_peak;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_CB1;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_CB2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_Dsa1;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_DsstKKst;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > >     yield_comb; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_dsstd0;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >     yield_dsdst0;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  yield_HORNS;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  yield_HORNS2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  yield_HILL;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  yield_HILL2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  yield_LITTLEHORNS;

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  Signal_total;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  Low_Mass_total;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  PR_total_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  frac;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  frac_HORNS;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  frac_HILL;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  frac_HILL2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  frac_LITTLEHORNS;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > >  frac_pi0;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > > >  frac_hel;
  
  //Signal Shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map< std::string, std::map< std::string, std::map< std::string, RooRealVar* > > > > > > > > mean_B;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  sigma;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooFormulaVar*> > > > > > > >  sigma2;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > >  avg_sigma;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > sig_frac;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > ncb;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > alpha;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > > sig_ratio;
  
  //Combinatoric
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > > > >  comb_slope;

  //PartReco BG RooHORNSDini
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_a;       //first horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_b;       //second horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_csi;     //slope parameter, <1 if peaks of unequal height
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_shift;   //shift along B mass - rigid shift of full shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_sigma;   //width of core resolution Gaussian
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_R;       //ratio of widths of two Gaussians
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS_f;       //between 0 and 1, fraction of total area in core Gauss
  //PartReco BG RooLITTLEHORNSDini
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_a;       //first horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_b;       //second horn peak mass
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_csi;     //slope parameter, <1 if peaks of unequal height
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_shift;   //shift along B mass - rigid shift of full shape
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_sigma;   //width of core resolution Gaussian
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_R;       //ratio of widths of two Gaussians
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_f;       //between 0 and 1, fraction of total area in core Gauss
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > LITTLEHORNS_g;       // no idea
  
  //PartReco BG RooHILLSDini
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_a;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_b;      
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_csi;    
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_shift;   
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_sigma; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_R;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL_f; 

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_a;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_b;      
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_csi;    
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_shift;   
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_sigma; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_R;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HILL2_f;

  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_a;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_b;      
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_csi;    
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_shift;   
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_sigma; 
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_R;       
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > HORNS2_f; 

  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindPrecision*> > > > B_yield;
  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindOffset*> > > > B_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooAbsReal*> > > > > > > > B_yield;
  //std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooUnblindUniform*> > > > B_yield;
  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > yield_double;


  std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > > fit_results;

};


#endif
