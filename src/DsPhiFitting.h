#ifndef DsPhiFitting_h
#define DsPhiFitting_h

#include <string>

#include "TPad.h"

#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooWorkspace.h"

#include "Base.h"
#include "DsPhiModel.h"
#include "PhiModel.h"
#include "DsModel.h"

class TTree;
class TNtuple;
class TPaveText;
class Parameters;
class TApplication;


class DsPhiFitting : public Base {
  public :
  DsPhiFitting(Parameters*,TApplication*);
  ~DsPhiFitting(){}

  void DefineModel();
  void DefineRooCategories();
  void AssignGlobalCategory();

  int  LoadDataSet();
  RooDataSet* MakeDataSet(std::string,std::string,std::string);
  std::string reducedFileName(std::string,std::string,std::string);
  RooDataSet *FinalDataSet(const std::string, const std::string, const std::string, const std::string, TTree*);
  //void IdMultipleCandidates(RooDataSet*);

  void PrintDataSet();

  void RunFullFit(bool);
  void RunManyFits();
  void RunEfficiency();
  void RunManyToys();
  void SetLHCbStyle(std::string);

  void OrderToys(int);
  void DisplayToys();

private:
  Parameters* par;
  std::string toys;
	
  DsPhiModel *model;
  DsPhiModel *model_gen;
  PhiModel   *model_phi;
  DsModel    *model_ds;
  RooDataSet* data;
  RooDataSet* cutdata;
  RooWorkspace *w;
  RooArgSet inputlist;
  RooArgSet fulllist;
  RooArgSet storelist;
  RooArgSet reducedlist;
  //RooRealVar eventNumber; //It doesn't like ULong64_t
  RooRealVar runNumber;
  RooRealVar mB;
  RooRealVar mD0;
  RooRealVar mDs;
  RooRealVar mPhi;
  //RooRealVar bach_dll;
  
  //Data BDT variables 
  RooRealVar bdtDs;
  RooRealVar bdtgDs;
  RooRealVar bdtbDs;
  RooRealVar bdtD0;
  RooRealVar bdtgD0;
  RooRealVar bdtbD0;
  RooRealVar bdtPhi;
  RooRealVar bdtgPhi;
  RooRealVar bdtbPhi;
  //MC BDT variables
  RooRealVar bdtMC;
  RooRealVar bdtgMC;
  RooRealVar bdtbMC;
  // PID variables
  RooRealVar D0_K0_pidk;
  RooRealVar D0_K1_pidk;
  RooRealVar D_K_pidk;
  RooRealVar D_K0_pidk;
  RooRealVar D_K1_pidk;
  RooRealVar D_P_pidk;
  RooRealVar D_P0_pidk;
  RooRealVar D_P1_pidk;
  RooRealVar D_P2_pidk;

  RooRealVar helicityD0;
  RooRealVar helicityPhi;
  RooRealVar mKK;
  RooRealVar bid;
  RooRealVar L0Hadron_TOS;
  RooRealVar L0_TIS;
  RooRealVar Hlt2IncPhi_TOS;
  RooRealVar Hlt2PhiIncPhi_TOS;
  // 2011 2012 Trigger Lines
  RooRealVar Hlt1TrackAllL0_TOS;
  RooRealVar Hlt2Topo2BodyBBDT_TOS;
  RooRealVar Hlt2Topo3BodyBBDT_TOS;
  RooRealVar Hlt2Topo4BodyBBDT_TOS;
  // 2015 Trigger Lines
  RooRealVar Hlt1TrackMVA_TOS;
  RooRealVar Hlt1TwoTrackMVA_TOS;
  RooRealVar Hlt2Topo2Body_TOS;
  RooRealVar Hlt2Topo3Body_TOS;
  RooRealVar Hlt2Topo4Body_TOS;
  
  // Veto Variables
  /// KKPi
  RooRealVar KKPi_D_Veto;
  RooRealVar KKPi_Lc_Veto;

  RooRealVar deltaMass;
  RooRealVar deltaMass2;
  /// KPiPi
  RooRealVar KPiPi_mPiKPi;
  RooRealVar KPiPi_mPiPiPi;
  RooRealVar KPiPi_mKKPi;

  //FDChi2
  RooRealVar D_FDCHI2;
  RooRealVar D0_FDCHI2;

  RooCategory type;
  RooCategory mode;
  RooCategory Bmode;
  RooCategory magnet;
  RooCategory charge;
  RooCategory bMassRegion; 
  RooCategory helBin;
  RooCategory DsBDTBin;
  RooCategory PhiBDTBin;

  float m_mD0_Mu;
  float m_mDs_Mu;
  float m_mPhi_Mu;
  double m_bdtCut;

  bool state;
  
  TNtuple *ntuple;
  float x[20];
	
  std::vector<std::string> typeList;
  std::vector<std::string> yearList;
  std::vector<std::string> modeList;
  std::vector<std::string> BmodeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> magnetList; 
  std::vector<std::string> HelBinList;
  std::vector<std::string> DsBDTBinList;
  std::vector<std::string> PhiBDTBinList;

 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,TPad*> > > > > > > canpad;
 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::string> > > > > > > > title;
 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::string> > > > > > > > bin_detail;
 
 
 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > > fit_results;
 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > fit_results_phi;
 std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > fit_results_ds;

 std::map<std::string,double> Initial_value;
};

#endif
