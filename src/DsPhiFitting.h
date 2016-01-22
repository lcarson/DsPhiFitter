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
  RooDataSet* MakeDataSet(std::string,std::string);
  std::string reducedFileName(std::string,std::string);
  RooDataSet *FinalDataSet(const std::string, const std::string, const std::string, TTree*);
  //void IdMultipleCandidates(RooDataSet*);

  void PrintDataSet();

  void RunFullFit(bool);
  //void RunManyFits();
  void RunManyToys();

  void OrderToys(int);
  void DisplayToys();

private:
  Parameters* par;
  std::string toys;
	
  DsPhiModel *model;
  DsPhiModel *model_gen;
  RooDataSet* data;
  RooWorkspace *w;
  RooArgSet inputlist;
  RooArgSet fulllist;
  RooArgSet storelist;
  RooArgSet reducedlist;
  //RooRealVar eventNumber; //It doesn't like ULong64_t
  RooRealVar runNumber;
  RooRealVar mB;
  RooRealVar mD;
  //RooRealVar bach_dll;
  RooRealVar bdt;
  RooRealVar bid;
  RooRealVar L0Hadron_TOS;
  RooRealVar L0_TIS;
  RooRealVar Hlt1TrackAllL0_TOS;
  RooRealVar Hlt2Topo2Body_TOS;
  RooRealVar Hlt2Topo3Body_TOS;
  RooRealVar Hlt2Topo4Body_TOS;

  RooCategory type;
  RooCategory mode;
  RooCategory magnet;
  RooCategory charge;
  RooCategory bMassRegion;

  float m_mD_Mu;
  double m_bdtCut;

  bool state;
  
  TNtuple *ntuple;
  float x[20];
	
  std::vector<std::string> typeList;
  std::vector<std::string> modeList;
  std::vector<std::string> chargeList;
  std::vector<std::string> magnetList;

  std::map<std::string,std::map<std::string,std::map<std::string,TPad*> > > canpad;
  std::map<std::string,std::map<std::string,std::map<std::string,std::string> > > title;
 

};

#endif
