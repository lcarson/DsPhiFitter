#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <math.h>

#include "TFile.h"
#include "TMath.h"

#include "RooRandom.h"
#include "RooAddPdf.h"
#include "RooRealConstant.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooThresholdCategory.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "RooArgusBG.h"
#include "RooWorkspace.h"
#include "RooUnblindUniform.h"
#include "RooUnblindPrecision.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooKeysPdf.h"
#include "RooEffProd.h"
//#include "RooLowMassPdf.h"
#include "RooMultiVarGaussian.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "PhiModel.h"
#include "Parameters.h"
//#include "RooDoubleCB.h"
#include "RooHORNSdini.h"
#include "RooHILLdini.h"
#include "TMatrixTSym.h"




PhiModel::PhiModel(Parameters* p, RooRealVar* pmPhi, RooRealVar* pmD0, bool nB)
: Base()
  //, misID(p)
  //, calibration(p)
  //, pid("pid","PID boolean")
, type("type","Type of data") 
, mode("mode","D_{s} decay mode")
, charge("charge","batchelor charge")
, magnet("magnet","magnet polarity")
, helBin("helBin","Helicity Angle bin")
, DsBDTBin("DsBDTBin","Ds BDT bin")
, PhiBDTBin("PhiBDTBin","Phi BDT bin")
, blindString("BlindLeadingTheBlind")
, needsBlinding(nB)
  //, pidList()
, typeList()
, modeList()
, BmodeList()
, chargeList()
, magnetList()
, HelBinList()
, DsBDTBinList()
, PhiBDTBinList()
  //, m_pidCut(-1000)
{
  par=p;
  mPhi=pmPhi;
  mD0=pmD0;
  particle_name[DsPhi] = Phi;
  particle_name[DsD0]  = D0;
  rand = new TRandom3(time(NULL));
  if(needsBlinding){
    RooRandom::randomGenerator()->SetSeed(20110602);
    std::cout<<"BLIND model envoked"<<std::endl;
  }
  //  m_pidCut=int(par->pidCut);
  DefineRooCategories();
  DefineModel();
}


void PhiModel::DefineRooCategories() //like D0h with "PIDcut" removed
{
  if(par->debug) std::cout<<"Running: PhiModel::DefineRooCategories()"<<std::endl;
  if(par->modes[Ds2PhiPi])     modeList.push_back(Ds2PhiPi);
  if(par->modes[Ds2KKPi])      modeList.push_back(Ds2KKPi);
  if(par->modes[Ds2PiPiPi])    modeList.push_back(Ds2PiPiPi);
  if(par->modes[Ds2KPiPi])     modeList.push_back(Ds2KPiPi);
  
  if(par->Bmodes[DsD0])        BmodeList.push_back(DsD0);
  if(par->Bmodes[DsPhi])       BmodeList.push_back(DsPhi);
 

  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){       
    mode.defineType( (*m).c_str()); 
  }

  for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){       
    Bmode.defineType( (*b).c_str()); 
  }
  

  if(par->sumOverCharges){
    chargeList.push_back(both);
  }else{
    chargeList.push_back(minus);
    chargeList.push_back(plus);
  }

  if(par->splitHel) {
    HelBinList.push_back(Helbin1);
    HelBinList.push_back(Helbin2);
  } else {
    HelBinList.push_back(both);
  }
  
  helBin.defineType(both.c_str());
  helBin.defineType(Helbin1.c_str());
  helBin.defineType(Helbin2.c_str());

  DsBDTBinList.push_back(DsBDTbin1);
  if(par->nDsBDTBins>1) DsBDTBinList.push_back(DsBDTbin2);

  PhiBDTBinList.push_back(PhiBDTbin1);
  if(par->nPhiBDTBins>1) PhiBDTBinList.push_back(PhiBDTbin2);
  
  DsBDTBin.defineType(DsBDTbin1.c_str());
  DsBDTBin.defineType(DsBDTbin2.c_str());

  PhiBDTBin.defineType(PhiBDTbin1.c_str());
  PhiBDTBin.defineType(PhiBDTbin2.c_str());

  charge.defineType(both.c_str());
  charge.defineType(plus.c_str());
  charge.defineType(minus.c_str());
  
  if(par->polarity==sensitive){
    magnetList.push_back(up);
    magnetList.push_back(dn);
  }else{
    magnetList.push_back(par->polarity);
  }  
  magnet.defineType(both.c_str());
  magnet.defineType(up.c_str());
  magnet.defineType(dn.c_str());

  cat = new RooCategory("cat","helBin/DsBDTBin/PhiBDTBin/Bmode/mode/charge/polarity");
  
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++)  {
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                std::stringstream str;
                str<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                cat->defineType(str.str().c_str());
                std::cout <<  " --> Creating Category : " <<str.str().c_str() <<std::endl;
              }  
            }
          }
        }
      }
    }
  }

}


void PhiModel::DefineModel()
{
  if(par->debug) std::cout<<"Running: PhiModel::DefineModel()"<<std::endl;

  std::map<std::string,std::string> mod;
  std::map<std::string,std::string> Bmod;
  mod[Ds2PhiPi]="#phi#pi"; mod[Ds2KKPi]="KK#pi"; mod[Ds2KPiPi]="K#pi#pi";    mod[Ds2KKPi]="#pi#pi#pi";
  Bmod[DsD0]="D_{s}D^{0}"; Bmod[DsPhi]="D_{s}#phi";

  // --------- Yields for DsD0, in an array (as line 447 of Model.C) ---------------
  if(par->debug) std::cout<<"Running: PhiModel::DefineModel() --> Making RooRealVars"<<std::endl;
  
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
        for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_%s_peak_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Peak D_{s}D^{0} (#phi#pi)",  1660, 0, 830000);
            yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_%s_peak_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield #phi Peak D_{s}D^{0} (KK#pi)",    1660, 0, 830000);
            yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_%s_peak_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield #phi Peak D_{s}D^{0} (#pi#pi#pi)", 350, 0, 175000);
            yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_%s_peak_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Peak D_{s}D^{0} (K#pi#pi)",   200, 0, 100000);

            yield_Phi_comb[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_%s_comb_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Comb. D_{s}D^{0} (#phi#pi)",   710, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_%s_comb_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield #phi Comb. D_{s}D^{0} (KK#pi)",     710, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_%s_comb_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield #phi Comb. D_{s}D^{0} (#pi#pi#pi)", 390, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_%s_comb_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",particle_name[DsD0].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Comb. D_{s}D^{0} (K#pi#pi)",   200, 0, 600000);

            // ---------------
            // DsPhi variables
            // ---------------
            yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_%s_peak_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Peak D_{s}#phi (#phi#pi)",  1660, 0, 83000);
            yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_%s_peak_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield #phi Peak D_{s}#phi (KK#pi)",    1660, 0, 83000);
            yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_%s_peak_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield #phi Peak D_{s}#phi (#pi#pi#pi)", 350, 0, 17500);
            yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_%s_peak_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Peak D_{s}#phi (K#pi#pi)",   200, 0, 10000);

            yield_Phi_comb[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_%s_comb_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Comb. D_{s}#phi (#phi#pi)",   710, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_%s_comb_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield #phi Comb. D_{s}#phi (KK#pi)",     710, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_%s_comb_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield #phi Comb. D_{s}#phi (#pi#pi#pi)", 390, 0, 600000);
            yield_Phi_comb[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_%s_comb_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield #phi Comb. D_{s}#phi (K#pi#pi)",   200, 0, 600000);


            // -------------------------------------------------------------------------

            if(par->debug) std::cout<<"Running: PhiModel::DefineModel() --> Setting up blinding"<<std::endl;
            // RooCategory used for blind/unblind switching.
            TString blind("blind"), unblind("unblind");
            RooCategory blindCatBu("blindCatBu","Bu blind state Category");
            blindCatBu.defineType(unblind, 0);
            blindCatBu.defineType(blind, 1);

            if (needsBlinding)
                blindCatBu.setLabel(blind);
            else
                blindCatBu.setLabel(unblind);

            
            if(needsBlinding){  
              Phi_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];

              //Phi_yield[DsD0][Ds2PhiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2PhiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "nsig Bu blind (#phi#pi)",  "nsigBuDsD0blindedPhiPi",  1000., *yield_Phi_peak[DsD0][Ds2PhiPi][*c][*a] );
              //Phi_yield[DsD0][Ds2KKPi][*c][*a]   = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "nsig Bu blind (KK#pi)",    "nsigBuDsD0blindedKKPi",   1000., *yield_Phi_peak[DsD0][Ds2KKPi][*c][*a]  );
              //Phi_yield[DsD0][Ds2PiPiPi][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "nsig Bu blind (#pi#pi#pi)","nsigBuDsD0blindedPiPiPi", 100., *yield_Phi_peak[DsD0][Ds2PiPiPi][*c][*a] );
              //Phi_yield[DsD0][Ds2KPiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "nsig Bu blind (K#pi#pi)",  "nsigBuDsD0blindedKPiPi",  100., *yield_Phi_peak[DsD0][Ds2KPiPi][*c][*a]  );
              
              // Only Blind DsPhi numbers...
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooUnblindUniform(Form("%s_nsig_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s", particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (#phi#pi)",  "nsigBuDsPhiblindePhiPi",   40., *yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a] );
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooUnblindUniform(Form("%s_nsig_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",  particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (KK#pi)",    "nsigBuDsPhiblindedKKPi",   40., *yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]  );
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooUnblindUniform(Form("%s_nsig_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (#pi#pi#pi)","nsigBuDsPhiblindedPiPiPi", 20., *yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]);
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooUnblindUniform(Form("%s_nsig_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s", particle_name[DsPhi].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (K#pi#pi)",  "nsigBuDsPhiblindedKPiPi",  10., *yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a] );
            } else {
              Phi_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = yield_Phi_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];

              Phi_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a];
              Phi_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = yield_Phi_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a];
            }
          } //end of loop over chargeList
        }  //end of loop over magnetList
      } 
    } 
  }



  // --------- Mean B mass ---------------

  double default_Phi=1019.461;
  double Phi_low    = 994.461;
  double Phi_high   = 1044.461;
  
  double default_D0 =1864.84;
  double D0_low     = 1839.84;
  double D0_high    = 1889.84;
  
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            bool isPhi = false;
            if((*b).c_str()==DsPhi) isPhi = true; 
            mean_Phi[*h][*ds][*ph][*b][*m][both][both] = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_%s"           ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][plus][up]   = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_plus_up_%s"   ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][minus][dn]  = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_minus_dn_%s"  ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][plus][dn]   = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_plus_dn_%s"   ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][minus][up]  = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_minus_up_%s"  ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][plus][both] = new RooRealVar(Form("mean_%s_%s_%s_%s_%s_plus_both_%s" ,particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][minus][both]= new RooRealVar(Form("mean_%s_%s_%s_%s_%s_minus_both_%s",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean Phi mass", (isPhi?default_Phi:default_D0), (isPhi?Phi_low:D0_low), (isPhi?Phi_high:D0_high));
            mean_Phi[*h][*ds][*ph][*b][*m][both][up]		= mean_Phi[*h][*ds][*ph][*b][*m][both][both];
            mean_Phi[*h][*ds][*ph][*b][*m][both][dn]		= mean_Phi[*h][*ds][*ph][*b][*m][both][both];
            mean_Phi[*h][*ds][*ph][*b][Ds2KPiPi][both][both]	= mean_Phi[*h][*ds][*ph][*b][Ds2KKPi][both][both];
            mean_Phi[*h][*ds][*ph][*b][Ds2PiPiPi][both][both] = mean_Phi[*h][*ds][*ph][*b][Ds2KKPi][both][both];
            mean_Phi[*h][*ds][*ph][*b][Ds2PhiPi][both][both]  = mean_Phi[*h][*ds][*ph][*b][Ds2KKPi][both][both];
          }
        }
      }
    }
  }

  if(par->debug) std::cout<<"Running: PhiModel::DefineModel() --> Making more RooRealVars"<<std::endl;
  // --------- Signal width ---------------
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            sigma_Phi[*h][*ds][*ph][*b][*m][both][both]  = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_%s",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())          ,Form("%s sigma",mod[*b].c_str()),10,  1,  20); //change D0 to D^0
            sigma_Phi[*h][*ds][*ph][*b][*m][plus][both]  = sigma_Phi[*h][*ds][*ph][*b][*m][both][both];
            sigma_Phi[*h][*ds][*ph][*b][*m][minus][both] = sigma_Phi[*h][*ds][*ph][*b][*m][both][both];
            sigma_Phi[*h][*ds][*ph][*b][*m][plus][up]    = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_%s_plus_up",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())  ,Form("%s sigma",mod[*b].c_str()),10,  1,  20);
            sigma_Phi[*h][*ds][*ph][*b][*m][plus][dn]    = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_%s_plus_dn",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())  ,Form("%s sigma",mod[*b].c_str()),10,  1,  20);
            sigma_Phi[*h][*ds][*ph][*b][*m][minus][up]   = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_%s_minus_up",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()) ,Form("%s sigma",mod[*b].c_str()),10,  1,  20);//sigma_Phi[*m][minus][up] = sigma_Phi[*m][plus][dn];
            sigma_Phi[*h][*ds][*ph][*b][*m][minus][dn]   = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_%s_minus_dn",particle_name[*b].c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()) ,Form("%s sigma",mod[*b].c_str()),10,  1,  20);//sigma_Phi[*m][minus][dn] = sigma_Phi[*m][plus][up];
            sigma_Phi[*h][*ds][*ph][*b][*m][both][up]    = sigma_Phi[*h][*ds][*ph][*b][*m][both][both]; 
            sigma_Phi[*h][*ds][*ph][*b][*m][both][dn]    = sigma_Phi[*h][*ds][*ph][*b][*m][both][both];            
            sigma_Phi[*h][*ds][*ph][DsPhi][*m][both][dn] = sigma_Phi[*h][*ds][*ph][DsD0][*m][both][both];
            sigma_Phi[*h][*ds][*ph][DsPhi][*m][both][up] = sigma_Phi[*h][*ds][*ph][DsD0][*m][both][both];
          }
        }
      }
    }
  }

  // --------- Comb. background slope ---------------
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            comb_a1[*h][*ds][*ph][*b][*m] = new RooRealVar(Form("comb_a1_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())   ,Form("%s comb. a1",mod[*b].c_str()),  0.4,   -1.0, 1.0);
            comb_a2[*h][*ds][*ph][*b][*m] = new RooRealVar(Form("comb_a2_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())   ,Form("%s comb. a2",mod[*b].c_str()), -0.004, -1.0, 1.0);
          }
        }
      }
    }
  }
  
  // --------- Fix RooRealVars to each other  -------------
  // ------------------------------------------------------
  // ------------------------------------------------------
   
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){           
 
            // Fix comb slope to be the same 
            comb_a1[*h][*ds][*ph][*b][*m] = comb_a1[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];
            comb_a2[*h][*ds][*ph][*b][*m] = comb_a2[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];

            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                
                // Use same mean value for all plots 
                mean_Phi[*h][*ds][*ph][*b][*m][*a][*c] = mean_Phi[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][both][both]; 

                // Allow sigma_Phi to vary between Ds decay modes
                sigma_Phi[*h][*ds][*ph][*b][*m][*a][*c] = sigma_Phi[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*m][*magnetList.begin()][*chargeList.begin()]; 

              }
            }
          }
        }
      }
    }
  }



  // --------- Define the simultaneous PDFs ---------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  
  if(par->debug) std::cout<<"Running: PhiModel::DefineModel() --> Setting up simultaneous PDF"<<std::endl;
  sim = new RooSimultaneous("model","Simultaneous model",*cat);

  RooAbsPdf *pdf_peak = 0;
  RooAbsPdf *pdf_comb = 0;
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                std::string tag=(*h)+underscore+(*ds)+underscore+(*ph)+underscore+(*b)+underscore+(*m)+underscore+(*c)+underscore+(*a);
                if(par->debug) std::cout<<"Running: PhiModel::DefineModel() --> Adding to sim pdf: "<< tag <<std::endl;
                
                //DsD0
                RooAbsReal* mub     = mean_Phi[*h][*ds][*ph][*b][*m][*c][*a];
                RooAbsReal* sig     = sigma_Phi[*h][*ds][*ph][*b][*m][*c][*a];
                pdf_peak  = new RooGaussian( Form("pdf_Phi_peak_%s",tag.c_str()), "", (((*b).c_str()==DsPhi)?*mPhi:*mD0), *mub, *sig );

                //comb bkg
                RooRealVar* comba1 = comb_a1[*h][*ds][*ph][*b][*m];
                RooRealVar* comba2 = comb_a2[*h][*ds][*ph][*b][*m];
                pdf_comb = new RooChebychev(Form("pdf_Phi_comb_%s",tag.c_str()),"",(((*b).c_str()==DsPhi)?*mPhi:*mD0),RooArgSet(*comba1,*comba2)) ;
              
                //List of PDFs and yields
                RooArgSet pdflist( *pdf_peak, *pdf_comb);
                RooArgSet nevents( *Phi_yield[*h][*ds][*ph][*b][*m][*c][*a], *yield_Phi_comb[*h][*ds][*ph][*b][*m][*c][*a]); 
                
                //Add to master PDF
                RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_Phi_%s",tag.c_str())       ,"",pdflist, nevents);
        	      
                std::stringstream str;
                str<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                
                sim->addPdf(*masterPdf,str.str().c_str());

              }
            }
          } 
        }//closing of for loops over b, m, a and c
      }
    }
  } // closing loop over h
} //end of funcn DefineModel()



RooArgSet* PhiModel::GetParameters()
{
  if(par->debug) std::cout<<"Running: PhiModel::GetParameters()"<<std::endl;
  RooArgSet* pars = sim->getParameters(RooArgSet(*mPhi,helBin,DsBDTBin,PhiBDTBin,Bmode,mode,charge,magnet));
  return pars;
}



void PhiModel::PrintResult()
{
  if(par->debug) std::cout<<"Running: PhiModel::PrintResults()"<<std::endl;
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              if(needsBlinding && (*b)==DsPhi){
                std::cout << "DsPhi Mode is BLIND" << std::endl;
              } else {
                if(par->sumOverCharges){
                  float y_peak =   yield_Phi_peak[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =   yield_Phi_comb[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  std::cout<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                  std::cout<< " Peak:   " <<y_peak<<std::endl;
                  std::cout<< " Comb:   "<<y_comb<<std::endl;
                  std::cout<< "---------------------------"<<std::endl;
                  std::cout<< "---------------------------"<<std::endl;
                }else{
                  float y_peak_minus =   yield_Phi_peak[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_comb_minus =   yield_Phi_comb[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_peak_plus =    yield_Phi_peak[*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                  float y_comb_plus =    yield_Phi_comb[*h][*ds][*ph][*b][*m][plus][*a]->getVal();

               
                  std::cout<<(*h)<<", "<<(*ds)<<", "<<(*ph)<<", "<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                  std::cout<<" Peak minus: "<<y_peak_minus<<" comb. minus: "<<y_comb_minus<<std::endl;
                  std::cout<<" Peak plus:  "<<y_peak_plus <<" comb. plus:  "<<y_comb_plus <<std::endl;
                  std::cout<<"--------------------------------------------------"<<std::endl;
                  std::cout<<" Peak total: "<<y_peak_minus+y_peak_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<std::endl;

                  std::cout<<" Peak minus: "<<y_peak_minus<<" comb. minus: "<<y_comb_minus<<std::endl;
                  std::cout<<" Peak plus:  "<<y_peak_plus<< " comb. plus:  "<<y_comb_plus<<std::endl;
                  std::cout<<"--------------------------------------------------"<<std::endl;
                  std::cout<<" Peak total: "<<y_peak_minus+y_peak_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
 //closes if(not needsBlinding)
  std::cout<<std::endl;
}

std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > PhiModel::GetResult()
{
  if(par->debug) std::cout<<"Running: PhiModel::GetResults()"<<std::endl;
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                if (needsBlinding && (*b)==DsPhi) {
                  std::cout << "DsPhi mode is BLIND, therefore no values returned " << std::endl;
                } else {


                  std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());
                  std::cout<< " Category: " << cat_str << std::endl;
                  fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Peak"]    = yield_Phi_peak[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  //fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Peakerr"] = yield_Phi_peak[*h][*ds][*ph][*b][*m][*c][*a]->getError();
                  fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Comb"]    = yield_Phi_comb[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  std::cout<< " Peak:   " <<fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Peak"]<<std::endl;
                  std::cout<< " Comb:   " <<fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Comb"]<<std::endl;
                  
                  //mPhi->setRange("signal", 5233.5, 5333.5) ;
                  //double combslope = comb_slope[*h][*ds][*ph][*b][*m]->getVal();
                  //RooRealVar* fixed_slope = new RooRealVar("fixed_slope","slope",  combslope);
                  //RooAbsPdf* pdf_comb = new RooExponential( "Exp", "",*mPhi,*fixed_slope);
                  //RooAbsReal* ibkg_sig = pdf_comb->createIntegral(*mPhi,RooFit::NormSet(*mPhi),RooFit::Range("signal"));  
                  //fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"] = (ibkg_sig->getVal()*yield_Phi_comb[*h][*ds][*ph][*b][*m][*c][*a]->getVal());
                  

                  //std::cout<< " Comb_Reduced:   " <<fit_results_Phi[*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"]<<std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
 //closes if(not needsBlinding)

return fit_results_Phi;
}
