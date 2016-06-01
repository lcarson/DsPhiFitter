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
#include "RooGenericPdf.h"
#include "RooKeysPdf.h"
#include "RooEffProd.h"
//#include "RooLowMassPdf.h"
#include "RooMultiVarGaussian.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "DsPhiModel.h"
#include "Parameters.h"
//#include "RooDoubleCB.h"
#include "RooHORNSdini.h"
#include "RooHILLdini.h"
#include "TMatrixTSym.h"


//**********************
//***----BLINDING----***
//**********************
Bool_t isBlind(kTRUE);


DsPhiModel::DsPhiModel(Parameters* p, RooRealVar* pmB, bool nB)
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
  mB=pmB;
  rand = new TRandom3(time(NULL));
  if(needsBlinding){
    RooRandom::randomGenerator()->SetSeed(20110602);
    std::cout<<"BLIND model envoked"<<std::endl;
  }
  //  m_pidCut=int(par->pidCut);
  DefineRooCategories();
  DefineModel();
}


void DsPhiModel::DefineRooCategories() //like D0h with "PIDcut" removed
{
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineRooCategories()"<<std::endl;
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


void DsPhiModel::DefineModel()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel()"<<std::endl;

  std::map<std::string,std::string> mod;
  std::map<std::string,std::string> Bmod;
  mod[Ds2PhiPi]="#phi#pi"; mod[Ds2KKPi]="KK#pi"; mod[Ds2KPiPi]="K#pi#pi";    mod[Ds2KKPi]="#pi#pi#pi";
  Bmod[DsD0]="D_{s}D^{0}"; Bmod[DsPhi]="D_{s}#phi";

  // --------- Yields for DsD0, in an array (as line 447 of Model.C) ---------------
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Making RooRealVars"<<std::endl;
  
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
        for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            
            yield_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_peak_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Peak D_{s}D^{0} (#phi#pi)",  1660, 0, 830000);
            yield_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_peak_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield Peak D_{s}D^{0} (KK#pi)",    1660, 0, 830000);
            yield_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_peak_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield Peak D_{s}D^{0} (#pi#pi#pi)", 350, 0, 175000);
            yield_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_peak_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Peak D_{s}D^{0} (K#pi#pi)",   200, 0, 100000);

            yield_comb[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_comb_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Comb. D_{s}D^{0} (#phi#pi)",   710, 0, 600000);
            yield_comb[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_comb_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield Comb. D_{s}D^{0} (KK#pi)",     710, 0, 600000);
            yield_comb[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_comb_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield Comb. D_{s}D^{0} (#pi#pi#pi)", 390, 0, 600000);
            yield_comb[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_comb_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Comb. D_{s}D^{0} (K#pi#pi)",   200, 0, 600000);

            frac[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  =   new RooRealVar(Form("frac_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),    "HORNS/total PR ratio D_{s}D^{0} (#phi#pi)",   0.5,  0, 1);
            frac[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   =   new RooRealVar(Form("frac_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),     "HORNS/total PR ratio D_{s}D^{0} (KK#pi)",     0.5,  0, 1);
            frac[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  =   new RooRealVar(Form("frac_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),    "HORNS/total PR ratio D_{s}D^{0} (K#pi#pi)",   0.5,  0, 1);
            frac[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] =   new RooRealVar(Form("frac_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "HORNS/total PR ratio D_{s}D^{0} (#pi#pi#pi)", 0.5,  0, 1);
            
            // Tom: Don't understand what these lines do...
            // Do these fractions need to be fixed to one another?
            //frac[Ds2KPiPi][*c][*a] = frac[Ds2KKPi][*c][*a];
            //frac[Ds2PiPiPi][*c][*a] = frac[Ds2KKPi][*c][*a];
            //frac[Ds2KKPi][plus][*a] = frac[Ds2KKPi][both][*a];
            //frac[Ds2KKPi][minus][*a] = frac[Ds2KKPi][both][*a];
            //frac[Ds2KKPi][both][up] = frac[Ds2KKPi][both][both];
            //frac[Ds2KKPi][both][dn] = frac[Ds2KKPi][both][both];

            
            PR_total_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_PR_total_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Total PR yield D_{s}D^{0} (#phi#pi)",  380, 0, 1000000);
            PR_total_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_PR_total_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Total PR yield D_{s}D^{0} (KK#pi)",    380, 0, 1000000);
            PR_total_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_PR_total_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Total PR yield D_{s}D^{0} (#pi#pi#pi)", 50, 0, 2000000 );
            PR_total_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_PR_total_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Total PR yield D_{s}D^{0} (K#pi#pi)",   30, 0, 1000000 );

            //HORNS:
            yield_XXst[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooFormulaVar(Form("yield_XXst_DsD0_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}D^{0}^{*} (#phi#pi)", "@0*@1",   RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a])  );
            yield_XXst[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooFormulaVar(Form("yield_XXst_DsD0_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield D_{s}D^{0}^{*} (KK#pi)", "@0*@1",     RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a],   *PR_total_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a])   );
            yield_XXst[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooFormulaVar(Form("yield_XXst_DsD0_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield D_{s}D^{0}^{*} (#pi#pi#pi)", "@0*@1", RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a], *PR_total_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a]) );
            yield_XXst[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooFormulaVar(Form("yield_XXst_DsD0_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}D^{0}^{*} (K#pi#pi)", "@0*@1",   RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a])  );
            //HILL:
            yield_XstX[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsD0_XstX_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}^{*}D^{0} (#phi#pi)", "(1-@0)*@1",   RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a])  );
            yield_XstX[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = new RooFormulaVar(Form("yield_DsD0_XstX_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield D_{s}^{*}D^{0} (KK#pi)", "(1-@0)*@1",     RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a],   *PR_total_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a])   );
            yield_XstX[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = new RooFormulaVar(Form("yield_DsD0_XstX_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield D_{s}^{*}D^{0} (#pi#pi#pi)", "(1-@0)*@1", RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a], *PR_total_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a]) );
            yield_XstX[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsD0_XstX_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}^{*}D^{0} (K#pi#pi)", "(1-@0)*@1",   RooArgList(*frac[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a])  );

            // ---------------
            // DsPhi variables
            // ---------------

            yield_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_peak_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Peak D_{s}#phi (#phi#pi)",  1660, 0, 83000);
            yield_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_peak_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield Peak D_{s}#phi (KK#pi)",    1660, 0, 83000);
            yield_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_peak_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield Peak D_{s}#phi (#pi#pi#pi)", 350, 0, 17500);
            yield_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_peak_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Peak D_{s}#phi (K#pi#pi)",   200, 0, 10000);

            yield_comb[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_comb_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Comb. D_{s}#phi (#phi#pi)",   710, 0, 600000);
            yield_comb[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_comb_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield Comb. D_{s}#phi (KK#pi)",     710, 0, 600000);
            yield_comb[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_comb_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield Comb. D_{s}#phi (#pi#pi#pi)", 390, 0, 600000);
            yield_comb[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_comb_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield Comb. D_{s}#phi (K#pi#pi)",   200, 0, 600000);

            frac[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  =   new RooRealVar(Form("frac_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),    "HORNS/total PR ratio D_{s}#phi (#phi#pi)",   0.5,  0, 1);
            frac[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   =   new RooRealVar(Form("frac_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),     "HORNS/total PR ratio D_{s}#phi (KK#pi)",     0.5,  0, 1);
            frac[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  =   new RooRealVar(Form("frac_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),    "HORNS/total PR ratio D_{s}#phi (K#pi#pi)",   0.5,  0, 1);
            frac[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] =   new RooRealVar(Form("frac_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "HORNS/total PR ratio D_{s}#phi (#pi#pi#pi)", 0.5,  0, 1);
            
            // Tom: Don't understand what these lines do...
            //frac[Ds2KPiPi][*c][*a] = frac[Ds2KKPi][*c][*a];
            //frac[Ds2PiPiPi][*c][*a] = frac[Ds2KKPi][*c][*a];
            //frac[Ds2KKPi][plus][*a] = frac[Ds2KKPi][both][*a];
            //frac[Ds2KKPi][minus][*a] = frac[Ds2KKPi][both][*a];
            //frac[Ds2KKPi][both][up] = frac[Ds2KKPi][both][both];
            //frac[Ds2KKPi][both][dn] = frac[Ds2KKPi][both][both];

            
            PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooRealVar(Form("yield_DsPhi_Ds2PhiPi_PR_total_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Total PR yield D_{s}#phi (#phi#pi)",   380, 0, 1000000);
            PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooRealVar(Form("yield_DsPhi_Ds2KKPi_PR_total_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Total PR yield D_{s}#phi (KK#pi)",    380, 0, 1000000);
            PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_DsPhi_Ds2PiPiPi_PR_total_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Total PR yield D_{s}#phi (#pi#pi#pi)", 50, 0, 200000 );
            PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooRealVar(Form("yield_DsPhi_Ds2KPiPi_PR_total_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Total PR yield D_{s}#phi (K#pi#pi)",   30, 0, 100000 );

            //HORNS:
            yield_XXst[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsPhi_XXst_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}#phi^{*} (#phi#pi)", "@0*@1",   RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a])  );
            yield_XXst[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooFormulaVar(Form("yield_DsPhi_XXst_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield D_{s}#phi^{*} (KK#pi)", "@0*@1",     RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a],   *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a])   );
            yield_XXst[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooFormulaVar(Form("yield_DsPhi_XXst_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield D_{s}#phi^{*} (#pi#pi#pi)", "@0*@1", RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a], *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]) );
            yield_XXst[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsPhi_XXst_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}#phi^{*} (K#pi#pi)", "@0*@1",   RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a])  );
            //HILL:
            yield_XstX[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsPhi_XstX_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}^{*}#phi (#phi#pi)", "(1-@0)*@1",   RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a])  );
            yield_XstX[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooFormulaVar(Form("yield_DsPhi_XstX_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "Yield D_{s}^{*}#phi (KK#pi)", "(1-@0)*@1",     RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a],   *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a])   );
            yield_XstX[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooFormulaVar(Form("yield_DsPhi_XstX_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "Yield D_{s}^{*}#phi (#pi#pi#pi)", "(1-@0)*@1", RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a], *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]) );
            yield_XstX[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooFormulaVar(Form("yield_DsPhi_XstX_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "Yield D_{s}^{*}#phi (K#pi#pi)", "(1-@0)*@1",   RooArgList(*frac[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a],  *PR_total_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a])  );




            // -------------------------------------------------------------------------

          if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Setting up blinding"<<std::endl;
          // RooCategory used for blind/unblind switching.
          TString blind("blind"), unblind("unblind");
          RooCategory blindCatBu("blindCatBu","Bu blind state Category");
          blindCatBu.defineType(unblind, 0);
          blindCatBu.defineType(blind, 1);
          //if(isBlind)
          if (needsBlinding)
              blindCatBu.setLabel(blind);
          else
              blindCatBu.setLabel(unblind);

          //RooUnblindUniform
          
          if(needsBlinding){
          //if(isBlind){      
            B_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = yield_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = yield_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];

            //B_yield[DsD0][Ds2PhiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2PhiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "nsig Bu blind (#phi#pi)",  "nsigBuDsD0blindedPhiPi",  1000., *yield_peak[DsD0][Ds2PhiPi][*c][*a] );
            //B_yield[DsD0][Ds2KKPi][*c][*a]   = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "nsig Bu blind (KK#pi)",    "nsigBuDsD0blindedKKPi",   1000., *yield_peak[DsD0][Ds2KKPi][*c][*a]  );
            //B_yield[DsD0][Ds2PiPiPi][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "nsig Bu blind (#pi#pi#pi)","nsigBuDsD0blindedPiPiPi", 100., *yield_peak[DsD0][Ds2PiPiPi][*c][*a] );
            //B_yield[DsD0][Ds2KPiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsD0_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "nsig Bu blind (K#pi#pi)",  "nsigBuDsD0blindedKPiPi",  100., *yield_peak[DsD0][Ds2KPiPi][*c][*a]  );
            
            // Only Blind DsPhi numbers...
            B_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2PhiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (#phi#pi)",  "nsigBuDsPhiblindePhiPi",   40., *yield_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a] );
            B_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2KKPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),   "nsig Bu blind (KK#pi)",    "nsigBuDsPhiblindedKKPi",   40., *yield_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]  );
            B_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2PiPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()), "nsig Bu blind (#pi#pi#pi)","nsigBuDsPhiblindedPiPiPi", 20., *yield_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]);
            B_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2KPiPi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str()),  "nsig Bu blind (K#pi#pi)",  "nsigBuDsPhiblindedKPiPi",  10., *yield_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a] );
          } else {
            B_yield[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]   = yield_peak[*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a] = yield_peak[*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];

            B_yield[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]   = yield_peak[*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a];
            B_yield[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]  = yield_peak[*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a];
            B_yield[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a] = yield_peak[*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a];
          }

          // -------------------------------------------------------------------------



          } //end of loop over chargeList
        }  //end of loop over magnetList
      } 
    } 
  }



  // --------- Mean B mass ---------------
  double default_mB=5279.29;
  
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            mean_B[*h][*ds][*ph][*b][*m][both][both] = new RooRealVar(Form("mean_B_%s_%s_%s_%s_%s"           ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][plus][up]   = new RooRealVar(Form("mean_B_%s_%s_%s_%s_plus_up_%s"   ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][minus][dn]  = new RooRealVar(Form("mean_B_%s_%s_%s_%s_minus_dn_%s"  ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][plus][dn]   = new RooRealVar(Form("mean_B_%s_%s_%s_%s_plus_dn_%s"   ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][minus][up]  = new RooRealVar(Form("mean_B_%s_%s_%s_%s_minus_up_%s"  ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][plus][both] = new RooRealVar(Form("mean_B_%s_%s_%s_%s_plus_both_%s" ,(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            mean_B[*h][*ds][*ph][*b][*m][minus][both]= new RooRealVar(Form("mean_B_%s_%s_%s_%s_minus_both_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
            //if(*m==d2kpi||*m==d2pik||*m==d2kk||*m==d2pipi||*m==d2kpipi0||*m==d2pikpi0){
            //  mean_B[*m][plus][both] = mean_B[*m][both][both];
            //  mean_B[*m][minus][both]= mean_B[*m][both][both];
            //}
            mean_B[*h][*ds][*ph][*b][*m][both][up]		= mean_B[*h][*ds][*ph][*b][*m][both][both];
            mean_B[*h][*ds][*ph][*b][*m][both][dn]		= mean_B[*h][*ds][*ph][*b][*m][both][both];
            mean_B[*h][*ds][*ph][*b][Ds2KPiPi][both][both]	= mean_B[*h][*ds][*ph][*b][Ds2KKPi][both][both];
            mean_B[*h][*ds][*ph][*b][Ds2PiPiPi][both][both] = mean_B[*h][*ds][*ph][*b][Ds2KKPi][both][both];
            mean_B[*h][*ds][*ph][*b][Ds2PhiPi][both][both]  = mean_B[*h][*ds][*ph][*b][Ds2KKPi][both][both];
          }
        }
      }
    }
  }

  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Making more RooRealVars"<<std::endl;
  // --------- Signal width ---------------
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            sigma[*h][*ds][*ph][*b][*m][both][both]  = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())          ,Form("%s sigma",mod[*b].c_str()),10,  5,  20); //change D0 to D^0
            sigma[*h][*ds][*ph][*b][*m][plus][both]  = sigma[*h][*ds][*ph][*b][*m][both][both];
            sigma[*h][*ds][*ph][*b][*m][minus][both] = sigma[*h][*ds][*ph][*b][*m][both][both];
            sigma[*h][*ds][*ph][*b][*m][plus][up]    = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_plus_up",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())  ,Form("%s sigma",mod[*b].c_str()),10,  5,  20);
            sigma[*h][*ds][*ph][*b][*m][plus][dn]    = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_plus_dn",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())  ,Form("%s sigma",mod[*b].c_str()),10,  5,  20);
            sigma[*h][*ds][*ph][*b][*m][minus][up]   = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_minus_up",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()) ,Form("%s sigma",mod[*b].c_str()),10,  5,  20);//sigma[*m][minus][up] = sigma[*m][plus][dn];
            sigma[*h][*ds][*ph][*b][*m][minus][dn]   = new RooRealVar(Form("sigma_%s_%s_%s_%s_%s_minus_dn",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()) ,Form("%s sigma",mod[*b].c_str()),10,  5,  20);//sigma[*m][minus][dn] = sigma[*m][plus][up];
            sigma[*h][*ds][*ph][*b][*m][both][up]    = sigma[*h][*ds][*ph][*b][*m][both][both]; 
            sigma[*h][*ds][*ph][*b][*m][both][dn]    = sigma[*h][*ds][*ph][*b][*m][both][both];
            
            sigma[*h][*ds][*ph][DsPhi][*m][both][dn] = sigma[*h][*ds][*ph][DsD0][*m][both][both];
            sigma[*h][*ds][*ph][DsPhi][*m][both][up] = sigma[*h][*ds][*ph][DsD0][*m][both][both];
          }
        }
      }
    }
  }

  // --------- Signal tail params ---------------
  
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            ncb[*h][*ds][*ph][*b][*m]     = new RooRealVar(Form("ncb_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())     , Form("%s nCB",mod[*b].c_str()),   1.000, 0,  100);
            alpha[*h][*ds][*ph][*b][*m]   = new RooRealVar(Form("alphaL_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())  , Form("%s alpha",mod[*b].c_str()), 2.2, 1, 10.0); //in Model.C it was 0.113, 0, 0.2  
          }
        }
        // Fixing parameters to KKPi
        alpha[*h][*ds][*ph][DsD0][Ds2PiPiPi] = alpha[*h][*ds][*ph][DsD0][Ds2KKPi]; 
        alpha[*h][*ds][*ph][DsD0][Ds2KPiPi]  = alpha[*h][*ds][*ph][DsD0][Ds2KKPi]; 
        alpha[*h][*ds][*ph][DsD0][Ds2PhiPi]  = alpha[*h][*ds][*ph][DsD0][Ds2KKPi]; 
        ncb[*h][*ds][*ph][DsD0][Ds2KPiPi]    = ncb[*h][*ds][*ph][DsD0][Ds2KKPi];
        ncb[*h][*ds][*ph][DsD0][Ds2PiPiPi]   = ncb[*h][*ds][*ph][DsD0][Ds2KKPi];
        ncb[*h][*ds][*ph][DsD0][Ds2PhiPi]    = ncb[*h][*ds][*ph][DsD0][Ds2KKPi];

        alpha[*h][*ds][*ph][DsPhi][Ds2PiPiPi] = alpha[*h][*ds][*ph][DsPhi][Ds2KKPi]; 
        alpha[*h][*ds][*ph][DsPhi][Ds2KPiPi]  = alpha[*h][*ds][*ph][DsPhi][Ds2KKPi]; 
        alpha[*h][*ds][*ph][DsPhi][Ds2PhiPi]  = alpha[*h][*ds][*ph][DsPhi][Ds2KKPi]; 
        ncb[*h][*ds][*ph][DsPhi][Ds2KPiPi]    = ncb[*h][*ds][*ph][DsPhi][Ds2KKPi];
        ncb[*h][*ds][*ph][DsPhi][Ds2PiPiPi]   = ncb[*h][*ds][*ph][DsPhi][Ds2KKPi];
        ncb[*h][*ds][*ph][DsPhi][Ds2PhiPi]    = ncb[*h][*ds][*ph][DsPhi][Ds2KKPi];

      }
    }
  }



  // --------- Comb. background slope ---------------
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            comb_slope[*h][*ds][*ph][*b][*m] = new RooRealVar(Form("comb_slope_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str())   ,Form("%s comb. slope",mod[*b].c_str()), -0.003, -1.0, -0.0000000001);
          }
        }
      }
    }
  }

  // --------- PartReco bkg: RooHORNSDini for DsD0*  ---------------
  for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            HORNS_a[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_a_%s_%s_%s%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HORNS a",      5051.8);
            HORNS_b[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_b_%s_%s_%s%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HORNS b",      5124.5);
            HORNS_csi[*h][*ds][*ph][*b][*m]    = new RooRealVar(Form("HORNS_csi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),   "HORNS csi",   0.72864);
            //HORNS_shift[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_shift_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HORNS shift", 0.59647);
            HORNS_shift[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_shift_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HORNS shift",  0.529,    0,    20);
            HORNS_sigma[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_sigma_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HORNS sigma",  8.6610);
            HORNS_R[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_R_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HORNS R",        15.5);
            HORNS_f[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_f_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HORNS f",     0.98133);
            
            //HORNS_a[*m]      = new RooRealVar(Form("HORNS_a_%s",(*m).c_str()),      "HORNS a",     5051.8, 5000,  5100);
            //HORNS_b[*m]      = new RooRealVar(Form("HORNS_b_%s",(*m).c_str()),      "HORNS b",     5124.5, 5100,  5180);
            //HORNS_csi[*m]    = new RooRealVar(Form("HORNS_csi_%s",(*m).c_str()),    "HORNS csi",   0.7286,    0,     1);
            //HORNS_shift[*m]  = new RooRealVar(Form("HORNS_shift_%s",(*m).c_str()),  "HORNS shift",  0.529,    0,    20);
            //HORNS_sigma[*m]  = new RooRealVar(Form("HORNS_sigma_%s",(*m).c_str()),  "HORNS sigma", 8.0954,    0,    30);
            //HORNS_R[*m]      = new RooRealVar(Form("HORNS_R_%s",(*m).c_str()),      "HORNS R",       0.68,    0,    50);
            //HORNS_f[*m]      = new RooRealVar(Form("HORNS_f_%s",(*m).c_str()),      "HORNS f",       0.98,    0,     1);
          }
        }


        // Fixing parameters to KKPi
        HORNS_shift[*h][*ds][*ph][DsD0][Ds2PiPiPi]  = HORNS_shift[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_shift[*h][*ds][*ph][DsD0][Ds2KPiPi]   = HORNS_shift[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_shift[*h][*ds][*ph][DsD0][Ds2PhiPi]   = HORNS_shift[*h][*ds][*ph][DsD0][Ds2KKPi];

        HORNS_sigma[*h][*ds][*ph][DsD0][Ds2PiPiPi]  = HORNS_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_sigma[*h][*ds][*ph][DsD0][Ds2KPiPi]   = HORNS_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_sigma[*h][*ds][*ph][DsD0][Ds2PhiPi]   = HORNS_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];

        HORNS_R[*h][*ds][*ph][DsD0][Ds2PiPiPi]      = HORNS_R[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_R[*h][*ds][*ph][DsD0][Ds2KPiPi]       = HORNS_R[*h][*ds][*ph][DsD0][Ds2KKPi];
        HORNS_R[*h][*ds][*ph][DsD0][Ds2PhiPi]       = HORNS_R[*h][*ds][*ph][DsD0][Ds2KKPi];

        HORNS_shift[*h][*ds][*ph][DsPhi][Ds2PiPiPi]  = HORNS_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_shift[*h][*ds][*ph][DsPhi][Ds2KPiPi]   = HORNS_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_shift[*h][*ds][*ph][DsPhi][Ds2PhiPi]   = HORNS_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];

        HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2PiPiPi]  = HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2KPiPi]   = HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2PhiPi]   = HORNS_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];

        HORNS_R[*h][*ds][*ph][DsPhi][Ds2PiPiPi]      = HORNS_R[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_R[*h][*ds][*ph][DsPhi][Ds2KPiPi]       = HORNS_R[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HORNS_R[*h][*ds][*ph][DsPhi][Ds2PhiPi]       = HORNS_R[*h][*ds][*ph][DsPhi][Ds2KKPi];

        // --------- PartReco bkg: RooHILLDini for Ds*D0 ---------------     
        for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            HILL_a[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_a_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HILL a",       4870.0);
            HILL_b[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_b_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HILL b",       5231.9);
            HILL_csi[*h][*ds][*ph][*b][*m]    = new RooRealVar(Form("HILL_csi_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),   "HILL csi",    0.80687);
            //HILL_shift[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_shift_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HILL shift",  -13.620);
            HILL_shift[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_shift_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HILL shift",   -14.2,   -100,   100);
            HILL_sigma[*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_sigma_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()), "HILL sigma",        8);
            HILL_R[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_R_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HILL R",          0.5);
            HILL_f[*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_f_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str()),     "HILL f",       0.8821);

            //HILL_a[*m]      = new RooRealVar(Form("HILL_a_%s",(*m).c_str()),     "HILL a",      4779.5,   4700,  4900);
            //HILL_b[*m]      = new RooRealVar(Form("HILL_b_%s",(*m).c_str()),     "HILL b",      5231.9,   5127,  5327);
            //HILL_csi[*m]    = new RooRealVar(Form("HILL_csi_%s",(*m).c_str()),   "HILL csi",     0.807,    -10,    10);
            //HILL_shift[*m]  = new RooRealVar(Form("HILL_shift_%s",(*m).c_str()), "HILL shift",   -14.2,   -100,   100);
            //HILL_sigma[*m]  = new RooRealVar(Form("HILL_sigma_%s",(*m).c_str()), "HILL sigma",   0.918,      0,    20);
            //HILL_R[*m]      = new RooRealVar(Form("HILL_R_%s",(*m).c_str()),     "HILL R",        11.3,      0,    20);
            //HILL_f[*m]      = new RooRealVar(Form("HILL_f_%s",(*m).c_str()),     "HILL f",       0.882,      0,     1);
          }  
        }

        // Fixing parameters to KKPi
      	HILL_a[*h][*ds][*ph][DsD0][Ds2PiPiPi]	    = HILL_a[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_a[*h][*ds][*ph][DsD0][Ds2KPiPi]      = HILL_a[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_a[*h][*ds][*ph][DsD0][Ds2PhiPi]      = HILL_a[*h][*ds][*ph][DsD0][Ds2KKPi];

      	HILL_shift[*h][*ds][*ph][DsD0][Ds2PiPiPi]	= HILL_shift[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_shift[*h][*ds][*ph][DsD0][Ds2KPiPi]  = HILL_shift[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_shift[*h][*ds][*ph][DsD0][Ds2PhiPi]  = HILL_shift[*h][*ds][*ph][DsD0][Ds2KKPi];

      	HILL_sigma[*h][*ds][*ph][DsD0][Ds2PiPiPi]	= HILL_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_sigma[*h][*ds][*ph][DsD0][Ds2KPiPi]  = HILL_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_sigma[*h][*ds][*ph][DsD0][Ds2PhiPi]  = HILL_sigma[*h][*ds][*ph][DsD0][Ds2KKPi];

      	HILL_R[*h][*ds][*ph][DsD0][Ds2PiPiPi]	    = HILL_R[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_R[*h][*ds][*ph][DsD0][Ds2KPiPi]      = HILL_R[*h][*ds][*ph][DsD0][Ds2KKPi];
        HILL_R[*h][*ds][*ph][DsD0][Ds2PhiPi]      = HILL_R[*h][*ds][*ph][DsD0][Ds2KKPi];

        HILL_a[*h][*ds][*ph][DsPhi][Ds2PiPiPi]     = HILL_a[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_a[*h][*ds][*ph][DsPhi][Ds2KPiPi]      = HILL_a[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_a[*h][*ds][*ph][DsPhi][Ds2PhiPi]      = HILL_a[*h][*ds][*ph][DsPhi][Ds2KKPi];

        HILL_shift[*h][*ds][*ph][DsPhi][Ds2PiPiPi] = HILL_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_shift[*h][*ds][*ph][DsPhi][Ds2KPiPi]  = HILL_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_shift[*h][*ds][*ph][DsPhi][Ds2PhiPi]  = HILL_shift[*h][*ds][*ph][DsPhi][Ds2KKPi];

        HILL_sigma[*h][*ds][*ph][DsPhi][Ds2PiPiPi] = HILL_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_sigma[*h][*ds][*ph][DsPhi][Ds2KPiPi]  = HILL_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_sigma[*h][*ds][*ph][DsPhi][Ds2PhiPi]  = HILL_sigma[*h][*ds][*ph][DsPhi][Ds2KKPi];

        HILL_R[*h][*ds][*ph][DsPhi][Ds2PiPiPi]     = HILL_R[*h][*ds][*ph][DsPhi][Ds2KKPi];
        HILL_R[*h][*ds][*ph][DsPhi][Ds2KPiPi]      = HILL_R[*h][*ds][*ph][DsPhi][Ds2KKPi];  
        HILL_R[*h][*ds][*ph][DsPhi][Ds2PhiPi]      = HILL_R[*h][*ds][*ph][DsPhi][Ds2KKPi];  

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
            // Use single value of tail parameters for all modes 
            ncb[*h][*ds][*ph][*b][*m]   = ncb[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];
            alpha[*h][*ds][*ph][*b][*m] = alpha[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];

            // Fix comb slope to be the same 
            comb_slope[*h][*ds][*ph][*b][*m] = comb_slope[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];

            // Want Horns and Hills to be different between D and B modes
            HILL_a[*h][*ds][*ph][*b][*m]     = HILL_a[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_b[*h][*ds][*ph][*b][*m]     = HILL_b[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_csi[*h][*ds][*ph][*b][*m]   = HILL_csi[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_shift[*h][*ds][*ph][*b][*m] = HILL_shift[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_sigma[*h][*ds][*ph][*b][*m] = HILL_sigma[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_R[*h][*ds][*ph][*b][*m]     = HILL_R[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HILL_f[*h][*ds][*ph][*b][*m]     = HILL_f[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];

            HORNS_shift[*h][*ds][*ph][*b][*m]  = HORNS_shift[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_sigma[*h][*ds][*ph][*b][*m]  = HORNS_sigma[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_R[*h][*ds][*ph][*b][*m]      = HORNS_R[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_csi[*h][*ds][*ph][*b][*m]    = HORNS_csi[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_f[*h][*ds][*ph][*b][*m]      = HORNS_f[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_a[*h][*ds][*ph][*b][*m]      = HORNS_a[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];
            HORNS_b[*h][*ds][*ph][*b][*m]      = HORNS_b[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m];


            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                
                // Use same mean value for all plots 
                mean_B[*h][*ds][*ph][*b][*m][*a][*c] = mean_B[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][both][both]; 

                // Allow sigma to vary between Ds decay modes
                sigma[*h][*ds][*ph][*b][*m][*a][*c] = sigma[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*m][*magnetList.begin()][*chargeList.begin()]; 

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
  
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Setting up simultaneous PDF"<<std::endl;
  sim = new RooSimultaneous("model","Simultaneous model",*cat);

  RooAbsPdf *pdf_peak = 0;
  RooAbsPdf *pdf_XXst = 0;
  RooAbsPdf *pdf_XstX = 0;
  RooAbsPdf *pdf_comb = 0;
  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                std::string tag=(*h)+underscore+(*ds)+underscore+(*ph)+underscore+(*b)+underscore+(*m)+underscore+(*c)+underscore+(*a);
                if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Adding to sim pdf: "<< tag <<std::endl;
                //DsD0
                RooAbsReal* mub     = mean_B[*h][*ds][*ph][*b][*m][*c][*a];
                RooAbsReal* sig_CB  = sigma[*h][*ds][*ph][*b][*m][*c][*a];
                RooAbsReal* n_CB    = ncb[*h][*ds][*ph][*b][*m];
                RooAbsReal* alp_CB  = alpha[*h][*ds][*ph][*b][*m];
                pdf_peak  = new RooCBShape( Form("pdf_peak_%s",tag.c_str()), "", *mB, *mub, *sig_CB, *alp_CB, *n_CB );

                //comb bkg
                RooRealVar* comb = comb_slope[*h][*ds][*ph][*b][*m];
                pdf_comb = new RooExponential( Form("pdf_comb_%s",tag.c_str()), "",*mB,*comb);
            
                //XX*
                RooRealVar* XXst_a     = HORNS_a[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_b     = HORNS_b[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_csi   = HORNS_csi[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_shift = HORNS_shift[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_sigma = HORNS_sigma[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_R     = HORNS_R[*h][*ds][*ph][*b][*m];
                RooRealVar* XXst_f     = HORNS_f[*h][*ds][*ph][*b][*m];
                pdf_XXst = new RooHORNSdini(Form("pdf_XXst_%s",tag.c_str()), "", *mB, *XXst_a, *XXst_b, *XXst_csi, *XXst_shift, *XXst_sigma, *XXst_R, *XXst_f);

                //X*X
                RooRealVar* XstX_a     = HILL_a[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_b     = HILL_b[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_csi   = HILL_csi[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_shift = HILL_shift[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_sigma = HILL_sigma[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_R     = HILL_R[*h][*ds][*ph][*b][*m];
                RooRealVar* XstX_f     = HILL_f[*h][*ds][*ph][*b][*m];
                pdf_XstX = new RooHILLdini(Form("pdf_XstX_%s",tag.c_str()), "", *mB, *XstX_a, *XstX_b, *XstX_csi, *XstX_shift, *XstX_sigma, *XstX_R, *XstX_f);

                //List of PDFs and yields
                //With PartReco
                RooArgSet pdflist( *pdf_peak, *pdf_comb, *pdf_XXst, *pdf_XstX );
                RooArgSet nevents(       *B_yield[*h][*ds][*ph][*b][*m][*c][*a], *yield_comb[*h][*ds][*ph][*b][*m][*c][*a], *yield_XXst[*h][*ds][*ph][*b][*m][*c][*a], *yield_XstX[*h][*ds][*ph][*b][*m][*c][*a] ); 
                
                //Add to master PDF
                RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_%s",tag.c_str())       ,"",pdflist, nevents);
        	      
                std::stringstream str;
                str<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                
                sim->addPdf(*masterPdf,str.str().c_str());

                /*
                //--------------------------------------------------------------------------
                //if(isBlind)
                if(needsBlinding)
                  sim->addPdf(*masterPdf_BLIND,str.str().c_str());
                  //sim->addPdf(*masterPdf,str.str().c_str());
                else
                  sim->addPdf(*masterPdf_BLIND,str.str().c_str());
                  //sim->addPdf(*masterPdf,str.str().c_str());
                //--------------------------------------------------------------------------
                */
              }
            }
          } 
        }//closing of for loops over b, m, a and c
      }
    }
  } // closing loop over h
} //end of funcn DefineModel()



RooArgSet* DsPhiModel::GetParameters()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::GetParameters()"<<std::endl;
  RooArgSet* pars = sim->getParameters(RooArgSet(*mB,helBin,DsBDTBin,PhiBDTBin,Bmode,mode,charge,magnet));
  return pars;
}



void DsPhiModel::PrintResult()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::PrintResults()"<<std::endl;
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
                  float y_peak =   yield_peak[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =   yield_comb[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_XXst =   yield_XXst[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_XstX =   yield_XstX[*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  std::cout<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                  std::cout<< " Peak:   " <<y_peak<<std::endl;
                  std::cout<< " Comb:   "<<y_comb<<std::endl;
                  std::cout<< " XX^{*}: "<<y_XXst<<std::endl;
                  std::cout<< " X^{*}X: "<<y_XstX<<std::endl;
                  std::cout<< "---------------------------"<<std::endl;
                  std::cout<< "---------------------------"<<std::endl;
                }else{
                  float y_peak_minus =   yield_peak[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_comb_minus =   yield_comb[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_XXst_minus =   yield_XXst[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_XstX_minus =   yield_XstX[*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                  float y_peak_plus =    yield_peak[*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                  float y_comb_plus =    yield_comb[*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                  float y_XXst_plus =    yield_XXst[*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                  float y_XstX_plus =    yield_XstX[*h][*ds][*ph][*b][*m][plus][*a]->getVal();

               
                  std::cout<<(*h)<<", "<<(*ds)<<", "<<(*ph)<<", "<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                  std::cout<<" Peak minus: "<<y_peak_minus<<" comb. minus: "<<y_comb_minus<<" XX^{*} minus: "<<y_XXst_minus<<" X^{*}X: minus: "<<y_XstX_minus<<std::endl;
                  std::cout<<" Peak plus:  "<<y_peak_plus <<" comb. plus:  "<<y_comb_plus <<" XX^{*} plus:  "<<y_XXst_plus <<" X^{*}X: plus:  "<<y_XstX_plus <<std::endl;
                  std::cout<<"--------------------------------------------------"<<std::endl;
                  std::cout<<" Peak total: "<<y_peak_minus+y_peak_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<" XX^{*} total: "<<y_XXst_minus+y_XXst_plus<<" X^{*}X: total: "<<y_XstX_minus+y_XstX_plus<<std::endl;

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

std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > DsPhiModel::GetResult()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::GetResults()"<<std::endl;
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
                  fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Peak"]    = yield_peak[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  //fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Peakerr"] = yield_peak[*h][*ds][*ph][*b][*m][*c][*a]->getError();
                  fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Comb"]    = yield_comb[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  fit_results[*h][*ds][*ph][*b][*m][*c][*a]["XXst"]    = yield_XXst[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  fit_results[*h][*ds][*ph][*b][*m][*c][*a]["XstX"]    = yield_XstX[*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                  std::cout<< " Peak:   " <<fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Peak"]<<std::endl;
                  std::cout<< " Comb:   " <<fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Comb"]<<std::endl;
                  std::cout<< " XXst:   " <<fit_results[*h][*ds][*ph][*b][*m][*c][*a]["XXst"]<<std::endl;
                  std::cout<< " XstX:   " <<fit_results[*h][*ds][*ph][*b][*m][*c][*a]["XstX"]<<std::endl;
                  
                  mB->setRange("signal", 5233.5, 5333.5) ;
                  //RooAbsPdf* sum = sim->getPdf(cat_str.c_str());
                  //RooAbsPdf* comb = sim->getPdf(Form("pdf_comb_%s",cat_str.c_str()));
                  //RooAbsPdf* comb = (RooAbsPdf*) sum->FindObject(Form("pdf_comb_%s",cat_str.c_str()));
                  //RooAbsReal* ibkg_sig = sum->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"),RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())));  
                  double combslope = comb_slope[*h][*ds][*ph][*b][*m]->getVal();
                  RooRealVar* fixed_slope = new RooRealVar("fixed_slope","slope",  combslope);
                  RooAbsPdf* pdf_comb = new RooExponential( "Exp", "",*mB,*fixed_slope);
                  RooAbsReal* ibkg_sig = pdf_comb->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"));  
                  fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"] = (ibkg_sig->getVal()*yield_comb[*h][*ds][*ph][*b][*m][*c][*a]->getVal());
                  

                  std::cout<< " Comb_Reduced:   " <<fit_results[*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"]<<std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
 //closes if(not needsBlinding)

return fit_results;
}
