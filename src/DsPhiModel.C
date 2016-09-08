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
#include "RooFFTConvPdf.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "DsPhiModel.h"
#include "Parameters.h"
//#include "RooDoubleCB.h"
#include "RooLITTLEHORNSdini.h"
#include "RooHORNSdini.h"
#include "RooHILLdini.h"
#include "TMatrixTSym.h"

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
  if(par->Bmodes[DsPhiSide])   BmodeList.push_back(DsPhiSide);
 

  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){       
    mode.defineType( (*m).c_str()); 
  }

  for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){       
    Bmode.defineType( (*b).c_str()); 
  }
  
  if(par->genToys||par->readToys){
    typeList.push_back(toy);
  } else {
    if(par->splitYears) {   
      if(par->dsetsReq[s21])   typeList.push_back(s21);
      if(par->dsetsReq[s21r1]) typeList.push_back(s21r1);
      if(par->dsetsReq[s24])   typeList.push_back(s24);
      if(par->dsetsReq[s26])   typeList.push_back(s26);
      if(0==typeList.size())   typeList.push_back("UNDEFTYPE"); 
    } else {
      typeList.push_back(both);
    }
  }

  type.defineType(toy.c_str());
  type.defineType(s21.c_str());
  type.defineType(s21r1.c_str());
  type.defineType(s24.c_str());
  type.defineType(s26.c_str());
  type.defineType(both.c_str());


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

  cat = new RooCategory("cat","type/helBin/DsBDTBin/PhiBDTBin/Bmode/mode/charge/polarity");
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
          for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++)  {
              for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                  std::stringstream str;
                  str<<(*t)<<underscore<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
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
}


void DsPhiModel::DefineModel()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel()"<<std::endl;

  //import PartReco PDFs
  std::string locationRooKeys = "/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/pdfs/";

  TFile f1((locationRooKeys + "kest1_Bs02Dsa1_KKst.root").c_str(), "read");
  RooWorkspace* workspace1 = (RooWorkspace*)f1.Get("workspace");
  RooKeysPdf* PartRecoPDF_Bs02Dsa1     = (RooKeysPdf*)workspace1->pdf("kest1");
  f1.Close();

  TFile f2((locationRooKeys + "kest1_Bs02DsstKKst.root").c_str(), "read");
  RooWorkspace* workspace2 = (RooWorkspace*)f2.Get("workspace");
  RooKeysPdf* PartRecoPDF_Bs02DsstKKst  = (RooKeysPdf*)workspace2->pdf("kest1");
  f2.Close();

  // Convolve Dsa1 shapes with gaussians to take into account MC and Data resolutions differences
  // Sigma values taken from sigma = sqrt((sigma_Data)^2 - (sigma_MC)^2) 
  // Values taken from smallest sigma (sigma1) in DsD0 mode 
  // Assumed MC/Data ratio is the same in DsPhi and DsD0

  RooRealVar* mg        = new RooRealVar("mg","mg",0) ;
  RooRealVar* sg_KKPi   = new RooRealVar("sg_KKPi",   "sg_KKPi",    4.39) ;
  RooRealVar* sg_KPiPi  = new RooRealVar("sg_KPiPi",  "sg_KPiPi",  14.09) ;
  RooRealVar* sg_PiPiPi = new RooRealVar("sg_PiPiPi", "sg_PiPiPi", 11.03) ;
  RooGaussian* gauss_KKPi   = new RooGaussian("gauss_KKPi",  "gauss_KKPi",  *mB,*mg,*sg_KKPi  ) ;
  RooGaussian* gauss_KPiPi  = new RooGaussian("gauss_KPiPi", "gauss_KPiPi", *mB,*mg,*sg_KPiPi ) ;
  RooGaussian* gauss_PiPiPi = new RooGaussian("gauss_PiPiPi","gauss_PiPiPi",*mB,*mg,*sg_PiPiPi) ;

  //std::map<std::string,RooFFTConvPdf*> PartRecoPDF_Bs02Dsa1_Convolved;
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2KKPi]   = new RooFFTConvPdf("kest1_adjust_KKPi",  "kest1_adjust_KKPi",  *mB,*PartRecoPDF_Bs02Dsa1, *gauss_KKPi  );
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2PhiPi]  = new RooFFTConvPdf("kest1_adjust_PhiPi", "kest1_adjust_PhiPi", *mB,*PartRecoPDF_Bs02Dsa1, *gauss_KKPi  );   
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2KPiPi]  = new RooFFTConvPdf("kest1_adjust_KPiPi", "kest1_adjust_KPiPi", *mB,*PartRecoPDF_Bs02Dsa1, *gauss_KPiPi );
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2PiPiPi] = new RooFFTConvPdf("kest1_adjust_PiPiPi","kest1_adjust_PiPiPi",*mB,*PartRecoPDF_Bs02Dsa1, *gauss_PiPiPi);

  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2KKPi]->setCacheObservables(  RooArgSet(*mg,*sg_KKPi));
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2PhiPi]->setCacheObservables( RooArgSet(*mg,*sg_KKPi));
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2KPiPi]->setCacheObservables( RooArgSet(*mg,*sg_KPiPi));
  //PartRecoPDF_Bs02Dsa1_Convolved[Ds2PiPiPi]->setCacheObservables(RooArgSet(*mg,*sg_PiPiPi));



  std::map<std::string,std::string> mod;
  std::map<std::string,std::string> Bmod;
  mod[Ds2PhiPi]="#phi#pi"; mod[Ds2KKPi]="KK#pi"; mod[Ds2KPiPi]="K#pi#pi";    mod[Ds2PiPiPi]="#pi#pi#pi";
  Bmod[DsD0]="D_{s}D^{0}"; Bmod[DsPhi]="D_{s}#phi"; Bmod[DsPhiSide]="D_{s}#phi Sideband";

  // ======================================================== 
  // =========== Values for Combinatoric ====================
  // ======================================================== 
  std::map<std::string,std::map<std::string,std::map<std::string,double>>> Yield_comb;
  Yield_comb[DsD0][Ds2PhiPi][both]  = 92; 
  Yield_comb[DsD0][Ds2KKPi][both]   = 348;
  Yield_comb[DsD0][Ds2PiPiPi][both] = 118;
  Yield_comb[DsD0][Ds2KPiPi][both]  = 164;

  Yield_comb[DsPhi][Ds2PhiPi][both]  = 123; 
  Yield_comb[DsPhi][Ds2KKPi][both]   = 280;
  Yield_comb[DsPhi][Ds2PiPiPi][both] = 146;
  Yield_comb[DsPhi][Ds2KPiPi][both]  = 182;
  
  Yield_comb[DsPhiSide][Ds2PhiPi][both]  = 22; 
  Yield_comb[DsPhiSide][Ds2KKPi][both]   = 73;
  Yield_comb[DsPhiSide][Ds2PiPiPi][both] = 20;
  Yield_comb[DsPhiSide][Ds2KPiPi][both]  = 37;


  Yield_comb[DsD0][Ds2KKPi][Helbin1]   = 157;
  Yield_comb[DsD0][Ds2KKPi][Helbin2]   = 62;
  Yield_comb[DsD0][Ds2KPiPi][Helbin1]  = 64;
  Yield_comb[DsD0][Ds2KPiPi][Helbin2]  = 51;
  Yield_comb[DsD0][Ds2PhiPi][Helbin1]  = 71; 
  Yield_comb[DsD0][Ds2PhiPi][Helbin2]  = 19; 
  Yield_comb[DsD0][Ds2PiPiPi][Helbin1] = 95;
  Yield_comb[DsD0][Ds2PiPiPi][Helbin2] = 39;
 
  Yield_comb[DsPhi][Ds2KKPi][Helbin1]   = 123;
  Yield_comb[DsPhi][Ds2KKPi][Helbin2]   = 102;
  Yield_comb[DsPhi][Ds2KPiPi][Helbin1]  = 68;
  Yield_comb[DsPhi][Ds2KPiPi][Helbin2]  = 37;
  Yield_comb[DsPhi][Ds2PhiPi][Helbin1]  = 113;
  Yield_comb[DsPhi][Ds2PhiPi][Helbin2]  = 67;
  Yield_comb[DsPhi][Ds2PiPiPi][Helbin1] = 140;
  Yield_comb[DsPhi][Ds2PiPiPi][Helbin2] = 82;

  
  Yield_comb[DsPhiSide][Ds2KKPi][Helbin1]   = 44;
  Yield_comb[DsPhiSide][Ds2KKPi][Helbin2]   = 42;
  Yield_comb[DsPhiSide][Ds2KPiPi][Helbin1]  = 49;
  Yield_comb[DsPhiSide][Ds2KPiPi][Helbin2]  = 36;
  Yield_comb[DsPhiSide][Ds2PhiPi][Helbin1]  = 26; 
  Yield_comb[DsPhiSide][Ds2PhiPi][Helbin2]  = 22;
  Yield_comb[DsPhiSide][Ds2PiPiPi][Helbin1] = 34;
  Yield_comb[DsPhiSide][Ds2PiPiPi][Helbin2] = 12;

  RooRealVar* global_comb_slope = new RooRealVar("global_comb_slope"  ,"Global Comb. Slope", -0.00258, -1.0, -0.0000000001);
  //global_comb_slope->setConstant();
  // ======================================================== 
  // =========== Values for Signal  =========================
  // ======================================================== 

  std::map<std::string,std::map<std::string,double>> Yield_CB_input, Fixed_CB_n, Fixed_CB_alpha, Fixed_CB_Sigma_ratio, Fixed_CB_fraction;
 
  Yield_CB_input[DsD0][Ds2KKPi]   = 1076;
  Yield_CB_input[DsD0][Ds2KPiPi]  = 162;
  Yield_CB_input[DsD0][Ds2PhiPi]  = 649; 
  Yield_CB_input[DsD0][Ds2PiPiPi] = 371;
 
  Yield_CB_input[DsPhi][Ds2KKPi]   = 25.6;
  Yield_CB_input[DsPhi][Ds2KPiPi]  = 2.5;
  Yield_CB_input[DsPhi][Ds2PhiPi]  = 14.2;
  Yield_CB_input[DsPhi][Ds2PiPiPi] = 16.2;

  //Yield_CB_input[DsPhiSide][Ds2PhiPi]  = 3.17; 
  //Yield_CB_input[DsPhiSide][Ds2KKPi]   = 5.70;
  //Yield_CB_input[DsPhiSide][Ds2PiPiPi] = 2.90;
  //Yield_CB_input[DsPhiSide][Ds2KPiPi]  = 0.65;

  // CB n -> Fixed to 1
  Fixed_CB_n[DsD0][Ds2PhiPi]      = 1.0; Fixed_CB_n[DsD0][Ds2KKPi]      = 1.0; Fixed_CB_n[DsD0][Ds2KPiPi]      = 1.0; Fixed_CB_n[DsD0][Ds2PiPiPi]      = 1.0;
  Fixed_CB_n[DsPhi][Ds2PhiPi]     = 1.0; Fixed_CB_n[DsPhi][Ds2KKPi]     = 1.0; Fixed_CB_n[DsPhi][Ds2KPiPi]     = 1.0; Fixed_CB_n[DsPhi][Ds2PiPiPi]     = 1.0;
  Fixed_CB_n[DsPhiSide][Ds2PhiPi] = 1.0; Fixed_CB_n[DsPhiSide][Ds2KKPi] = 1.0; Fixed_CB_n[DsPhiSide][Ds2KPiPi] = 1.0; Fixed_CB_n[DsPhiSide][Ds2PiPiPi] = 1.0;

  // CB alpha -> From MC 
  Fixed_CB_alpha[DsD0][Ds2PhiPi]      = 3.06956; Fixed_CB_alpha[DsD0][Ds2KKPi]      = 3.06956; Fixed_CB_alpha[DsD0][Ds2KPiPi]      = 3.56031; Fixed_CB_alpha[DsD0][Ds2PiPiPi]      = 3.41357;
  Fixed_CB_alpha[DsPhi][Ds2PhiPi]     = 3.10875; Fixed_CB_alpha[DsPhi][Ds2KKPi]     = 3.10875; Fixed_CB_alpha[DsPhi][Ds2KPiPi]     = 3.01317; Fixed_CB_alpha[DsPhi][Ds2PiPiPi]     = 3.10875;
  Fixed_CB_alpha[DsPhiSide][Ds2PhiPi] = 3.10875; Fixed_CB_alpha[DsPhiSide][Ds2KKPi] = 3.10875; Fixed_CB_alpha[DsPhiSide][Ds2KPiPi] = 3.01317; Fixed_CB_alpha[DsPhiSide][Ds2PiPiPi] = 3.10875;

  // 2 CB sigma ratios -> From MC 
  Fixed_CB_Sigma_ratio[DsD0][Ds2PhiPi]      = 0.525039; Fixed_CB_Sigma_ratio[DsD0][Ds2KKPi]      = 0.525039; Fixed_CB_Sigma_ratio[DsD0][Ds2KPiPi]      = 0.548737; Fixed_CB_Sigma_ratio[DsD0][Ds2PiPiPi]      = 0.548367;
  Fixed_CB_Sigma_ratio[DsPhi][Ds2PhiPi]     = 0.500433; Fixed_CB_Sigma_ratio[DsPhi][Ds2KKPi]     = 0.500433; Fixed_CB_Sigma_ratio[DsPhi][Ds2KPiPi]     = 0.459650; Fixed_CB_Sigma_ratio[DsPhi][Ds2PiPiPi]     = 0.503421;
  Fixed_CB_Sigma_ratio[DsPhiSide][Ds2PhiPi] = 0.500433; Fixed_CB_Sigma_ratio[DsPhiSide][Ds2KKPi] = 0.500433; Fixed_CB_Sigma_ratio[DsPhiSide][Ds2KPiPi] = 0.459650; Fixed_CB_Sigma_ratio[DsPhiSide][Ds2PiPiPi] = 0.503421;

  // 2 CB sigma ratios -> From MC 
  Fixed_CB_fraction[DsD0][Ds2PhiPi]      = 0.806492; Fixed_CB_fraction[DsD0][Ds2KKPi]      = 0.806492; Fixed_CB_fraction[DsD0][Ds2KPiPi]      = 0.806492; Fixed_CB_fraction[DsD0][Ds2PiPiPi]      = 0.750798;
  Fixed_CB_fraction[DsPhi][Ds2PhiPi]     = 0.818999; Fixed_CB_fraction[DsPhi][Ds2KKPi]     = 0.818999; Fixed_CB_fraction[DsPhi][Ds2KPiPi]     = 0.861946; Fixed_CB_fraction[DsPhi][Ds2PiPiPi]     = 0.797389;
  Fixed_CB_fraction[DsPhiSide][Ds2PhiPi] = 0.818999; Fixed_CB_fraction[DsPhiSide][Ds2KKPi] = 0.818999; Fixed_CB_fraction[DsPhiSide][Ds2KPiPi] = 0.861946; Fixed_CB_fraction[DsPhiSide][Ds2PiPiPi] = 0.797389;

  std::map<std::string,double> Fixed_Norm_Sigma_ratio;
  // Ratio of sigmas from DsD0 to DsPhi (small sigma)
  // Ratio = Sigma(DsPhi)/Sigma(DsD0)
  Fixed_Norm_Sigma_ratio[Ds2PhiPi] = 1.33580;  Fixed_Norm_Sigma_ratio[Ds2KKPi] = 1.33580;  Fixed_Norm_Sigma_ratio[Ds2KPiPi] = 1.42894;  Fixed_Norm_Sigma_ratio[Ds2PiPiPi] = 1.30976;

  // Initial Sigma values
  std::map<std::string,double> Inital_Sigma_values;
  Inital_Sigma_values[Ds2KKPi]   =  8.1;
  Inital_Sigma_values[Ds2KPiPi]  =  9.2;
  Inital_Sigma_values[Ds2PhiPi]  =  7.9;
  Inital_Sigma_values[Ds2PiPiPi] =  9.0; 


  // ======================================================== 
  // =========== Values for backgrounds      ================
  // ======================================================== 
  std::map<std::string,std::map<std::string,double>> Yield_Dsa1_input;
 
  Yield_Dsa1_input[DsPhi][Ds2PhiPi]  =  19;
  Yield_Dsa1_input[DsPhi][Ds2KKPi]   = 110; 
  Yield_Dsa1_input[DsPhi][Ds2PiPiPi] =  51; 
  Yield_Dsa1_input[DsPhi][Ds2KPiPi]  =  13; 

  std::map<std::string,std::map<std::string,double>> Yield_PR_total_input;
 
  Yield_PR_total_input[DsD0][Ds2KKPi]   = 2186; 
  Yield_PR_total_input[DsD0][Ds2KPiPi]  =  330;
  Yield_PR_total_input[DsD0][Ds2PhiPi]  = 1395; 
  Yield_PR_total_input[DsD0][Ds2PiPiPi] =  699;

  Yield_PR_total_input[DsPhi][Ds2KKPi]   =   127; 
  Yield_PR_total_input[DsPhi][Ds2KPiPi]  =   0.1;
  Yield_PR_total_input[DsPhi][Ds2PhiPi]  =    16; 
  Yield_PR_total_input[DsPhi][Ds2PiPiPi] =   115; 

  //Fraction of Ds*D0 initial value
  std::map<std::string,double> Fraction_DsstD0_initial;
  Fraction_DsstD0_initial[Ds2PhiPi]  = 0.37;//0.60;
  Fraction_DsstD0_initial[Ds2KKPi]   = 0.37; //0.68;
  Fraction_DsstD0_initial[Ds2PiPiPi] = 0.49;//0.75;
  Fraction_DsstD0_initial[Ds2KPiPi]  = 0.37;//0.56;
  RooRealVar* Fraction_DsstD0        = new RooRealVar("Fraction_DsstD0", "" , 0.49, 0 ,1); 
  RooRealVar* Fraction_DsstD0_Helbin1= new RooRealVar("Fraction_DsstD0_Helbin1", "" , 0.429, 0 ,1); 
  RooRealVar* Fraction_DsstD0_Helbin2= new RooRealVar("Fraction_DsstD0_Helbin2", "" , 0.362, 0 ,1); 

  //Fraction_DsstD0_Helbin1->setConstant();
  //Fraction_DsstD0_Helbin2->setConstant();


  RooRealVar* DsPhi_BG_fraction      = new RooRealVar("Dsa1_to_DsstPhi_fraction", "" , 0.15 , 0 , 1);
  //DsPhi_BG_fraction->setConstant();
  // ======================================================== 
  // =========== Fractions in different cats ================
  // ======================================================== 

  //Signal fractions 
  // Fraction of signal decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* Signal_DsPhi_Fraction   = new RooRealVar("Signal_DsPhi_Fraction",   "", 0.958294693);
  // Fraction of signal decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Signal_Helbin1_Fraction = new RooRealVar("Signal_Helbin1_Fraction", "", 0.784945192);

  //Dsa1 fractions 
  // Fraction of Dsa1 decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* Dsa1_DsPhi_Fraction     = new RooRealVar("Dsa1_DsPhi_Fraction",     "", 0.667117117);
  // Fraction of Dsa1 decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Dsa1_Helbin1_Fraction   = new RooRealVar("Dsa1_Helbin1_Fraction",   "", 0.641091721);

  //Ds*Phi fractions 
  // Fraction of Ds*Phi decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* DsstPhi_DsPhi_Fraction = new RooRealVar("DsstPhi_DsPhi_Fraction", "", 0.958294693);

  // ======================================================== 
  // =========== Global variables      ======================
  // ========================================================
  RooRealVar* global_csi_hill   = new RooRealVar("global_csi_hill",   "Global csi HILL",    1);
  RooRealVar* global_csi   = new RooRealVar("global_csi",   "Global csi",     1.12, 0.1,    10);
  RooRealVar* global_shift = new RooRealVar("global_shift", "Global shift",  -2.33, -10,    10);
  //global_shift->setConstant();
  //global_csi->setConstant();

  //RooRealVar* global_csi   = new RooRealVar("global_csi",   "Global csi",     0.76);
  //RooRealVar* global_shift = new RooRealVar("global_shift", "Global shift",  -4.15);

  // Fraction of D* -> D pi0 vs. D* -> D pi0 gamma adjusted for fraction in mass window
  RooRealVar* fraction_Dst0_D0pi0_DsD0_adj = new RooRealVar("fraction_Dst0_D0pi0_DsD0", "",  0.66 * (0.838/0.713));
  RooRealVar* fraction_Dsst_Dspi0_DsD0_adj = new RooRealVar("fraction_Dsst_Dspi0_DsD0_adj", "",  0.06 * (0.841/0.732));

  RooRealVar* fraction_Dsst_Dspi0_DsPhi_010_adj = new RooRealVar("fraction_Dsst_Dspi0_DsPhi_adj", "",  0.06 * (0.531/0.610));

  RooRealVar* fraction_Dsst_Dspi0_DsPhi_101_adj = new RooRealVar("fraction_Dsst_Dspi0_DsPhi_adj", "",  0.06 * (0.779/0.544));
  //RooRealVar* fraction_hel        = new RooRealVar("fraction_hel",        "",  0.83, 0, 1);
  RooRealVar* fraction_hel        = new RooRealVar("fraction_hel",        "",  0.5);


  // ======================================================== 
  // =========== Fixing ratios between Ds modes =============
  // =========== in Hel split mode to improve   =============
  // =========== stability                      =============
  // ========================================================

  // Fraction of DsD0 peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsD0_peak_fraction = new RooRealVar("splitHel_DsD0_peak_fraction", "", 0.535, 0.0 ,1.0 );
  //splitHel_DsD0_peak_fraction->setConstant();

  // Fraction of DsD0 PR peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsD0_PR_peak_fraction = new RooRealVar("splitHel_DsD0_PR_peak_fraction", "", 0.531, 0.0 ,1.0 );
  //splitHel_DsD0_PR_peak_fraction->setConstant();
  
  // Fraction of Ds*Phi peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsstPhi_peak_fraction = new RooRealVar("splitHel_DsstPhi_peak_fraction", "", 0.920, 0.0 ,1.0 );
  //splitHel_DsstPhi_peak_fraction->setConstant();
  // Fraction of Dsa1 peak in two helicity bins to be the same for all Ds modes
  //RooRealVar* splitHel_Dsa1_peak_fraction   = new RooRealVar("splitHel_Dsa1_peak_fraction", "", 0.51, 0.0 ,1.0 );
  RooRealVar* splitHel_Dsa1_peak_fraction   = new RooRealVar("splitHel_Dsa1_peak_fraction", "", 0.641091721 );


  //RooRealVar* Dsa1_Fraction     = new RooRealVar("Dsa1_fraction",   "Dsa1_fraction",    0.5);
  //RooRealVar* DsstKKst_Fraction = new RooRealVar("DsstKKst_fraction",   "DsstKKst_fraction",  0.5, 0.0, 1.0);
  //
  
  // ========================================================
  // =========== Blinding numbers               =============
  // ========================================================
  std::map<std::string,std::map<std::string,std::map<std::string,double>>> Blind_number; 
  Blind_number[DsPhi][Ds2PhiPi][both]  = 15;
  Blind_number[DsPhi][Ds2KKPi][both]   = 25;
  Blind_number[DsPhi][Ds2PiPiPi][both] = 11;
  Blind_number[DsPhi][Ds2KPiPi][both]  = 5;
  
  Blind_number[DsPhi][Ds2PhiPi][Helbin1]  = 7;
  Blind_number[DsPhi][Ds2KKPi][Helbin1]   = 10;
  Blind_number[DsPhi][Ds2PiPiPi][Helbin1] = 5;
  Blind_number[DsPhi][Ds2KPiPi][Helbin1]  = 3;
  
  Blind_number[DsPhi][Ds2PhiPi][Helbin2]  = 7;
  Blind_number[DsPhi][Ds2KKPi][Helbin2]   = 10;
  Blind_number[DsPhi][Ds2PiPiPi][Helbin2] = 5;
  Blind_number[DsPhi][Ds2KPiPi][Helbin2]  = 3;

  Blind_number[DsPhiSide][Ds2PhiPi][both]  = 1.5;
  Blind_number[DsPhiSide][Ds2KKPi][both]   = 2.5;
  Blind_number[DsPhiSide][Ds2PiPiPi][both] = 1.1;
  Blind_number[DsPhiSide][Ds2KPiPi][both]  = 1.0;
  
  Blind_number[DsPhiSide][Ds2PhiPi][Helbin1]  = 1.0;
  Blind_number[DsPhiSide][Ds2KKPi][Helbin1]   = 1.0;
  Blind_number[DsPhiSide][Ds2PiPiPi][Helbin1] = 1.0;
  Blind_number[DsPhiSide][Ds2KPiPi][Helbin1]  = 1.0;
  
  Blind_number[DsPhiSide][Ds2PhiPi][Helbin2]  = 1.0;
  Blind_number[DsPhiSide][Ds2KKPi][Helbin2]   = 1.0;
  Blind_number[DsPhiSide][Ds2PiPiPi][Helbin2] = 1.0;
  Blind_number[DsPhiSide][Ds2KPiPi][Helbin2]  = 1.0;


  // --------- Yields for DsD0, in an array (as line 447 of Model.C) ---------------
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Making RooRealVars"<<std::endl;
 
  // ------------------------------------------------------
  // --------- Setup yields for signals       -------------
  // ------------------------------------------------------
 
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){  
        for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){

              // ------------------------------------------------------
              // --------- Setup yields for signals       -------------
              // ------------------------------------------------------

              // Fix yields for not helicity split: both
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              
              yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",      DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), Yield_CB_input[DsD0][*m],    0, 100000);
              //((RooRealVar*)yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
              Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]   = new RooRealVar(   Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -2, 100000); 
              
              //if(needsBlinding){
              //  Signal_total_unblind[*t][both][*ds][*ph][DsPhi][*m][*c][*a]= new RooRealVar(       Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -2, 100000); 
              //  Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]        = new RooUnblindUniform(Form("B_nsig_total_%s",cat_name.c_str()),  "nsig total",  "nsigBuDsPhiblindeTOTAL",   Blind_number[*m][both], *Signal_total_unblind[*t][both][*ds][*ph][DsPhi][*m][*c][*a] );
              //}else{
              //    Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]      = new RooRealVar(       Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -2, 100000); 
              //}
              
              yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_DsPhi_Fraction, *Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(1-@0)*@1" , RooArgList(*Signal_DsPhi_Fraction, *Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              // Blinding
              if(needsBlinding){
                B_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DsPhi_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhi",       Blind_number[DsPhi][*m][both],     *yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a] );
                B_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhiSide_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhiSide",   Blind_number[DsPhiSide][*m][both], *yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] );
              } else {
                B_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a];
                B_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a];
              }  
              // Fix yields for helicity split state: Helbin1
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              yield_peak[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsD0].c_str(),       mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*splitHel_DsD0_peak_fraction, *yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),      mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_peak[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhiSide].c_str(),  mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              // Blinding
              if(needsBlinding){
                B_yield[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DsPhi_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhi",       Blind_number[DsPhi][*m][Helbin1],     *yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a] );
                B_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhiSide_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhiSide",   Blind_number[DsPhiSide][*m][Helbin1], *yield_peak[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] );
              } else {
                B_yield[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
                B_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = yield_peak[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a];
              }
              
              // Fix yields for helicity split state: Helbin2
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              yield_peak[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*splitHel_DsD0_peak_fraction, *yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              yield_peak[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_peak[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              if(needsBlinding){
                B_yield[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DsPhi_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhi",       Blind_number[DsPhi][*m][Helbin2],     *yield_peak[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a] );
                B_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhiSide_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhiSide",   Blind_number[DsPhiSide][*m][Helbin2], *yield_peak[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] );
              } else {
                B_yield[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = yield_peak[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a];
                B_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = yield_peak[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a];
              }   




              /*
              yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), Yield_CB_input[DsD0][*m],    0, 100000);

              yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m], -2, 100000);
              yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*@1", RooArgList(*Signal_DsPhi_Fraction,*yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a])); 
              
              // Fix yields for helicity split state: Helbin1
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              
              //yield_peak[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), Signal_Helbin1_Fraction->getVal()*Yield_CB_input[DsD0][*m],    0, 100000);
              
              yield_peak[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), "@0*@1" ,         RooArgList(*splitHel_DsD0_peak_fraction, *yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              //yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Signal_Helbin1_Fraction->getVal()*Yield_CB_input[DsPhi][*m], -2, 100000);
              //yield_peak[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*@1", RooArgList(*Signal_DsPhi_Fraction,*yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              
              yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "@0*@1", RooArgList(*Signal_Helbin1_Fraction,*yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_peak[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "@0*@1", RooArgList(*Signal_Helbin1_Fraction,*yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
              
              // Fix yields for helicity split state: Helbin1
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              //yield_peak[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h2",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), (1-Signal_Helbin1_Fraction->getVal())*Yield_CB_input[DsD0][*m],   0, 100000);
              
              yield_peak[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), "(1-@0)*@1" ,                 RooArgList(*splitHel_DsD0_peak_fraction, *yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              //yield_peak[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h2",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "((1-@0)/@0)*@1",             RooArgList(*Signal_Helbin1_Fraction,*yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              //yield_peak[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s h2",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*((1-@1)/@1)*@2", RooArgList(*Signal_DsPhi_Fraction,  *Signal_Helbin1_Fraction,*yield_peak[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              yield_peak[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s h2",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*Signal_Helbin1_Fraction,*yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a])); 
              yield_peak[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s h1",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*Signal_Helbin1_Fraction,*yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
              
              */

              // ------------------------------------------------------
              // --------- Setup yields for backgrounds   -------------
              // ------------------------------------------------------

              
              // Both
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());              
               
              Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsD0][*m],  0,  1000000);
              Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsPhi][*m], 0,  1000000);
              //((RooRealVar*)Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
              //((RooRealVar*)Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a])->setConstant();


              PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "@0",             RooArgList(*Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a] ));//;
              PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",          RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));//;
              PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1", RooArgList(*DsstPhi_DsPhi_Fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1" ,    RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "((1-@0)/@0)*@1",RooArgList(*Dsa1_DsPhi_Fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a])); 
              
              // Helbin 1
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              PR_total_yield[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "@0*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));
              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
              
              
              // Helbin 2
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              PR_total_yield[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));
              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "(1-@0)*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
              





              /*

              PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  Yield_PR_total_input[DsD0][*m],  0,  1000000);
              
              Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsPhi][*m], 0,  1000000);
              
              PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",  RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]) ); //Yield_PR_total_input[DsPhi][*m], 0,  1000000);
              //PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsPhi][*m], 0,  1000000);
              PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1", RooArgList(*DsstPhi_DsPhi_Fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              
              yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a]      = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),   "(1-@0)*@1" , RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              //yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a]      = new RooRealVar(   Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),   Dsa1_DsPhi_Fraction->getVal()*Yield_Dsa1_input[DsPhi][*m], 0, 1000000);
              yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]  = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*@1", RooArgList(*Dsa1_DsPhi_Fraction,  *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a])); 
                 
              

              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              //yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h1",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),   Yield_Dsa1_input[DsPhi][*m], 0, 1000000);
              
              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s",   Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "(1-@0)*@1" ,     RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h1",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*@1", RooArgList(*Dsa1_DsPhi_Fraction,  *yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              //yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h1",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_Dsa1_peak_fraction,*yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] );
              
              PR_total_yield[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "@0*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              

              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h2",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "((1-@0)/@0)*@1",             RooArgList(*Dsa1_Helbin1_Fraction,  *yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhiSide.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h2",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()), "((1-@0)/@0)*((1-@1)/@1)*@2", RooArgList(*Dsa1_DsPhi_Fraction,    *Dsa1_Helbin1_Fraction,*yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])); 
              
              //yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s h1",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_Dsa1_peak_fraction,*yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] );
              PR_total_yield[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",     Bmod[DsD0].c_str(),      mod[(*m)].c_str()),  "(1-@0)*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",     Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",     Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              */

            }
          }
        }
      }
    }
  }

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
          
              for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
                for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
                  sig_frac[*t][*h][*ds][*ph][*b][*m]           = new RooRealVar(   Form("Sigma_frac_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("%s %s Sigma Fraction",Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_fraction[*b][*m]); 
                  
                  yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a] = new RooRealVar(   Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), Yield_comb[*b][*m][*h], 0, 600000);
                  //((RooRealVar*)yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a])->setConstant();
                }
              } 

       
              for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
                
                //frac[*t][*h][*ds][*ph][DsD0][*m][*c][*a]                 = new RooRealVar(   Form("frac_%s_%s_%s",          DsD0.c_str(),      (*m).c_str(),cat_name.c_str()), Form("D_{s}*D^{0} frac %s", mod[(*m)].c_str()), Fraction_DsstD0_initial[*m],  0, 1);
                frac[*t][both][*ds][*ph][DsD0][*m][*c][*a]               = Fraction_DsstD0;
                frac[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]            = Fraction_DsstD0_Helbin1;
                frac[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]            = Fraction_DsstD0_Helbin2;
                //PR_total_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a]       = new RooRealVar(   Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(),cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  Yield_PR_total_input[DsD0][*m],  0,  1000000);
                //PR_total_yield[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]      = new RooRealVar(   Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(),cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsPhi][*m], 0,  1000000);
                //PR_total_yield[*t][*h][*ds][*ph][DsPhiSide][*m][*c][*a]  = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(),cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1", RooArgList(*DsstPhi_DsPhi_Fraction,*PR_total_yield[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]));
                
                yield_dsstd0[*t][*h][*ds][*ph][DsD0][*m][*c][*a]     = new RooFormulaVar(Form("yield_DsstD0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D_{s}* D^{0} %s",mod[(*m)].c_str()), "@0*@1",     RooArgList(*frac[*t][*h][*ds][*ph][DsD0][*m][*c][*a],  *PR_total_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a])  );
                yield_dsdst0[*t][*h][*ds][*ph][DsD0][*m][*c][*a]     = new RooFormulaVar(Form("yield_DsDst0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D_{s} D*^{0} %s",mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*frac[*t][*h][*ds][*ph][DsD0][*m][*c][*a],  *PR_total_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a])  );
              

                frac_HORNS[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HORNS_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*@1",        RooArgList(*fraction_Dsst_Dspi0_DsPhi_010_adj,*fraction_hel));
                frac_HILL[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]        = new RooFormulaVar(Form("fraction_HILL_%s_%s",       (*m).c_str(),cat_name.c_str()), "", "(1-@0)*@1",    RooArgList(*fraction_Dsst_Dspi0_DsPhi_010_adj,*fraction_hel));
                frac_HILL2[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HILL2_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*(1-@1)",    RooArgList(*fraction_Dsst_Dspi0_DsPhi_101_adj,*fraction_hel));
                frac_LITTLEHORNS[*t][*h][*ds][*ph][DsPhi][*m][*c][*a] = new RooFormulaVar(Form("fraction_LITTLEHORNS_%s_%s",(*m).c_str(),cat_name.c_str()), "", "(1-@0)*(1-@1)",RooArgList(*fraction_Dsst_Dspi0_DsPhi_101_adj,*fraction_hel));
              }

              // -------------------------------------------------------------------------

              if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Setting up blinding"<<std::endl;
              // RooCategory used for blind/unblind switching.
              TString blind("blind"), unblind("unblind");
              RooCategory blindCatBu("blindCatBu","Bu blind state Category");
              blindCatBu.defineType(unblind, 0);
              blindCatBu.defineType(blind, 1);
              if (needsBlinding)
                blindCatBu.setLabel(blind);
              else
                blindCatBu.setLabel(unblind);
              /*
              if(needsBlinding){

                B_yield[*t][*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]       = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]        = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]       = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a]      = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];
                
                // Only Blind DsPhi numbers...
                double split_ratio = 1.0;
                if(par->splitHel) split_ratio = 2.0;

                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]      = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2PhiPi_%s",cat_name.c_str()),  "nsig Bu blind DsPhi (#phi#pi)",  "nsigBuDsPhiblindePhiPi",   15.0/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a] );
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]       = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2KKPi_%s",cat_name.c_str()),   "nsig Bu blind DsPhi (KK#pi)",    "nsigBuDsPhiblindedKKPi",   25.0/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]  );
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2PiPiPi_%s",cat_name.c_str()), "nsig Bu blind DsPhi (#pi#pi#pi)","nsigBuDsPhiblindedPiPiPi", 15.0/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]);
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]      = new RooUnblindUniform(Form("B_nsig_DsPhi_Ds2KPiPi_%s",cat_name.c_str()),  "nsig Bu blind DsPhi (K#pi#pi)",  "nsigBuDsPhiblindedKPiPi",   5.0/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a] );
                
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2PhiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsPhiSide_Ds2PhiPi_%s",cat_name.c_str()),  "nsig Bu blind DsPhiSide (#phi#pi)",  "nsigBuDsPhiSideblindePhiPi",   1.50/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2PhiPi][*c][*a] );
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2KKPi][*c][*a]   = new RooUnblindUniform(Form("B_nsig_DsPhiSide_Ds2KKPi_%s",cat_name.c_str()),   "nsig Bu blind DsPhiSide (KK#pi)",    "nsigBuDsPhiSideblindedKKPi",   2.50/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2KKPi][*c][*a]  );
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2PiPiPi][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhiSide_Ds2PiPiPi_%s",cat_name.c_str()), "nsig Bu blind DsPhiSide (#pi#pi#pi)","nsigBuDsPhiSideblindedPiPiPi", 1.50/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2PiPiPi][*c][*a]);
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2KPiPi][*c][*a]  = new RooUnblindUniform(Form("B_nsig_DsPhiSide_Ds2KPiPi_%s",cat_name.c_str()),  "nsig Bu blind DsPhiSide (K#pi#pi)",  "nsigBuDsPhiSideblindedKPiPi",  0.50/split_ratio, *yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2KPiPi][*c][*a] );
              
              } else {
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a]       = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2PhiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a]        = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2KKPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a]       = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2KPiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a]      = yield_peak[*t][*h][*ds][*ph][DsD0][Ds2PiPiPi][*c][*a];

                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a]      = yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2PhiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a]       = yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2KKPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a]      = yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2KPiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a]     = yield_peak[*t][*h][*ds][*ph][DsPhi][Ds2PiPiPi][*c][*a];

                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2PhiPi][*c][*a]  = yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2PhiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2KKPi][*c][*a]   = yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2KKPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2KPiPi][*c][*a]  = yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2KPiPi][*c][*a];
                B_yield[*t][*h][*ds][*ph][DsPhiSide][Ds2PiPiPi][*c][*a] = yield_peak[*t][*h][*ds][*ph][DsPhiSide][Ds2PiPiPi][*c][*a];
              }
              */

            } //end of loop over chargeList
          }  //end of loop over magnetList
        } 
      } 
    }
  }


  // --------- Mean B mass ---------------
  double default_mB=5279.29;

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str());
              mean_B[*t][*h][*ds][*ph][*b][*m][both][both] = new RooRealVar(Form("mean_B_%s"           ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][up]   = new RooRealVar(Form("mean_B_%s_plus_up"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][dn]  = new RooRealVar(Form("mean_B_%s_minus_dn"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][dn]   = new RooRealVar(Form("mean_B_%s_plus_dn"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][up]  = new RooRealVar(Form("mean_B_%s_minus_up"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][both] = new RooRealVar(Form("mean_B_%s_plus_both" ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][both]= new RooRealVar(Form("mean_B_%s_minus_both",cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              //mean_B[*t][*h][*ds][*ph][*b][*m][both][both]->setConstant();

              mean_B[*t][*h][*ds][*ph][*b][*m][both][up]		      = mean_B[*t][*h][*ds][*ph][*b][*m][both][both];
              mean_B[*t][*h][*ds][*ph][*b][*m][both][dn]		      = mean_B[*t][*h][*ds][*ph][*b][*m][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2KPiPi][both][both]	= mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2PiPiPi][both][both] = mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2PhiPi][both][both]  = mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
            }
          }
        }
      }
    }
  }
  
  // --------- Signal width ---------------
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
        for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              //for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
                
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());

              sigma[*t][both][*ds][*ph][DsD0][*m][*c][*a]   = new RooRealVar(   Form("sigma_DsD0_%s", cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
              sigma[*t][both][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma_DsPhi_%s",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma[*t][both][*ds][*ph][DsD0][*m][*c][*a])); 
              sigma[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma[*t][both][*ds][*ph][DsPhi][*m][*c][*a];
              //sig_ratio[*t][both][*ds][*ph][*b][*m]         = new RooRealVar(Form("Sigma_ratio_%s", cat_name.c_str()), Form("%s %s Sigma ratio ",mod[*m].c_str(),Bmod[*b].c_str()), Fixed_CB_Sigma_ratio[*b][*m]);
      
              sigma2[*t][both][*ds][*ph][DsD0][*m][*c][*a]   = new RooFormulaVar(Form("sigma2_DsD0_%s", cat_name.c_str()), Form("Sigma 2 DsD0 %s",(*m).c_str()),  Form("@0/%f",Fixed_CB_Sigma_ratio[DsD0][*m] ), RooArgList(*sigma[*t][both][*ds][*ph][DsD0][*m][*c][*a]   ));
              sigma2[*t][both][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma2_DsPhi_%s",cat_name.c_str()), Form("Sigma 2 DsPhi %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DsPhi][*m]), RooArgList(*sigma[*t][both][*ds][*ph][DsPhi][*m][*c][*a]  ));
              sigma2[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma2[*t][both][*ds][*ph][DsPhi][*m][*c][*a];
              
              avg_sigma[*t][both][*ds][*ph][DsD0][*m][both][both]      = new RooFormulaVar(Form("avg_sigma_DsD0_%s",  cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsD0].c_str(), mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsD0][*m], Fixed_CB_fraction[DsD0][*m]  ), RooArgList(*sigma[*t][both][*ds][*ph][DsD0][*m][*c][*a], *sigma2[*t][both][*ds][*ph][DsD0][*m][*c][*a]  ));
              avg_sigma[*t][both][*ds][*ph][DsPhi][*m][both][both]     = new RooFormulaVar(Form("avg_sigma_DsPhi_%s", cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsPhi].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsPhi][*m],Fixed_CB_fraction[DsPhi][*m] ), RooArgList(*sigma[*t][both][*ds][*ph][DsPhi][*m][*c][*a],*sigma2[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));
              avg_sigma[*t][both][*ds][*ph][DsPhiSide][*m][both][both] = avg_sigma[*t][both][*ds][*ph][DsPhi][*m][both][both];
              
              cat_name = Form("%s_%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());

              sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   = new RooRealVar(   Form("sigma_DsD0_%s", cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
              //((RooRealVar*)sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
              sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma_DsPhi_%s",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a])); 
              sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
              //sig_ratio[*t][Helbin1][*ds][*ph][*b][*m]         = new RooRealVar(Form("Sigma_ratio_%s", cat_name.c_str()), Form("%s %s Sigma ratio ",mod[*m].c_str(),Bmod[*b].c_str()), Fixed_CB_Sigma_ratio[*b][*m]);
      
              sigma2[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   = new RooFormulaVar(Form("sigma2_DsD0_%s", cat_name.c_str()), Form("Sigma 2 DsD0 %s",(*m).c_str()),  Form("@0/%f",Fixed_CB_Sigma_ratio[DsD0][*m] ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   ));
              sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma2_DsPhi_%s",cat_name.c_str()), Form("Sigma 2 DsPhi %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DsPhi][*m]), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  ));
              sigma2[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
              
              avg_sigma[*t][Helbin1][*ds][*ph][DsD0][*m][both][both]      = new RooFormulaVar(Form("avg_sigma_DsD0_%s",  cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsD0].c_str(), mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsD0][*m], Fixed_CB_fraction[DsD0][*m]  ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a], *sigma2[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]  ));
              avg_sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][both][both]     = new RooFormulaVar(Form("avg_sigma_DsPhi_%s", cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsPhi].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsPhi][*m],Fixed_CB_fraction[DsPhi][*m] ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a] ));
              avg_sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][both][both] = avg_sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][both][both];
            

              sigma[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]              = sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a];
              sigma[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]             = sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
              sigma[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]         = sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a];

              sigma2[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]             = sigma2[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a];        
              sigma2[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]            = sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];         
              sigma2[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]        = sigma2[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a];      
              
              avg_sigma[*t][Helbin2][*ds][*ph][DsD0][*m][both][both]      = avg_sigma[*t][Helbin1][*ds][*ph][DsD0][*m][both][both]; 
              avg_sigma[*t][Helbin2][*ds][*ph][DsPhi][*m][both][both]     = avg_sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][both][both]; 
              avg_sigma[*t][Helbin2][*ds][*ph][DsPhiSide][*m][both][both] = avg_sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][both][both];  
 
              //sigma[*t][*h][*ds][*ph][*b][*m][both][both]  = new RooRealVar(Form("sigma_%s",         cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),10,  5,  20); //change D0 to D^0
              //sigma[*t][*h][*ds][*ph][*b][*m][plus][up]    = new RooRealVar(Form("sigma_%s_plus_up", cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),10,  5,  20);
              //sigma[*t][*h][*ds][*ph][*b][*m][plus][dn]    = new RooRealVar(Form("sigma_%s_plus_dn", cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),10,  5,  20);
              //sigma[*t][*h][*ds][*ph][*b][*m][minus][up]   = new RooRealVar(Form("sigma_%s_minus_up",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),10,  5,  20);//sigma[*m][minus][up] = sigma[*m][plus][dn];
              //sigma[*t][*h][*ds][*ph][*b][*m][minus][dn]   = new RooRealVar(Form("sigma_%s_minus_dn",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),10,  5,  20);//sigma[*m][minus][dn] = sigma[*m][plus][up];
              /*
              sigma[*t][*h][*ds][*ph][*b][*m][plus][both]  = sigma[*t][*h][*ds][*ph][*b][*m][both][both];
              sigma[*t][*h][*ds][*ph][*b][*m][minus][both] = sigma[*t][*h][*ds][*ph][*b][*m][both][both];             
              sigma[*t][*h][*ds][*ph][*b][*m][both][up]    = sigma[*t][*h][*ds][*ph][*b][*m][both][both]; 
              sigma[*t][*h][*ds][*ph][*b][*m][both][dn]    = sigma[*t][*h][*ds][*ph][*b][*m][both][both];
              
              sigma[*t][*h][*ds][*ph][DsPhi][*m][both][dn] = sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
              sigma[*t][*h][*ds][*ph][DsPhi][*m][both][up] = sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
              
      
              sigma2[*t][*h][*ds][*ph][*b][*m][both][both]  = new RooFormulaVar(Form("sigma2_%s",         cat_name.c_str()), Form("Sigma 2 %s %s", (*b).c_str(),(*m).c_str()), "@1/@0", RooArgList(*sig_ratio[*t][*h][*ds][*ph][*b][*m],  *sigma[*t][*h][*ds][*ph][*b][*m][both][both]   ));
              sigma2[*t][*h][*ds][*ph][*b][*m][plus][up]    = new RooFormulaVar(Form("sigma2_%s_plus_up", cat_name.c_str()) ,Form("Sigma 2 %s %s", (*b).c_str(),(*m).c_str()), "@1/@0", RooArgList(*sig_ratio[*t][*h][*ds][*ph][*b][*m],  *sigma[*t][*h][*ds][*ph][*b][*m][plus][up]     ));
              sigma2[*t][*h][*ds][*ph][*b][*m][plus][dn]    = new RooFormulaVar(Form("sigma2_%s_plus_dn", cat_name.c_str()) ,Form("Sigma 2 %s %s", (*b).c_str(),(*m).c_str()), "@1/@0", RooArgList(*sig_ratio[*t][*h][*ds][*ph][*b][*m],  *sigma[*t][*h][*ds][*ph][*b][*m][plus][dn]     ));
              sigma2[*t][*h][*ds][*ph][*b][*m][minus][up]   = new RooFormulaVar(Form("sigma2_%s_minus_up",cat_name.c_str()) ,Form("Sigma 2 %s %s", (*b).c_str(),(*m).c_str()), "@1/@0", RooArgList(*sig_ratio[*t][*h][*ds][*ph][*b][*m],  *sigma[*t][*h][*ds][*ph][*b][*m][minus][up]    ));//sigma2[*m][minus][up] = sigma2[*m][plus][dn];
              sigma2[*t][*h][*ds][*ph][*b][*m][minus][dn]   = new RooFormulaVar(Form("sigma2_%s_minus_dn",cat_name.c_str()) ,Form("Sigma 2 %s %s", (*b).c_str(),(*m).c_str()), "@1/@0", RooArgList(*sig_ratio[*t][*h][*ds][*ph][*b][*m],  *sigma[*t][*h][*ds][*ph][*b][*m][minus][dn]    ));//sigma2[*m][minus][dn] = sigma2[*m][plus][up];

              sigma2[*t][*h][*ds][*ph][*b][*m][plus][both]  = sigma2[*t][*h][*ds][*ph][*b][*m][both][both];
              sigma2[*t][*h][*ds][*ph][*b][*m][minus][both] = sigma2[*t][*h][*ds][*ph][*b][*m][both][both];             
              sigma2[*t][*h][*ds][*ph][*b][*m][both][up]    = sigma2[*t][*h][*ds][*ph][*b][*m][both][both]; 
              sigma2[*t][*h][*ds][*ph][*b][*m][both][dn]    = sigma2[*t][*h][*ds][*ph][*b][*m][both][both];
              
              sigma2[*t][*h][*ds][*ph][DsPhi][*m][both][dn] = sigma2[*t][*h][*ds][*ph][DsD0][*m][both][both];
              sigma2[*t][*h][*ds][*ph][DsPhi][*m][both][up] = sigma2[*t][*h][*ds][*ph][DsD0][*m][both][both];
              */
              //}
            }
          }
        }
      }
    }
  }
  // --------- Signal tail params ---------------
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str());
              
              ncb[*t][*h][*ds][*ph][*b][*m]         = new RooRealVar(Form("ncb_%s",        cat_name.c_str()), Form("%s %s  nCB",          Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_n[*b][*m]);
              alpha[*t][*h][*ds][*ph][*b][*m]       = new RooRealVar(Form("alphaL_%s",     cat_name.c_str()), Form("%s %s alpha",         Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_alpha[*b][*m]); 
            }
          }
        }
      }
    }
  }



  // --------- Comb. background slope ---------------
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
              std::string cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str());
              comb_slope[*t][*h][*ds][*ph][*b][*m] = global_comb_slope;
              //comb_slope[*t][*h][*ds][*ph][*b][*m] = new RooRealVar(Form("comb_slope_%s",cat_name.c_str())   ,Form("%s comb. slope", mod[*b].c_str()), -0.00112, -1.0, -0.0000000001);
            }
          }
        }
      }
    }
  }


  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
            std::string cat_name = Form("%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*m).c_str());
            
            // ------------------------------------------
            // PartReco shapes for DsD0 
            // ------------------------------------------
            
            // Ds*D0 pi0
            HORNS_a[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS_Ds*D0_a_%s",     cat_name.c_str()), "HORNS a",      5051.3  );
            HORNS_b[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS_Ds*D0_b_%s",     cat_name.c_str()), "HORNS b",      5132.9  );
            //HORNS_csi[*t][*h][*ds][*ph][DsD0][*m]    = new RooRealVar(Form("HORNS_Ds*D0_csi_%s",   cat_name.c_str()), "HORNS csi",    0.72864 );
            //HORNS_sigma[*t][*h][*ds][*ph][DsD0][*m]  = new RooRealVar(Form("HORNS_Ds*D0_sigma_%s", cat_name.c_str()), "HORNS sigma",  8.6610  );
            HORNS_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HORNS_R[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS_Ds*D0_R_%s",     cat_name.c_str()), "HORNS R",      1.0    );
            HORNS_f[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS_Ds*D0_f_%s",     cat_name.c_str()), "HORNS f",      1.0     );
            HORNS_csi[*t][*h][*ds][*ph][DsD0][*m]    = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsD0][*m]  = global_shift;
            
            // Ds*D0 gamma
            HILL_a[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL_Ds*D0_a_%s",      cat_name.c_str()), "HILL a",      4976.730   );
            HILL_b[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL_Ds*D0_b_%s",      cat_name.c_str()), "HILL b",      5213.055   );
            //HILL_csi[*t][*h][*ds][*ph][DsD0][*m]     = new RooRealVar(Form("HILL_Ds*D0_csi_%s",    cat_name.c_str()), "HILL csi",    0.80687  );
            //HILL_sigma[*t][*h][*ds][*ph][DsD0][*m]   = new RooRealVar(Form("HILL_Ds*D0_sigma_%s",  cat_name.c_str()), "HILL sigma",  8         );
            HILL_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HILL_R[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL_Ds*D0_R_%s",      cat_name.c_str()), "HILL R",      1.0       );
            HILL_f[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL_Ds*D0_f_%s",      cat_name.c_str()), "HILL f",      1.0       );              // PartReco shapes for DsD0 
            HILL_csi[*t][*h][*ds][*ph][DsD0][*m]     = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsD0][*m]   = global_shift;
            
            // DsD*0 pi0
            HORNS2_a[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS2_DsD*0_a_%s",     cat_name.c_str()), "HORNS2 a",      5051.478  );
            HORNS2_b[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS2_DsD*0_b_%s",     cat_name.c_str()), "HORNS2 b",      5128.577  );
            //HORNS2_sigma[*t][*h][*ds][*ph][DsD0][*m]  = new RooRealVar(Form("HORNS2_DsD*0_sigma_%s", cat_name.c_str()), "HORNS2 sigma",  8.6610  );
            HORNS2_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HORNS2_R[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS2_DsD*0_R_%s",     cat_name.c_str()), "HORNS2 R",      1.0    );
            HORNS2_f[*t][*h][*ds][*ph][DsD0][*m]      = new RooRealVar(Form("HORNS2_DsD*0_f_%s",     cat_name.c_str()), "HORNS2 f",      1.0     );
            HORNS2_csi[*t][*h][*ds][*ph][DsD0][*m]    = global_csi;
            HORNS2_shift[*t][*h][*ds][*ph][DsD0][*m]  = global_shift;
            
            // DsD*0 gamma
            HILL2_a[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL2_DsD*0_a_%s",      cat_name.c_str()), "HILL2 a",      4970.140   );
            HILL2_b[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL2_DsD*0_b_%s",      cat_name.c_str()), "HILL2 b",      5216.106   );
            //HILL2_sigma[*t][*h][*ds][*ph][DsD0][*m]   = new RooRealVar(Form("HILL2_DsD*0_sigma_%s",  cat_name.c_str()), "HILL2 sigma",  8         );
            HILL2_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HILL2_R[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL2_DsD*0_R_%s",      cat_name.c_str()), "HILL2 R",      1.0       );
            HILL2_f[*t][*h][*ds][*ph][DsD0][*m]       = new RooRealVar(Form("HILL2_DsD*0_f_%s",      cat_name.c_str()), "HILL2 f",      1.0       );  
            HILL2_csi[*t][*h][*ds][*ph][DsD0][*m]     = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsD0][*m]   = global_shift;            // PartReco shapes for DsD0 

            // ------------------------------------------
            // PartReco shapes for Ds*Phi
            // ------------------------------------------

            // missed pi0 components 
            // pi0 010
            HORNS_a[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsPhi_a_%s",     cat_name.c_str()), "HORNS a",      5026.751  );
            HORNS_b[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsPhi_b_%s",     cat_name.c_str()), "HORNS b",      5124.790  );
            //HORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = new RooRealVar(Form("HORNS_DsPhi_csi_%s",   cat_name.c_str()), "HORNS csi",    0.72864 );
            //HORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = new RooRealVar(Form("HORNS_DsPhi_sigma_%s", cat_name.c_str()), "HORNS sigma",  8.6610  );
            HORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both];
            HORNS_R[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsPhi_R_%s",     cat_name.c_str()), "HORNS R",      1.0    );
            HORNS_f[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsPhi_f_%s",     cat_name.c_str()), "HORNS f",      1.0     );
            HORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;
            
            //pi0 101
            HILL2_a[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HILL2_DsPhi_a_%s",      cat_name.c_str()), "HILL2 a",      5026.751   );
            HILL2_b[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HILL2_DsPhi_b_%s",      cat_name.c_str()), "HILL2 b",      5124.790   );
            //HILL2_csi[*t][*h][*ds][*ph][DsPhi][*m]    = new RooRealVar(Form("HILL2_DsPhi_csi_%s",    cat_name.c_str()), "HILL2 csi",    0.80687  );
            //HILL2_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = new RooRealVar(Form("HILL2_DsPhi_sigma_%s",  cat_name.c_str()), "HILL2 sigma",  8         );
            HILL2_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both];
            HILL2_R[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HILL2_DsPhi_R_%s",      cat_name.c_str()), "HILL2 R",      1.0       );
            HILL2_f[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HILL2_DsPhi_f_%s",      cat_name.c_str()), "HILL2 f",      1.0       );
            HILL2_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;
            
            // Missed gamma components 
            //gamma 010
            HILL_a[*t][*h][*ds][*ph][DsPhi][*m]       = new RooRealVar(Form("HILL_DsPhi_a_%s",      cat_name.c_str()), "HILL a",      4936.387   );
            HILL_b[*t][*h][*ds][*ph][DsPhi][*m]       = new RooRealVar(Form("HILL_DsPhi_b_%s",      cat_name.c_str()), "HILL b",      5220.646   );
            //HILL_csi[*t][*h][*ds][*ph][DsPhi][*m]     = new RooRealVar(Form("HILL_DsPhi_csi_%s",    cat_name.c_str()), "HILL csi",    0.80687  );]  
            //HILL_sigma[*t][*h][*ds][*ph][DsPhi][*m]   = new RooRealVar(Form("HILL_DsPhi_sigma_%s",  cat_name.c_str()), "HILL_sigma",  8         );  
            HILL_sigma[*t][*h][*ds][*ph][DsPhi][*m]   = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both]; 
            HILL_R[*t][*h][*ds][*ph][DsPhi][*m]       = new RooRealVar(Form("HILL_DsPhi_R_%s",      cat_name.c_str()), "HILL R",      1.0       );
            HILL_f[*t][*h][*ds][*ph][DsPhi][*m]       = new RooRealVar(Form("HILL_DsPhi_f_%s",      cat_name.c_str()), "HILL f",      1.0       );
            HILL_csi[*t][*h][*ds][*ph][DsPhi][*m]     = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsPhi][*m]   = global_shift;

            //gamma 101
            LITTLEHORNS_a[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhi_a_%s",     cat_name.c_str()), "LITTLEHORNS a",      4936.387  );
            LITTLEHORNS_b[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhi_b_%s",     cat_name.c_str()), "LITTLEHORNS b",      5220.646  );
            //LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = new RooRealVar(Form("LITTLEHORNS_DsPhi_csi_%s",   cat_name.c_str()), "LITTLEHORNS csi",    0.72864 );
            //LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = new RooRealVar(Form("LITTLEHORNS_DsPhi_sigma_%s", cat_name.c_str()), "LITTLEHORNS sigma",  8.6610  );
            LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both]; 
            LITTLEHORNS_R[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhi_R_%s",     cat_name.c_str()), "LITTLEHORNS R",      1.0    );
            LITTLEHORNS_f[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhi_f_%s",     cat_name.c_str()), "LITTLEHORNS f",      1.0     );
            LITTLEHORNS_g[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhi_g_%s",     cat_name.c_str()), "LITTLEHORNS g",      0.0     );
            LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi;
            LITTLEHORNS_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;

            // ------------------------------------------
            // PartReco shapes for Ds*Phi Sideband
            // ------------------------------------------

            // missed pi0 components 
            // pi0 010
            HORNS_a[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsPhiSide_a_%s",     cat_name.c_str()), "HORNS a",      5026.751  );
            HORNS_b[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsPhiSide_b_%s",     cat_name.c_str()), "HORNS b",      5124.790  );
            //HORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = new RooRealVar(Form("HORNS_DsPhiSide_csi_%s",   cat_name.c_str()), "HORNS csi",    0.72864 );
            //HORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = new RooRealVar(Form("HORNS_DsPhiSide_sigma_%s", cat_name.c_str()), "HORNS sigma",  8.6610  );
            HORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both];
            HORNS_R[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsPhiSide_R_%s",     cat_name.c_str()), "HORNS R",      1.0    );
            HORNS_f[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsPhiSide_f_%s",     cat_name.c_str()), "HORNS f",      1.0     );
            HORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsPhiSide][*m]  = global_shift;
            
            //pi0 101
            HILL2_a[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HILL2_DsPhiSide_a_%s",      cat_name.c_str()), "HILL2 a",      5026.751   );
            HILL2_b[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HILL2_DsPhiSide_b_%s",      cat_name.c_str()), "HILL2 b",      5124.790   );
            //HILL2_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = new RooRealVar(Form("HILL2_DsPhiSide_csi_%s",    cat_name.c_str()), "HILL2 csi",    0.80687  );
            //HILL2_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = new RooRealVar(Form("HILL2_DsPhiSide_sigma_%s",  cat_name.c_str()), "HILL2 sigma",  8         );
            HILL2_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both];
            HILL2_R[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HILL2_DsPhiSide_R_%s",      cat_name.c_str()), "HILL2 R",      1.0       );
            HILL2_f[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HILL2_DsPhiSide_f_%s",      cat_name.c_str()), "HILL2 f",      1.0       );
            HILL2_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsPhiSide][*m]  = global_shift;
            
            // Missed gamma components 
            //gamma 010
            HILL_a[*t][*h][*ds][*ph][DsPhiSide][*m]       = new RooRealVar(Form("HILL_DsPhiSide_a_%s",      cat_name.c_str()), "HILL a",      4936.387   );
            HILL_b[*t][*h][*ds][*ph][DsPhiSide][*m]       = new RooRealVar(Form("HILL_DsPhiSide_b_%s",      cat_name.c_str()), "HILL b",      5220.646   );
            //HILL_csi[*t][*h][*ds][*ph][DsPhiSide][*m]     = new RooRealVar(Form("HILL_DsPhiSide_csi_%s",    cat_name.c_str()), "HILL csi",    0.80687  );]  
            //HILL_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]   = new RooRealVar(Form("HILL_DsPhiSide_sigma_%s",  cat_name.c_str()), "HILL_sigma",  8         );  
            HILL_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]   = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both]; 
            HILL_R[*t][*h][*ds][*ph][DsPhiSide][*m]       = new RooRealVar(Form("HILL_DsPhiSide_R_%s",      cat_name.c_str()), "HILL R",      1.0       );
            HILL_f[*t][*h][*ds][*ph][DsPhiSide][*m]       = new RooRealVar(Form("HILL_DsPhiSide_f_%s",      cat_name.c_str()), "HILL f",      1.0       );
            HILL_csi[*t][*h][*ds][*ph][DsPhiSide][*m]     = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsPhiSide][*m]   = global_shift;

            //gamma 101
            LITTLEHORNS_a[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_a_%s",     cat_name.c_str()), "LITTLEHORNS a",      4936.387  );
            LITTLEHORNS_b[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_b_%s",     cat_name.c_str()), "LITTLEHORNS b",      5220.646  );
            //LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_csi_%s",   cat_name.c_str()), "LITTLEHORNS csi",    0.72864 );
            //LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_sigma_%s", cat_name.c_str()), "LITTLEHORNS sigma",  8.6610  );
            LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both]; 
            LITTLEHORNS_R[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_R_%s",     cat_name.c_str()), "LITTLEHORNS R",      1.0    );
            LITTLEHORNS_f[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_f_%s",     cat_name.c_str()), "LITTLEHORNS f",      1.0     );
            LITTLEHORNS_g[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("LITTLEHORNS_DsPhiSide_g_%s",     cat_name.c_str()), "LITTLEHORNS g",      0.0     );
            LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = global_csi;
            LITTLEHORNS_shift[*t][*h][*ds][*ph][DsPhiSide][*m]  = global_shift;

            //HORNS_a[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_a_%s",     cat_name.c_str()), "HORNS a",      5051.8  );
            //HORNS_b[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_b_%s",     cat_name.c_str()), "HORNS b",      5124.5  );
            //HORNS_csi[*t][*h][*ds][*ph][*b][*m]    = new RooRealVar(Form("HORNS_csi_%s",   cat_name.c_str()), "HORNS csi",    0.72864 );
            //HORNS_shift[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_shift_%s",cat_name.c_str()), "HORNS shift", 0.59647);
            //HORNS_shift[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_shift_%s",cat_name.c_str()), "HORNS shift",   0.529,    -10,    10);
            //HORNS_sigma[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HORNS_sigma_%s",cat_name.c_str()), "HORNS sigma",   8.6610  );
            //HORNS_R[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_R_%s",    cat_name.c_str()), "HORNS R",       15.5    );
            //HORNS_f[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HORNS_f_%s",    cat_name.c_str()), "HORNS f",       0.98133 );
            
            //HORNS_a[*m]      = new RooRealVar(Form("HORNS_a_%s",(*m).c_str()),      "HORNS a",     5051.8, 5000,  5100);
            //HORNS_b[*m]      = new RooRealVar(Form("HORNS_b_%s",(*m).c_str()),      "HORNS b",     5124.5, 5100,  5180);
            //HORNS_csi[*m]    = new RooRealVar(Form("HORNS_csi_%s",(*m).c_str()),    "HORNS csi",   0.7286,    0,     1);
            //HORNS_shift[*m]  = new RooRealVar(Form("HORNS_shift_%s",(*m).c_str()),  "HORNS shift",  0.529,    0,    20);
            //HORNS_sigma[*m]  = new RooRealVar(Form("HORNS_sigma_%s",(*m).c_str()),  "HORNS sigma", 8.0954,    0,    30);
            //HORNS_R[*m]      = new RooRealVar(Form("HORNS_R_%s",(*m).c_str()),      "HORNS R",       0.68,    0,    50);
            //HORNS_f[*m]      = new RooRealVar(Form("HORNS_f_%s",(*m).c_str()),      "HORNS f",       0.98,    0,     1);

            //HILL_a[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_a_%s",    cat_name.c_str()), "HILL a",      4870.0    );
            //HILL_b[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_b_%s",    cat_name.c_str()), "HILL b",      5231.9    );
            //HILL_csi[*t][*h][*ds][*ph][*b][*m]    = new RooRealVar(Form("HILL_csi_%s",  cat_name.c_str()), "HILL csi",       0.80687);
            //HILL_shift[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_shift_%s",cat_name.c_str()), "HILL shift",  -13.620);
            //HILL_shift[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_shift_%s",cat_name.c_str()), "HILL shift",   -14.2,   -100,   100);
            //HILL_sigma[*t][*h][*ds][*ph][*b][*m]  = new RooRealVar(Form("HILL_sigma_%s",cat_name.c_str()), "HILL sigma",     8      );
            //HILL_R[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_R_%s",    cat_name.c_str()), "HILL R",         0.5    );
            //HILL_f[*t][*h][*ds][*ph][*b][*m]      = new RooRealVar(Form("HILL_f_%s",    cat_name.c_str()), "HILL f",         0.8821 );

            //HILL_a[*m]      = new RooRealVar(Form("HILL_a_%s",(*m).c_str()),     "HILL a",      4779.5,   4700,  4900);
            //HILL_b[*m]      = new RooRealVar(Form("HILL_b_%s",(*m).c_str()),     "HILL b",      5231.9,   5127,  5327);
            //HILL_csi[*m]    = new RooRealVar(Form("HILL_csi_%s",(*m).c_str()),   "HILL csi",     0.807,    -10,    10);
            //HILL_shift[*m]  = new RooRealVar(Form("HILL_shift_%s",(*m).c_str()), "HILL shift",   -14.2,   -100,   100);
            //HILL_sigma[*m]  = new RooRealVar(Form("HILL_sigma_%s",(*m).c_str()), "HILL sigma",   0.918,      0,    20);
            //HILL_R[*m]      = new RooRealVar(Form("HILL_R_%s",(*m).c_str()),     "HILL R",        11.3,      0,    20);
            //HILL_f[*m]      = new RooRealVar(Form("HILL_f_%s",(*m).c_str()),     "HILL f",       0.882,      0,     1);
   
            /*
            //HILL_a[*m]      = new RooRealVar(Form("HILL_a_%s",(*m).c_str()),     "HILL a",      4779.5,   4700,  4900);
            //HILL_b[*m]      = new RooRealVar(Form("HILL_b_%s",(*m).c_str()),     "HILL b",      5231.9,   5127,  5327);
            //HILL_csi[*m]    = new RooRealVar(Form("HILL_csi_%s",(*m).c_str()),   "HILL csi",     0.807,    -10,    10);
            //HILL_shift[*m]  = new RooRealVar(Form("HILL_shift_%s",(*m).c_str()), "HILL shift",   -14.2,   -100,   100);
            //HILL_sigma[*m]  = new RooRealVar(Form("HILL_sigma_%s",(*m).c_str()), "HILL sigma",   0.918,      0,    20);
            //HILL_R[*m]      = new RooRealVar(Form("HILL_R_%s",(*m).c_str()),     "HILL R",        11.3,      0,    20);
            //HILL_f[*m]      = new RooRealVar(Form("HILL_f_%s",(*m).c_str()),     "HILL f",       0.882,      0,     1);
            */
          } 
        }
      }
    }
  }
  


  // ------------------------------------------------------
  // --------- Fix RooRealVars to each other  -------------
  // ------------------------------------------------------
 
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){           

              // Fix comb slope to be the same 
              comb_slope[*t][*h][*ds][*ph][*b][*m] = comb_slope[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];

              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                  
                  // Use same mean value for all plots 
                  mean_B[*t][*h][*ds][*ph][*b][*m][*c][*a] = mean_B[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()]; 

                  // Allow sigma to vary between Ds decay modes
                  //sigma[*t][*h][*ds][*ph][*b][*m][*c][*a] = sigma[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m][*chargeList.begin()][*magnetList.begin()]; 
                }
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
  RooAbsPdf *pdf_peak_cb1 = 0;
  RooAbsPdf *pdf_peak_cb2 = 0;
  RooAbsPdf *pdf_PartReco = 0;
  RooAbsPdf *pdf_HORNS = 0;
  RooAbsPdf *pdf_HORNS2 = 0;
  RooAbsPdf *pdf_LITTLEHORNS = 0;
  RooAbsPdf *pdf_HILL = 0;
  RooAbsPdf *pdf_HILL2 = 0;
  RooAbsPdf *pdf_comb = 0;
  RooAbsPdf *pdf_Dsa1 = 0;
  RooAbsPdf *pdf_DsstD0 = 0;
  RooAbsPdf *pdf_DsDst0 = 0;

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
          for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
            for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                  std::string tag=(*t)+underscore+(*h)+underscore+(*ds)+underscore+(*ph)+underscore+(*b)+underscore+(*m)+underscore+(*c)+underscore+(*a);
                  //if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Adding to sim pdf: "<< tag <<std::endl;
                  //Both b modes

                  // Double CB shape from MC
                  RooAbsReal* mub     = mean_B[*t][*h][*ds][*ph][*b][*m][*c][*a];
                  RooAbsReal* sig_CB  = sigma[*t][*h][*ds][*ph][*b][*m][*c][*a];
                  RooAbsReal* n_CB    = ncb[*t][*h][*ds][*ph][*b][*m];
                  RooAbsReal* alp_CB  = alpha[*t][*h][*ds][*ph][*b][*m];
                  pdf_peak_cb1  = new RooCBShape( Form("pdf_peak_cb1_%s",tag.c_str()), "", *mB, *mub, *sig_CB,  *alp_CB, *n_CB );

                  RooAbsReal* sig2_CB = sigma2[*t][*h][*ds][*ph][*b][*m][*c][*a];
                  pdf_peak_cb2  = new RooCBShape( Form("pdf_peak_cb2_%s",tag.c_str()), "", *mB, *mub, *sig2_CB, *alp_CB, *n_CB );
                  
                  pdf_peak      = new RooAddPdf(  Form("pdf_peak_%s",    tag.c_str()), "", RooArgSet(*pdf_peak_cb1,*pdf_peak_cb2),*sig_frac[*t][*h][*ds][*ph][*b][*m]);
                  
                  //comb bkg
                  RooRealVar* comb = comb_slope[*t][*h][*ds][*ph][*b][*m];
                  pdf_comb = new RooExponential( Form("pdf_comb_%s",tag.c_str()), "",*mB,*comb);
  
                  if( (*b).c_str()==DsD0){
                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Adding DsD0 mode... " <<std::endl;
                    // Ds*D0 pi0
                    RooAbsReal* HORNS_var_a     = HORNS_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_b     = HORNS_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_csi   = HORNS_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_shift = HORNS_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_sigma = HORNS_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_R     = HORNS_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_f     = HORNS_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HORNS = new RooHORNSdini(Form("pdf_HORNS_%s",tag.c_str()), "", *mB, *HORNS_var_a, *HORNS_var_b, *HORNS_var_csi, *HORNS_var_shift, *HORNS_var_sigma, *HORNS_var_R, *HORNS_var_f);
                
                    // Ds*D0 gamma
                    RooAbsReal* HILL_var_a     = HILL_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_b     = HILL_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_csi   = HILL_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_shift = HILL_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_sigma = HILL_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_R     = HILL_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_f     = HILL_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HILL = new RooHILLdini(Form("pdf_HILL_%s",tag.c_str()), "", *mB, *HILL_var_a, *HILL_var_b, *HILL_var_csi, *HILL_var_shift, *HILL_var_sigma, *HILL_var_R, *HILL_var_f);
 
                    // Ds*D0 Shape
                    pdf_DsstD0 = new RooAddPdf(  Form("pdf_DsstD0_%s",tag.c_str()), "", RooArgSet(*pdf_HORNS,*pdf_HILL),*fraction_Dsst_Dspi0_DsD0_adj);
                  
                    //DsD*0 pi0
                    RooAbsReal* HORNS2_var_a     = HORNS2_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_b     = HORNS2_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_csi   = HORNS2_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_shift = HORNS2_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_sigma = HORNS2_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_R     = HORNS2_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS2_var_f     = HORNS2_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HORNS2 = new RooHORNSdini(Form("pdf_HORNS2_%s",tag.c_str()), "", *mB, *HORNS2_var_a, *HORNS2_var_b, *HORNS2_var_csi, *HORNS2_var_shift, *HORNS2_var_sigma, *HORNS2_var_R, *HORNS2_var_f);
                
                    //DsD*0 gamma
                    RooAbsReal* HILL2_var_a     = HILL2_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_b     = HILL2_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_csi   = HILL2_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_shift = HILL2_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_sigma = HILL2_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_R     = HILL2_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_f     = HILL2_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HILL2 = new RooHILLdini(Form("pdf_HILL2_%s",tag.c_str()), "", *mB, *HILL2_var_a, *HILL2_var_b, *HILL2_var_csi, *HILL2_var_shift, *HILL2_var_sigma, *HILL2_var_R, *HILL2_var_f);
                    
                    // Ds*D0 Shape
                    pdf_DsDst0 = new RooAddPdf(  Form("pdf_DsDst0_%s",tag.c_str()), "", RooArgSet(*pdf_HORNS2,*pdf_HILL2),*fraction_Dst0_D0pi0_DsD0_adj);

                    RooArgSet pdflist(  *pdf_comb
                                        ,*pdf_peak 
                                        ,*pdf_DsstD0 
                                        ,*pdf_DsDst0 
                                        );
                    RooArgSet nevents(  *yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*B_yield[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_dsstd0[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_dsdst0[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ); 

                    //Add to master PDF
                    RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_%s",tag.c_str())       ,"",pdflist, nevents);
                    
                    std::stringstream str;
                    str<<(*t)<<underscore<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                    
                    sim->addPdf(*masterPdf,str.str().c_str());
                    //sim->ls();
                    //sim->Print();
                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Added to sim pdf: "<< str.str().c_str() <<std::endl;
                    
                  } else if((*b).c_str()==DsPhi || (*b).c_str()==DsPhiSide){
                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Adding DsPhi mode... " <<std::endl; 
                    
                    //RooHORNSdini
                    RooAbsReal* HORNS_var_a     = HORNS_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_b     = HORNS_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_csi   = HORNS_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_shift = HORNS_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_sigma = HORNS_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_R     = HORNS_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_var_f     = HORNS_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HORNS = new RooHORNSdini(Form("pdf_HORNS_%s",tag.c_str()), "", *mB, *HORNS_var_a, *HORNS_var_b, *HORNS_var_csi, *HORNS_var_shift, *HORNS_var_sigma, *HORNS_var_R, *HORNS_var_f);
                
                    //RooHILLdini
                    RooAbsReal* HILL_var_a     = HILL_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_b     = HILL_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_csi   = HILL_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_shift = HILL_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_sigma = HILL_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_R     = HILL_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_var_f     = HILL_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_HILL = new RooHILLdini(Form("pdf_HILL_%s",tag.c_str()), "", *mB, *HILL_var_a, *HILL_var_b, *HILL_var_csi, *HILL_var_shift, *HILL_var_sigma, *HILL_var_R, *HILL_var_f);
 
                    // RooHILLdini
                    RooAbsReal* HILL2_var_a     = HILL2_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_b     = HILL2_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_csi   = HILL2_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_shift = HILL2_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_sigma = HILL2_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_R     = HILL2_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL2_var_f     = HILL2_f[*t][*h][*ds][*ph][*b][*m];
                    
                    pdf_HILL2 = new RooHILLdini(Form("pdf_HILL2_%s",tag.c_str()), "", *mB, *HILL2_var_a, *HILL2_var_b, *HILL2_var_csi, *HILL2_var_shift, *HILL2_var_sigma, *HILL2_var_R, *HILL2_var_f);
                    
                    //RooLITTLEHORNSdini
                    RooAbsReal* LITTLEHORNS_var_a     = LITTLEHORNS_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_b     = LITTLEHORNS_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_csi   = LITTLEHORNS_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_shift = LITTLEHORNS_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_sigma = LITTLEHORNS_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_R     = LITTLEHORNS_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_f     = LITTLEHORNS_f[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* LITTLEHORNS_var_g     = LITTLEHORNS_g[*t][*h][*ds][*ph][*b][*m];
                    pdf_LITTLEHORNS = new RooLITTLEHORNSdini(Form("pdf_LITTLEHORNS_%s",tag.c_str()), "", *mB, *LITTLEHORNS_var_a, *LITTLEHORNS_var_b, *LITTLEHORNS_var_csi, *LITTLEHORNS_var_shift, *LITTLEHORNS_var_sigma, *LITTLEHORNS_var_R, *LITTLEHORNS_var_f, *LITTLEHORNS_var_g);

                    std::cout << "Making PartReco RooAddPdf" << std::endl;
                    pdf_PartReco = new RooAddPdf(Form("pdf_PartReco_%s", tag.c_str()),
                                                 "",
                                                 RooArgSet(*pdf_HORNS,
                                                           *pdf_HILL,
                                                           *pdf_HILL2,
                                                           *pdf_LITTLEHORNS), 
                                                 RooArgSet(*frac_HORNS[*t][*h][*ds][*ph][DsPhi][*m][*c][*a],
                                                           *frac_HILL[*t][*h][*ds][*ph][DsPhi][*m][*c][*a],
                                                           *frac_HILL2[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]) );
                    
                    std::cout << "Made PartReco RooAddPdf" << std::endl;

                    //Partially reconstructed Bs0 -> Ds (a1-> K Kst0)
                    pdf_Dsa1     = new RooKeysPdf(*PartRecoPDF_Bs02Dsa1, Form("pdf_Dsa1_%s",tag.c_str()) );
                    //RooFFTConvPdf* temp = PartRecoPDF_Bs02Dsa1_Convolved[*m];f
                    //pdf_Dsa1     = new RooFFTConvPdf(*(PartRecoPDF_Bs02Dsa1_Convolved[*m]), Form("pdf_Dsa1_%s",tag.c_str()) );
                    //Partially reconstructed Bs0 -> Ds* K Kst0
                    //pdf_DsstKKst = new RooKeysPdf(*PartRecoPDF_Bs02DsstKKst, Form("pdf_DsstKKst_%s",tag.c_str()) );
                    std::cout << "Made RooKeysPdf" << std::endl;

                    RooArgSet pdflist( *pdf_peak, 
                                        *pdf_comb,
                                        *pdf_PartReco, 
                                        *pdf_Dsa1 
                                       );
                                       // *pdf_DsstKKst );
                    std::cout << "Made pdf list" << std::endl;

                    RooArgSet nevents(  *B_yield[*t][*h][*ds][*ph][*b][*m][*c][*a], 
                                        *yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a],
                                        *PR_total_yield[*t][*h][*ds][*ph][*b][*m][*c][*a],
                                        *yield_Dsa1[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        );
                                        // *yield_DsstKKst[*t][*h][*ds][*ph][*b][*m][*c][*a] );
                    std::cout << "Making Master PDF" << std::endl;
                    RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_%s",tag.c_str()),"",pdflist, nevents);
                    
                    std::stringstream str;
                    str<<(*t)<<underscore<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                    
                    sim->addPdf(*masterPdf,str.str().c_str());
                    //sim->ls();
                    //sim->Print();
                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Added to sim pdf: "<< str.str().c_str() <<std::endl;
                    
                  } /*else if((*b).c_str()==DsPhiSide){

                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Adding DsPhiSide mode... " <<std::endl;

                    //Partially reconstructed Bs0 -> Ds (a1-> K Kst)
                    pdf_Dsa1 = new RooKeysPdf(*PartRecoPDF_Bs02Dsa1, Form("pdf_Dsa1_%s",tag.c_str()) );
                    //Partially reconstructed Bs0 -> Ds* K Kst0
                    //pdf_DsstKKst = new RooKeysPdf(*PartRecoPDF_Bs02DsstKKst, Form("pdf_DsstKKst_%s",tag.c_str()) );

                    RooArgSet pdflist( *pdf_peak, 
                                       *pdf_comb,  
                                       *pdf_Dsa1 );
                                       // *pdf_DsstKKst  );

                    RooArgSet nevents(  *B_yield[*t][*h][*ds][*ph][*b][*m][*c][*a], 
                                        *yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a], 
                                        *yield_Dsa1[*t][*h][*ds][*ph][*b][*m][*c][*a] );
                                        // *yield_DsstKKst[*t][*h][*ds][*ph][*b][*m][*c][*a] );

                    RooAddPdf* masterPdf = new RooAddPdf(Form("masterPdf_%s",tag.c_str())       ,"",pdflist, nevents);
                    
                    std::stringstream str;
                    str<<(*t)<<underscore<<(*h)<<underscore<<(*ds)<<underscore<<(*ph)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                    
                    sim->addPdf(*masterPdf,str.str().c_str());
                    //sim->ls();
                    //sim->Print();
                    if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Added to sim pdf: "<< str.str().c_str() <<std::endl; 
                  
                  }  //exit mode if*/
                } //close c
              } // close a
            } //close m
          } //close b
        } // close ph
      } // close ds 
    } // close  h
  } // close t
}//end of funcn DefineModel()

RooArgSet* DsPhiModel::GetParameters()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::GetParameters()"<<std::endl;
  RooArgSet* pars = sim->getParameters(RooArgSet(*mB,type,helBin,DsBDTBin,PhiBDTBin,Bmode,mode,charge,magnet));
  return pars;
}



void DsPhiModel::PrintResult()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::PrintResults()"<<std::endl;
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
          for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                if((*b)==DsPhi){
                  if (needsBlinding) std::cout << "DsPhi Mode is BLIND: Printing blinded yield" << std::endl;
                  float y_peak =  B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_dsa1 =  yield_Dsa1[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_partreco =   PR_total_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  
                  std::cout<< "\n---------------------------"<<std::endl;
                  std::cout<<(*b)<<", "<<(*m)<<", magnet: "<<(*a)<<", year: "<<(*t)<<"HelBin: " << (*h)<<std::endl;
                  std::cout<< " B Yield: "<<y_peak<<std::endl;
                  //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                  std::cout<< " Comb:    "<<y_comb<<std::endl;
                  std::cout<< " Dsa1:    "<<y_dsa1<<std::endl;
                  std::cout<< " Ds*Phi:  "<<y_partreco<<std::endl;
                  std::cout<< "---------------------------\n"<<std::endl;
                  //B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->printValue();
                } else if((*b)==DsPhiSide){
                  if (needsBlinding) std::cout << "DsPhi Mode is BLIND: Printing blinded yield" << std::endl;
                  float y_peak =  B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_dsa1 =  yield_Dsa1[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  
                  std::cout<< "\n---------------------------"<<std::endl;
                  std::cout<<(*b)<<", "<<(*m)<<", magnet: "<<(*a)<<", year: "<<(*t)<<", HelBin: " << (*h)<<std::endl;
                  std::cout<< " B Yield: "<<y_peak<<std::endl;
                  //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                  std::cout<< " Comb:    "<<y_comb<<std::endl;
                  std::cout<< " Dsa1:    "<<y_dsa1<<std::endl;

                  std::cout<< "---------------------------\n"<<std::endl;

                } else {
                  if(par->sumOverCharges){
                    float y_peak   =  B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_comb   =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_dsdst0 =  yield_dsdst0[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_dsstd0 =  yield_dsstd0[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    
                    std::cout<< "\n---------------------------"<<std::endl;
                    std::cout<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<", year-"<<(*t)<<" :"<<", HelBin: " << (*h)<<std::endl;
                    std::cout<< " B Yield: "<<y_peak<<std::endl;
                    //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                    std::cout<< " Comb:   "<<y_comb<<std::endl;
                    std::cout<< " DsD*0:  "<<y_dsdst0<<std::endl;
                    std::cout<< " Ds*D0:  "<<y_dsstd0<<std::endl;
                    std::cout<< "---------------------------\n"<<std::endl;
                  }else{
                    float y_peak_minus   =  yield_peak[*t][*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                    float y_comb_minus   =  yield_comb[*t][*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                    float y_dsstd0_minus =  yield_dsstd0[*t][*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                    float y_dsdst0_minus =  yield_dsdst0[*t][*h][*ds][*ph][*b][*m][minus][*a]->getVal();
                    float y_peak_plus    =  yield_peak[*t][*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                    float y_comb_plus    =  yield_comb[*t][*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                    float y_dsstd0_plus  =  yield_dsstd0[*t][*h][*ds][*ph][*b][*m][plus][*a]->getVal();
                    float y_dsdst0_plus  =  yield_dsdst0[*t][*h][*ds][*ph][*b][*m][plus][*a]->getVal();

                 
                    std::cout<<(*t)<<", "<<(*h)<<", "<<(*ds)<<", "<<(*ph)<<", "<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                    std::cout<<" Peak minus: "<<y_peak_minus<<" comb. minus: "<<y_comb_minus<<" XX^{*} minus: "<<y_dsstd0_minus<<" X^{*}X: minus: "<<y_dsdst0_minus<<std::endl;
                    std::cout<<" Peak plus:  "<<y_peak_plus <<" comb. plus:  "<<y_comb_plus <<" XX^{*} plus:  "<<y_dsstd0_plus <<" X^{*}X: plus:  "<<y_dsdst0_plus <<std::endl;
                    std::cout<<"--------------------------------------------------"<<std::endl;
                    std::cout<<" Peak total: "<<y_peak_minus+y_peak_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<" XX^{*} total: "<<y_dsstd0_minus+y_dsstd0_plus<<" X^{*}X: total: "<<y_dsdst0_minus+y_dsstd0_plus<<std::endl;

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
  }
 //closes if(not needsBlinding)
  std::cout<<std::endl;
}

std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > > > DsPhiModel::GetResult()
{
  if(par->debug) std::cout<<"Running: DsPhiModel::GetResults()"<<std::endl;
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
          for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
                  if (needsBlinding && ((*b)==DsPhi||(*b)==DsPhiSide) ) {
                    std::cout << "DsPhi mode is BLIND, therefore no values returned " << std::endl;
                  } else {


                    std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());
                    std::cout<< " Category: " << cat_str << std::endl;
                    fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Peak"]    = yield_peak[*t][*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                    //fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Peakerr"] = yield_peak[*t][*h][*ds][*ph][*b][*m][*c][*a]->getError();
                    fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Comb"]    = yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                    fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["XXst"]    = yield_dsdst0[*t][*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                    fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["XstX"]    = yield_dsstd0[*t][*h][*ds][*ph][*b][*m][*c][*a]->getVal();
                    std::cout<< " Peak:   " <<fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Peak"]<<std::endl;
                    std::cout<< " Comb:   " <<fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Comb"]<<std::endl;
                    std::cout<< " XXst:   " <<fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["XXst"]<<std::endl;
                    std::cout<< " XstX:   " <<fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["XstX"]<<std::endl;
                    
                    mB->setRange("signal", 5233.5, 5333.5) ;
                    //RooAbsPdf* sum = sim->getPdf(cat_str.c_str());
                    //RooAbsPdf* comb = sim->getPdf(Form("pdf_comb_%s",cat_str.c_str()));
                    //RooAbsPdf* comb = (RooAbsPdf*) sum->FindObject(Form("pdf_comb_%s",cat_str.c_str()));
                    //RooAbsReal* ibkg_sig = sum->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"),RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())));  
                    double combslope = comb_slope[*t][*h][*ds][*ph][*b][*m]->getVal();
                    RooRealVar* fixed_slope = new RooRealVar("fixed_slope","slope",  combslope);
                    RooAbsPdf* pdf_comb = new RooExponential( "Exp", "",*mB,*fixed_slope);
                    RooAbsReal* ibkg_sig = pdf_comb->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"));  
                    fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"] = (ibkg_sig->getVal()*yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a]->getVal());
                    

                    std::cout<< " Comb_Reduced:   " <<fit_results[*t][*h][*ds][*ph][*b][*m][*c][*a]["Comb_Reduced"]<<std::endl;
                  }
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
