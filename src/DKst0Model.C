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

#include "DKst0Model.h"
#include "Parameters.h"
//#include "RooDoubleCB.h"
#include "RooLITTLEHORNSdini.h"
#include "RooHORNSdini.h"
#include "RooHILLdini.h"
#include "TMatrixTSym.h"

DKst0Model::DKst0Model(Parameters* p, RooRealVar* pmB, bool nB, bool vA, bool gen)
: Base()
  //, misID(p)
  //, calibration(p)
  //, pid("pid","PID boolean")
, type("type","Type of data") 
, mode("mode","D_{s} decay mode")
, charge("charge","batchelor charge")
, magnet("magnet","magnet polarity")
, helBin("helBin","Helicity Angle bin")
, blindString("BlindLeadingTheBlind")
, needsBlinding(nB)
, varyAllowed(vA)
, genModel(gen)
  //, pidList()
, typeList()
, modeList()
, BmodeList()
, chargeList()
, magnetList()
, HelBinList()
  //, m_pidCut(-1000)
{
  par=p;
  mB=pmB;
  UInt_t seed = 20110602;
  if(par->useSeed) seed = (UInt_t) par->seed;
  std::cout <<  "====================== " << std::endl;
  std::cout <<  "Using Seed: " << seed << std::endl;
  std::cout <<  "====================== " << std::endl;
  RooRandom::randomGenerator()->SetSeed(seed); 
  //rand = new TRandom3(time(NULL));
  rand = new TRandom3(seed);

  if(needsBlinding){
    std::cout<<"BLIND model envoked"<<std::endl;
  }

  DefineRooCategories();
  DefineModel();
}


void DKst0Model::DefineRooCategories() //like D0h with "PIDcut" removed
{
  if(par->debug) std::cout<<"Running: DKst0Model::DefineRooCategories()"<<std::endl;
  if(par->modes[D2PiKPi])     modeList.push_back(D2PiKPi);
  if(par->modes[D2KKPi])      modeList.push_back(D2KKPi);
  if(par->modes[D2PiPiPi])    modeList.push_back(D2PiPiPi);
  
  if(par->Bmodes[DD0])        BmodeList.push_back(DD0);
  if(par->Bmodes[DKst0])      BmodeList.push_back(DKst0);
  if(par->Bmodes[DKst0Side])  BmodeList.push_back(DKst0Side);
 

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
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++)  {
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              std::stringstream str;
              str<<(*t)<<underscore<<(*h)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
              cat->defineType(str.str().c_str());
              std::cout <<  " --> Creating Category : " <<str.str().c_str() <<std::endl;
            }  
          }
        }
      }
    }
  }
}



void DKst0Model::DefineModel()
{
  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel()"<<std::endl;
  
  // ======================================================== 
  // =========== Switches for quick/significance ============
  // ======================================================== 
  
  // Make fit run super quick --> for testing
  bool allConst_exYield = false;
  
  // Fix yields to zero --> to calculate sigma wilks 
  
  bool NoBR = false;
  
  bool NoYield_all = false;  
  std::map<std::string,bool> NoYield; 
  NoYield[D2PiKPi]  = false;
  NoYield[D2KKPi]   = false;  
  NoYield[D2PiPiPi] = false; 

  if(NoYield_all){ 
    NoYield[D2PiKPi]  = true; 
    NoYield[D2KKPi]   = true; 
    NoYield[D2PiPiPi] = true; 
  }


  RooRealVar* global_shift      = new RooRealVar("global_shift",      "Global shift",      -1.162, -10.0,    10.0);
  RooRealVar* global_shift_test = new RooRealVar("global_shift_test", "Global shift test", -1.162, -10.0,    10.0);

  // ======================================================== 
  // =========== import DsKK*0 PDF       ====================
  // ======================================================== 
  
  std::string locationRooKeys = "pdfs/";
  
  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing Dsa1 shape"<<std::endl;
  
  TFile f1((locationRooKeys + "kest1_Bs2DKPiPi_4800_5900.root").c_str(), "read"); 
  RooWorkspace* workspace1 = (RooWorkspace*)f1.Get("workspace_Bs2DKPiPi");
  if (!workspace1) std::cout<<"Running: DKst0Model::DefineModel() --> Can't find Dsa1 workspace"<<std::endl;   
  RooKeysPdf* PartRecoPDF_Bs2DKPiPi     = (RooKeysPdf*)workspace1->pdf("kest1_Bs2DKPiPi");
  if(!PartRecoPDF_Bs2DKPiPi) std::cout<<"Running: DKst0Model::DefineModel() --> Can't find DKPiPi pdf"<<std::endl; 
  f1.Close();
  
  RooRealVar* sg_test       = new RooRealVar("sg_test",   "sg_test",    1  ) ;
  RooGaussian* gauss_test   = new RooGaussian("gauss_test",  "gauss_test",  *mB,*global_shift_test,*sg_test  ) ;
  RooFFTConvPdf* PartRecoPDF_Bs2DKPiPi_Conv   = new RooFFTConvPdf("kest1_Bs2DKPiPi_conv",  "kest1_Bs2DKPiPi_conv",  *mB,*PartRecoPDF_Bs2DKPiPi, *gauss_test  );
  

  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing DsstKKst shape"<<std::endl;
  TFile f2((locationRooKeys + "kest1_Bu2DstKPi_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace2 = (RooWorkspace*)f2.Get("workspace_Bu2DstKPi");
  RooKeysPdf* PartRecoPDF_Bu2DstKPi  = (RooKeysPdf*)workspace2->pdf("kest1_Bu2DstKPi");
  f2.Close();
  RooFFTConvPdf* PartRecoPDF_Bu2DstKPi_Conv   = new RooFFTConvPdf("kest1_Bu2DstKPi_conv",  "kest1_Bu2DstKPi_conv",  *mB,*PartRecoPDF_Bu2DstKPi, *gauss_test  );


  // Import Gen MC shapes 

  RooRealVar* gen_res       = new RooRealVar("gen_res",   "gen_res",    10  ) ;
  RooGaussian* gauss_res   = new RooGaussian("gauss_res",  "gauss_res",  *mB,*global_shift,*gen_res  ) ;

  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing Gen shape 1"<<std::endl;
  TFile f3((locationRooKeys + "GEN_kest1_B02Dst0Kst0_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace3 = (RooWorkspace*)f3.Get("workspace_B02Dst0Kst0");
  RooKeysPdf* PartRecoPDF_B02Dst0Kst0  = (RooKeysPdf*)workspace3->pdf("kest1_B02Dst0Kst0");
  f3.Close();

  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing Gen shape 2"<<std::endl;
  TFile f4((locationRooKeys + "GEN_kest1_B02Dst0Kst0_2_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace4 = (RooWorkspace*)f4.Get("workspace_B02Dst0Kst0_2");
  RooKeysPdf* PartRecoPDF_B02Dst0Kst0_2  = (RooKeysPdf*)workspace4->pdf("kest1_B02Dst0Kst0_2");
  f4.Close();
 
  RooRealVar* gen_frac = new RooRealVar("gen_frac","gen_frac", 0.732);
  //RooRealVar* gen_frac = new RooRealVar("gen_frac","gen_frac", 0.732, 0.0 ,1.0);
  RooAbsPdf*  pdf_gen  = new RooAddPdf(  "pdf_Gen", "", RooArgSet(*PartRecoPDF_Bu2DstKPi_Conv,*PartRecoPDF_B02Dst0Kst0_2),*gen_frac);  

  RooFFTConvPdf* PartRecoPDF_B02Dst0Kst0_Conv   = new RooFFTConvPdf("kest1_B02Dst0Kst0_conv",  "kest1_B02Dst0Kst0_conv",  *mB,*pdf_gen, *gauss_res  );

  


  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing Gen shape 3"<<std::endl;
  TFile f5((locationRooKeys + "GEN_kest1_Bs2DKst0pi_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace5 = (RooWorkspace*)f5.Get("workspace_Bs2DKst0pi");
  RooKeysPdf* PartRecoPDF_Bs2DKst0pi  = (RooKeysPdf*)workspace5->pdf("kest1_Bs2DKst0pi");
  f5.Close();

  RooFFTConvPdf* PartRecoPDF_Bs2DKst0pi_Conv   = new RooFFTConvPdf("kest1_Bs2DKst0pi_conv",  "kest1_Bs2DKst0pi_conv",  *mB,*PartRecoPDF_Bs2DKst0pi, *gauss_test  );


  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Importing Gen shape 4"<<std::endl;
  TFile f6((locationRooKeys + "GEN_kest1_Bs2DKst0pi_2_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace6 = (RooWorkspace*)f6.Get("workspace_Bs2DKst0pi_2");
  RooKeysPdf* PartRecoPDF_Bs2DKst0pi_2  = (RooKeysPdf*)workspace6->pdf("kest1_Bs2DKst0pi_2");
  f6.Close();

  RooFFTConvPdf* PartRecoPDF_Bs2DKst0pi_2_Conv   = new RooFFTConvPdf("kest1_Bs2DKst0pi_2_conv",  "kest1_Bs2DKst0pi_2_conv",  *mB,*PartRecoPDF_Bs2DKst0pi_2, *gauss_test  );



  RooRealVar* frac_Dsa1_DsstKKst = new RooRealVar("frac_Dsa1_DsstKKst","frac_Dsa1_DsstKKst", 0.816, 0.0, 1.0);

  //RooRealVar* frac_Dsa1_DsstKKst = new RooRealVar("frac_Dsa1_DsstKKst","frac_Dsa1_DsstKKst",0.997633244);
  // Convolve Dsa1 shapes with gaussians to take into account MC and Data resolutions differences
  // Sigma values taken from sigma = sqrt((sigma_Data)^2 - (sigma_MC)^2) 
  // Values taken from smallest sigma (sigma1) in DsD0 mode 
  // Assumed MC/Data ratio is the same in DKst0 and DsD0
  /*
  RooRealVar* mg        = new RooRealVar("mg","mg",0 + (varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,1.0):0) ) ;
  RooRealVar* sg_KKPi   = new RooRealVar("sg_KKPi",   "sg_KKPi",    0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg_PiKPi  = new RooRealVar("sg_PiKPi",  "sg_PiKPi",   0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg_PiPiPi = new RooRealVar("sg_PiPiPi", "sg_PiPiPi",  0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooGaussian* gauss_KKPi   = new RooGaussian("gauss_KKPi",  "gauss_KKPi",  *mB,*mg,*sg_KKPi  ) ;
  RooGaussian* gauss_PiKPi  = new RooGaussian("gauss_PiKPi", "gauss_PiKPi", *mB,*mg,*sg_PiKPi ) ;
  RooGaussian* gauss_PiPiPi = new RooGaussian("gauss_PiPiPi","gauss_PiPiPi",*mB,*mg,*sg_PiPiPi) ;

  RooRealVar* mg2        = new RooRealVar("mg2","mg2",0 + (varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,1.0):0) ) ;
  RooRealVar* sg2_KKPi   = new RooRealVar("sg2_KKPi",   "sg2_KKPi",    0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg2_PiKPi  = new RooRealVar("sg2_PiKPi",  "sg2_PiKPi",   0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg2_PiPiPi = new RooRealVar("sg2_PiPiPi", "sg2_PiPiPi",  0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooGaussian* gauss2_KKPi   = new RooGaussian("gauss2_KKPi",  "gauss2_KKPi",  *mB,*mg2,*sg2_KKPi  ) ;
  RooGaussian* gauss2_PiKPi  = new RooGaussian("gauss2_PiKPi", "gauss2_PiKPi", *mB,*mg2,*sg2_PiKPi ) ;
  RooGaussian* gauss2_PiPiPi = new RooGaussian("gauss2_PiPiPi","gauss2_PiPiPi",*mB,*mg2,*sg2_PiPiPi) ;


  std::map<std::string,RooFFTConvPdf*> PartRecoPDF_Bs02Dsa1_Convolved;
  PartRecoPDF_Bs2DKPiPi_Convolved[D2KKPi]   = new RooFFTConvPdf("kest1_adjust_KKPi",  "kest1_adjust_KKPi",  *mB,*PartRecoPDF_Bs2DKPiPi_Conv, *gauss_KKPi  );   
  PartRecoPDF_Bs2DKPiPi_Convolved[D2PiKPi]  = new RooFFTConvPdf("kest1_adjust_PiKPi", "kest1_adjust_PiKPi", *mB,*PartRecoPDF_Bs2DKPiPi_Conv, *gauss_PiKPi );
  PartRecoPDF_Bs2DKPiPi_Convolved[D2PiPiPi] = new RooFFTConvPdf("kest1_adjust_PiPiPi","kest1_adjust_PiPiPi",*mB,*PartRecoPDF_Bs2DKPiPi_Conv, *gauss_PiPiPi);


  std::map<std::string,RooFFTConvPdf*> PartRecoPDF_Bu2DstKPi_Convolved;
  PartRecoPDF_Bu2DstKPi_Convolved[D2KKPi]   = new RooFFTConvPdf("kest1_adjust_KKPi",  "kest1_adjust_KKPi",  *mB,*PartRecoPDF_Bu2DstKPi_Conv, *gauss2_KKPi  );  
  PartRecoPDF_Bu2DstKPi_Convolved[D2PiKPi]  = new RooFFTConvPdf("kest1_adjust_PiKPi", "kest1_adjust_PiKPi", *mB,*PartRecoPDF_Bu2DstKPi_Conv, *gauss2_PiKPi );
  PartRecoPDF_Bu2DstKPi_Convolved[D2PiPiPi] = new RooFFTConvPdf("kest1_adjust_PiPiPi","kest1_adjust_PiPiPi",*mB,*PartRecoPDF_Bu2DstKPi_Conv, *gauss2_PiPiPi);

  */

  std::map<std::string,std::string> mod;
  std::map<std::string,std::string> Bmod;
  mod[D2PiKPi]="#piK#pi";  mod[D2KKPi]="KK#pi";       mod[D2PiPiPi]="#pi#pi#pi";
  Bmod[DD0]="D_{d}D^{0}";  Bmod[DKst0]="D_{d}K*^{0}"; Bmod[DKst0Side]="D_{d}K*^{0} Sideband"; 

  // ======================================================== 
  // =========== Values for Combinatoric ====================
  // ======================================================== 
  std::map<std::string,std::map<std::string,std::map<std::string,double>>> Yield_comb;
  Yield_comb[DD0][D2PiKPi][both]        = 920; 
  Yield_comb[DKst0][D2PiKPi][both]      = 1230;
  Yield_comb[DKst0Side][D2PiKPi][both]  = 1230;

  Yield_comb[DD0][D2PiKPi][Helbin1]   = 1764;
  Yield_comb[DD0][D2PiKPi][Helbin2]   =  604;
 
  Yield_comb[DKst0][D2PiKPi][Helbin1]   = 561;
  Yield_comb[DKst0][D2PiKPi][Helbin2]   = 293;
  
  Yield_comb[DKst0Side][D2PiKPi][Helbin1]   = 583;
  Yield_comb[DKst0Side][D2PiKPi][Helbin2]   = 296;


  RooRealVar* global_comb_slope = new RooRealVar("global_comb_slope"  ,"Global Comb. Slope", -0.00183, -1.0, -0.0000000001);
  if(allConst_exYield) global_comb_slope->setConstant();
  
  std::map<std::string,std::map<std::string,RooRealVar*>> comb_slope_sys;
  comb_slope_sys[DD0][D2PiKPi]  = new RooRealVar("comb_slope_DD0_D2PiKPi"   ,"", -0.00258, -1.0, -0.0000000001);
  comb_slope_sys[DD0][D2KKPi]   = new RooRealVar("comb_slope_DD0_D2KKPi"    ,"", -0.00258, -1.0, -0.0000000001);
  comb_slope_sys[DD0][D2PiPiPi] = new RooRealVar("comb_slope_DD0_D2PiPiPi"  ,"", -0.00258, -1.0, -0.0000000001);

  comb_slope_sys[DKst0][D2PiKPi]  = comb_slope_sys[DD0][D2PiKPi] ;
  comb_slope_sys[DKst0][D2KKPi]   = comb_slope_sys[DD0][D2KKPi]  ;
  comb_slope_sys[DKst0][D2PiPiPi] = comb_slope_sys[DD0][D2PiPiPi];
  
  comb_slope_sys[DKst0Side][D2PiKPi]  = comb_slope_sys[DD0][D2PiKPi] ;
  comb_slope_sys[DKst0Side][D2KKPi]   = comb_slope_sys[DD0][D2KKPi]  ;
  comb_slope_sys[DKst0Side][D2PiPiPi] = comb_slope_sys[DD0][D2PiPiPi];


  // ======================================================== 
  // =========== Values for Signal  =========================
  // ======================================================== 

  std::map<std::string,std::map<std::string,double>> Yield_CB_input, Fixed_CB_n, Fixed_CB_alpha, Fixed_CB_Sigma_ratio, Fixed_CB_fraction;
 
  Yield_CB_input[DD0][D2PiKPi]  = 1624;
  Yield_CB_input[DD0][D2KKPi]   = 600;
  Yield_CB_input[DD0][D2PiPiPi] = 300;
 
  Yield_CB_input[DKst0][D2PiKPi]  = 50;
  Yield_CB_input[DKst0][D2KKPi]   = 40;
  Yield_CB_input[DKst0][D2PiPiPi] = 30;


  double factor = 1.0;
  if(par->variation[fixedSig_double]) factor = 2.0;
  // CB n -> Fixed to 1  
  Fixed_CB_n[DD0][D2PiKPi]            = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DD0][D2KKPi]             = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DD0][D2PiPiPi]           = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  Fixed_CB_n[DKst0][D2PiKPi]          = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DKst0][D2KKPi]           = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DKst0][D2PiPiPi]         = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  Fixed_CB_n[DKst0Side][D2PiKPi]      = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DKst0Side][D2KKPi]       = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DKst0Side][D2PiPiPi]     = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  // CB alpha -> From MC 
  Fixed_CB_alpha[DD0][D2PiKPi]        = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.13*factor):0);
  Fixed_CB_alpha[DD0][D2KKPi]         = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.13*factor):0);
  Fixed_CB_alpha[DD0][D2PiPiPi]       = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.13*factor):0);

  Fixed_CB_alpha[DKst0][D2PiKPi]      = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.25*factor):0);
  Fixed_CB_alpha[DKst0][D2KKPi]       = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.18*factor):0);
  Fixed_CB_alpha[DKst0][D2PiPiPi]     = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.26*factor):0);

  Fixed_CB_alpha[DKst0Side][D2PiKPi]  = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.25*factor):0);
  Fixed_CB_alpha[DKst0Side][D2KKPi]   = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.18*factor):0);
  Fixed_CB_alpha[DKst0Side][D2PiPiPi] = 3.25053 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.26*factor):0);

  // 2 CB sigma ratios -> From MC 
  Fixed_CB_Sigma_ratio[DD0][D2PiKPi]       = 0.450676 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.013):0);
  Fixed_CB_Sigma_ratio[DD0][D2KKPi]        = 0.450676 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.013):0);
  Fixed_CB_Sigma_ratio[DD0][D2PiPiPi]      = 0.450676 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.013):0);

  Fixed_CB_Sigma_ratio[DKst0][D2KKPi]      = 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);
  Fixed_CB_Sigma_ratio[DKst0][D2PiKPi]     = 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);
  Fixed_CB_Sigma_ratio[DKst0][D2PiPiPi]    = 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);

  Fixed_CB_Sigma_ratio[DKst0Side][D2KKPi]  = 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);
  Fixed_CB_Sigma_ratio[DKst0Side][D2PiKPi] = 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);
  Fixed_CB_Sigma_ratio[DKst0Side][D2PiPiPi]= 0.491303 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.021):0);


  // 2 CB sigma fractions -> From MC 
  Fixed_CB_fraction[DD0][D2PiKPi]        = 0.855871 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.010):0);
  Fixed_CB_fraction[DD0][D2KKPi]         = 0.855871 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.010):0);
  Fixed_CB_fraction[DD0][D2PiPiPi]       = 0.855871 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.010):0);

  Fixed_CB_fraction[DKst0][D2PiKPi]      = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);
  Fixed_CB_fraction[DKst0][D2KKPi]       = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);
  Fixed_CB_fraction[DKst0][D2PiPiPi]     = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);

  Fixed_CB_fraction[DKst0Side][D2PiKPi]  = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);
  Fixed_CB_fraction[DKst0Side][D2KKPi]   = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);
  Fixed_CB_fraction[DKst0Side][D2PiPiPi] = 0.815666 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.020):0);


  std::map<std::string,double> Fixed_Norm_Sigma_ratio;
  // Ratio of sigmas from DD0 to DKst0 (small sigma)
  // Ratio = Sigma(DKst0)/Sigma(DD0)
  Fixed_Norm_Sigma_ratio[D2PiKPi]        = 1.40958 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.172):0);
  Fixed_Norm_Sigma_ratio[D2KKPi]         = 1.40958 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.172):0);
  Fixed_Norm_Sigma_ratio[D2PiPiPi]       = 1.40958 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.172):0);

  // Initial Sigma values
  std::map<std::string,double> Inital_Sigma_values;
  Inital_Sigma_values[D2PiKPi]  =  8.2;
  Inital_Sigma_values[D2KKPi]   =  7.8;
  Inital_Sigma_values[D2PiPiPi] =  8.6; 


  // ======================================================== 
  // =========== Values for backgrounds      ================
  // ======================================================== 
  std::map<std::string,std::map<std::string,double>> Yield_Dsa1_input;
 
  Yield_Dsa1_input[DKst0][D2KKPi]   = 110; 
  Yield_Dsa1_input[DKst0][D2PiPiPi] =  51; 
  Yield_Dsa1_input[DKst0][D2PiKPi]  =  13; 

  std::map<std::string,std::map<std::string,double>> Yield_PR_total_input;
 
  Yield_PR_total_input[DD0][D2KKPi]         = 2177; 
  Yield_PR_total_input[DD0][D2PiKPi]        = 3897;
  Yield_PR_total_input[DD0][D2PiPiPi]       =  693;

  Yield_PR_total_input[DKst0][D2KKPi]       =   131; 
  Yield_PR_total_input[DKst0][D2PiKPi]      =   711;
  Yield_PR_total_input[DKst0][D2PiPiPi]     =    55; 

  Yield_PR_total_input[DKst0Side][D2KKPi]   =   131; 
  Yield_PR_total_input[DKst0Side][D2PiKPi]  =   600;
  Yield_PR_total_input[DKst0Side][D2PiPiPi] =    55; 

  //Fraction of D*D0 initial value
  std::map<std::string,double> Fraction_DstD0_initial;
  Fraction_DstD0_initial[D2KKPi]   = 0.37; //0.68;
  Fraction_DstD0_initial[D2PiPiPi] = 0.49;//0.75;
  Fraction_DstD0_initial[D2PiKPi]  = 0.37;//0.56;
  RooRealVar* Fraction_DstD0        = new RooRealVar("Fraction_DstD0",         "" , 0.49,  0 ,1); 
  RooRealVar* Fraction_DstD0_Helbin1= new RooRealVar("Fraction_DstD0_Helbin1", "" , 0.252, 0 ,1); 
  RooRealVar* Fraction_DstD0_Helbin2= new RooRealVar("Fraction_DstD0_Helbin2", "" , 0.126, 0 ,1); 

  if(allConst_exYield) Fraction_DstD0_Helbin1->setConstant();
  if(allConst_exYield) Fraction_DstD0_Helbin2->setConstant(); 


  // ======================================================== 
  // =========== Global variables      ======================
  // ========================================================
  
  RooRealVar* global_csi_hill   = new RooRealVar("global_csi_hill",   "Global csi HILL",    1.0             );
  RooRealVar* global_csi        = new RooRealVar("global_csi",        "Global csi",         0.970,  0,    2 );
  if(allConst_exYield) global_shift->setConstant();
  if(allConst_exYield) global_csi->setConstant();

  // Fraction of D* -> D pi0 vs. D* -> D pi0 gamma adjusted for fraction in mass window
  RooRealVar* fraction_Dst0_D0pi0_DD0_adj     = new RooRealVar("fraction_Dst0_D0pi0_DD0",         "",  0.619                * (0.838/0.713) + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,0.03):0));
  RooRealVar* fraction_Dst_Dpi0_DD0_adj       = new RooRealVar("fraction_Dst_Dpi0_DD0_adj",       "",  0.307*0.841/(0.3017*0.841+0.016*0.732) + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,0.007):0));

  //RooRealVar* fraction_Dst_Dpi0_DKst0_010_adj   = new RooRealVar("fraction_Dst_Dpi0_DKst0_adj",   "",  0.307*0.531/(0.3017*0.531+0.016*0.610) + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,0.007):0));
  //RooRealVar* fraction_Dst_Dpi0_DKst0_101_adj   = new RooRealVar("fraction_Dst_Dpi0_DKst0_adj",   "",  0.307*0.779/(0.3017*0.779+0.016*0.544) + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,0.007):0));
 
  RooRealVar* fraction_Dst_Dpi0_DKst0_010_adj   = new RooRealVar("fraction_Dst_Dpi0_DKst0_adj",   "",  0.307/(0.3017+0.016) + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,0.007):0));
  RooRealVar* fraction_Dst_Dpi0_DKst0_101_adj   = new RooRealVar("fraction_Dst_Dpi0_DKst0_adj",   "",  0.307/(0.3017+0.016) + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,0.007):0));
  
  // ======================================================== 
  // ===========    D*Kst0 shape parameters    ==============
  // ========================================================
  RooRealVar* DstKst0_long_Helbin1_frac = new RooRealVar("DstKst0_long_Helbin1_frac", "",  0.936 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.1):0));
  RooRealVar* DstKst0_tran_Helbin1_frac = new RooRealVar("DstKst0_tran_Helbin1_frac", "",  0.432 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.1):0));

  RooRealVar* fraction_hel              = new RooRealVar("fraction_hel",   "Polariztion Fraction",  0.5 + (varyAllowed&&par->variation[fixedBG_hel]?rand->Gaus(0,0.2):0) , 0,1 );
  if(fraction_hel->getVal() > 1 ) fraction_hel->setVal(1.0);
  if(fraction_hel->getVal() < 0 ) fraction_hel->setVal(0.0);

  if(allConst_exYield) fraction_hel->setConstant();
  //fraction_hel->setConstant();

  RooFormulaVar* DstKst0_Helbin1_frac  = new RooFormulaVar("DstKst0_Helbin1_frac", "",  "@0*@1 + (1-@0)*@2",         RooArgList(*fraction_hel,*DstKst0_long_Helbin1_frac,*DstKst0_tran_Helbin1_frac));
  

  // ======================================================== 
  // =========== Fixing ratios between Ds modes =============
  // =========== in Hel split mode to improve   =============
  // =========== stability                      =============
  // ========================================================

  // Fraction of DD0 peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DD0_peak_fraction    = new RooRealVar("splitHel_DD0_peak_fraction",    "", 0.580, 0.0 ,1.0 );
  if(allConst_exYield) splitHel_DD0_peak_fraction->setConstant();

  // Fraction of DD0 PR peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DD0_PR_peak_fraction = new RooRealVar("splitHel_DD0_PR_peak_fraction", "", 0.579, 0.0 ,1.0 );
  if(allConst_exYield) splitHel_DD0_PR_peak_fraction->setConstant();
  
  // Fraction of Ds*Phi peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DstKst0_peak_fraction = new RooRealVar("splitHel_DstKst0_peak_fraction", "", 0.920, 0.0 ,1.0 );
  splitHel_DstKst0_peak_fraction->setConstant();


  // ========================================================
  // =========== Blinding numbers               =============
  // ========================================================
  std::map<std::string,std::map<std::string,std::map<std::string,double>>> Blind_number; 
  Blind_number[DKst0][D2KKPi][both]          = 20;
  Blind_number[DKst0][D2PiPiPi][both]        = 20;
  Blind_number[DKst0][D2PiKPi][both]         = 20;
  
  Blind_number[DKst0][D2KKPi][Helbin1]       = 20;
  Blind_number[DKst0][D2PiPiPi][Helbin1]     = 20;
  Blind_number[DKst0][D2PiKPi][Helbin1]      = 20;
  
  Blind_number[DKst0][D2KKPi][Helbin2]       = 20;
  Blind_number[DKst0][D2PiPiPi][Helbin2]     = 20;
  Blind_number[DKst0][D2PiKPi][Helbin2]      = 20;

  Blind_number[DKst0Side][D2KKPi][both]      = 20;
  Blind_number[DKst0Side][D2PiPiPi][both]    = 20;
  Blind_number[DKst0Side][D2PiKPi][both]     = 20;
  
  Blind_number[DKst0Side][D2KKPi][Helbin1]   = 20;
  Blind_number[DKst0Side][D2PiPiPi][Helbin1] = 20;
  Blind_number[DKst0Side][D2PiKPi][Helbin1]  = 20;
  
  Blind_number[DKst0Side][D2KKPi][Helbin2]   = 20;
  Blind_number[DKst0Side][D2PiPiPi][Helbin2] = 20;
  Blind_number[DKst0Side][D2PiKPi][Helbin2]  = 20;

 
  // ======================================================== 
  // =========== Create Branching fraction ==================
  // ======================================================== 
  

  //============= Efficiency ratios weighted for years 
  //== eff taken from ANA note
  //== weighting taken from average of number of events in DsD0 vs. DsPhi + DsPhiSide roodatasets
  //=========================================================
  //== Ratio is eff(DsD0)/eff(DsPhi)

  // ratio * weight 
  //====================|==== 2011 ===|==== 2012 ====|==== 2015 ====|=== 2016 ===|

  eff_ratio[D2PiKPi]  =  1.243*0.12  +  1.265*0.28  +  1.284*0.09  +  1.278*0.52;

  eff_ratio[D2KKPi]   =  0.749*0.20  +  0.769*0.49  +  0.723*0.15  +  0.740*0.16;
  eff_ratio[D2PiPiPi] =  0.881*0.16  +  0.895*0.41  +  0.895*0.23  +  0.912*0.21;

  // ratio * weight 
  //========================|==== 2011 ===|==== 2012 ====|==== 2015 ====|=== 2016 ===|

  eff_ratio_err[D2PiKPi]  =  0.022*0.18  +  0.023*0.48  +  0.029*0.16  +  0.030*0.16;

  eff_ratio_err[D2KKPi]   =  0.011*0.20  +  0.012*0.49  +  0.015*0.15  +  0.016*0.16;
  eff_ratio_err[D2PiPiPi] =  0.014*0.16  +  0.015*0.41  +  0.017*0.23  +  0.017*0.21;
 


  eff_ratio_rrv[D2KKPi]   = new RooRealVar("eff_ratio_D2KKPi",  "", eff_ratio[D2KKPi]  );
  eff_ratio_rrv[D2PiKPi]  = new RooRealVar("eff_ratio_D2PiKPi", "", eff_ratio[D2PiKPi] );
  eff_ratio_rrv[D2PiPiPi] = new RooRealVar("eff_ratio_D2PiPiPi","", eff_ratio[D2PiPiPi]);



 
  bool fitSingleBR = par->fitBr;
  bool fitFourBr = par->fitFourBr; 
  
  // Single Branching fraction
  //Branching_fraction_all = new RooRealVar("Branching_fraction",       "",  eff_ratio_rrv[Ds2KKPi]->getVal()*(Yield_CB_input[DsPhi][Ds2KKPi]/Yield_CB_input[DsD0][Ds2KKPi]) ,0.0,1.0 );
    RooRealVar* Correction_factor = new RooRealVar("Correction_factor","",224.01);

  Branching_fraction_all = new RooRealVar("Branching_fraction",       "Branching Fraction (#times10^{-7})",  1.00, -1000.0, 10000.0 );
  if(NoBR&&!genModel){
    Branching_fraction_all->setVal(0.0);
    Branching_fraction_all->setConstant();
  }

  if(par->doSensitivity&&genModel){ 
    Branching_fraction_all->setVal(par->sensitivityBR);
  } else if(par->doSensitivity&&!genModel){
    Branching_fraction_all->setVal(0.0);
    Branching_fraction_all->setConstant();  
  }

  // Four seperate Branching fractions 
  Branching_fraction[D2KKPi]    = new RooRealVar("Branching_fraction_D2KKPi",   "",  1.0, -1000.0, 100000.0 );
  Branching_fraction[D2PiKPi]   = new RooRealVar("Branching_fraction_D2PiKPi",  "",  1.0, -1000.0, 100000.0 );
  Branching_fraction[D2PiPiPi]  = new RooRealVar("Branching_fraction_D2PiPiPi", "",  1.0, -1000.0, 100000.0 );
 

  // ======================================================== 
  // =========== Fractions in different cats ================
  // ======================================================== 


  //=============== Signal fractions ================= 
  // Fraction of signal decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* Signal_DKst0_Fraction   = new RooRealVar("Signal_DKst0_Fraction",   "", 0.79852);
  // Fraction of signal decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Signal_Helbin1_Fraction = new RooRealVar("Signal_Helbin1_Fraction", "", 0.92863);

  


  //=============== Non-resonant fractions ===========
  // Fraction of DKPi decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DKPi_DKst0_Fraction   = new RooRealVar("DKPi_DKst0_Fraction",   "", 0.34767);
  // Fraction of DKPi decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* DKPi_Helbin1_Fraction = new RooRealVar("DKPi_Helbin1_Fraction", "", 0.79422);


  //=============== Ds*Phi fractions =================
  // Fraction of D*Kst0 decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DstKst0_DKst0_Fraction  = new RooRealVar("DstKst0_DKst0_Fraction",  "", 0.79272 + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,0.01):0));

  

  //=============== Background fractions =============
  // Fraction of DKPiPi decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DKPiPi_DKst0_Fraction   = new RooRealVar("DKPiPi_DKst0_Fraction",   "", 0.33464);
  // Fraction of DKPiPi decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* DKPiPi_Helbin1_Fraction = new RooRealVar("DKPiPi_Helbin1_Fraction", "", 0.55435);

  // Fraction of DstKPi decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DstKPi_DKst0_Fraction   = new RooRealVar("DstKPi_DKst0_Fraction",   "", 0.29489);
  // Fraction of DstKPi decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* DstKPi_Helbin1_Fraction = new RooRealVar("DstKPi_Helbin1_Fraction", "", 0.80158);

  // Fraction of Dst0Kst0 decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* Dst0Kst0_DKst0_Fraction   = new RooRealVar("Dst0Kst0_DKst0_Fraction",   "", 0.71788);
  // Fraction of Dst0Kst0 decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Dst0Kst0_Helbin1_Fraction = new RooRealVar("Dst0Kst0_Helbin1_Fraction", "", 0.93319);
  
  // Fraction of Dst0Kst0 decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* Dst0Kst0_2_DKst0_Fraction   = new RooRealVar("Dst0Kst0_2_DKst0_Fraction",   "", 0.35118);
  // Fraction of Dst0Kst0 decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Dst0Kst0_2_Helbin1_Fraction = new RooRealVar("Dst0Kst0_2_Helbin1_Fraction", "", 0.85142);
  
  
  // Fraction of Dst0Kst0 decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* Dst0Kst0_Avg_DKst0_Fraction   = new RooRealVar("Dst0Kst0_Avg_DKst0_Fraction",   "", 0.59448);
  // Fraction of Dst0Kst0 decays in Helbin1 / (Helbin1+ Helbin2)
  //RooRealVar* Dst0Kst0_Avg_Helbin1_Fraction = new RooRealVar("Dst0Kst0_Avg_Helbin1_Fraction", "", 0.91173);
  RooRealVar* Dst0Kst0_Avg_Helbin1_Fraction = new RooRealVar("Dst0Kst0_Avg_Helbin1_Fraction", "", 0.771353); // , 0.0,1.0);


  // Fraction of DKst0pi decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DKst0pi_DKst0_Fraction   = new RooRealVar("DKst0pi_DKst0_Fraction",   "", 0.79852);
  // Fraction of DKst0pi decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* DKst0pi_Helbin1_Fraction = new RooRealVar("DKst0pi_Helbin1_Fraction", "", 0.52511); //, 0, 1);
  
  
  // Fraction of DKst0pi decays in DKst0 / (DKst0Side + DKst0)
  RooRealVar* DKst0pi_2_DKst0_Fraction   = new RooRealVar("DKst0pi_2_DKst0_Fraction",   "", 0.41380); //, 0, 1);
  // Fraction of DKst0pi decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* DKst0pi_2_Helbin1_Fraction = new RooRealVar("DKst0pi_2_Helbin1_Fraction", "", 0.52511);// , 0, 1);
  
  
  
  if(par->Bmodes[DKst0Side]&&!(par->Bmodes[DKst0Side])){
    Signal_DKst0_Fraction->setVal(1.0);
    DstKst0_DKst0_Fraction->setVal(1.0);
    DKPiPi_DKst0_Fraction->setVal(1.0);
    DstKPi_DKst0_Fraction->setVal(1.0);
    Dst0Kst0_DKst0_Fraction->setVal(1.0);
    Dst0Kst0_2_DKst0_Fraction->setVal(1.0);
    Dst0Kst0_Avg_DKst0_Fraction->setVal(1.0);
    DKst0pi_DKst0_Fraction->setVal(1.0);
    DKst0pi_2_DKst0_Fraction->setVal(1.0); 
  }



  // --------- Yields for DD0, in an array (as line 447 of Model.C) ---------------
  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Making RooRealVars"<<std::endl;
 
  // ------------------------------------------------------
  // --------- Setup yields for signals       -------------
  // ------------------------------------------------------
 

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){    
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
        for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
          std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*c).c_str(),(*a).c_str());
          yield_peak[*t][both][DD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",      DD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DD0].c_str(),      mod[(*m)].c_str()), Yield_CB_input[DD0][*m],    0, 100000);    
          if(allConst_exYield) ((RooRealVar*)yield_peak[*t][both][DD0][*m][*c][*a])->setConstant();
          
          if(fitFourBr){
            Signal_total[*t][both][DKst0][*m][*c][*a]   = new RooFormulaVar(Form("yield_peak_total_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),"(@0*@1)/(@2*@3)" , RooArgList(*Branching_fraction[*m],*yield_peak[*t][both][DD0][*m][*c][*a], *eff_ratio_rrv[*m],*Correction_factor));
          }else if(fitSingleBR){
            Signal_total[*t][both][DKst0][*m][*c][*a]   = new RooFormulaVar(Form("yield_peak_total_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),"(@0*@1)/(@2*@3)" , RooArgList(*Branching_fraction_all,*yield_peak[*t][both][DD0][*m][*c][*a], *eff_ratio_rrv[*m],*Correction_factor)); 
          } else {
            Signal_total[*t][both][DKst0][*m][*c][*a]   = new RooRealVar(   Form("yield_peak_total_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DKst0][*m],  0, 100000); 
            if(NoYield[*m]){
              ((RooRealVar*)Signal_total[*t][both][DKst0][*m][*c][*a])->setVal(0.0);
              ((RooRealVar*)Signal_total[*t][both][DKst0][*m][*c][*a])->setConstant();
            }
          }

          yield_peak[*t][both][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),"@0*@1"     , RooArgList( *Signal_DKst0_Fraction, *Signal_total[*t][both][DKst0][*m][*c][*a]));
          yield_peak[*t][both][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0Side.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),"(1-@0)*@1" , RooArgList( *Signal_DKst0_Fraction, *Signal_total[*t][both][DKst0][*m][*c][*a]));
          
          // Blinding
          if(needsBlinding){
            B_yield[*t][both][DD0][*m][*c][*a]      = yield_peak[*t][both][DD0][*m][*c][*a];
            B_yield[*t][both][DKst0][*m][*c][*a]    = new RooUnblindUniform(Form("B_nsig_DKst0_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),  Form("nsigBuDKst0blindePhi%s",    cat_name.c_str()),   Blind_number[DKst0][*m][both],     *yield_peak[*t][both][DKst0][*m][*c][*a] );
            B_yield[*t][both][DKst0Side][*m][*c][*a]= new RooUnblindUniform(Form("B_nsig_DKst0Side_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),  Form("nsigBuDKst0SideblindePhi%s",cat_name.c_str()),   Blind_number[DKst0Side][*m][both], *yield_peak[*t][both][DKst0Side][*m][*c][*a] );
          } else {
            B_yield[*t][both][DD0][*m][*c][*a]      = yield_peak[*t][both][DD0][*m][*c][*a];
            B_yield[*t][both][DKst0][*m][*c][*a]    = yield_peak[*t][both][DKst0][*m][*c][*a];
            B_yield[*t][both][DKst0Side][*m][*c][*a]= yield_peak[*t][both][DKst0Side][*m][*c][*a];
          }

          // Fix yields for helicity split state: Helbin1
          cat_name = Form("%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*c).c_str(),(*a).c_str());
          yield_peak[*t][Helbin1][DD0][*m][*c][*a]       = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DD0].c_str(),        mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*splitHel_DD0_peak_fraction,  *yield_peak[*t][both][DD0][*m][*c][*a]));
          yield_peak[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),      mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][DKst0][*m][*c][*a]));
          yield_peak[*t][Helbin1][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0Side.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0Side].c_str(),  mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][DKst0Side][*m][*c][*a]));
          // Blinding
          if(needsBlinding){
            B_yield[*t][Helbin1][DD0][*m][*c][*a]       = yield_peak[*t][Helbin1][DD0][*m][*c][*a];
            B_yield[*t][Helbin1][DKst0][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DKst0_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),  Form("nsigBuDKst0blindePhi%s",    cat_name.c_str()),   Blind_number[DKst0][*m][Helbin1],     *yield_peak[*t][Helbin1][DKst0][*m][*c][*a] );
            B_yield[*t][Helbin1][DKst0Side][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DKst0Side_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),  Form("nsigBuDKst0SideblindePhi%s",cat_name.c_str()),   Blind_number[DKst0Side][*m][Helbin1], *yield_peak[*t][Helbin1][DKst0Side][*m][*c][*a] );
          } else {
            B_yield[*t][Helbin1][DD0][*m][*c][*a]       = yield_peak[*t][Helbin1][DD0][*m][*c][*a];
            B_yield[*t][Helbin1][DKst0][*m][*c][*a]     = yield_peak[*t][Helbin1][DKst0][*m][*c][*a];
            B_yield[*t][Helbin1][DKst0Side][*m][*c][*a] = yield_peak[*t][Helbin1][DKst0Side][*m][*c][*a];
          }
          
          // Fix yields for helicity split state: Helbin2
          cat_name = Form("%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*c).c_str(),(*a).c_str());
          yield_peak[*t][Helbin2][DD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DD0].c_str(),       mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*splitHel_DD0_peak_fraction,  *yield_peak[*t][both][DD0][*m][*c][*a]));
          yield_peak[*t][Helbin2][DKst0][*m][*c][*a]    = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][DKst0][*m][*c][*a]));
          yield_peak[*t][Helbin2][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("yield_peak_%s_%s_%s",DKst0Side.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),"(1-@0)*@1" ,     RooArgList(*Signal_Helbin1_Fraction,     *yield_peak[*t][both][DKst0Side][*m][*c][*a]));
          
          if(needsBlinding){
            B_yield[*t][Helbin2][DD0][*m][*c][*a]       = yield_peak[*t][Helbin2][DD0][*m][*c][*a];
            B_yield[*t][Helbin2][DKst0][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DKst0_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),  Form("nsigBuDKst0blindePhi%s",    cat_name.c_str()),   Blind_number[DKst0][*m][Helbin2],     *yield_peak[*t][Helbin2][DKst0][*m][*c][*a] );
            B_yield[*t][Helbin2][DKst0Side][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DKst0Side_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),  Form("nsigBuDKst0SideblindePhi%s",cat_name.c_str()),   Blind_number[DKst0Side][*m][Helbin2], *yield_peak[*t][Helbin2][DKst0Side][*m][*c][*a] );
          } else {
            B_yield[*t][Helbin2][DD0][*m][*c][*a]      = yield_peak[*t][Helbin2][DD0][*m][*c][*a];
            B_yield[*t][Helbin2][DKst0][*m][*c][*a]    = yield_peak[*t][Helbin2][DKst0][*m][*c][*a];
            B_yield[*t][Helbin2][DKst0Side][*m][*c][*a]= yield_peak[*t][Helbin2][DKst0Side][*m][*c][*a];
          }   

          // ------------------------------------------------------
          // --------- Setup yields for backgrounds   -------------
          // ------------------------------------------------------

          
          // Both
          cat_name = Form("%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*c).c_str(),(*a).c_str());              
           
          Low_Mass_total[*t][both][DD0][*m][*c][*a]              = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DD0.c_str(),       (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DD0].c_str(),      mod[(*m)].c_str()),  Yield_PR_total_input[DD0][*m],   0,  1000000);
          Low_Mass_total[*t][both][DKst0][*m][*c][*a]            = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DKst0][*m], 0,  1000000);
          if(allConst_exYield) ((RooRealVar*)Low_Mass_total[*t][both][DD0][*m][*c][*a])->setConstant();
          if(allConst_exYield) ((RooRealVar*)Low_Mass_total[*t][both][DKst0][*m][*c][*a])->setConstant();

 
          PR_total_yield[*t][both][DD0][*m][*c][*a]              = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DD0.c_str(),       (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DD0].c_str(),      mod[(*m)].c_str()),  "@0",          RooArgList(*Low_Mass_total[*t][both][DD0][*m][*c][*a] ));//;
          PR_total_yield[*t][both][DKst0][*m][*c][*a]            = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1",       RooArgList(*DstKst0_DKst0_Fraction, *Low_Mass_total[*t][both][DKst0][*m][*c][*a] ));//;
          PR_total_yield[*t][both][DKst0Side][*m][*c][*a]        = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1",   RooArgList(*DstKst0_DKst0_Fraction, *Low_Mass_total[*t][both][DKst0][*m][*c][*a] ));//;
          
          DKPi_total_yield[*t][both][DKst0][*m][*c][*a]          = new RooRealVar(   Form("DKPi_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  10, 0,  1000000);
          DKPi_total_yield[*t][both][DKst0Side][*m][*c][*a]      = new RooFormulaVar(Form("DKPi_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DKPi_DKst0_Fraction,*DKPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          
          if(needsBlinding){
            DKPi_unblind_yield[*t][both][DKst0][*m][*c][*a]      = new RooUnblindUniform(Form("DKPi_unblind_DKst0_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind DKPi Yield %s %s",Bmod[DKst0].c_str(),     mod[(*m)].c_str()),  Form("nsigDKpiblinde%s",    cat_name.c_str()),   Blind_number[DKst0][*m][both],     *DKPi_total_yield[*t][both][DKst0][*m][*c][*a] );
            DKPi_unblind_yield[*t][both][DKst0Side][*m][*c][*a]  = new RooUnblindUniform(Form("DKPi_unblind_DKst0Side_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind DKPi Yield %s %s",Bmod[DKst0Side].c_str(), mod[(*m)].c_str()),  Form("nsigDKpiSideblinde%s",cat_name.c_str()),   Blind_number[DKst0Side][*m][both], *DKPi_total_yield[*t][both][DKst0Side][*m][*c][*a] );           
          } else {
            DKPi_unblind_yield[*t][both][DKst0][*m][*c][*a]      = DKPi_total_yield[*t][both][DKst0][*m][*c][*a];
            DKPi_unblind_yield[*t][both][DKst0Side][*m][*c][*a]  = DKPi_total_yield[*t][both][DKst0Side][*m][*c][*a];
          }

          DKPiPi_total_yield[*t][both][DKst0][*m][*c][*a]        = new RooRealVar(   Form("DKPiPi_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  1567, 0,  1000000);
          DKPiPi_total_yield[*t][both][DKst0Side][*m][*c][*a]    = new RooFormulaVar(Form("DKPiPi_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DKPiPi_DKst0_Fraction,*DKPiPi_total_yield[*t][both][DKst0][*m][*c][*a]));

          DstKPi_total_yield[*t][both][DKst0][*m][*c][*a]        = new RooRealVar(   Form("DstKPi_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  1, 0,  1000000);
          DstKPi_total_yield[*t][both][DKst0Side][*m][*c][*a]    = new RooFormulaVar(Form("DstKPi_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DstKPi_DKst0_Fraction,*DstKPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          
          Dst0Kst0_total_yield[*t][both][DKst0][*m][*c][*a]      = new RooRealVar(   Form("Dst0Kst0_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  2012, 0,  1000000);
          Dst0Kst0_total_yield[*t][both][DKst0Side][*m][*c][*a]  = new RooFormulaVar(Form("Dst0Kst0_total_%s_%s_%s",DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*Dst0Kst0_Avg_DKst0_Fraction,*Dst0Kst0_total_yield[*t][both][DKst0][*m][*c][*a]));
          
          DKst0pi_total_yield[*t][both][DKst0][*m][*c][*a]       = new RooRealVar(   Form("DKst0pi_total_%s_%s_%s", DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0pi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  2012, 0,  1000000);
          DKst0pi_total_yield[*t][both][DKst0Side][*m][*c][*a]   = new RooFormulaVar(Form("DKst0pi_total_%s_%s_%s", DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0pi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DKst0pi_DKst0_Fraction,*DKst0pi_total_yield[*t][both][DKst0][*m][*c][*a]));
          
          DKst0pi_2_total_yield[*t][both][DKst0][*m][*c][*a]     = new RooRealVar(   Form("DKst0pi_2_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0pi_2 yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  2012, 0,  1000000);
          DKst0pi_2_total_yield[*t][both][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DKst0pi_2_total_%s_%s_%s",DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0pi_2 yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DKst0pi_2_DKst0_Fraction,*DKst0pi_2_total_yield[*t][both][DKst0][*m][*c][*a]));

          
          // Helbin 1
          cat_name = Form("%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*c).c_str(),(*a).c_str());
          PR_total_yield[*t][Helbin1][DD0][*m][*c][*a]           = new RooFormulaVar(Form("yield_PR_H1_total_%s_%s_%s",DD0.c_str(),       (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DD0].c_str(),      mod[(*m)].c_str()),  "@0*@1" ,RooArgList(*splitHel_DD0_PR_peak_fraction,*PR_total_yield[*t][both][DD0][*m][*c][*a]));
          PR_total_yield[*t][Helbin1][DKst0][*m][*c][*a]         = new RooFormulaVar(Form("yield_PR_H1_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DstKst0_Helbin1_frac,         *PR_total_yield[*t][both][DKst0][*m][*c][*a]));
          PR_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_H1_total_%s_%s_%s",DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DstKst0_Helbin1_frac,         *PR_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DKPi_unblind_yield[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DKPi_H1_unblind_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s", Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKPi_Helbin1_Fraction,        *DKPi_unblind_yield[*t][both][DKst0][*m][*c][*a]));
          DKPi_unblind_yield[*t][Helbin1][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DKPi_H1_unblind_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s", Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKPi_Helbin1_Fraction,        *DKPi_unblind_yield[*t][both][DKst0Side][*m][*c][*a]));
                    
          DKPiPi_total_yield[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DKPiPi_H1_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKPiPi_Helbin1_Fraction,    *DKPiPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKPiPi_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DKPiPi_H1_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKPiPi_Helbin1_Fraction,    *DKPiPi_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DstKPi_total_yield[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DstKPi_H1_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DstKPi_Helbin1_Fraction,    *DstKPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DstKPi_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DstKPi_H1_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DstKPi_Helbin1_Fraction,    *DstKPi_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          Dst0Kst0_total_yield[*t][Helbin1][DKst0][*m][*c][*a]   = new RooFormulaVar(Form("Dst0Kst0_H1_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*Dst0Kst0_Avg_Helbin1_Fraction,    *Dst0Kst0_total_yield[*t][both][DKst0][*m][*c][*a]));
          Dst0Kst0_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("Dst0Kst0_H1_total_%s_%s_%s",DKst0Side.c_str(),(*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*Dst0Kst0_Avg_Helbin1_Fraction,    *Dst0Kst0_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DKst0pi_total_yield[*t][Helbin1][DKst0][*m][*c][*a]    = new RooFormulaVar(Form("DKst0pi_H1_total_%s_%s_%s", DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKst0pi_Helbin1_Fraction,    *DKst0pi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKst0pi_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("DKst0pi_H1_total_%s_%s_%s", DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKst0pi_Helbin1_Fraction,    *DKst0pi_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DKst0pi_2_total_yield[*t][Helbin1][DKst0][*m][*c][*a]    = new RooFormulaVar(Form("DKst0pi_2_H1_total_%s_%s_%s", DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0pi_2 yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKst0pi_2_Helbin1_Fraction,    *DKst0pi_2_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKst0pi_2_total_yield[*t][Helbin1][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("DKst0pi_2_H1_total_%s_%s_%s", DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0Pi_2 yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DKst0pi_2_Helbin1_Fraction,    *DKst0pi_2_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          

          // Helbin 2
          cat_name = Form("%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*c).c_str(),(*a).c_str());
          PR_total_yield[*t][Helbin2][DD0][*m][*c][*a]           = new RooFormulaVar(Form("yield_PR_H2_total_%s_%s_%s",DD0.c_str(),       (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DD0].c_str(),      mod[(*m)].c_str()),   "(1-@0)*@1" ,RooArgList(*splitHel_DD0_PR_peak_fraction,*PR_total_yield[*t][both][DD0][*m][*c][*a]));
          PR_total_yield[*t][Helbin2][DKst0][*m][*c][*a]         = new RooFormulaVar(Form("yield_PR_H2_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0].c_str(),    mod[(*m)].c_str()),   "(1-@0)*@1", RooArgList(*DstKst0_Helbin1_frac,         *PR_total_yield[*t][both][DKst0][*m][*c][*a]));            
          PR_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_H2_total_%s_%s_%s",DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),   "(1-@0)*@1", RooArgList(*DstKst0_Helbin1_frac,         *PR_total_yield[*t][both][DKst0Side][*m][*c][*a]));            
          
          DKPi_unblind_yield[*t][Helbin2][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DKPi_H2_unblind_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s", Bmod[DKst0].c_str(),    mod[(*m)].c_str()),   "(1-@0)*@1", RooArgList(*DKPi_Helbin1_Fraction,        *DKPi_unblind_yield[*t][both][DKst0][*m][*c][*a]));
          DKPi_unblind_yield[*t][Helbin2][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DKPi_H2_unblind_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPi yield %s %s", Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),   "(1-@0)*@1", RooArgList(*DKPi_Helbin1_Fraction,        *DKPi_unblind_yield[*t][both][DKst0Side][*m][*c][*a]));
 
          DKPiPi_total_yield[*t][Helbin2][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DKPiPi_H2_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DKPiPi_Helbin1_Fraction,      *DKPiPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKPiPi_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DKPiPi_H2_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKPiPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DKPiPi_Helbin1_Fraction,      *DKPiPi_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DstKPi_total_yield[*t][Helbin2][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("DstKPi_H2_total_%s_%s_%s",  DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DstKPi_Helbin1_Fraction,      *DstKPi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DstKPi_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a] = new RooFormulaVar(Form("DstKPi_H2_total_%s_%s_%s",  DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DstKPi_Helbin1_Fraction,      *DstKPi_total_yield[*t][both][DKst0Side][*m][*c][*a])); 
            
          Dst0Kst0_total_yield[*t][Helbin2][DKst0][*m][*c][*a]   = new RooFormulaVar(Form("Dst0Kst0_H2_total_%s_%s_%s",DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*Dst0Kst0_Avg_Helbin1_Fraction,*Dst0Kst0_total_yield[*t][both][DKst0][*m][*c][*a]));
          Dst0Kst0_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("Dst0Kst0_H2_total_%s_%s_%s",DKst0Side.c_str(),(*m).c_str(), cat_name.c_str()), Form("DstKPi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*Dst0Kst0_Avg_Helbin1_Fraction,*Dst0Kst0_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DKst0pi_total_yield[*t][Helbin2][DKst0][*m][*c][*a]    = new RooFormulaVar(Form("DKst0pi_H2_total_%s_%s_%s", DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*DKst0pi_Helbin1_Fraction,     *DKst0pi_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKst0pi_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("DKst0pi_H2_total_%s_%s_%s", DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*DKst0pi_Helbin1_Fraction,     *DKst0pi_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
          DKst0pi_2_total_yield[*t][Helbin2][DKst0][*m][*c][*a]    = new RooFormulaVar(Form("DKst0pi_2_H2_total_%s_%s_%s", DKst0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0].c_str(),    mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*DKst0pi_2_Helbin1_Fraction,    *DKst0pi_2_total_yield[*t][both][DKst0][*m][*c][*a]));
          DKst0pi_2_total_yield[*t][Helbin2][DKst0Side][*m][*c][*a]= new RooFormulaVar(Form("DKst0pi_2_H2_total_%s_%s_%s", DKst0Side.c_str(), (*m).c_str(), cat_name.c_str()), Form("DKst0Pi yield %s %s",Bmod[DKst0Side].c_str(),mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*DKst0pi_2_Helbin1_Fraction,    *DKst0pi_2_total_yield[*t][both][DKst0Side][*m][*c][*a]));
          
        }
      }
    }
  }




  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
      for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
          for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
            for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 

              std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*c).c_str(),(*a).c_str());
              sig_frac[*t][*h][*b][*m]           = new RooRealVar(   Form("Sigma_frac_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("%s %s Sigma Fraction",Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_fraction[*b][*m]); 
              yield_comb[*t][*h][*b][*m][*c][*a] = new RooRealVar(   Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), Yield_comb[*b][*m][*h], 0, 600000);
              if(allConst_exYield) ((RooRealVar*)yield_comb[*t][*h][*b][*m][*c][*a])->setConstant();
            }
          }
        } 

        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){
            std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*c).c_str(),(*a).c_str());

            frac[*t][both][DD0][*m][*c][*a]               = Fraction_DstD0;
            frac[*t][Helbin1][DD0][*m][*c][*a]            = Fraction_DstD0_Helbin1;
            frac[*t][Helbin2][DD0][*m][*c][*a]            = Fraction_DstD0_Helbin2;

            yield_dstd0[*t][*h][DD0][*m][*c][*a]     = new RooFormulaVar(Form("yield_DstD0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D* D^{0} %s",mod[(*m)].c_str()), "@0*@1",     RooArgList(*frac[*t][*h][DD0][*m][*c][*a],  *PR_total_yield[*t][*h][DD0][*m][*c][*a])  );
            yield_ddst0[*t][*h][DD0][*m][*c][*a]     = new RooFormulaVar(Form("yield_DDst0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D D*^{0} %s",mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*frac[*t][*h][DD0][*m][*c][*a],  *PR_total_yield[*t][*h][DD0][*m][*c][*a])  );
          

            frac_HORNS[*t][*h][DKst0][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HORNS_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*@1",        RooArgList(*fraction_Dst_Dpi0_DKst0_010_adj,*fraction_hel));
            frac_HILL[*t][*h][DKst0][*m][*c][*a]        = new RooFormulaVar(Form("fraction_HILL_%s_%s",       (*m).c_str(),cat_name.c_str()), "", "(1-@0)*@1",    RooArgList(*fraction_Dst_Dpi0_DKst0_010_adj,*fraction_hel));
            frac_HILL2[*t][*h][DKst0][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HILL2_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*(1-@1)",    RooArgList(*fraction_Dst_Dpi0_DKst0_101_adj,*fraction_hel));
            frac_LITTLEHORNS[*t][*h][DKst0][*m][*c][*a] = new RooFormulaVar(Form("fraction_LITTLEHORNS_%s_%s",(*m).c_str(),cat_name.c_str()), "", "(1-@0)*(1-@1)",RooArgList(*fraction_Dst_Dpi0_DKst0_101_adj,*fraction_hel));  
          }
        }
      } 
    }
  }

  // --------- Mean B mass ---------------
  double default_mB=5279.29;

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*b).c_str(),(*m).c_str());
          mean_B[*t][*h][*b][*m][both][both] = new RooRealVar(Form("mean_B_%s"           ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][plus][up]   = new RooRealVar(Form("mean_B_%s_plus_up"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][minus][dn]  = new RooRealVar(Form("mean_B_%s_minus_dn"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][plus][dn]   = new RooRealVar(Form("mean_B_%s_plus_dn"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][minus][up]  = new RooRealVar(Form("mean_B_%s_minus_up"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][plus][both] = new RooRealVar(Form("mean_B_%s_plus_both" ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          mean_B[*t][*h][*b][*m][minus][both]= new RooRealVar(Form("mean_B_%s_minus_both",cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
          if(allConst_exYield) mean_B[*t][*h][*b][*m][both][both]->setConstant();

          //mean_B[*t][*h][*b][*m][both][up]          = mean_B[*t][*h][*b][*m][both][both];
          //mean_B[*t][*h][*b][*m][both][dn]          = mean_B[*t][*h][*b][*m][both][both];
          //mean_B[*t][*h][*b][D2PiKPi][both][both]   = mean_B[*t][*h][*b][Ds2KKPi][both][both];
          //mean_B[*t][*h][*b][D2PiPiPi][both][both]  = mean_B[*t][*h][*b][Ds2KKPi][both][both];
        }
      }
    }
  }

  // --------- Signal width ---------------
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
        for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){

                
          std::string cat_name = Form("%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());

          sigma[*t][both][DD0][*m][*c][*a]       = new RooRealVar(   Form("sigma_DD0_%s",  cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
          sigma[*t][both][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("sigma_DKst0_%s",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma[*t][both][DD0][*m][*c][*a])); 
          sigma[*t][both][DKst0Side][*m][*c][*a] = sigma[*t][both][DKst0][*m][*c][*a];

          sigma2[*t][both][DD0][*m][*c][*a]       = new RooFormulaVar(Form("sigma2_DD0_%s",  cat_name.c_str()), Form("Sigma 2 DD0 %s",  (*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DD0][*m] ),  RooArgList(*sigma[*t][both][DD0][*m][*c][*a]    ));
          sigma2[*t][both][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("sigma2_DKst0_%s",cat_name.c_str()), Form("Sigma 2 DKst0 %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DKst0][*m]), RooArgList(*sigma[*t][both][DKst0][*m][*c][*a]  ));
          sigma2[*t][both][DKst0Side][*m][*c][*a] = sigma2[*t][both][DKst0][*m][*c][*a];

          avg_sigma[*t][both][DD0][*m][both][both]       = new RooFormulaVar(Form("avg_sigma_DD0_%s",   cat_name.c_str()), Form("Average sigma %s %s",Bmod[DD0].c_str(),  mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DD0][*m],  Fixed_CB_fraction[DD0][*m]   ), RooArgList(*sigma[*t][both][DD0][*m][*c][*a],  *sigma2[*t][both][DD0][*m][*c][*a]   ));
          avg_sigma[*t][both][DKst0][*m][both][both]     = new RooFormulaVar(Form("avg_sigma_DKst0_%s", cat_name.c_str()), Form("Average sigma %s %s",Bmod[DKst0].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DKst0][*m],Fixed_CB_fraction[DKst0][*m] ), RooArgList(*sigma[*t][both][DKst0][*m][*c][*a],*sigma2[*t][both][DKst0][*m][*c][*a] ));
          avg_sigma[*t][both][DKst0Side][*m][both][both] = avg_sigma[*t][both][DKst0][*m][both][both];

          cat_name = Form("%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());

          sigma[*t][Helbin1][DD0][*m][*c][*a]       = new RooRealVar(   Form("sigma_DD0_%s",  cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
          if(allConst_exYield) ((RooRealVar*)sigma[*t][Helbin1][DD0][*m][*c][*a])->setConstant();
          sigma[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("sigma_DKst0_%s",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma[*t][Helbin1][DD0][*m][*c][*a])); 
          sigma[*t][Helbin1][DKst0Side][*m][*c][*a] = sigma[*t][Helbin1][DKst0][*m][*c][*a] ; 

          //sig_ratio[*t][Helbin1][*b][*m]         = new RooRealVar(Form("Sigma_ratio_%s", cat_name.c_str()), Form("%s %s Sigma ratio ",mod[*m].c_str(),Bmod[*b].c_str()), Fixed_CB_Sigma_ratio[*b][*m]);
  
          sigma2[*t][Helbin1][DD0][*m][*c][*a]       = new RooFormulaVar(Form("sigma2_DD0_%s",  cat_name.c_str()), Form("Sigma 2 DD0 %s",  (*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DD0][*m] ),  RooArgList(*sigma[*t][Helbin1][DD0][*m][*c][*a]    ));
          sigma2[*t][Helbin1][DKst0][*m][*c][*a]     = new RooFormulaVar(Form("sigma2_DKst0_%s",cat_name.c_str()), Form("Sigma 2 DKst0 %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DKst0][*m]), RooArgList(*sigma[*t][Helbin1][DKst0][*m][*c][*a]  ));
          sigma2[*t][Helbin1][DKst0Side][*m][*c][*a] = sigma2[*t][Helbin1][DKst0][*m][*c][*a];

          avg_sigma[*t][Helbin1][DD0][*m][both][both]       = new RooFormulaVar(Form("avg_sigma_DD0_%s",   cat_name.c_str()), Form("Average sigma %s %s",Bmod[DD0].c_str(),  mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DD0][*m],  Fixed_CB_fraction[DD0][*m]  ),  RooArgList(*sigma[*t][Helbin1][DD0][*m][*c][*a],  *sigma2[*t][Helbin1][DD0][*m][*c][*a]   ));
          avg_sigma[*t][Helbin1][DKst0][*m][both][both]     = new RooFormulaVar(Form("avg_sigma_DKst0_%s", cat_name.c_str()), Form("Average sigma %s %s",Bmod[DKst0].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DKst0][*m],Fixed_CB_fraction[DKst0][*m] ), RooArgList(*sigma[*t][Helbin1][DKst0][*m][*c][*a],*sigma2[*t][Helbin1][DKst0][*m][*c][*a] ));
          avg_sigma[*t][Helbin1][DKst0Side][*m][both][both] = avg_sigma[*t][Helbin1][DKst0][*m][both][both];
          

          sigma[*t][Helbin2][DD0][*m][*c][*a]               = sigma[*t][Helbin1][DD0][*m][*c][*a];
          sigma[*t][Helbin2][DKst0][*m][*c][*a]             = sigma[*t][Helbin1][DKst0][*m][*c][*a];
          sigma[*t][Helbin2][DKst0Side][*m][*c][*a]         = sigma[*t][Helbin1][DKst0Side][*m][*c][*a];

          sigma2[*t][Helbin2][DD0][*m][*c][*a]              = sigma2[*t][Helbin1][DD0][*m][*c][*a];        
          sigma2[*t][Helbin2][DKst0][*m][*c][*a]            = sigma2[*t][Helbin1][DKst0][*m][*c][*a];      
          sigma2[*t][Helbin2][DKst0Side][*m][*c][*a]        = sigma2[*t][Helbin1][DKst0Side][*m][*c][*a];          
          
          avg_sigma[*t][Helbin2][DD0][*m][both][both]       = avg_sigma[*t][Helbin1][DD0][*m][both][both]; 
          avg_sigma[*t][Helbin2][DKst0][*m][both][both]     = avg_sigma[*t][Helbin1][DKst0][*m][both][both];   
          avg_sigma[*t][Helbin2][DKst0Side][*m][both][both] = avg_sigma[*t][Helbin1][DKst0Side][*m][both][both];   
        }
      }
    }
  }

  // --------- Signal tail params ---------------
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*b).c_str(),(*m).c_str());
          
          ncb[*t][*h][*b][*m]         = new RooRealVar(Form("ncb_%s",        cat_name.c_str()), Form("%s %s  nCB",          Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_n[*b][*m]);
          alpha[*t][*h][*b][*m]       = new RooRealVar(Form("alphaL_%s",     cat_name.c_str()), Form("%s %s alpha",         Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_alpha[*b][*m]); 
        }
      }
    }
  }



  // --------- Comb. background slope ---------------
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){ 
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          std::string cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*b).c_str(),(*m).c_str());
          comb_slope[*t][*h][*b][*m] = global_comb_slope;
          if(varyAllowed&&par->variation[fixedBG_slope]) comb_slope[*t][*h][*b][*m] = comb_slope_sys[*b][*m];
        }
      }
    }
  }


  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
        std::string cat_name = Form("%s_%s_%s",(*t).c_str(),(*h).c_str(),(*m).c_str());
        
        // ------------------------------------------
        // PartReco shapes for DsD0 
        // ------------------------------------------
        
        // D*D0 pi0
        HORNS_a[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS_D*D0_a_%s",     cat_name.c_str()), "HORNS a",     5054.0  + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0) );
        HORNS_b[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS_D*D0_b_%s",     cat_name.c_str()), "HORNS b",     5124.8  + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0) );
        HORNS_sigma[*t][*h][DD0][*m]  = avg_sigma[*t][*h][DD0][*m][both][both];
        HORNS_R[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS_D*D0_R_%s",     cat_name.c_str()), "HORNS R",     1.0  );
        HORNS_f[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS_D*D0_f_%s",     cat_name.c_str()), "HORNS f",     1.0 );
        HORNS_csi[*t][*h][DD0][*m]    = global_csi;
        HORNS_shift[*t][*h][DD0][*m]  = global_shift;
        
        // D*D0 gamma
        HILL_a[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL_D*D0_a_%s",      cat_name.c_str()), "HILL a",      4967.013 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0)  );
        HILL_b[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL_D*D0_b_%s",      cat_name.c_str()), "HILL b",      5217.935 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0)  );
        HILL_sigma[*t][*h][DD0][*m]   = avg_sigma[*t][*h][DD0][*m][both][both];
        HILL_R[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL_D*D0_R_%s",      cat_name.c_str()), "HILL R",      1.0    );
        HILL_f[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL_D*D0_f_%s",      cat_name.c_str()), "HILL f",      1.0    ); 
        HILL_csi[*t][*h][DD0][*m]     = global_csi_hill;
        HILL_shift[*t][*h][DD0][*m]   = global_shift;
        
        // DD*0 pi0
        HORNS2_a[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS2_DD*0_a_%s",     cat_name.c_str()), "HORNS2 a",  5047.614 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0) );
        HORNS2_b[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS2_DD*0_b_%s",     cat_name.c_str()), "HORNS2 b",  5127.168 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0) );
        HORNS2_sigma[*t][*h][DD0][*m]  = avg_sigma[*t][*h][DD0][*m][both][both];
        HORNS2_R[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS2_DD*0_R_%s",     cat_name.c_str()), "HORNS2 R",  1.0   );
        HORNS2_f[*t][*h][DD0][*m]      = new RooRealVar(Form("HORNS2_DD*0_f_%s",     cat_name.c_str()), "HORNS2 f",  1.0   );
        HORNS2_csi[*t][*h][DD0][*m]    = global_csi;
        HORNS2_shift[*t][*h][DD0][*m]  = global_shift;
        
        // DD*0 gamma
        HILL2_a[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL2_DD*0_a_%s",      cat_name.c_str()), "HILL2 a",   4963.587 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0)  );
        HILL2_b[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL2_DD*0_b_%s",      cat_name.c_str()), "HILL2 b",   5217.383 + (varyAllowed&&par->variation[fixedBG_DD0]?rand->Gaus(0,1):0)  );
        HILL2_sigma[*t][*h][DD0][*m]   = avg_sigma[*t][*h][DD0][*m][both][both];
        HILL2_R[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL2_DD*0_R_%s",      cat_name.c_str()), "HILL2 R",   1.0   );
        HILL2_f[*t][*h][DD0][*m]       = new RooRealVar(Form("HILL2_DD*0_f_%s",      cat_name.c_str()), "HILL2 f",   1.0   );  
        HILL2_csi[*t][*h][DD0][*m]     = global_csi_hill;
        HILL2_shift[*t][*h][DD0][*m]   = global_shift;            // PartReco shapes for DD0 

        // ------------------------------------------
        // PartReco shapes for D*KSt0
        // ------------------------------------------

        // missed pi0 components 
        // pi0 010
        HORNS_a[*t][*h][DKst0][*m]      = new RooRealVar(Form("HORNS_DKst0_a_%s",     cat_name.c_str()), "HORNS a",    5028.156  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HORNS_b[*t][*h][DKst0][*m]      = new RooRealVar(Form("HORNS_DKst0_b_%s",     cat_name.c_str()), "HORNS b",    5113.659  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HORNS_sigma[*t][*h][DKst0][*m]  = avg_sigma[*t][*h][DKst0][*m][both][both];
        HORNS_R[*t][*h][DKst0][*m]      = new RooRealVar(Form("HORNS_DKst0_R_%s",     cat_name.c_str()), "HORNS R",    1.0    );
        HORNS_f[*t][*h][DKst0][*m]      = new RooRealVar(Form("HORNS_DKst0_f_%s",     cat_name.c_str()), "HORNS f",    1.0    );
        HORNS_csi[*t][*h][DKst0][*m]    = global_csi;
        HORNS_shift[*t][*h][DKst0][*m]  = global_shift;
        
        //pi0 101
        HILL2_a[*t][*h][DKst0][*m]      = new RooRealVar(Form("HILL2_DKst0_a_%s",      cat_name.c_str()), "HILL2 a",   5028.156  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL2_b[*t][*h][DKst0][*m]      = new RooRealVar(Form("HILL2_DKst0_b_%s",      cat_name.c_str()), "HILL2 b",   5113.659  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL2_sigma[*t][*h][DKst0][*m]  = avg_sigma[*t][*h][DKst0][*m][both][both];
        HILL2_R[*t][*h][DKst0][*m]      = new RooRealVar(Form("HILL2_DKst0_R_%s",      cat_name.c_str()), "HILL2 R",   1.0    );
        HILL2_f[*t][*h][DKst0][*m]      = new RooRealVar(Form("HILL2_DKst0_f_%s",      cat_name.c_str()), "HILL2 f",   1.0    );
        HILL2_csi[*t][*h][DKst0][*m]    = global_csi_hill;
        HILL2_shift[*t][*h][DKst0][*m]  = global_shift;
        
        // Missed gamma components 
        //gamma 010
        HILL_a[*t][*h][DKst0][*m]       = new RooRealVar(Form("HILL_DKst0_a_%s",      cat_name.c_str()), "HILL a",     4922.512  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL_b[*t][*h][DKst0][*m]       = new RooRealVar(Form("HILL_DKst0_b_%s",      cat_name.c_str()), "HILL b",     5225.377  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL_sigma[*t][*h][DKst0][*m]   = avg_sigma[*t][*h][DKst0][*m][both][both]; 
        HILL_R[*t][*h][DKst0][*m]       = new RooRealVar(Form("HILL_DKst0_R_%s",      cat_name.c_str()), "HILL R",      1.0   );
        HILL_f[*t][*h][DKst0][*m]       = new RooRealVar(Form("HILL_DKst0_f_%s",      cat_name.c_str()), "HILL f",      1.0   );
        HILL_csi[*t][*h][DKst0][*m]     = global_csi_hill;
        HILL_shift[*t][*h][DKst0][*m]   = global_shift;

        //gamma 101
        LITTLEHORNS_a[*t][*h][DKst0][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0_a_%s",     cat_name.c_str()), "LITTLEHORNS a",  4922.512  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        LITTLEHORNS_b[*t][*h][DKst0][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0_b_%s",     cat_name.c_str()), "LITTLEHORNS b",  5225.377  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        LITTLEHORNS_sigma[*t][*h][DKst0][*m]  = avg_sigma[*t][*h][DKst0][*m][both][both]; 
        LITTLEHORNS_R[*t][*h][DKst0][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0_R_%s",     cat_name.c_str()), "LITTLEHORNS R",  1.0    );
        LITTLEHORNS_f[*t][*h][DKst0][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0_f_%s",     cat_name.c_str()), "LITTLEHORNS f",  1.0    );
        LITTLEHORNS_g[*t][*h][DKst0][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0_g_%s",     cat_name.c_str()), "LITTLEHORNS g",  0.0    );
        LITTLEHORNS_csi[*t][*h][DKst0][*m]    = global_csi;
        LITTLEHORNS_shift[*t][*h][DKst0][*m]  = global_shift;

        // ------------------------------------------
        // PartReco shapes for D*KSt0 Side
        // ------------------------------------------

        // missed pi0 components 
        // pi0 010
        HORNS_a[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HORNS_DKst0Side_a_%s",     cat_name.c_str()), "HORNS a",    5028.156  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HORNS_b[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HORNS_DKst0Side_b_%s",     cat_name.c_str()), "HORNS b",    5113.659  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HORNS_sigma[*t][*h][DKst0Side][*m]  = avg_sigma[*t][*h][DKst0Side][*m][both][both];
        HORNS_R[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HORNS_DKst0Side_R_%s",     cat_name.c_str()), "HORNS R",    1.0    );
        HORNS_f[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HORNS_DKst0Side_f_%s",     cat_name.c_str()), "HORNS f",    1.0    );
        HORNS_csi[*t][*h][DKst0Side][*m]    = global_csi;
        HORNS_shift[*t][*h][DKst0Side][*m]  = global_shift;
        
        //pi0 101
        HILL2_a[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HILL2_DKst0Side_a_%s",      cat_name.c_str()), "HILL2 a",   5028.156  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL2_b[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HILL2_DKst0Side_b_%s",      cat_name.c_str()), "HILL2 b",   5113.659  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL2_sigma[*t][*h][DKst0Side][*m]  = avg_sigma[*t][*h][DKst0Side][*m][both][both];
        HILL2_R[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HILL2_DKst0Side_R_%s",      cat_name.c_str()), "HILL2 R",   1.0    );
        HILL2_f[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("HILL2_DKst0Side_f_%s",      cat_name.c_str()), "HILL2 f",   1.0    );
        HILL2_csi[*t][*h][DKst0Side][*m]    = global_csi_hill;
        HILL2_shift[*t][*h][DKst0Side][*m]  = global_shift;
        
        // Missed gamma components 
        //gamma 010
        HILL_a[*t][*h][DKst0Side][*m]       = new RooRealVar(Form("HILL_DKst0Side_a_%s",      cat_name.c_str()), "HILL a",     4922.512  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL_b[*t][*h][DKst0Side][*m]       = new RooRealVar(Form("HILL_DKst0Side_b_%s",      cat_name.c_str()), "HILL b",     5225.377  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        HILL_sigma[*t][*h][DKst0Side][*m]   = avg_sigma[*t][*h][DKst0Side][*m][both][both]; 
        HILL_R[*t][*h][DKst0Side][*m]       = new RooRealVar(Form("HILL_DKst0Side_R_%s",      cat_name.c_str()), "HILL R",      1.0   );
        HILL_f[*t][*h][DKst0Side][*m]       = new RooRealVar(Form("HILL_DKst0Side_f_%s",      cat_name.c_str()), "HILL f",      1.0   );
        HILL_csi[*t][*h][DKst0Side][*m]     = global_csi_hill;
        HILL_shift[*t][*h][DKst0Side][*m]   = global_shift;

        //gamma 101
        LITTLEHORNS_a[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0Side_a_%s",     cat_name.c_str()), "LITTLEHORNS a",  4922.512  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        LITTLEHORNS_b[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0Side_b_%s",     cat_name.c_str()), "LITTLEHORNS b",  5225.377  + (varyAllowed&&par->variation[fixedBG_DKst0]?rand->Gaus(0,1):0)  );
        LITTLEHORNS_sigma[*t][*h][DKst0Side][*m]  = avg_sigma[*t][*h][DKst0Side][*m][both][both]; 
        LITTLEHORNS_R[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0Side_R_%s",     cat_name.c_str()), "LITTLEHORNS R",  1.0    );
        LITTLEHORNS_f[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0Side_f_%s",     cat_name.c_str()), "LITTLEHORNS f",  1.0    );
        LITTLEHORNS_g[*t][*h][DKst0Side][*m]      = new RooRealVar(Form("LITTLEHORNS_DKst0Side_g_%s",     cat_name.c_str()), "LITTLEHORNS g",  0.0    );
        LITTLEHORNS_csi[*t][*h][DKst0Side][*m]    = global_csi;
        LITTLEHORNS_shift[*t][*h][DKst0Side][*m]  = global_shift;

      } 
    }
  }

  


  // ------------------------------------------------------
  // --------- Fix RooRealVars to each other  -------------
  // ------------------------------------------------------
 
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){           
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              
              // Use same mean value for all plots 
              mean_B[*t][*h][*b][*m][*c][*a] = mean_B[*typeList.begin()][*HelBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()]; 

              // Allow sigma to vary between Ds decay modes
              //sigma[*t][*h][*b][*m][*c][*a] = sigma[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*b][*m][*chargeList.begin()][*magnetList.begin()]; 
            }
          }
        }
      }
    }
  }




  // --------- Define the simultaneous PDFs ---------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  
  if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Setting up simultaneous PDF"<<std::endl;
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
  RooAbsPdf *pdf_DstD0 = 0;
  RooAbsPdf *pdf_DDst0 = 0;
  RooAbsPdf *pdf_Bs2DKPiPi = 0;
  RooAbsPdf *pdf_B2DKPi = 0;
  RooAbsPdf *pdf_Bu2DstKPi= 0;
  RooAbsPdf *pdf_B02Dst0Kst0= 0;
  RooAbsPdf *pdf_Bs2DKst0pi= 0;
  RooAbsPdf *pdf_Bs2DKst0pi_2= 0;

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              
              std::string tag=(*t)+underscore+(*h)+underscore+(*b)+underscore+(*m)+underscore+(*c)+underscore+(*a);
              //if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Adding to sim pdf: "<< tag <<std::endl;
              //Both b modes

              // Double CB shape from MC
              RooAbsReal* mub     = mean_B[*t][*h][*b][*m][*c][*a];
              RooAbsReal* sig_CB  = sigma[*t][*h][*b][*m][*c][*a];
              RooAbsReal* n_CB    = ncb[*t][*h][*b][*m];
              RooAbsReal* alp_CB  = alpha[*t][*h][*b][*m];
              pdf_peak_cb1  = new RooCBShape( Form("pdf_peak_cb1_%s",tag.c_str()), "", *mB, *mub, *sig_CB,  *alp_CB, *n_CB );

              RooAbsReal* sig2_CB = sigma2[*t][*h][*b][*m][*c][*a];
              pdf_peak_cb2  = new RooCBShape( Form("pdf_peak_cb2_%s",tag.c_str()), "", *mB, *mub, *sig2_CB, *alp_CB, *n_CB );
              
              pdf_peak      = new RooAddPdf(  Form("pdf_peak_%s",    tag.c_str()), "", RooArgSet(*pdf_peak_cb1,*pdf_peak_cb2),*sig_frac[*t][*h][*b][*m]);
              
              

              //comb bkg
              RooRealVar* comb = comb_slope[*t][*h][*b][*m];
              pdf_comb = new RooExponential( Form("pdf_comb_%s",tag.c_str()), "",*mB,*comb);

              if( (*b).c_str()==DD0){
                if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Adding DD0 mode... " <<std::endl;
                // D*D0 pi0
                RooAbsReal* HORNS_var_a     = HORNS_a[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_b     = HORNS_b[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_csi   = HORNS_csi[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_shift = HORNS_shift[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_sigma = HORNS_sigma[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_R     = HORNS_R[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_f     = HORNS_f[*t][*h][*b][*m];
                pdf_HORNS = new RooHORNSdini(Form("pdf_HORNS_%s",tag.c_str()), "", *mB, *HORNS_var_a, *HORNS_var_b, *HORNS_var_csi, *HORNS_var_shift, *HORNS_var_sigma, *HORNS_var_R, *HORNS_var_f);
            
                // D*D0 gamma
                RooAbsReal* HILL_var_a     = HILL_a[*t][*h][*b][*m];
                RooAbsReal* HILL_var_b     = HILL_b[*t][*h][*b][*m];
                RooAbsReal* HILL_var_csi   = HILL_csi[*t][*h][*b][*m];
                RooAbsReal* HILL_var_shift = HILL_shift[*t][*h][*b][*m];
                RooAbsReal* HILL_var_sigma = HILL_sigma[*t][*h][*b][*m];
                RooAbsReal* HILL_var_R     = HILL_R[*t][*h][*b][*m];
                RooAbsReal* HILL_var_f     = HILL_f[*t][*h][*b][*m];
                pdf_HILL = new RooHILLdini(Form("pdf_HILL_%s",tag.c_str()), "", *mB, *HILL_var_a, *HILL_var_b, *HILL_var_csi, *HILL_var_shift, *HILL_var_sigma, *HILL_var_R, *HILL_var_f);

                // D*D0 Shape
                pdf_DstD0 = new RooAddPdf(  Form("pdf_DstD0_%s",tag.c_str()), "", RooArgSet(*pdf_HORNS,*pdf_HILL),*fraction_Dst_Dpi0_DD0_adj);
              
                //DsD*0 pi0
                RooAbsReal* HORNS2_var_a     = HORNS2_a[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_b     = HORNS2_b[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_csi   = HORNS2_csi[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_shift = HORNS2_shift[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_sigma = HORNS2_sigma[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_R     = HORNS2_R[*t][*h][*b][*m];
                RooAbsReal* HORNS2_var_f     = HORNS2_f[*t][*h][*b][*m];
                pdf_HORNS2 = new RooHORNSdini(Form("pdf_HORNS2_%s",tag.c_str()), "", *mB, *HORNS2_var_a, *HORNS2_var_b, *HORNS2_var_csi, *HORNS2_var_shift, *HORNS2_var_sigma, *HORNS2_var_R, *HORNS2_var_f);
            
                //DD*0 gamma
                RooAbsReal* HILL2_var_a     = HILL2_a[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_b     = HILL2_b[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_csi   = HILL2_csi[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_shift = HILL2_shift[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_sigma = HILL2_sigma[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_R     = HILL2_R[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_f     = HILL2_f[*t][*h][*b][*m];
                pdf_HILL2 = new RooHILLdini(Form("pdf_HILL2_%s",tag.c_str()), "", *mB, *HILL2_var_a, *HILL2_var_b, *HILL2_var_csi, *HILL2_var_shift, *HILL2_var_sigma, *HILL2_var_R, *HILL2_var_f);
                
                // D*D0 Shape
                pdf_DDst0 = new RooAddPdf(  Form("pdf_DDst0_%s",tag.c_str()), "", RooArgSet(*pdf_HORNS2,*pdf_HILL2),*fraction_Dst0_D0pi0_DD0_adj);

                RooArgSet pdflist(  *pdf_comb
                                    ,*pdf_peak 
                                    ,*pdf_DstD0 
                                    ,*pdf_DDst0 
                                    );
                RooArgSet nevents(  *yield_comb[*t][*h][*b][*m][*c][*a] 
                                    ,*B_yield[*t][*h][*b][*m][*c][*a] 
                                    ,*yield_dstd0[*t][*h][*b][*m][*c][*a] 
                                    ,*yield_ddst0[*t][*h][*b][*m][*c][*a] 
                                    ); 

                //Add to master PDF
                RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_%s",tag.c_str())       ,"",pdflist, nevents);
                
                std::stringstream str;
                str<<(*t)<<underscore<<(*h)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                
                sim->addPdf(*masterPdf,str.str().c_str());
                //sim->ls();
                //sim->Print();
                if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Added to sim pdf: "<< str.str().c_str() <<std::endl;
                
              } else if((*b).c_str()==DKst0||(*b).c_str()==DKst0Side){
                if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Adding "<<(*b)<<" mode... " <<std::endl; 
                
                //RooHORNSdini
                RooAbsReal* HORNS_var_a     = HORNS_a[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_b     = HORNS_b[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_csi   = HORNS_csi[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_shift = HORNS_shift[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_sigma = HORNS_sigma[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_R     = HORNS_R[*t][*h][*b][*m];
                RooAbsReal* HORNS_var_f     = HORNS_f[*t][*h][*b][*m];
                pdf_HORNS = new RooHORNSdini(Form("pdf_HORNS_%s",tag.c_str()), "", *mB, *HORNS_var_a, *HORNS_var_b, *HORNS_var_csi, *HORNS_var_shift, *HORNS_var_sigma, *HORNS_var_R, *HORNS_var_f);
            
                //RooHILLdini
                RooAbsReal* HILL_var_a     = HILL_a[*t][*h][*b][*m];
                RooAbsReal* HILL_var_b     = HILL_b[*t][*h][*b][*m];
                RooAbsReal* HILL_var_csi   = HILL_csi[*t][*h][*b][*m];
                RooAbsReal* HILL_var_shift = HILL_shift[*t][*h][*b][*m];
                RooAbsReal* HILL_var_sigma = HILL_sigma[*t][*h][*b][*m];
                RooAbsReal* HILL_var_R     = HILL_R[*t][*h][*b][*m];
                RooAbsReal* HILL_var_f     = HILL_f[*t][*h][*b][*m];
                pdf_HILL = new RooHILLdini(Form("pdf_HILL_%s",tag.c_str()), "", *mB, *HILL_var_a, *HILL_var_b, *HILL_var_csi, *HILL_var_shift, *HILL_var_sigma, *HILL_var_R, *HILL_var_f);

                // RooHILLdini
                RooAbsReal* HILL2_var_a     = HILL2_a[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_b     = HILL2_b[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_csi   = HILL2_csi[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_shift = HILL2_shift[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_sigma = HILL2_sigma[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_R     = HILL2_R[*t][*h][*b][*m];
                RooAbsReal* HILL2_var_f     = HILL2_f[*t][*h][*b][*m];
                
                pdf_HILL2 = new RooHILLdini(Form("pdf_HILL2_%s",tag.c_str()), "", *mB, *HILL2_var_a, *HILL2_var_b, *HILL2_var_csi, *HILL2_var_shift, *HILL2_var_sigma, *HILL2_var_R, *HILL2_var_f);
                
                //RooLITTLEHORNSdini
                RooAbsReal* LITTLEHORNS_var_a     = LITTLEHORNS_a[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_b     = LITTLEHORNS_b[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_csi   = LITTLEHORNS_csi[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_shift = LITTLEHORNS_shift[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_sigma = LITTLEHORNS_sigma[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_R     = LITTLEHORNS_R[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_f     = LITTLEHORNS_f[*t][*h][*b][*m];
                RooAbsReal* LITTLEHORNS_var_g     = LITTLEHORNS_g[*t][*h][*b][*m];
                pdf_LITTLEHORNS = new RooLITTLEHORNSdini(Form("pdf_LITTLEHORNS_%s",tag.c_str()), "", *mB, *LITTLEHORNS_var_a, *LITTLEHORNS_var_b, *LITTLEHORNS_var_csi, *LITTLEHORNS_var_shift, *LITTLEHORNS_var_sigma, *LITTLEHORNS_var_R, *LITTLEHORNS_var_f, *LITTLEHORNS_var_g);

                //std::cout << "Making PartReco RooAddPdf" << std::endl;
                pdf_PartReco = new RooAddPdf(Form("pdf_PartReco_%s", tag.c_str()),
                                             "",
                                             RooArgSet(*pdf_HORNS,
                                                       *pdf_HILL,
                                                       *pdf_HILL2,
                                                       *pdf_LITTLEHORNS), 
                                             RooArgSet(*frac_HORNS[*t][*h][DKst0][*m][*c][*a],
                                                       *frac_HILL[*t][*h][DKst0][*m][*c][*a],
                                                       *frac_HILL2[*t][*h][DKst0][*m][*c][*a]) );
              
                //std::cout << "Made PartReco RooAddPdf" << std::endl;
                pdf_Bs2DKPiPi     = new RooFFTConvPdf(*PartRecoPDF_Bs2DKPiPi_Conv,     Form("pdf_Bs2DKPiPi_%s",tag.c_str())    );
                pdf_Bu2DstKPi     = new RooFFTConvPdf(*PartRecoPDF_Bu2DstKPi_Conv,     Form("pdf_Bu2DstKPi_%s",tag.c_str())    );
                pdf_B02Dst0Kst0   = new RooFFTConvPdf(*PartRecoPDF_B02Dst0Kst0_Conv,   Form("pdf_B02Dst0Kst0_%s",tag.c_str())  );
                pdf_Bs2DKst0pi    = new RooFFTConvPdf(*PartRecoPDF_Bs2DKst0pi_Conv,    Form("pdf_Bs2DKst0pi_%s",tag.c_str())   );
                pdf_Bs2DKst0pi_2  = new RooFFTConvPdf(*PartRecoPDF_Bs2DKst0pi_2_Conv,  Form("pdf_Bs2DKst0pi_2_%s",tag.c_str()) );

                //pdf_B2DKPi        = new RooAddPdf(    (RooAddPdf)*pdf_peak,            Form("pdf_B2DKPi_%s",tag.c_str())       );
                pdf_B2DKPi        = new RooAddPdf(  Form("pdf_B2DKPi_%s",    tag.c_str()), "", RooArgSet(*pdf_peak_cb1,*pdf_peak_cb2),*sig_frac[*t][*h][*b][*m]);

                RooArgSet pdflist( *pdf_peak
                                  ,*pdf_comb
                                  ,*pdf_PartReco
                                  //,*pdf_B2DKPi
                                  ,*pdf_Bs2DKPiPi
                                  ,*pdf_Bu2DstKPi
                                  ,*pdf_B02Dst0Kst0
                                  ,*pdf_Bs2DKst0pi
                                  ,*pdf_Bs2DKst0pi_2
                                   );


                                   // *pdf_DsstKKst );
                //std::cout << "Made pdf list" << std::endl;

                RooArgSet nevents(  *B_yield[*t][*h][*b][*m][*c][*a] 
                                   ,*yield_comb[*t][*h][*b][*m][*c][*a]
                                   ,*PR_total_yield[*t][*h][*b][*m][*c][*a]
                                   //,*DKPi_unblind_yield[*t][*h][*b][*m][*c][*a]
                                   ,*DKPiPi_total_yield[*t][*h][*b][*m][*c][*a]
                                   ,*DstKPi_total_yield[*t][*h][*b][*m][*c][*a] 
                                   ,*Dst0Kst0_total_yield[*t][*h][*b][*m][*c][*a]
                                   ,*DKst0pi_total_yield[*t][*h][*b][*m][*c][*a]
                                   ,*DKst0pi_2_total_yield[*t][*h][*b][*m][*c][*a]
                                    );
                                    // *yield_DsstKKst[*t][*h][*b][*m][*c][*a] );
                //std::cout << "Making Master PDF" << std::endl;
                RooAddPdf* masterPdf       = new RooAddPdf(Form("masterPdf_%s",tag.c_str()),"",pdflist, nevents);
                
                std::stringstream str;
                str<<(*t)<<underscore<<(*h)<<underscore<<(*b)<<underscore<<(*m)<<underscore<<(*c)<<underscore<<(*a);
                
                sim->addPdf(*masterPdf,str.str().c_str());
                //sim->ls();
                //sim->Print();
                if(par->debug) std::cout<<"Running: DKst0Model::DefineModel() --> Added to sim pdf: "<< str.str().c_str() <<std::endl;
                
              } 
            } //close c
          } // close a
        } //close m
      } //close b 
    } // close  h
  } // close t
}//end of funcn DefineModel()

RooArgSet* DKst0Model::GetParameters()
{
  if(par->debug) std::cout<<"Running: DKst0Model::GetParameters()"<<std::endl;
  RooArgSet* pars = sim->getParameters(RooArgSet(*mB,type,helBin,Bmode,mode,charge,magnet));
  return pars;
}



void DKst0Model::PrintResult()
{
  if(par->debug) std::cout<<"Running: DKst0Model::PrintResults()"<<std::endl;
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            if((*b)==DKst0||(*b)==DKst0Side){
              if (needsBlinding) std::cout << "DKst0 Mode is BLIND: Printing blinded yield" << std::endl;
              float y_peak =  B_yield[*t][*h][*b][*m][both][*a]->getVal();
              //float y_yp   =  yield_peak[*t][*h][*b][*m][both][*a]->getVal();
              float y_comb =  yield_comb[*t][*h][*b][*m][both][*a]->getVal();
              float y_partreco =   PR_total_yield[*t][*h][*b][*m][both][*a]->getVal();
              
              std::cout<< "\n---------------------------"<<std::endl;
              std::cout<<(*b)<<", "<<(*m)<<", magnet: "<<(*a)<<", year: "<<(*t)<<"HelBin: " << (*h)<<std::endl;
              //std::cout<< " B Yield: "<<y_peak<<std::endl;
              //std::cout<< " Yield peak: "<<y_yp<<std::endl;
              std::cout<< " Comb:    "<<y_comb<<std::endl;
              std::cout<< " D*Kst0:  "<<y_partreco<<std::endl;
              std::cout<< "---------------------------\n"<<std::endl;
              //B_yield[*t][*h][*b][*m][both][*a]->printValue();
            }  else {
              if(par->sumOverCharges){
                float y_peak   =  B_yield[*t][*h][*b][*m][both][*a]->getVal();
                //float y_yp   =  yield_peak[*t][*h][*b][*m][both][*a]->getVal();
                float y_comb   =  yield_comb[*t][*h][*b][*m][both][*a]->getVal();
                float y_ddst0  =  yield_ddst0[*t][*h][*b][*m][both][*a]->getVal();
                float y_dstd0  =  yield_dstd0[*t][*h][*b][*m][both][*a]->getVal();
                
                std::cout<< "\n---------------------------"<<std::endl;
                std::cout<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<", year-"<<(*t)<<" :"<<", HelBin: " << (*h)<<std::endl;
                std::cout<< " B Yield: "<<y_peak<<std::endl;
                //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                std::cout<< " Comb:   "<<y_comb<<std::endl;
                std::cout<< " DD*0:  "<<y_ddst0<<std::endl;
                std::cout<< " D*D0:  "<<y_dstd0<<std::endl;
                std::cout<< "---------------------------\n"<<std::endl;
              }else{
                float y_peak_minus   =  yield_peak[*t][*h][*b][*m][minus][*a]->getVal();
                float y_comb_minus   =  yield_comb[*t][*h][*b][*m][minus][*a]->getVal();
                float y_dstd0_minus  =  yield_dstd0[*t][*h][*b][*m][minus][*a]->getVal();
                float y_ddst0_minus  =  yield_ddst0[*t][*h][*b][*m][minus][*a]->getVal();
                float y_peak_plus    =  yield_peak[*t][*h][*b][*m][plus][*a]->getVal();
                float y_comb_plus    =  yield_comb[*t][*h][*b][*m][plus][*a]->getVal();
                float y_dstd0_plus   =  yield_dstd0[*t][*h][*b][*m][plus][*a]->getVal();
                float y_ddst0_plus   =  yield_ddst0[*t][*h][*b][*m][plus][*a]->getVal();

             
                std::cout<<(*t)<<", "<<(*h)<<", "<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
                std::cout<<" Peak minus: "<<y_peak_minus<<" comb. minus: "<<y_comb_minus<<" XX^{*} minus: "<<y_dstd0_minus<<" X^{*}X: minus: "<<y_ddst0_minus<<std::endl;
                std::cout<<" Peak plus:  "<<y_peak_plus <<" comb. plus:  "<<y_comb_plus <<" XX^{*} plus:  "<<y_dstd0_plus <<" X^{*}X: plus:  "<<y_ddst0_plus <<std::endl;
                std::cout<<"--------------------------------------------------"<<std::endl;
                std::cout<<" Peak total: "<<y_peak_minus+y_peak_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<" XX^{*} total: "<<y_dstd0_minus+y_dstd0_plus<<" X^{*}X: total: "<<y_ddst0_minus+y_dstd0_plus<<std::endl;

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

  

  std::cout<<std::endl;
}

std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,double> > > > > > > DKst0Model::GetResult()
{
  if(par->debug) std::cout<<"Running: DKst0Model::GetResults()"<<std::endl;
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              if (needsBlinding && ((*b)==DKst0||(*b)==DKst0Side) ) {
                std::cout << "DKst0 mode is BLIND, therefore no values returned " << std::endl;
              } else {


                std::string cat_str = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*b).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());
                std::cout<< " Category: " << cat_str << std::endl;
                fit_results[*t][*h][*b][*m][*c][*a]["Peak"]    = yield_peak[*t][*h][*b][*m][*c][*a]->getVal();
                //fit_results[*t][*h][*b][*m][*c][*a]["Peakerr"] = yield_peak[*t][*h][*b][*m][*c][*a]->getError();
                fit_results[*t][*h][*b][*m][*c][*a]["Comb"]    = yield_comb[*t][*h][*b][*m][*c][*a]->getVal();
                fit_results[*t][*h][*b][*m][*c][*a]["XXst"]    = yield_ddst0[*t][*h][*b][*m][*c][*a]->getVal();
                fit_results[*t][*h][*b][*m][*c][*a]["XstX"]    = yield_dstd0[*t][*h][*b][*m][*c][*a]->getVal();
                std::cout<< " Peak:   " <<fit_results[*t][*h][*b][*m][*c][*a]["Peak"]<<std::endl;
                std::cout<< " Comb:   " <<fit_results[*t][*h][*b][*m][*c][*a]["Comb"]<<std::endl;
                std::cout<< " XXst:   " <<fit_results[*t][*h][*b][*m][*c][*a]["XXst"]<<std::endl;
                std::cout<< " XstX:   " <<fit_results[*t][*h][*b][*m][*c][*a]["XstX"]<<std::endl;
                
                mB->setRange("signal", 5233.5, 5333.5) ;
                //RooAbsPdf* sum = sim->getPdf(cat_str.c_str());
                //RooAbsPdf* comb = sim->getPdf(Form("pdf_comb_%s",cat_str.c_str()));
                //RooAbsPdf* comb = (RooAbsPdf*) sum->FindObject(Form("pdf_comb_%s",cat_str.c_str()));
                //RooAbsReal* ibkg_sig = sum->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"),RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())));  
                double combslope = comb_slope[*t][*h][*b][*m]->getVal();
                RooRealVar* fixed_slope = new RooRealVar("fixed_slope","slope",  combslope);
                RooAbsPdf* pdf_comb = new RooExponential( "Exp", "",*mB,*fixed_slope);
                RooAbsReal* ibkg_sig = pdf_comb->createIntegral(*mB,RooFit::NormSet(*mB),RooFit::Range("signal"));  
                fit_results[*t][*h][*b][*m][*c][*a]["Comb_Reduced"] = (ibkg_sig->getVal()*yield_comb[*t][*h][*b][*m][*c][*a]->getVal());
                

                std::cout<< " Comb_Reduced:   " <<fit_results[*t][*h][*b][*m][*c][*a]["Comb_Reduced"]<<std::endl;
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
