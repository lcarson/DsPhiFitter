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

DsPhiModel::DsPhiModel(Parameters* p, RooRealVar* pmB, bool nB, bool vA, bool gen)
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
, varyAllowed(vA)
, genModel(gen)
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
                  
                  //if(*b == DsD0){
                  //  mB->setRange(Form("fitRange_%s",str.str().c_str()),5060,5900);
                  //  std::cout <<  "   --> Set fit range : 5060,5900" << std::endl;
                  //} else {   
                  //  mB->setRange(Form("fitRange_%s",str.str().c_str()),4820,5900);
                  //  std::cout <<  "   --> Set fit range : 48200,5900" << std::endl;
                  //}
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
  
  // ======================================================== 
  // =========== Switches for quick/significance ============
  // ======================================================== 
  
  // Make fit run super quick --> for testing
  bool allConst_exYield =true;
  
  // Fix yields to zero --> to calculate sigma wilks 
  
  bool NoBR = false;
  
  bool NoYield_all = false;  
  std::map<std::string,bool> NoYield;
  NoYield[Ds2PhiPi]  = false; 
  NoYield[Ds2KKPi]   = false; 
  NoYield[Ds2KPiPi]  = false; 
  NoYield[Ds2PiPiPi] = false; 

  if(NoYield_all){
    NoYield[Ds2PhiPi]  = true; 
    NoYield[Ds2KKPi]   = true; 
    NoYield[Ds2KPiPi]  = true; 
    NoYield[Ds2PiPiPi] = true; 
  }

  // Mean
  double default_mB=5279.29;

  RooRealVar* global_mean = new RooRealVar("global_mean","global_mean", default_mB, 5270, 5290);


  // ======================================================== 
  // =========== import DsKK*0 PDF       ====================
  // ======================================================== 
  
  std::string locationRooKeys = "pdfs/";
  //TFile f1((locationRooKeys + "kest1_Bs02Dsa1_KKst.root").c_str(), "read");
  
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Importing Dsa1 shape"<<std::endl;
  TFile f1((locationRooKeys + "kest1_Bs2Dsa1_4600_5900.root").c_str(), "read"); 
  RooWorkspace* workspace1 = (RooWorkspace*)f1.Get("workspace_Bs2Dsa1");
  if (!workspace1) std::cout<<"Running: DsPhiModel::DefineModel() --> Can't find Dsa1 workspace"<<std::endl; 
  RooKeysPdf* PartRecoPDF_Bs02Dsa1     = (RooKeysPdf*)workspace1->pdf("kest1_Bs2Dsa1");
  if(!PartRecoPDF_Bs02Dsa1) std::cout<<"Running: DsPhiModel::DefineModel() --> Can't find Dsa1 pdf"<<std::endl; 
  f1.Close();
  
  RooRealVar* global_shift      = new RooRealVar("global_shift",      "Global shift",      -1.81, -10.0,    10.0);
  RooRealVar* global_shift_test = new RooRealVar("global_shift_test", "Global shift test", -3.01, -100.0,    100.0);

  RooRealVar* sg_test      = new RooRealVar("sg_test",   "sg_test",    1  ) ;
  RooGaussian* gauss_test  = new RooGaussian("gauss_test",  "gauss_test",  *mB,*global_shift,*sg_test  ) ;
  RooFFTConvPdf* PartRecoPDF_Bs02Dsa1_Conv   = new RooFFTConvPdf("kest1_Dasa1_conv",  "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1, *gauss_test  );
  
  // ======================================================== 
  // =========== import DsKK*0 PDF       ====================
  // ======================================================== 
  
  TFile f_in_H1( (locationRooKeys + "kest1_Bs2Dsa1_4600_5900_in_H1.root" ).c_str(), "read"); 
  TFile f_in_H2( (locationRooKeys + "kest1_Bs2Dsa1_4600_5900_in_H2.root" ).c_str(), "read"); 
  TFile f_out_H1((locationRooKeys + "kest1_Bs2Dsa1_4600_5900_out_H1.root").c_str(), "read"); 
  TFile f_out_H2((locationRooKeys + "kest1_Bs2Dsa1_4600_5900_out_H2.root").c_str(), "read"); 

  RooWorkspace* workspace_in_H1  = (RooWorkspace*)f_in_H1.Get("workspace_Bs2Dsa1");
  RooWorkspace* workspace_in_H2  = (RooWorkspace*)f_in_H2.Get("workspace_Bs2Dsa1");
  RooWorkspace* workspace_out_H1 = (RooWorkspace*)f_out_H1.Get("workspace_Bs2Dsa1");
  RooWorkspace* workspace_out_H2 = (RooWorkspace*)f_out_H2.Get("workspace_Bs2Dsa1");

  RooKeysPdf* PartRecoPDF_Bs02Dsa1_in_H1     = (RooKeysPdf*)workspace_in_H1->pdf("kest1_Bs2Dsa1");
  RooKeysPdf* PartRecoPDF_Bs02Dsa1_in_H2     = (RooKeysPdf*)workspace_in_H2->pdf("kest1_Bs2Dsa1");
  RooKeysPdf* PartRecoPDF_Bs02Dsa1_out_H1    = (RooKeysPdf*)workspace_out_H1->pdf("kest1_Bs2Dsa1");
  RooKeysPdf* PartRecoPDF_Bs02Dsa1_out_H2    = (RooKeysPdf*)workspace_out_H2->pdf("kest1_Bs2Dsa1");

  f_in_H1.Close();
  f_in_H2.Close();
  f_out_H1.Close();
  f_out_H2.Close();

  std::map<std::string,std::map<std::string,RooFFTConvPdf*>> PartRecoPDF_Bs02Dsa1_Conv_split;

  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhi][Helbin1]     = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhi_Helbin1",      "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1_in_H1, *gauss_test  );
  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhi][Helbin2]     = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhi_Helbin2",      "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1_in_H2, *gauss_test  );
  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhiSide][Helbin1] = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhiSide_Helbin1",  "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1_out_H1, *gauss_test  );
  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhiSide][Helbin2] = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhiSide_Helbin2",  "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1_out_H2, *gauss_test  );

  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhi][both]     = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhi_both",      "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1, *gauss_test  );
  PartRecoPDF_Bs02Dsa1_Conv_split[DsPhiSide][both] = new RooFFTConvPdf("kest1_Dasa1_conv_DsPhiSide_both",  "kest1_Dsa1_conv",  *mB,*PartRecoPDF_Bs02Dsa1, *gauss_test  );




  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Importing DsstKKst shape"<<std::endl;
  //TFile f2((locationRooKeys + "kest1_Bs02DsstKKst.root").c_str(), "read");
  TFile f2((locationRooKeys + "kest1_Bs2DsstKKst_4600_5900.root").c_str(), "read");
  RooWorkspace* workspace2 = (RooWorkspace*)f2.Get("workspace_Bs2DsstKKst");
  RooKeysPdf* PartRecoPDF_Bs02DsstKKst  = (RooKeysPdf*)workspace2->pdf("kest1_Bs2DsstKKst");
  f2.Close();
  RooFFTConvPdf* PartRecoPDF_Bs02DsstKKst_Conv   = new RooFFTConvPdf("kest1_DsstKKst_conv",  "kest1_DsstKKst_conv",  *mB,*PartRecoPDF_Bs02DsstKKst, *gauss_test  );

  //TFile f3((locationRooKeys + "kest1_Bs2DsDs_4800_5900.root").c_str(), "read");
  //RooWorkspace* workspace3 = (RooWorkspace*)f3.Get("workspace");
  //RooKeysPdf* PartRecoPDF_Bs02DsDs  = (RooKeysPdf*)workspace3->pdf("kest1");
  //f3.Close();

  RooRealVar* frac_Dsa1_DsstKKst = new RooRealVar("frac_Dsa1_DsstKKst","frac_Dsa1_DsstKKst", 0.774, 0.0, 1.0);
  if(allConst_exYield) frac_Dsa1_DsstKKst->setConstant();

  
  // Import other shapes 

  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Importing B02DsD shape"<<std::endl;
  
  TFile f3((locationRooKeys + "kest1_B2DsD_4800_5900.root").c_str(), "read");
  RooWorkspace* workspace3 = (RooWorkspace*)f3.Get("workspace_B2DsD");
  RooKeysPdf* PartRecoPDF_B2DsD  = (RooKeysPdf*)workspace3->pdf("kest1_B2DsD");
  f3.Close();
  RooFFTConvPdf* PartRecoPDF_B2DsD_Conv   = new RooFFTConvPdf("kest1_B2DsD_conv",  "kest1_B2DsD_conv",  *mB,*PartRecoPDF_B2DsD, *gauss_test  );


  // Import DD' shapes

  TFile f4((locationRooKeys + "kest1_Bs2DsDs_4600_5900.root").c_str(), "read");
  RooWorkspace* workspace4 = (RooWorkspace*)f4.Get("workspace_Bs2DsDs");
  RooKeysPdf* PartRecoPDF_Bs2DsDs  = (RooKeysPdf*)workspace4->pdf("kest1_Bs2DsDs");
  f4.Close();
  RooFFTConvPdf* PartRecoPDF_Bs2DsDs_Conv   = new RooFFTConvPdf("kest1_Bs2DsDs_conv",  "kest1_Bs2DsDs_conv",  *mB,*PartRecoPDF_Bs2DsDs, *gauss_test  );

  TFile f5((locationRooKeys + "kest1_Bs2DsstDs_4600_5900.root").c_str(), "read");
  RooWorkspace* workspace5 = (RooWorkspace*)f5.Get("workspace_Bs2DsstDs");
  RooKeysPdf* PartRecoPDF_Bs2DsstDs  = (RooKeysPdf*)workspace5->pdf("kest1_Bs2DsstDs");
  f5.Close();
  RooFFTConvPdf* PartRecoPDF_Bs2DsstDs_Conv   = new RooFFTConvPdf("kest1_Bs2DsstDs_conv",  "kest1_Bs2DsstDs_conv",  *mB,*PartRecoPDF_Bs2DsstDs, *gauss_test  );


  // Shift DsDs to make DsD

  RooRealVar* offset                       = new RooRealVar("offset","offset",40.0 ) ; 
  RooFormulaVar* global_shift_minus_offset = new RooFormulaVar("global_shift_minus_offset", "" , "@0-@1",RooArgList(*global_shift,*offset));
  RooGaussian* gauss_shift_test            = new RooGaussian("gauss_shift_test",  "gauss_shift_test",  *mB,*global_shift_minus_offset,*sg_test  ) ;

  RooFFTConvPdf* PartRecoPDF_B02DsD_Conv   = new RooFFTConvPdf("kest1_B02DsD_conv",  "kest1_B02DsD_conv",  *mB,*PartRecoPDF_Bs2DsDs, *gauss_shift_test  );

  //RooRealVar* frac_Dsa1_DsstKKst = new RooRealVar("frac_Dsa1_DsstKKst","frac_Dsa1_DsstKKst",0.997633244);
  // Convolve Dsa1 shapes with gaussians to take into account MC and Data resolutions differences
  // Sigma values taken from sigma = sqrt((sigma_Data)^2 - (sigma_MC)^2) 
  // Values taken from smallest sigma (sigma1) in DsD0 mode 
  // Assumed MC/Data ratio is the same in DsPhi and DsD0

  RooRealVar* mg             = new RooRealVar("mg","mg",0 + (varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,1.0):0) ) ;
  RooRealVar* sg_KKPi        = new RooRealVar("sg_KKPi",   "sg_KKPi",    0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg_KPiPi       = new RooRealVar("sg_KPiPi",  "sg_KPiPi",   0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg_PiPiPi      = new RooRealVar("sg_PiPiPi", "sg_PiPiPi",  0 + fabs(varyAllowed&&par->variation[fixedBG_Dsa1_smear]?rand->Gaus(0,10.0):0)) ;
  RooGaussian* gauss_KKPi    = new RooGaussian("gauss_KKPi",  "gauss_KKPi",  *mB,*mg,*sg_KKPi  ) ;
  RooGaussian* gauss_KPiPi   = new RooGaussian("gauss_KPiPi", "gauss_KPiPi", *mB,*mg,*sg_KPiPi ) ;
  RooGaussian* gauss_PiPiPi  = new RooGaussian("gauss_PiPiPi","gauss_PiPiPi",*mB,*mg,*sg_PiPiPi) ;

  RooRealVar* mg2            = new RooRealVar("mg2","mg2",0 + (varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,1.0):0) ) ;
  RooRealVar* sg2_KKPi       = new RooRealVar("sg2_KKPi",   "sg2_KKPi",    0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg2_KPiPi      = new RooRealVar("sg2_KPiPi",  "sg2_KPiPi",   0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooRealVar* sg2_PiPiPi     = new RooRealVar("sg2_PiPiPi", "sg2_PiPiPi",  0 + fabs(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]?rand->Gaus(0,10.0):0)) ;
  RooGaussian* gauss2_KKPi   = new RooGaussian("gauss2_KKPi",  "gauss2_KKPi",  *mB,*mg2,*sg2_KKPi  ) ;
  RooGaussian* gauss2_KPiPi  = new RooGaussian("gauss2_KPiPi", "gauss2_KPiPi", *mB,*mg2,*sg2_KPiPi ) ;
  RooGaussian* gauss2_PiPiPi = new RooGaussian("gauss2_PiPiPi","gauss2_PiPiPi",*mB,*mg2,*sg2_PiPiPi) ;

  std::map<std::string,RooFFTConvPdf*> PartRecoPDF_Bs02Dsa1_Convolved;
  PartRecoPDF_Bs02Dsa1_Convolved[Ds2KKPi]   = new RooFFTConvPdf("kest1_adjust_KKPi",  "kest1_adjust_KKPi",  *mB,*PartRecoPDF_Bs02Dsa1_Conv, *gauss_KKPi  );
  PartRecoPDF_Bs02Dsa1_Convolved[Ds2PhiPi]  = new RooFFTConvPdf("kest1_adjust_PhiPi", "kest1_adjust_PhiPi", *mB,*PartRecoPDF_Bs02Dsa1_Conv, *gauss_KKPi  );   
  PartRecoPDF_Bs02Dsa1_Convolved[Ds2KPiPi]  = new RooFFTConvPdf("kest1_adjust_KPiPi", "kest1_adjust_KPiPi", *mB,*PartRecoPDF_Bs02Dsa1_Conv, *gauss_KPiPi );
  PartRecoPDF_Bs02Dsa1_Convolved[Ds2PiPiPi] = new RooFFTConvPdf("kest1_adjust_PiPiPi","kest1_adjust_PiPiPi",*mB,*PartRecoPDF_Bs02Dsa1_Conv, *gauss_PiPiPi);


  std::map<std::string,RooFFTConvPdf*> PartRecoPDF_Bs02DsstKKst_Convolved;
  PartRecoPDF_Bs02DsstKKst_Convolved[Ds2KKPi]   = new RooFFTConvPdf("kest1_adjust_KKPi",  "kest1_adjust_KKPi",  *mB,*PartRecoPDF_Bs02DsstKKst_Conv, *gauss2_KKPi  );
  PartRecoPDF_Bs02DsstKKst_Convolved[Ds2PhiPi]  = new RooFFTConvPdf("kest1_adjust_PhiPi", "kest1_adjust_PhiPi", *mB,*PartRecoPDF_Bs02DsstKKst_Conv, *gauss2_KKPi  );   
  PartRecoPDF_Bs02DsstKKst_Convolved[Ds2KPiPi]  = new RooFFTConvPdf("kest1_adjust_KPiPi", "kest1_adjust_KPiPi", *mB,*PartRecoPDF_Bs02DsstKKst_Conv, *gauss2_KPiPi );
  PartRecoPDF_Bs02DsstKKst_Convolved[Ds2PiPiPi] = new RooFFTConvPdf("kest1_adjust_PiPiPi","kest1_adjust_PiPiPi",*mB,*PartRecoPDF_Bs02DsstKKst_Conv, *gauss2_PiPiPi);

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
  
  Yield_comb[DsD0][Ds2KKPi][both]   = 154;
  Yield_comb[DsD0][Ds2KPiPi][both]  = 58;
  Yield_comb[DsD0][Ds2PhiPi][both]  = 56; 
  Yield_comb[DsD0][Ds2PiPiPi][both] = 124;
 
  Yield_comb[DsPhiSide][Ds2KKPi][both]   = 86;
  Yield_comb[DsPhiSide][Ds2KPiPi][both]  = 120;
  Yield_comb[DsPhiSide][Ds2PhiPi][both]  = 54;
  Yield_comb[DsPhiSide][Ds2PiPiPi][both] = 111;
   
  Yield_comb[DsPhi][Ds2KKPi][both]   = 203;
  Yield_comb[DsPhi][Ds2KPiPi][both]  = 41;
  Yield_comb[DsPhi][Ds2PhiPi][both]  = 121;
  Yield_comb[DsPhi][Ds2PiPiPi][both] = 126;


  Yield_comb[DsD0][Ds2KKPi][Helbin1]   = 229;
  Yield_comb[DsD0][Ds2KKPi][Helbin2]   = 83;
  Yield_comb[DsD0][Ds2KPiPi][Helbin1]  = 83;
  Yield_comb[DsD0][Ds2KPiPi][Helbin2]  = 35;
  Yield_comb[DsD0][Ds2PhiPi][Helbin1]  = 75; 
  Yield_comb[DsD0][Ds2PhiPi][Helbin2]  = 37; 
  Yield_comb[DsD0][Ds2PiPiPi][Helbin1] = 129;
  Yield_comb[DsD0][Ds2PiPiPi][Helbin2] = 116;
  
  Yield_comb[DsPhiSide][Ds2KKPi][Helbin1]   =  56;
  Yield_comb[DsPhiSide][Ds2KKPi][Helbin2]   =  57;
  Yield_comb[DsPhiSide][Ds2KPiPi][Helbin1]  = 101;
  Yield_comb[DsPhiSide][Ds2KPiPi][Helbin2]  =  58;
  Yield_comb[DsPhiSide][Ds2PhiPi][Helbin1]  =  16; 
  Yield_comb[DsPhiSide][Ds2PhiPi][Helbin2]  =  37;
  Yield_comb[DsPhiSide][Ds2PiPiPi][Helbin1] =  61;
  Yield_comb[DsPhiSide][Ds2PiPiPi][Helbin2] =  54;

  Yield_comb[DsPhi][Ds2KKPi][Helbin1]   = 100;
  Yield_comb[DsPhi][Ds2KKPi][Helbin2]   =  76;
  Yield_comb[DsPhi][Ds2KPiPi][Helbin1]  =  12;
  Yield_comb[DsPhi][Ds2KPiPi][Helbin2]  =  32;
  Yield_comb[DsPhi][Ds2PhiPi][Helbin1]  =  40;
  Yield_comb[DsPhi][Ds2PhiPi][Helbin2]  =  62;
  Yield_comb[DsPhi][Ds2PiPiPi][Helbin1] =  46;
  Yield_comb[DsPhi][Ds2PiPiPi][Helbin2] =  28;


  RooRealVar* global_comb_slope = new RooRealVar("global_comb_slope"  ,"Global Comb. Slope", -0.00351, -1.0, -0.0000000001);
  //if(allConst_exYield) global_comb_slope->setConstant();
  
  std::map<std::string,std::map<std::string,RooRealVar*>> comb_slope_sys;
  comb_slope_sys[DsD0][Ds2PhiPi]  = new RooRealVar("comb_slope_DsD0_Ds2PhiPi"   ,"", -0.00286, -1.0, -0.0000000001);
  comb_slope_sys[DsD0][Ds2KKPi]   = new RooRealVar("comb_slope_DsD0_Ds2KKPi"    ,"", -0.00286, -1.0, -0.0000000001);
  comb_slope_sys[DsD0][Ds2KPiPi]  = new RooRealVar("comb_slope_DsD0_Ds2KPiPi"   ,"", -0.00286, -1.0, -0.0000000001);
  comb_slope_sys[DsD0][Ds2PiPiPi] = new RooRealVar("comb_slope_DsD0_Ds2PiPiPi"  ,"", -0.00286, -1.0, -0.0000000001);

  comb_slope_sys[DsPhi][Ds2PhiPi]  = comb_slope_sys[DsD0][Ds2PhiPi] ;
  comb_slope_sys[DsPhi][Ds2KKPi]   = comb_slope_sys[DsD0][Ds2KKPi]  ;
  comb_slope_sys[DsPhi][Ds2KPiPi]  = comb_slope_sys[DsD0][Ds2KPiPi] ;
  comb_slope_sys[DsPhi][Ds2PiPiPi] = comb_slope_sys[DsD0][Ds2PiPiPi];

  comb_slope_sys[DsPhiSide][Ds2PhiPi]  = comb_slope_sys[DsPhi][Ds2PhiPi] ;
  comb_slope_sys[DsPhiSide][Ds2KKPi]   = comb_slope_sys[DsPhi][Ds2KKPi]  ;
  comb_slope_sys[DsPhiSide][Ds2KPiPi]  = comb_slope_sys[DsPhi][Ds2KPiPi] ;
  comb_slope_sys[DsPhiSide][Ds2PiPiPi] = comb_slope_sys[DsPhi][Ds2PiPiPi];



  std::map<std::string,std::map<std::string,RooRealVar*>> comb_slope_test;
  comb_slope_test[DsD0][Helbin1]   = new RooRealVar("comb_slope_DsD0"   ,"", -0.00286, -1.0, -0.0000000001);
  comb_slope_test[DsD0][Helbin2]   = comb_slope_test[DsD0][Helbin1];
  comb_slope_test[DsD0][both]      = comb_slope_test[DsD0][Helbin1];


  comb_slope_test[DsPhi][both]      = new RooRealVar("comb_slope_DsPhi"   ,   "", -0.00143, -1.0, -0.0000000001);
  comb_slope_test[DsPhiSide][both]  = new RooRealVar("comb_slope_DsPhiSide",  "", -0.00183, -1.0, -0.0000000001);

  comb_slope_test[DsPhi][Helbin1]      = new RooRealVar("comb_slope_DsPhi_H1"   ,   "", -0.001025, -1.0, -0.0000000001);
  comb_slope_test[DsPhi][Helbin2]      = new RooRealVar("comb_slope_DsPhi_H2"  ,    "", -0.002676, -1.0, -0.0000000001);
  comb_slope_test[DsPhiSide][Helbin1]  = new RooRealVar("comb_slope_DsPhiSide_H1",  "", -0.001812, -1.0, -0.0000000001);
  comb_slope_test[DsPhiSide][Helbin2]  = new RooRealVar("comb_slope_DsPhiSide_H2",  "", -0.001882, -1.0, -0.0000000001);

  RooRealVar* Comb_Helbin1_Fraction   = new RooRealVar("Comb_Helbin1_Fraction",  ""  , 0.6);

  RooRealVar* comb_flat_fraction = new RooRealVar("comb_flat_fraction", "" , 0.5, 0.0, 1.0);
  RooRealVar* comb_line_offset   = new RooRealVar("comb_line_offset",    "" , 0.59, 0.0,100);
  RooRealVar* comb_line_slope    = new RooRealVar("comb_line_slope",     "" , -0.000099,-0.001,0.0) ;
  //RooRealVar* comb_line_offset = new RooRealVar("a0","a0",0.4,-1,1);
  //RooRealVar* comb_line_slope  = new RooRealVar("a1","a1",-0.004,-1,1);

  //comb_slope_sys[DsPhi][Ds2PhiPi]  = new RooRealVar("comb_slope_DsPhi_Ds2PhiPi"   ,"", -0.00258, -1.0, -0.0000000001);
  //comb_slope_sys[DsPhi][Ds2KKPi]   = new RooRealVar("comb_slope_DsPhi_Ds2KKPi"    ,"", -0.00258, -1.0, -0.0000000001);
  //comb_slope_sys[DsPhi][Ds2KPiPi]  = new RooRealVar("comb_slope_DsPhi_Ds2KPiPi"   ,"", -0.00258, -1.0, -0.0000000001);
  //comb_slope_sys[DsPhi][Ds2PiPiPi] = new RooRealVar("comb_slope_DsPhi_Ds2PiPiPi"  ,"", -0.00258, -1.0, -0.0000000001);


  // Set fraction of combinatoric in the two Hel bins to be the same for all modes 
  //RooRealVar* splitHel_comb_fraction = new RooRealVar("splitHel_comb_fraction", "", 0.60, 0.0 ,1.0 );


  // ======================================================== 
  // =========== Values for Signal  =========================
  // ======================================================== 

  std::map<std::string,std::map<std::string,double>> Yield_CB_input, Fixed_CB_n, Fixed_CB_alpha, Fixed_CB_Sigma_ratio, Fixed_CB_fraction;
 
  Yield_CB_input[DsD0][Ds2KKPi]   = 1363;
  Yield_CB_input[DsD0][Ds2KPiPi]  = 166;
  Yield_CB_input[DsD0][Ds2PhiPi]  = 848; 
  Yield_CB_input[DsD0][Ds2PiPiPi] = 386;
 
  Yield_CB_input[DsPhi][Ds2KKPi]   = 45.8;
  Yield_CB_input[DsPhi][Ds2KPiPi]  = 3.8;
  Yield_CB_input[DsPhi][Ds2PhiPi]  = 28.8;
  Yield_CB_input[DsPhi][Ds2PiPiPi] = 11.0;


  double factor = 1.0;
  if(par->variation[fixedSig_double]) factor = 2.0;
  // CB n -> Fixed to 1  
  Fixed_CB_n[DsD0][Ds2PhiPi]       = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsD0][Ds2KKPi]        = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsD0][Ds2KPiPi]       = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsD0][Ds2PiPiPi]      = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  Fixed_CB_n[DsPhi][Ds2PhiPi]      = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhi][Ds2KKPi]       = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhi][Ds2KPiPi]      = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhi][Ds2PiPiPi]     = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  Fixed_CB_n[DsPhiSide][Ds2PhiPi]  = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhiSide][Ds2KKPi]   = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhiSide][Ds2KPiPi]  = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0); 
  Fixed_CB_n[DsPhiSide][Ds2PiPiPi] = 1.0 + (varyAllowed&&par->variation[fixedSig_n]?rand->Gaus(0,0.1*factor):0);

  // CB alpha -> From MC 
  Fixed_CB_alpha[DsD0][Ds2PhiPi]       = 2.91078 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0556065*factor):0);
  Fixed_CB_alpha[DsD0][Ds2KKPi]        = 2.91078 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0556065*factor):0);
  Fixed_CB_alpha[DsD0][Ds2KPiPi]       = 3.36188 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.260552*factor):0);
  Fixed_CB_alpha[DsD0][Ds2PiPiPi]      = 3.53538 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.254103*factor):0);

  Fixed_CB_alpha[DsPhi][Ds2PhiPi]      = 2.7604  + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0693368*factor):0);
  Fixed_CB_alpha[DsPhi][Ds2KKPi]       = 2.7604  + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0693368*factor):0);
  Fixed_CB_alpha[DsPhi][Ds2KPiPi]      = 3.06495 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.163778*factor):0);
  Fixed_CB_alpha[DsPhi][Ds2PiPiPi]     = 3.71842 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.235817*factor):0);

  Fixed_CB_alpha[DsPhiSide][Ds2PhiPi]  = 2.7604  + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0693368*factor):0);
  Fixed_CB_alpha[DsPhiSide][Ds2KKPi]   = 2.7604  + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.0693368*factor):0);
  Fixed_CB_alpha[DsPhiSide][Ds2KPiPi]  = 3.06495 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.163778*factor):0);
  Fixed_CB_alpha[DsPhiSide][Ds2PiPiPi] = 3.71842 + (varyAllowed&&par->variation[fixedSig_alpha]?rand->Gaus(0,0.26*factor):0);

  // 2 CB sigma ratios -> From MC 
  Fixed_CB_Sigma_ratio[DsD0][Ds2PhiPi]       = 0.427848 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0098877):0);
  Fixed_CB_Sigma_ratio[DsD0][Ds2KKPi]        = 0.427848 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0098877):0);
  Fixed_CB_Sigma_ratio[DsD0][Ds2KPiPi]       = 0.418909 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0119772):0);
  Fixed_CB_Sigma_ratio[DsD0][Ds2PiPiPi]      = 0.403931 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0081709):0);

  Fixed_CB_Sigma_ratio[DsPhi][Ds2PhiPi]      = 0.492471 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0103285):0);
  Fixed_CB_Sigma_ratio[DsPhi][Ds2KKPi]       = 0.492471 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0103285):0);
  Fixed_CB_Sigma_ratio[DsPhi][Ds2KPiPi]      = 0.467814 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0130645):0);
  Fixed_CB_Sigma_ratio[DsPhi][Ds2PiPiPi]     = 0.462447 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.00826237):0);

  Fixed_CB_Sigma_ratio[DsPhiSide][Ds2PhiPi]  = 0.492471 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0103285):0);
  Fixed_CB_Sigma_ratio[DsPhiSide][Ds2KKPi]   = 0.492471 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0103285):0);
  Fixed_CB_Sigma_ratio[DsPhiSide][Ds2KPiPi]  = 0.467814 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.0130645):0);
  Fixed_CB_Sigma_ratio[DsPhiSide][Ds2PiPiPi] = 0.462447 + (varyAllowed&&par->variation[fixedSig_Sigmaratio]?rand->Gaus(0,0.00826237):0);

  // 2 CB sigma fractions -> From MC 
  Fixed_CB_fraction[DsD0][Ds2PhiPi]       = 0.881908 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.00890296):0);
  Fixed_CB_fraction[DsD0][Ds2KKPi]        = 0.881908 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.00890296):0);
  Fixed_CB_fraction[DsD0][Ds2KPiPi]       = 0.878286 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.00908944):0);
  Fixed_CB_fraction[DsD0][Ds2PiPiPi]      = 0.885753 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.00702964):0);

  Fixed_CB_fraction[DsPhi][Ds2PhiPi]      = 0.807626 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0148011):0);
  Fixed_CB_fraction[DsPhi][Ds2KKPi]       = 0.807626 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0148011):0);
  Fixed_CB_fraction[DsPhi][Ds2KPiPi]      = 0.843349 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0136344):0);
  Fixed_CB_fraction[DsPhi][Ds2PiPiPi]     = 0.80677  + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0134398):0);

  Fixed_CB_fraction[DsPhiSide][Ds2PhiPi]  = 0.807626 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0148011):0);
  Fixed_CB_fraction[DsPhiSide][Ds2KKPi]   = 0.807626 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0148011):0);
  Fixed_CB_fraction[DsPhiSide][Ds2KPiPi]  = 0.843349 + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0136344):0);
  Fixed_CB_fraction[DsPhiSide][Ds2PiPiPi] = 0.80677  + (varyAllowed&&par->variation[fixedSig_Sigmafrac]?rand->Gaus(0,0.0134398):0);


  std::map<std::string,double> Fixed_Norm_Sigma_ratio;
  // Ratio of sigmas from DsD0 to DsPhi (small sigma)
  // Ratio = Sigma(DsPhi)/Sigma(DsD0)
  Fixed_Norm_Sigma_ratio[Ds2PhiPi]  = 1.27 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.02):0);
  Fixed_Norm_Sigma_ratio[Ds2KKPi]   = 1.27 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.02):0);
  Fixed_Norm_Sigma_ratio[Ds2KPiPi]  = 1.31 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.02):0);
  Fixed_Norm_Sigma_ratio[Ds2PiPiPi] = 1.26 + (varyAllowed&&par->variation[fixedSig_NormSigma]?rand->Gaus(0,0.02):0);

  // Initial Sigma values
  std::map<std::string,double> Inital_Sigma_values;
  Inital_Sigma_values[Ds2KKPi]   =  7.4;
  Inital_Sigma_values[Ds2KPiPi]  =  9.4;
  Inital_Sigma_values[Ds2PhiPi]  =  7.4;
  Inital_Sigma_values[Ds2PiPiPi] =  8.4; 


  // ======================================================== 
  // =========== Values for backgrounds      ================
  // ======================================================== 
  std::map<std::string,std::map<std::string,double>> Yield_Dsa1_input;
 
  Yield_Dsa1_input[DsPhi][Ds2PhiPi]  =  19;
  Yield_Dsa1_input[DsPhi][Ds2KKPi]   = 110; 
  Yield_Dsa1_input[DsPhi][Ds2PiPiPi] =  51; 
  Yield_Dsa1_input[DsPhi][Ds2KPiPi]  =  13; 

  std::map<std::string,std::map<std::string,double>> Yield_PR_total_input;
 
  Yield_PR_total_input[DsD0][Ds2KKPi]   = 3250; 
  Yield_PR_total_input[DsD0][Ds2KPiPi]  =  419;
  Yield_PR_total_input[DsD0][Ds2PhiPi]  = 1997; 
  Yield_PR_total_input[DsD0][Ds2PiPiPi] =  924;

  Yield_PR_total_input[DsPhi][Ds2KKPi]   =   101; 
  Yield_PR_total_input[DsPhi][Ds2KPiPi]  =    3.05;
  Yield_PR_total_input[DsPhi][Ds2PhiPi]  =    45; 
  Yield_PR_total_input[DsPhi][Ds2PiPiPi] =    18; 

  //Fraction of Ds*D0 initial value
  std::map<std::string,double> Fraction_DsstD0_initial;
  Fraction_DsstD0_initial[Ds2PhiPi]  = 0.37;//0.60;
  Fraction_DsstD0_initial[Ds2KKPi]   = 0.37; //0.68;
  Fraction_DsstD0_initial[Ds2PiPiPi] = 0.49;//0.75;
  Fraction_DsstD0_initial[Ds2KPiPi]  = 0.37;//0.56;
  RooRealVar* Fraction_DsstD0        = new RooRealVar("Fraction_DsstD0",         "" , 0.49,  0 ,1); 
  RooRealVar* Fraction_DsstD0_Helbin1= new RooRealVar("Fraction_DsstD0_Helbin1", "" , 0.286, 0 ,1); 
  RooRealVar* Fraction_DsstD0_Helbin2= new RooRealVar("Fraction_DsstD0_Helbin2", "" , 0.343, 0 ,1); 

  if(allConst_exYield) Fraction_DsstD0_Helbin1->setConstant();
  if(allConst_exYield) Fraction_DsstD0_Helbin2->setConstant();

  RooRealVar* DsPhi_BG_fraction      = new RooRealVar("Dsa1_to_DsstPhi_fraction", "" , 0.115 , 0 , 1);
  if(allConst_exYield) DsPhi_BG_fraction->setConstant();
  if(varyAllowed&&par->variation[fixedBG_noDsstPhi]) {
    DsPhi_BG_fraction->setVal(0.0);
    DsPhi_BG_fraction->setConstant();
  }
  if(varyAllowed&&par->variation[fixedBG_noDsa1]){
    DsPhi_BG_fraction->setVal(1.0);   
    DsPhi_BG_fraction->setConstant();
  }   

  
  std::map<std::string,double> Yield_DD_input;
 
  Yield_DD_input[Ds2KKPi]   = 326; 
  Yield_DD_input[Ds2KPiPi]  =  16;
  Yield_DD_input[Ds2PhiPi]  = 133; 
  Yield_DD_input[Ds2PiPiPi] =  48;
  
  std::map<std::string,double> Yield_DsstDst0_input;
 
  Yield_DsstDst0_input[Ds2KKPi]   = 3440; 
  Yield_DsstDst0_input[Ds2KPiPi]  =  427;
  Yield_DsstDst0_input[Ds2PhiPi]  = 2087; 
  Yield_DsstDst0_input[Ds2PiPiPi] =  959;


  RooRealVar* DsstDst0_mean      = new RooRealVar("DsstDst0_mean"    ,"", 5000, 4600,  5100 );
  RooRealVar* DsstDst0_sigma1    = new RooRealVar("DsstDst0_sigma1"  ,"", 10.0, 2.0,  100.0 ); 
  RooRealVar* DsstDst0_sigma2    = new RooRealVar("DsstDst0_sigma2"  ,"", 10.0, 2.0,  100.0);
  RooRealVar* DsstDst0_alpha     = new RooRealVar("DsstDst0_alpha"   ,"", -1.0, -30.0 , 0.0);
  RooRealVar* DsstDst0_n         = new RooRealVar("DsstDst0_n"       ,"", 1.0);
  RooRealVar* DsstDst0_fraction  = new RooRealVar("DsstDst0_fraction","", 0.5 , 0.0 ,1.0 ); 

  // ======================================================== 
  // =========== Fractions in different cats ================
  // ======================================================== 

  //Signal fractions 
  // Fraction of signal decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* Signal_DsPhi_Fraction   = new RooRealVar("Signal_DsPhi_Fraction",   "", 0.959 + (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.0039):0));
  // Fraction of signal decays in Helbin1 / (Helbin1+ Helbin2)
  RooRealVar* Signal_Helbin1_Fraction = new RooRealVar("Signal_Helbin1_Fraction", "", 0.929 + (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.0040):0));

  //Dsa1 fractions 
  // Fraction of Dsa1 decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* Dsa1_DsPhi_Fraction     = new RooRealVar("Dsa1_DsPhi_Fraction",     "", 0.625 + (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.007):0));
  // Fraction of Dsa1 decays in Helbin1 / (Helbin1+ Helbin2)
  //RooRealVar* Dsa1_Helbin1_Fraction   = new RooRealVar("Dsa1_Helbin1_Fraction",   "", 0.641091721 + (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.1):0));

  //Ds*Phi fractions 
  // Fraction of Ds*Phi decays in DsPhi / (DsPhiSide + DsPhi)
  RooRealVar* DsstPhi_DsPhi_Fraction = new RooRealVar("DsstPhi_DsPhi_Fraction", "", 0.959 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.0039):0));

  // if inner phi bin is <10 MeV
  if(true){
    Signal_DsPhi_Fraction->setVal( 0.882 + (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.0039):0));
    Dsa1_DsPhi_Fraction->setVal(   0.311 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.0070):0));
    DsstPhi_DsPhi_Fraction->setVal(0.882 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.0039):0));

  }

  RooRealVar* DsstKK_DsPhi_Fraction     = new RooRealVar("DsstKK_DsPhi_Fraction",     "", 0.317 , 0.0 ,1.0 );
  RooRealVar* DsstKK_Helbin1_Fraction   = new RooRealVar("DsstKK_Helbin1_Fraction",   "", 0.5   , 0.0 ,1.0 );
  
  //RooRealVar* DsD_DsPhi_Fraction     = new RooRealVar("DsD_DsPhi_Fraction",     "", 0.199  );
  //RooRealVar* DsD_Helbin1_Fraction   = new RooRealVar("DsD_Helbin1_Fraction",   "", 0.541 );
 
  RooRealVar* DsD_DsPhi_H1_Fraction     = new RooRealVar("DsD_DsPhi_H1_Fraction" ,     "" , 0.06);
  RooRealVar* DsD_DsPhi_H2_Fraction     = new RooRealVar("DsD_DsPhi_H2_Fraction" ,     "" , 0.09);
  RooRealVar* DsD_DsPhiSide_H1_Fraction = new RooRealVar("DsD_DsPhiSide_H1_Fraction" , "" , 0.49);
  RooRealVar* DsD_DsPhiSide_H2_Fraction = new RooRealVar("DsD_DsPhiSide_H2_Fraction" , "" , 0.33);
  
  RooRealVar *Ratio_DsD_to_Dsa1         = new RooRealVar("Ratio_DsD_to_Dsa1",          "",  0.1 ,0.0 ,10.0); 

  // Fractions for DD' backgrounds in extended range

  RooRealVar* DD_DsPhi_H1_Fraction     = new RooRealVar("DD_DsPhi_H1_Fraction" ,     "" , 0.70);
  RooRealVar* DD_DsPhi_H2_Fraction     = new RooRealVar("DD_DsPhi_H2_Fraction" ,     "" , 0.08);
  RooRealVar* DD_DsPhiSide_H1_Fraction = new RooRealVar("DD_DsPhiSide_H1_Fraction" , "" , 0.13);
  RooRealVar* DD_DsPhiSide_H2_Fraction = new RooRealVar("DD_DsPhiSide_H2_Fraction" , "" , 0.07);

  RooRealVar* DD_DsDs_Fraction         = new RooRealVar("DD_DsDs_Fraction",          "" , 0.45);
  RooRealVar* DD_DsstDs_Fraction       = new RooRealVar("DD_DsstDs_Fraction",        "" , 0.33);
  RooRealVar* DD_DsD_Fraction          = new RooRealVar("DD_DsD_Fraction",           "" , 0.21);

  RooRealVar *Ratio_DD_to_Dsa1         = new RooRealVar("Ratio_DD_to_Dsa1",          "",  1.006 ,0.0 ,10.0);
  if(allConst_exYield) Ratio_DD_to_Dsa1->setConstant(); 
  //Range 4900 to 5900
  if(true){
    DD_DsPhi_H1_Fraction->setVal(    0.72 + (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.02):0));
    DD_DsPhi_H2_Fraction->setVal(    0.06 + (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsPhiSide_H1_Fraction->setVal(0.14 + (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.02):0));
    DD_DsPhiSide_H2_Fraction->setVal(0.06 + (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));

    DD_DsDs_Fraction->setVal(  0.666  + (varyAllowed&&par->variation[fixedBG_DD_Fractions]?rand->Gaus(0,0.01):0));
    DD_DsstDs_Fraction->setVal(0.039  + (varyAllowed&&par->variation[fixedBG_DD_Fractions]?rand->Gaus(0,0.01):0));
    DD_DsD_Fraction->setVal(   0.295  + (varyAllowed&&par->variation[fixedBG_DD_Fractions]?rand->Gaus(0,0.01):0));
  }

  // If outer phi bin in <80 MeV


  if(false){
    
    Signal_DsPhi_Fraction->setVal(0.871);
    DsstPhi_DsPhi_Fraction->setVal(0.871);

    Dsa1_DsPhi_Fraction->setVal(0.195);

    DD_DsPhi_H1_Fraction->setVal(0.67);
    DD_DsPhi_H2_Fraction->setVal(0.05);
    DD_DsPhiSide_H1_Fraction->setVal(0.19);
    DD_DsPhiSide_H2_Fraction->setVal(0.08);

    DD_DsDs_Fraction->setVal(0.669);
    DD_DsstDs_Fraction->setVal(0.036);
    DD_DsD_Fraction->setVal(0.295);

  }


  // ======================================================== 
  // =========== Fixing ratios between Ds modes =============
  // =========== in Hel split mode to improve   =============
  // =========== stability                      =============
  // ========================================================

  // Fraction of DsD0 peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsD0_peak_fraction    = new RooRealVar("splitHel_DsD0_peak_fraction",    "", 0.594, 0.0 ,1.0 );
  if(allConst_exYield) splitHel_DsD0_peak_fraction->setConstant();

  // Fraction of DsD0 PR peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsD0_PR_peak_fraction = new RooRealVar("splitHel_DsD0_PR_peak_fraction", "", 0.595, 0.0 ,1.0 );
  if(allConst_exYield) splitHel_DsD0_PR_peak_fraction->setConstant();
  
  // Fraction of Ds*Phi peak in two helicity bins to be the same for all Ds modes
  RooRealVar* splitHel_DsstPhi_peak_fraction = new RooRealVar("splitHel_DsstPhi_peak_fraction", "", 0.920, 0.0 ,1.0 );
  splitHel_DsstPhi_peak_fraction->setConstant();
  // Fraction of Dsa1 peak in two helicity bins to be the same for all Ds modes
  //RooRealVar* splitHel_Dsa1_peak_fraction   = new RooRealVar("splitHel_Dsa1_peak_fraction", "", 0.51, 0.0 ,1.0 );
  //RooRealVar* splitHel_Dsa1_peak_fraction    = new RooRealVar("splitHel_Dsa1_peak_fraction",    "", 0.641091721 + (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.1):0));
  RooRealVar* splitHel_Dsa1_peak_fraction    = new RooRealVar("splitHel_Dsa1_peak_fraction",    "", 0.662 + (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.1):0));




  // ======================================================== 
  // =========== Global variables      ======================
  // ========================================================
  
  RooRealVar* global_csi_hill   = new RooRealVar("global_csi_hill",   "Global csi HILL",    1.0             );
  RooRealVar* global_csi        = new RooRealVar("global_csi",        "Global csi",         0.85,   0,    2 );
  RooRealVar* global_R          = new RooRealVar("global_R",          "Global R",           1.0             );
  RooRealVar* global_f          = new RooRealVar("global_f",          "Global f",           1.0             );
  RooRealVar* global_G          = new RooRealVar("global_G",          "Global G",           1.0             );

  if(allConst_exYield) global_shift->setConstant();
  if(allConst_exYield) global_csi->setConstant();

  // Set a and b endpoints for Roo*Dini shapes

  // Ds*D0 pi0
  RooRealVar* DsstarD0_pi0_a     = new RooRealVar("DsstarD0_pi0_a",     "DsstarD0 pi0 a",     5051.351  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsstarD0_pi0_b     = new RooRealVar("DsstarD0_pi0_b",     "DsstarD0 pi0 b",     5132.858  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  // Ds*D0 gamma
  RooRealVar* DsstarD0_gamma_a   = new RooRealVar("DsstarD0_gamma_a",   "DsstarD0 gamma a",   4976.730  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsstarD0_gamma_b   = new RooRealVar("DsstarD0_gamma_b",   "DsstarD0 gamma b",   5213.055  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  // DsD*0 pi0
  RooRealVar* DsDstar0_pi0_a     = new RooRealVar("DsDstar0_pi0_a",     "DsDstar0 pi0 a",     5051.478  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsDstar0_pi0_b     = new RooRealVar("DsDstar0_pi0_b",     "DsDstar0 pi0 b",     5128.577  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );   
  // DsD*0 gamma
  RooRealVar* DsDstar0_gamma_a   = new RooRealVar("DsDstar0_gamma_a",   "DsDstar0 gamma a",   4970.140  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsDstar0_gamma_b   = new RooRealVar("DsDstar0_gamma_b",   "DsDstar0 gamma b",   5216.106  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );          

  // Ds*Phi pi0
  RooRealVar* DsstarPhi_pi0_a     = new RooRealVar("DsstarPhi_pi0_a",     "DsstarPhi pi0 a",     5026.751    + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsstarPhi_pi0_b     = new RooRealVar("DsstarPhi_pi0_b",     "DsstarPhi pi0 b",     5124.790    + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  //Ds*Phi gamma
  RooRealVar* DsstarPhi_gamma_a   = new RooRealVar("DsstarPhi_gamma_a",   "DsstarPhi gamma a",   4936.387    + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsstarPhi_gamma_b   = new RooRealVar("DsstarPhi_gamma_b",   "DsstarPhi gamma b",   5220.646    + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) ); 

  

  // Ds*D*0 gamma gamma
  RooRealVar* DsstarDstar0_gamma_gamma_a   = new RooRealVar("DsstarDstar0_gamma_gamma_a",   "DsstarDstar0 gamma a",   4667  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );
  RooRealVar* DsstarDstar0_gamma_gamma_b   = new RooRealVar("DsstarDstar0_gamma_gamma_b",   "DsstarDstar0 gamma b",   5150  + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,1):0) );          

  RooRealVar* DsstarDstar0_gamma_pi0_a     = new RooRealVar("DsstarDstar0_gamma_pi0_a",     "DsstarDstar0 pi0 a",     4748  + (varyAllowed&&par->variation[fixedBG_DsstDst0_endpoints]?rand->Gaus(0,10):0) );
  RooRealVar* DsstarDstar0_gamma_pi0_b     = new RooRealVar("DsstarDstar0_gamma_pi0_b",     "DsstarDstar0 pi0 b",     5062  + (varyAllowed&&par->variation[fixedBG_DsstDst0_endpoints]?rand->Gaus(0,10):0) );          

  RooRealVar *Ratio_DsstDst0_to_Lowmass = new RooRealVar("Ratio_DsstDst0_to_Lowmass","Ratio_DsstDst0_to_Lowmass",0.639 ,0.0 ,10.0);
  if(allConst_exYield) Ratio_DsstDst0_to_Lowmass->setConstant();

  // Fraction of D* -> D pi0 vs. D* -> D pi0 gamma adjusted for fraction in mass window
  RooRealVar* fraction_Dst0_D0pi0_DsD0_adj      = new RooRealVar("fraction_Dst0_D0pi0_DsD0",      "",  0.66 * (0.838/0.713) + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.03):0));
  RooRealVar* fraction_Dsst_Dspi0_DsD0_adj      = new RooRealVar("fraction_Dsst_Dspi0_DsD0_adj",  "",  0.06 * (0.841/0.732) + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.007):0));

  RooRealVar* fraction_Dsst_Dspi0_DsPhi_010_adj = new RooRealVar("fraction_Dsst_Dspi0_DsPhi_adj", "",  0.06 * (0.531/0.610) + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.007):0));
  RooRealVar* fraction_Dsst_Dspi0_DsPhi_101_adj = new RooRealVar("fraction_Dsst_Dspi0_DsPhi_adj", "",  0.06 * (0.779/0.544) + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.007):0));
  
  // When using extended range, whole shape in mass range
  if(true){
    fraction_Dst0_D0pi0_DsD0_adj->setVal(     0.66 + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.03):0));
    fraction_Dsst_Dspi0_DsD0_adj->setVal(     0.06 + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.007):0));
    fraction_Dsst_Dspi0_DsPhi_010_adj->setVal(0.06 + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.007):0));
    fraction_Dsst_Dspi0_DsPhi_101_adj->setVal(0.06 + (varyAllowed&&par->variation[fixedBG_DsD0]?rand->Gaus(0,0.007):0));
  }

  // ======================================================== 
  // ===========    Ds*Phi shape parameters    ==============
  // ========================================================
  RooRealVar* DsstPhi_long_Helbin1_frac = new RooRealVar("DsstPhi_long_Helbin1_frac", "",  0.936 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.1):0));
  RooRealVar* DsstPhi_tran_Helbin1_frac = new RooRealVar("DsstPhi_tran_Helbin1_frac", "",  0.432 + (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.1):0));

  RooRealVar* fraction_hel              = new RooRealVar("fraction_hel",   "Polariztion Fraction",  0.5 + (varyAllowed&&par->variation[fixedBG_hel]?rand->Gaus(0,0.2):0) , 0,1 );
  if(fraction_hel->getVal() > 1 ) fraction_hel->setVal(1.0);
  if(fraction_hel->getVal() < 0 ) fraction_hel->setVal(0.0);

  if(allConst_exYield) fraction_hel->setConstant();
  fraction_hel->setConstant();

  RooFormulaVar* DsstPhi_Helbin1_frac  = new RooFormulaVar("DsstPhi_Helbin1_frac", "",  "@0*@1 + (1-@0)*@2",         RooArgList(*fraction_hel,*DsstPhi_long_Helbin1_frac,*DsstPhi_tran_Helbin1_frac));
  

  // ==================================================================
  // =========== Modifiy BG ratios depending on the range =============
  // ==================================================================

  if( mB->getMin() == 4900) {
    Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
    Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
    DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

    Dsa1_DsPhi_Fraction->setVal(        0.323858+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
    splitHel_Dsa1_peak_fraction->setVal(0.660636+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123017):0));

    DD_DsPhi_H1_Fraction->setVal(       0.718579+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235058):0));
    DD_DsPhi_H2_Fraction->setVal(       0.0819672+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143387):0));
    DD_DsPhiSide_H1_Fraction->setVal(   0.128415+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174873):0));
    DD_DsPhiSide_H2_Fraction->setVal(   0.0710383+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0134278):0));

    DD_DsDs_Fraction->setVal(   0.666475+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsstDs_Fraction->setVal( 0.0402676+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsDs_Fraction->setVal(   0.293258+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4895) {
    Signal_Helbin1_Fraction->setVal(    0.931745+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419522):0));
    Signal_DsPhi_Fraction->setVal(      0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));
    DsstPhi_DsPhi_Fraction->setVal(     0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));

    Dsa1_DsPhi_Fraction->setVal(        0.324062+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397661):0));
    splitHel_Dsa1_peak_fraction->setVal(0.660502+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.01228):0));

    DD_DsPhi_H1_Fraction->setVal(       0.718085+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0232035):0));
    DD_DsPhi_H2_Fraction->setVal(       0.0851064+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143904):0));
    DD_DsPhiSide_H1_Fraction->setVal(   0.12766+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0172098):0));
    DD_DsPhiSide_H2_Fraction->setVal(   0.0691489+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.013084):0));

    DD_DsDs_Fraction->setVal(   0.651473+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsstDs_Fraction->setVal( 0.0536401+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsDs_Fraction->setVal(   0.294887+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4890) {
    Signal_Helbin1_Fraction->setVal(    0.931745+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419522):0));
    Signal_DsPhi_Fraction->setVal(      0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));
    DsstPhi_DsPhi_Fraction->setVal(     0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));

    Dsa1_DsPhi_Fraction->setVal(        0.323397+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397661):0));
    splitHel_Dsa1_peak_fraction->setVal(0.660841+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0122715):0));

    DD_DsPhi_H1_Fraction->setVal(       0.717617+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0229125):0));
    DD_DsPhi_H2_Fraction->setVal(       0.0854922+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0142319):0));
    DD_DsPhiSide_H1_Fraction->setVal(   0.126943+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0169446):0));
    DD_DsPhiSide_H2_Fraction->setVal(   0.0699482+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0129822):0));

    DD_DsDs_Fraction->setVal(   0.64109+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsstDs_Fraction->setVal( 0.0661085+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
    DD_DsDs_Fraction->setVal(   0.292801+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4885) {
      Signal_Helbin1_Fraction->setVal(    0.931745+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419522):0));
      Signal_DsPhi_Fraction->setVal(      0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883192+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397661):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323177+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397661):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660775+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0122471):0));

      DD_DsPhi_H1_Fraction->setVal(       0.714646+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0226929):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0883838+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0142641):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.128788+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0168326):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0681818+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0126664):0));

      DD_DsDs_Fraction->setVal(   0.641203+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0644504+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.294347+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4905) {
      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323334+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.661172+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123202):0));

      DD_DsPhi_H1_Fraction->setVal(       0.724432+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0238145):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0795455+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0144224):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.122159+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174542):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0738636+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0139406):0));

      DD_DsDs_Fraction->setVal(   0.683864+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0257769+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.290359+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));

  } else if( mB->getMin() == 4910) {
      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.322785+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.661339+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123426):0));

      DD_DsPhi_H1_Fraction->setVal(       0.726471+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0241753):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0764706+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0144123):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.120588+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0176608):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0764706+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0144123):0));

      DD_DsDs_Fraction->setVal(   0.693355+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0180381+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.288607+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4915) {
      Signal_Helbin1_Fraction->setVal(    0.931749+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419532):0));
      Signal_DsPhi_Fraction->setVal(      0.883184+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397665):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883184+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397665):0));

      Dsa1_DsPhi_Fraction->setVal(        0.322681+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397665):0));
      splitHel_Dsa1_peak_fraction->setVal(0.661765+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123614):0));

      DD_DsPhi_H1_Fraction->setVal(       0.724771+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0246987):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0764526+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0146944):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.119266+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0179228):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0795107+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0149606):0));

      DD_DsDs_Fraction->setVal(   0.702279+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.00949829+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.288222+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  } else if( mB->getMin() == 4902) {

      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323941+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660551+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123022):0));

      DD_DsPhi_H1_Fraction->setVal(       0.720994+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235732):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0801105+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0142678):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.127072+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0175049):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0718232+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0135704):0));

      DD_DsDs_Fraction->setVal(   0.677097+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0330891+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.289814+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  } else if( mB->getMin() == 4898) {

      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323614+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660567+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0122919):0));

      DD_DsPhi_H1_Fraction->setVal(       0.719346+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0234542):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0817439+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143013):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.128065+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174431):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0708447+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0133926):0));

      DD_DsDs_Fraction->setVal(   0.662672+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0399288+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.297399+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4899) {
      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323866+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660652+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0122977):0));

      DD_DsPhi_H1_Fraction->setVal(       0.718579+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235058):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0819672+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143387):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.128415+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174873):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0710383+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0134278):0));

      DD_DsDs_Fraction->setVal(   0.664997+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0401783+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.294824+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4901) {

      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323941+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660551+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123022):0));

      DD_DsPhi_H1_Fraction->setVal(       0.71978+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235396):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0824176+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0144139):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.126374+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174157):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0714286+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0134987):0));

      DD_DsDs_Fraction->setVal(   0.671998+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0408244+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.287177+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4900.5) {
      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323941+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660551+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0123022):0));

      DD_DsPhi_H1_Fraction->setVal(       0.717808+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235575):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0821918+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143762):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.128767+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0175317):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0712329+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0134632):0));

      DD_DsDs_Fraction->setVal(   0.668843+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0405214+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.290635+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  
  } else if( mB->getMin() == 4899.5) {
      Signal_Helbin1_Fraction->setVal(    0.931739+ (varyAllowed&&par->variation[fixedSig_BinRatios]?rand->Gaus(0,0.00419523):0));
      Signal_DsPhi_Fraction->setVal(      0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));
      DsstPhi_DsPhi_Fraction->setVal(     0.883201+ (varyAllowed&&par->variation[fixedBG_DsPhi]?rand->Gaus(0,0.00397669):0));

      Dsa1_DsPhi_Fraction->setVal(        0.323866+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.00397669):0));
      splitHel_Dsa1_peak_fraction->setVal(0.660652+ (varyAllowed&&par->variation[fixedBG_Dsa1]?rand->Gaus(0,0.0122977):0));

      DD_DsPhi_H1_Fraction->setVal(       0.718579+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0235058):0));
      DD_DsPhi_H2_Fraction->setVal(       0.0819672+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0143387):0));
      DD_DsPhiSide_H1_Fraction->setVal(   0.128415+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0174873):0));
      DD_DsPhiSide_H2_Fraction->setVal(   0.0710383+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.0134278):0));

      DD_DsDs_Fraction->setVal(   0.666475+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsstDs_Fraction->setVal( 0.0402676+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
      DD_DsDs_Fraction->setVal(   0.293258+ (varyAllowed&&par->variation[fixedBG_DD_Ratios]?rand->Gaus(0,0.01):0));
  }








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

 
  // ======================================================== 
  // =========== Create Branching fraction ==================
  // ======================================================== 
  

  //============= Efficiency ratios weighted for years 
  //== eff taken from ANA note
  //== weighting taken from number of DsD0 events in each year
  //=========================================================
  //== Ratio is eff(DsD0)/eff(DsPhi)

  // ratio * weight 
  //====================|==== 2011 ===|==== 2012 ====|==== 2015 ====|=== 2016 ===|
  eff_ratio[Ds2PhiPi]  =  0.770*0.13  +  0.795*0.32  +  0.710*0.07  +  0.739*0.48;
  eff_ratio[Ds2KKPi]   =  0.770*0.12  +  0.795*0.32  +  0.710*0.09  +  0.739*0.48;
  eff_ratio[Ds2KPiPi]  =  1.139*0.15  +  1.165*0.26  +  1.101*0.11  +  1.141*0.47;
  eff_ratio[Ds2PiPiPi] =  0.916*0.12  +  0.943*0.32  +  0.858*0.11  +  0.892*0.44;
  
  //eff_ratio[Ds2PhiPi]  =  0.773*0.13  +  0.799*0.32  +  0.666*0.07  +  0.719*0.48;
  //eff_ratio[Ds2KKPi]   =  0.773*0.12  +  0.799*0.32  +  0.666*0.09  +  0.719*0.48;
  //eff_ratio[Ds2KPiPi]  =  1.143*0.15  +  1.171*0.26  +  1.094*0.11  +  1.110*0.47;
  //eff_ratio[Ds2PiPiPi] =  0.916*0.12  +  0.945*0.32  +  0.833*0.11  +  0.871*0.44;

  // ratio * weight 
  //========================|==== 2011 ===|==== 2012 ====|==== 2015 ====|=== 2016 ===|
  eff_ratio_err[Ds2PhiPi]  =  0.011*0.20  +  0.012*0.51  +  0.015*0.13  +  0.016*0.15;
  eff_ratio_err[Ds2KKPi]   =  0.011*0.20  +  0.012*0.49  +  0.015*0.15  +  0.016*0.16;
  eff_ratio_err[Ds2KPiPi]  =  0.022*0.18  +  0.023*0.48  +  0.029*0.16  +  0.030*0.16;
  eff_ratio_err[Ds2PiPiPi] =  0.014*0.16  +  0.015*0.41  +  0.017*0.23  +  0.017*0.21;
 


  eff_ratio_rrv[Ds2PhiPi]  = new RooRealVar("eff_ratio_Ds2PhiPi", "", eff_ratio[Ds2PhiPi] );
  eff_ratio_rrv[Ds2KKPi]   = new RooRealVar("eff_ratio_Ds2KKPi",  "", eff_ratio[Ds2KKPi]  );
  eff_ratio_rrv[Ds2KPiPi]  = new RooRealVar("eff_ratio_Ds2KPiPi", "", eff_ratio[Ds2KPiPi] );
  eff_ratio_rrv[Ds2PiPiPi] = new RooRealVar("eff_ratio_Ds2PiPiPi","", eff_ratio[Ds2PiPiPi]);

  //eff_ratio_rrv[Ds2PhiPi]->setError(eff_ratio_err[Ds2PhiPi]);
  //eff_ratio_rrv[Ds2KKPi]->setError(eff_ratio_err[Ds2KKPi]);
  //eff_ratio_rrv[Ds2KPiPi]->setError(eff_ratio_err[Ds2KPiPi]);
  //eff_ratio_rrv[Ds2PiPiPi]->setError(eff_ratio_err[Ds2PiPiPi]);


 
  bool fitSingleBR = par->fitBr;
  bool fitFourBr = par->fitFourBr; 
  
  // Single Branching fraction
  //Branching_fraction_all = new RooRealVar("Branching_fraction",       "",  eff_ratio_rrv[Ds2KKPi]->getVal()*(Yield_CB_input[DsPhi][Ds2KKPi]/Yield_CB_input[DsD0][Ds2KKPi]) ,0.0,1.0 );
  
  // This factor is used to get the branching fraction in units of 10^{-7}

  ///////-----------------------------------------------------------------------
  ////   Correction factor = 10^7 * Br(B->DsD0) * ( Br(D0->KK)/Br(phi -> KK) )
  ////                     = 10^7 * 9.0e-3      * ( 4.01e-3   / 0.489        ) 
  ////                     = 738.0368098 
  ///////-----------------------------------------------------------------------
  RooRealVar* Correction_factor = new RooRealVar("Correction_factor","",738.0368098);
  
  Branching_fraction_all = new RooRealVar("Branching_fraction",       "Branching Fraction (#times10^{-7})",  10.0, -100.0, 10000.0 ); //*0.01355
  //Branching_fraction_all->setVal(3.0);
  
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
  Branching_fraction[Ds2PhiPi]   = new RooRealVar("Branching_fraction_Ds2PhiPi",  "",  18.7 , -100.0, 10000.0 );
  Branching_fraction[Ds2KKPi]    = new RooRealVar("Branching_fraction_Ds2KKPi",   "",  18.7 , -100.0, 10000.0 );
  Branching_fraction[Ds2KPiPi]   = new RooRealVar("Branching_fraction_Ds2KPiPi",  "",  18.7 , -100.0, 10000.0 );
  Branching_fraction[Ds2PiPiPi]  = new RooRealVar("Branching_fraction_Ds2PiPiPi", "",  18.7 , -100.0, 10000.0 );
 
 
 
  // ======================================================== 
  // =========== Create Asymmetry ===========================
  // ======================================================== 

  Asymmetry_all         = new RooRealVar("Asymmetry_all",       "Asymmetry_all",  0.0, -1.0, 1.0 );

  Asymmetry[Ds2PhiPi]   = new RooRealVar("Asymmetry_Ds2PhiPi",  "Asymmetry_Ds2PhiPi",  0.0 ,-1.0, 1.0 );
  Asymmetry[Ds2KKPi]    = new RooRealVar("Asymmetry_Ds2KKPi",   "Asymmetry_Ds2KKPi",   0.0 ,-1.0, 1.0 );
  Asymmetry[Ds2KPiPi]   = new RooRealVar("Asymmetry_Ds2KPiPi",  "Asymmetry_Ds2KPiPi",  0.0 ,-1.0, 1.0 );
  Asymmetry[Ds2PiPiPi]  = new RooRealVar("Asymmetry_Ds2PiPiPi", "Asymmetry_Ds2PiPiPi", 0.0 ,-1.0, 1.0 );
 
  if(needsBlinding){
    Asymmetry_all_unblinded         = new RooUnblindUniform("Asymmetry_UnBlinded",           "Asymmetry Unblind ",            "UnblindAsymmetryString",         0.2, *Asymmetry_all  );
    Asymmetry_unblinded[Ds2PhiPi]   = new RooUnblindUniform("Asymmetry_UnBlinded_Ds2PhiPi",  "Asymmetry UnBlinded Ds2PhiPi",  "UnblindAsymmetryStringDs2PhiPi", 0.2, *Asymmetry[Ds2PhiPi]);
    Asymmetry_unblinded[Ds2KKPi]    = new RooUnblindUniform("Asymmetry_UnBlinded_Ds2KKPi",   "Asymmetry UnBlinded Ds2KKPi",   "UnblindAsymmetryStringDs2KKPi",  0.2, *Asymmetry[Ds2KKPi]);
    Asymmetry_unblinded[Ds2KPiPi]   = new RooUnblindUniform("Asymmetry_UnBlinded_Ds2KPiPi",  "Asymmetry UnBlinded Ds2KPiPi",  "UnblindAsymmetryStringDs2KPiPi", 0.2, *Asymmetry[Ds2KPiPi]);
    Asymmetry_unblinded[Ds2PiPiPi]  = new RooUnblindUniform("Asymmetry_UnBlinded_Ds2PiPiPi", "Asymmetry UnBlinded Ds2PiPiPi", "UnblindAsymmetryStringDs2PiPiPi",0.2, *Asymmetry[Ds2PiPiPi]);           
  } else {
    Asymmetry_all_unblinded         = Asymmetry_all;      
    Asymmetry_unblinded[Ds2PhiPi]   = Asymmetry[Ds2PhiPi];
    Asymmetry_unblinded[Ds2KKPi]    = Asymmetry[Ds2KKPi];  
    Asymmetry_unblinded[Ds2KPiPi]   = Asymmetry[Ds2KPiPi];   
    Asymmetry_unblinded[Ds2PiPiPi]  = Asymmetry[Ds2PiPiPi];  
  }

  RooRealVar* Aprod = new RooRealVar("Aprod","Aprod", -0.005);

  std::map<std::string,RooRealVar*> Adet;

  Adet[Ds2KKPi]    = new RooRealVar("Adet_Ds2KKPi",  "Adet_Ds2KKPi",   0.0008);          // Api
  Adet[Ds2PhiPi]   = new RooRealVar("Adet_Ds2PhiPi", "Adet_Ds2PhiPi",  0.0008);          // Api
  Adet[Ds2PiPiPi]  = new RooRealVar("Adet_Ds2PiPiPi","Adet_Ds2PiPiPi", 0.0008);          // Api
  Adet[Ds2KPiPi]   = new RooRealVar("Adet_Ds2KPiPi", "Adet_Ds2KPiPi", -0.0106 + 0.0008); // AKpi + Api

  // --------- Yields for DsD0, in an array (as line 447 of Model.C) ---------------
  if(par->debug) std::cout<<"Running: DsPhiModel::DefineModel() --> Making RooRealVars"<<std::endl;
 
  // ------------------------------------------------------
  // --------- Setup yields for signals       -------------
  // ------------------------------------------------------
 std::string cat_name;


  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){  
        for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            




            //for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){

            // ------------------------------------------------------
            // --------- Setup yields for signals       -------------
            // ------------------------------------------------------

            // Fix yields for not helicity split: both
            for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("yield_peak_%s_%s_%s",      DsD0.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsD0].c_str(),      mod[(*m)].c_str()), Yield_CB_input[DsD0][*m],    0, 100000);
              if(allConst_exYield) ((RooRealVar*)yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
            }



            cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),both.c_str(),(*a).c_str());
            
            if(par->sumOverCharges){
              if(fitFourBr){
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a]   = new RooFormulaVar(Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(@0*@1)/(@2*@3)" , RooArgList(*Branching_fraction[*m],*yield_peak[*t][both][*ds][*ph][DsD0][*m][both][*a], *eff_ratio_rrv[*m],*Correction_factor));
              }else if(fitSingleBR){
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a]   = new RooFormulaVar(Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"(@0*@1)/(@2*@3)" , RooArgList(*Branching_fraction_all,*yield_peak[*t][both][*ds][*ph][DsD0][*m][both][*a], *eff_ratio_rrv[*m],*Correction_factor)); 
              } else {
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a]   = new RooRealVar(   Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -20, 100000); 
                if(NoYield[*m]){
                  ((RooRealVar*)Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a])->setVal(0.0);
                  ((RooRealVar*)Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a])->setConstant();
                }
              }

            } else {
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),both.c_str(),(*a).c_str());
              Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a]   = new RooRealVar(   Form("yield_peak_total_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield total %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), 2*Yield_CB_input[DsPhi][*m],  -20, 100000); 
                
              if(fitFourBr){
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][plus][*a]   = new RooFormulaVar(Form("yield_peak_plus_%s_%s_%s", DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield plus  %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "@0*(1-@1)/2-@2-@3" , RooArgList(*Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a],*Asymmetry_unblinded[*m],*Aprod,*Adet[*m])); 
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][minus][*a]  = new RooFormulaVar(Form("yield_peak_minus_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield minus %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "@0*(1+@1)/2-@2-@3" , RooArgList(*Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a],*Asymmetry_unblinded[*m],*Aprod,*Adet[*m])); 
              
              } else {
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][plus][*a]   = new RooFormulaVar(Form("yield_peak_plus_%s_%s_%s", DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield plus  %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "@0*(1-@1)/2-@2-@3" , RooArgList(*Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a],*Asymmetry_all_unblinded,*Aprod,*Adet[*m])); 
                Signal_total[*t][both][*ds][*ph][DsPhi][*m][minus][*a]  = new RooFormulaVar(Form("yield_peak_minus_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield minus %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), "@0*(1+@1)/2-@2-@3" , RooArgList(*Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a],*Asymmetry_all_unblinded,*Aprod,*Adet[*m])); 
                //Signal_total[*t][both][*ds][*ph][DsPhi][*m][plus][*a]   = new RooRealVar(   Form("yield_peak_plus_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield total %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -20, 100000); 
                //Signal_total[*t][both][*ds][*ph][DsPhi][*m][minus][*a]  = new RooRealVar(   Form("yield_peak_minus_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield total %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()), Yield_CB_input[DsPhi][*m],  -20, 100000); 
              

              }
            }

            for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());

              yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhi.c_str(),    (*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),"@0*@1" ,     RooArgList(*Signal_DsPhi_Fraction, *Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_peak_%s_%s_%s",DsPhiSide.c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield %s %s",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),"(1-@0)*@1" , RooArgList(*Signal_DsPhi_Fraction, *Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              // Blinding
              if(needsBlinding){
                B_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a];
                B_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooUnblindUniform(Form("B_nsig_DsPhi_%s_%s",    (*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "nsigBuDsPhiblindePhi",       Blind_number[DsPhi][*m][both],     *yield_peak[*t][both][*ds][*ph][DsPhi][*m][*c][*a] );
                B_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooUnblindUniform(Form("B_nsig_DsPhiSide_%s_%s",(*m).c_str(),cat_name.c_str()),  Form("Unblind Yield %s %s",Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "nsigBuDsPhiblindePhiSide",   Blind_number[DsPhiSide][*m][both], *yield_peak[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] );
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
            //}

              // Setup Branching fractions
              //Branching_fraction[*m]    = new RooFormulaVar(Form("Branching_fraction_%s",(*m).c_str()), "", "(@0/@1)", RooArgList(*Signal_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a],*yield_peak[*t][both][*ds][*ph][DsD0][*m][*c][*a]) ); //,*eff_ratio_rrv[*m]


              // ------------------------------------------------------
              // --------- Setup yields for backgrounds   -------------
              // ------------------------------------------------------

              
              // Both
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());              
               
              Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  Yield_PR_total_input[DsD0][*m],  0,  1000000);
              Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooRealVar(   Form("low_mass_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  Yield_PR_total_input[DsPhi][*m], 0,  1000000);
              if(allConst_exYield) ((RooRealVar*)Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
              if(allConst_exYield) ((RooRealVar*)Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a])->setConstant();


              PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "@0",             RooArgList(*Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a] ));//;
              PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",          RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));//;
              PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "((1-@0)/@0)*@1", RooArgList(*DsstPhi_DsPhi_Fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              //yield_DsstarDstar0[*t][both][*ds][*ph][DsD0][*m][*c][*a]  = new RooRealVar(   Form("yield_DsstarDstar0_%s_%s_%s",DsD0.c_str(),  (*m).c_str(), cat_name.c_str()), Form("Ds*D*0 yield %s %s",Bmod[DsD0].c_str(),    mod[(*m)].c_str()),  Yield_DsstDst0_input[*m],  0,  1000000);
              yield_DsstarDstar0[*t][both][*ds][*ph][DsD0][*m][*c][*a]  = new RooFormulaVar(Form("yield_DsstarDstar0_%s_%s_%s",DsD0.c_str(),  (*m).c_str(), cat_name.c_str()), Form("Ds*D*0 yield %s %s",Bmod[DsD0].c_str(),    mod[(*m)].c_str()),  "@0*@1",          RooArgList(*Ratio_DsstDst0_to_Lowmass,*Low_Mass_total[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              //if(allConst_exYield) ((RooRealVar*)yield_DsstarDstar0[*t][both][*ds][*ph][DsD0][*m][*c][*a])->setConstant();

              yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1" ,    RooArgList(*DsPhi_BG_fraction,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "((1-@0)/@0)*@1",RooArgList(*Dsa1_DsPhi_Fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a])); 
              
              yield_dsstKK[*t][both][*ds][*ph][DsPhi][*m][*c][*a]       = new RooRealVar(   Form("yield_dsstKK_%s_%s_%s",  DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  150, 0,  1000000);
              yield_dsstKK[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]   = new RooFormulaVar(Form("yield_dsstKK_%s_%s_%s",  DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DsstKK_DsPhi_Fraction,*yield_dsstKK[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              
              //yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]          = new RooRealVar(   Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  150, 0,  1000000);
              //yield_DsD[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "((1-@0)/@0)*@1",   RooArgList(*DsD_DsPhi_Fraction,*yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]          = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1",               RooArgList(*Ratio_DsD_to_Dsa1, *Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_DsD[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*((@1+@2)/(@3+@4))",RooArgList(*yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a],*DsD_DsPhiSide_H1_Fraction,*DsD_DsPhiSide_H2_Fraction,*DsD_DsPhi_H1_Fraction,*DsD_DsPhi_H2_Fraction));
              

              //yield_DD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]           = new RooRealVar(   Form("yield_DD_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  Yield_DD_input[*m] , 0,  1000000);
              yield_DD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]           = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1",               RooArgList(*Ratio_DD_to_Dsa1, *Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_DD[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]       = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*((@1+@2)/(@3+@4))",RooArgList(*yield_DD[*t][both][*ds][*ph][DsPhi][*m][*c][*a],*DD_DsPhiSide_H1_Fraction,*DD_DsPhiSide_H2_Fraction,*DD_DsPhi_H1_Fraction,*DD_DsPhi_H2_Fraction));
              //if(allConst_exYield) ((RooRealVar*)yield_DD[*t][both][*ds][*ph][DsPhi][*m][*c][*a])->setConstant();

              // Helbin 1
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              PR_total_yield[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "@0*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              //PR_total_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              //PR_total_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1", RooArgList(*DsstPhi_Helbin1_frac,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1", RooArgList(*DsstPhi_Helbin1_frac,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
               
              yield_DsstarDstar0[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]  = new RooFormulaVar(Form("yield_DsstarDstar0_%s_%s_%s",DsD0.c_str(),  (*m).c_str(), cat_name.c_str()), Form("Ds*D*0 yield %s %s",Bmod[DsD0].c_str(),    mod[(*m)].c_str()),  "@0*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*yield_DsstarDstar0[*t][both][*ds][*ph][DsD0][*m][*c][*a]));


              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));
              yield_Dsa1[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
              
              yield_dsstKK[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("yield_dsstKK_%s_%s_%s",  DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s",Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",   RooArgList(*DsstKK_Helbin1_Fraction,    *yield_dsstKK[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_dsstKK[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]   = new RooFormulaVar(Form("yield_dsstKK_%s_%s_%s",  DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s",Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1",   RooArgList(*DsstKK_Helbin1_Fraction,    *yield_dsstKK[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              //yield_DsD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]          = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",   RooArgList(*DsD_Helbin1_Fraction,    *yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              //yield_DsD[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*@1",   RooArgList(*DsD_Helbin1_Fraction,    *yield_DsD[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              yield_DsD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]          = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*@1",     RooArgList(*Ratio_DsD_to_Dsa1,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_DsD[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s",  Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DsD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DsD_DsPhiSide_H1_Fraction,*DsD_DsPhi_H1_Fraction));
              
              //yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]           = new RooRealVar(   Form("yield_DD_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  Yield_DD_input[*m] , 0,  1000000);
              yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]           = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*@1",     RooArgList(*Ratio_DD_to_Dsa1,*Low_Mass_total[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_DD[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a]       = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DD_DsPhiSide_H1_Fraction,*DD_DsPhi_H1_Fraction));
              //if(allConst_exYield) ((RooRealVar*)yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a])->setConstant();

              // Helbin 2
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin2.c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());
              PR_total_yield[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsD0.c_str(),      (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsD0].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsD0][*m][*c][*a]));
              //PR_total_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              //PR_total_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*splitHel_DsstPhi_peak_fraction,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DsstPhi_Helbin1_frac,*PR_total_yield[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              PR_total_yield[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a] = new RooFormulaVar(Form("yield_PR_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("PR yield %s %s",   Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1", RooArgList(*DsstPhi_Helbin1_frac,*PR_total_yield[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));
              
              yield_DsstarDstar0[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]  = new RooFormulaVar(Form("yield_DsstarDstar0_%s_%s_%s",DsD0.c_str(),  (*m).c_str(), cat_name.c_str()), Form("Ds*D*0 yield %s %s",Bmod[DsD0].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1" ,RooArgList(*splitHel_DsD0_PR_peak_fraction,*yield_DsstarDstar0[*t][both][*ds][*ph][DsD0][*m][*c][*a]));

              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]         = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "(1-@0)*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhi][*m][*c][*a] ));
              yield_Dsa1[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]     = new RooFormulaVar(Form("yield_Dsa1_%s_%s_%s",    DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("Dsa1 yield %s %s", Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "(1-@0)*@1",RooArgList(*splitHel_Dsa1_peak_fraction, *yield_Dsa1[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a])); 
         
              yield_dsstKK[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("yield_dsstKK_%s_%s_%s",  DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s",Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1",   RooArgList(*DsstKK_Helbin1_Fraction,    *yield_dsstKK[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              yield_dsstKK[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]   = new RooFormulaVar(Form("yield_dsstKK_%s_%s_%s",  DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsstKK yield %s %s",Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1",   RooArgList(*DsstKK_Helbin1_Fraction,    *yield_dsstKK[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));     
              
              //yield_DsD[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]          = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s", Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "(1-@0)*@1",   RooArgList(*DsD_Helbin1_Fraction,    *yield_DsD[*t][both][*ds][*ph][DsPhi][*m][*c][*a]));
              //yield_DsD[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s", Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "(1-@0)*@1",   RooArgList(*DsD_Helbin1_Fraction,    *yield_DsD[*t][both][*ds][*ph][DsPhiSide][*m][*c][*a]));     
              yield_DsD[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]          = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s", Bmod[DsPhi].c_str(),    mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DsD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DsD_DsPhi_H2_Fraction,    *DsD_DsPhi_H1_Fraction));
              yield_DsD[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsD_%s_%s_%s",     DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DsD yield %s %s", Bmod[DsPhiSide].c_str(),mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DsD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DsD_DsPhiSide_H2_Fraction,*DsD_DsPhi_H1_Fraction));
              
              yield_DD[*t][Helbin2][*ds][*ph][DsPhi][*m][*c][*a]           = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhi.c_str(),     (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhi].c_str(),     mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DD_DsPhi_H2_Fraction,    *DD_DsPhi_H1_Fraction));
              yield_DD[*t][Helbin2][*ds][*ph][DsPhiSide][*m][*c][*a]       = new RooFormulaVar(Form("yield_DD_total_%s_%s_%s",DsPhiSide.c_str(), (*m).c_str(), cat_name.c_str()), Form("DD yield %s %s",  Bmod[DsPhiSide].c_str(), mod[(*m)].c_str()),  "@0*(@1/@2)",RooArgList(*yield_DD[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*DD_DsPhiSide_H2_Fraction,*DD_DsPhi_H1_Fraction));
              
            }
          }
        }
      }
    }
  }

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
        for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
          for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
              

              for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 
                cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(both).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());                 
                yield_comb[*t][both][*ds][*ph][*b][*m][*c][*a] = new RooRealVar(   Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), Yield_comb[*b][*m][both], 0, 600000);                

                if(false){
                  cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(Helbin1).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());  
                  yield_comb[*t][Helbin1][*ds][*ph][*b][*m][*c][*a]    = new RooFormulaVar(Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), "@0*@1",    RooArgList(*Comb_Helbin1_Fraction,*yield_comb[*t][both][*ds][*ph][*b][*m][*c][*a]));
              
                  cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(Helbin2).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());  
                  yield_comb[*t][Helbin2][*ds][*ph][*b][*m][*c][*a]    = new RooFormulaVar(Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), "(1-@0)*@1",RooArgList(*Comb_Helbin1_Fraction,*yield_comb[*t][both][*ds][*ph][*b][*m][*c][*a]));
                
                } else {  
                  cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(Helbin1).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());  
                  yield_comb[*t][Helbin1][*ds][*ph][*b][*m][*c][*a] = new RooRealVar(   Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), Yield_comb[*b][*m][Helbin1], 0, 600000);
                  
                  cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(Helbin2).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());  
                  yield_comb[*t][Helbin2][*ds][*ph][*b][*m][*c][*a] = new RooRealVar(   Form("yield_comb_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("Yield Comb. %s %s",   Bmod[*b].c_str(),mod[*m].c_str()), Yield_comb[*b][*m][Helbin2], 0, 600000);         

                }

              if(allConst_exYield) ((RooRealVar*)yield_comb[*t][both][*ds][*ph][*b][*m][*c][*a])->setConstant();
              if(allConst_exYield) ((RooRealVar*)yield_comb[*t][Helbin1][*ds][*ph][*b][*m][*c][*a])->setConstant();
              if(allConst_exYield) ((RooRealVar*)yield_comb[*t][Helbin2][*ds][*ph][*b][*m][*c][*a])->setConstant();
              }


              for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){
                
                for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){ 

                  cat_name = Form("%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str());
                  sig_frac[*t][*h][*ds][*ph][*b][*m]             = new RooRealVar(   Form("Sigma_frac_%s_%s_%s",(*b).c_str(),(*m).c_str(), cat_name.c_str()), Form("%s %s Sigma Fraction",Bmod[*b].c_str(),mod[*m].c_str()), Fixed_CB_fraction[*b][*m]); 
                }

                cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*c).c_str(),(*a).c_str());

                frac[*t][both][*ds][*ph][DsD0][*m][*c][*a]            = Fraction_DsstD0;
                frac[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]         = Fraction_DsstD0_Helbin1;
                frac[*t][Helbin2][*ds][*ph][DsD0][*m][*c][*a]         = Fraction_DsstD0_Helbin2;
                
                yield_dsstd0[*t][*h][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsstD0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D_{s}* D^{0} %s",mod[(*m)].c_str()), "@0*@1",     RooArgList(*frac[*t][*h][*ds][*ph][DsD0][*m][*c][*a],  *PR_total_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a])  );
                yield_dsdst0[*t][*h][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("yield_DsDst0_%s_%s", (*m).c_str(), cat_name.c_str()), Form("Yield D_{s} D*^{0} %s",mod[(*m)].c_str()), "(1-@0)*@1", RooArgList(*frac[*t][*h][*ds][*ph][DsD0][*m][*c][*a],  *PR_total_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a])  );
              

                frac_HORNS[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HORNS_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*@1",        RooArgList(*fraction_Dsst_Dspi0_DsPhi_010_adj,*fraction_hel));
                frac_HILL[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]        = new RooFormulaVar(Form("fraction_HILL_%s_%s",       (*m).c_str(),cat_name.c_str()), "", "(1-@0)*@1",    RooArgList(*fraction_Dsst_Dspi0_DsPhi_010_adj,*fraction_hel));
                frac_HILL2[*t][*h][*ds][*ph][DsPhi][*m][*c][*a]       = new RooFormulaVar(Form("fraction_HILL2_%s_%s",      (*m).c_str(),cat_name.c_str()), "", "@0*(1-@1)",    RooArgList(*fraction_Dsst_Dspi0_DsPhi_101_adj,*fraction_hel));
                frac_LITTLEHORNS[*t][*h][*ds][*ph][DsPhi][*m][*c][*a] = new RooFormulaVar(Form("fraction_LITTLEHORNS_%s_%s",(*m).c_str(),cat_name.c_str()), "", "(1-@0)*(1-@1)",RooArgList(*fraction_Dsst_Dspi0_DsPhi_101_adj,*fraction_hel));
              }
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

          } //end of loop over chargeList
        }  //end of loop over magnetList  
      } 
    }
  }


  // --------- Mean B mass ---------------
  //double default_mB=5279.29;

  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){   
    for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
      for(std::vector<std::string>::iterator ds=allDsBDTbinList.begin();ds!=allDsBDTbinList.end();ds++){ 
        for(std::vector<std::string>::iterator ph=allPhiBDTbinList.begin();ph!=allPhiBDTbinList.end();ph++){
          for(std::vector<std::string>::iterator b=allBmodeList.begin();b!=allBmodeList.end();b++){  
            for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
              cat_name = Form("%s_%s_%s_%s_%s_%s",(*t).c_str(),(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str());
                
              for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){  
                for(std::vector<std::string>::iterator a=allmagnetList.begin();a!=allmagnetList.end();a++){
                  mean_B[*t][*h][*ds][*ph][*b][*m][*c][*a] = global_mean;
                  if(allConst_exYield) mean_B[*t][*h][*ds][*ph][*b][*m][*c][*a]->setConstant();
                }
              }

              /*
              mean_B[*t][*h][*ds][*ph][*b][*m][both][both] = new RooRealVar(Form("mean_B_%s"           ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][up]   = new RooRealVar(Form("mean_B_%s_plus_up"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][dn]  = new RooRealVar(Form("mean_B_%s_minus_dn"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][dn]   = new RooRealVar(Form("mean_B_%s_plus_dn"   ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][up]  = new RooRealVar(Form("mean_B_%s_minus_up"  ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][plus][both] = new RooRealVar(Form("mean_B_%s_plus_both" ,cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              mean_B[*t][*h][*ds][*ph][*b][*m][minus][both]= new RooRealVar(Form("mean_B_%s_minus_both",cat_name.c_str()), "Mean B mass", default_mB, 5270, 5290);
              if(allConst_exYield) mean_B[*t][*h][*ds][*ph][*b][*m][both][both]->setConstant();
              if(allConst_exYield) mean_B[*t][*h][*ds][*ph][*b][*m][minus][both]->setConstant();

              mean_B[*t][*h][*ds][*ph][*b][*m][both][up]		      = mean_B[*t][*h][*ds][*ph][*b][*m][both][both];
              mean_B[*t][*h][*ds][*ph][*b][*m][both][dn]		      = mean_B[*t][*h][*ds][*ph][*b][*m][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2KPiPi][both][both]	= mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2PiPiPi][both][both] = mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
              mean_B[*t][*h][*ds][*ph][*b][Ds2PhiPi][both][both]  = mean_B[*t][*h][*ds][*ph][*b][Ds2KKPi][both][both];
              */
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
          for(std::vector<std::string>::iterator a=allmagnetList.begin();a!=allmagnetList.end();a++){
            for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){
              //for(std::vector<std::string>::iterator h=allHelbinList.begin();h!=allHelbinList.end();h++){ 
                
              cat_name = Form("%s_%s_%s_%s_%s_%s_%s",(*t).c_str(),both.c_str(),(*ds).c_str(),(*ph).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());

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
              
              cat_name = Form("%s_%s_%s_%s_%s_%s_%s",(*t).c_str(),Helbin1.c_str(),(*ds).c_str(),(*ph).c_str(),(*m).c_str(),(both).c_str(),(*a).c_str());

              sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   = new RooRealVar(   Form("sigma_DsD0_%s", cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
              if(allConst_exYield) ((RooRealVar*)sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a])->setConstant();
              sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma_DsPhi_%s",cat_name.c_str()) ,Form("%s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a])); 
              sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
              //sig_ratio[*t][Helbin1][*ds][*ph][*b][*m]         = new RooRealVar(Form("Sigma_ratio_%s", cat_name.c_str()), Form("%s %s Sigma ratio ",mod[*m].c_str(),Bmod[*b].c_str()), Fixed_CB_Sigma_ratio[*b][*m]);
      
              sigma2[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   = new RooFormulaVar(Form("sigma2_DsD0_%s", cat_name.c_str()), Form("Sigma 2 DsD0 %s",(*m).c_str()),  Form("@0/%f",Fixed_CB_Sigma_ratio[DsD0][*m] ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]   ));
              sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  = new RooFormulaVar(Form("sigma2_DsPhi_%s",cat_name.c_str()), Form("Sigma 2 DsPhi %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DsPhi][*m]), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]  ));
              sigma2[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
              
              avg_sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]      = new RooFormulaVar(Form("avg_sigma_DsD0_%s",  cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsD0].c_str(), mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsD0][*m], Fixed_CB_fraction[DsD0][*m]  ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a], *sigma2[*t][Helbin1][*ds][*ph][DsD0][*m][*c][*a]  ));
              avg_sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a]     = new RooFormulaVar(Form("avg_sigma_DsPhi_%s", cat_name.c_str()), Form("Average sigma %s %s",Bmod[DsPhi].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsPhi][*m],Fixed_CB_fraction[DsPhi][*m] ), RooArgList(*sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a],*sigma2[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a] ));
              avg_sigma[*t][Helbin1][*ds][*ph][DsPhiSide][*m][*c][*a] = avg_sigma[*t][Helbin1][*ds][*ph][DsPhi][*m][*c][*a];
            

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


  // --------- Signal width simple ---------------

  for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
    sigma_simple[DsD0][*m]          = new RooRealVar(   Form("sigma_DsD0_%s", (*m).c_str()) ,Form("DsD0  %s sigma",mod[*m].c_str()),Inital_Sigma_values[*m],  2,  30);
    sigma_simple[DsPhi][*m]         = new RooFormulaVar(Form("sigma_DsPhi_%s",(*m).c_str()) ,Form("DsPhi %s sigma",mod[*m].c_str()),Form("%f*@0",Fixed_Norm_Sigma_ratio[*m]) ,RooArgList(*sigma_simple[DsD0][*m])); 
    sigma_simple[DsPhiSide][*m]     = sigma_simple[DsPhi][*m];
    
    sigma2_simple[DsD0][*m]         = new RooFormulaVar(Form("sigma2_DsD0_%s", (*m).c_str()), Form("Sigma 2 DsD0  %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DsD0][*m] ), RooArgList(*sigma_simple[DsD0][*m]   ));
    sigma2_simple[DsPhi][*m]        = new RooFormulaVar(Form("sigma2_DsPhi_%s",(*m).c_str()), Form("Sigma 2 DsPhi %s",(*m).c_str()), Form("@0/%f",Fixed_CB_Sigma_ratio[DsPhi][*m]), RooArgList(*sigma_simple[DsPhi][*m]  ));
    sigma2_simple[DsPhiSide][*m]    = sigma2_simple[DsPhi][*m];

    avg_sigma_simple[DsD0][*m]      = new RooFormulaVar(Form("avg_sigma_DsD0_%s",  (*m).c_str()), Form("Average sigma %s %s",Bmod[DsD0].c_str(), mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsD0][*m], Fixed_CB_fraction[DsD0][*m]  ), RooArgList(*sigma_simple[DsD0][*m], *sigma2_simple[DsD0][*m]  ));
    avg_sigma_simple[DsPhi][*m]     = new RooFormulaVar(Form("avg_sigma_DsPhi_%s", (*m).c_str()), Form("Average sigma %s %s",Bmod[DsPhi].c_str(),mod[*m].c_str()),  Form("(@0*%f + @1*(1-%f))",Fixed_CB_fraction[DsPhi][*m],Fixed_CB_fraction[DsPhi][*m] ), RooArgList(*sigma_simple[DsPhi][*m],*sigma2_simple[DsPhi][*m] ));
    avg_sigma_simple[DsPhiSide][*m] = avg_sigma_simple[DsPhi][*m];
  
    sig_frac_simple[DsD0][*m]       = new RooRealVar(   Form("Sigma_frac_%s_%s",DsD0.c_str(),     (*m).c_str()), Form("%s %s Sigma Fraction",Bmod[DsD0].c_str(),     mod[*m].c_str()), Fixed_CB_fraction[DsD0][*m]); 
    sig_frac_simple[DsPhi][*m]      = new RooRealVar(   Form("Sigma_frac_%s_%s",DsPhi.c_str(),    (*m).c_str()), Form("%s %s Sigma Fraction",Bmod[DsPhi].c_str(),    mod[*m].c_str()), Fixed_CB_fraction[DsPhi][*m]); 
    sig_frac_simple[DsPhiSide][*m]  = new RooRealVar(   Form("Sigma_frac_%s_%s",DsPhiSide.c_str(),(*m).c_str()), Form("%s %s Sigma Fraction",Bmod[DsPhiSide].c_str(),mod[*m].c_str()), Fixed_CB_fraction[DsPhiSide][*m]); 

    if(allConst_exYield) ((RooRealVar*)sigma_simple[DsD0][*m])->setConstant();

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
              //comb_slope[*t][*h][*ds][*ph][*b][*m] = comb_slope_test[*b][*h];

              
              if(varyAllowed&&par->variation[fixedBG_slope]) comb_slope[*t][*h][*ds][*ph][*b][*m] = comb_slope_sys[*b][*m];
              
              //comb_slope[*t][*h][*ds][*ph][*b][*m] = new RooReacomb_slope[*t][*h][*ds][*ph][*b][*m] = lVar(Form("comb_slope_%s",cat_name.c_str())   ,Form("%s comb. slope", mod[*b].c_str()), -0.00112, -1.0, -0.0000000001);
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
            HORNS_a[*t][*h][*ds][*ph][DsD0][*m]      = DsstarD0_pi0_a;
            HORNS_b[*t][*h][*ds][*ph][DsD0][*m]      = DsstarD0_pi0_b;
            //HORNS_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HORNS_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma_simple[DsD0][*m];
            HORNS_R[*t][*h][*ds][*ph][DsD0][*m]      = global_R;
            HORNS_f[*t][*h][*ds][*ph][DsD0][*m]      = global_f;
            HORNS_csi[*t][*h][*ds][*ph][DsD0][*m]    = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsD0][*m]  = global_shift;
            
            // Ds*D0 gamma
            HILL_a[*t][*h][*ds][*ph][DsD0][*m]       = DsstarD0_gamma_a;
            HILL_b[*t][*h][*ds][*ph][DsD0][*m]       = DsstarD0_gamma_b;
            HILL_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma_simple[DsD0][*m];
            //HILL_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HILL_R[*t][*h][*ds][*ph][DsD0][*m]       = global_R;
            HILL_f[*t][*h][*ds][*ph][DsD0][*m]       = global_f; 
            HILL_csi[*t][*h][*ds][*ph][DsD0][*m]     = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsD0][*m]   = global_shift;
            
            // DsD*0 pi0
            HORNS2_a[*t][*h][*ds][*ph][DsD0][*m]      = DsDstar0_pi0_a;
            HORNS2_b[*t][*h][*ds][*ph][DsD0][*m]      = DsDstar0_pi0_b;
            //HORNS2_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HORNS2_sigma[*t][*h][*ds][*ph][DsD0][*m]  = avg_sigma_simple[DsD0][*m];
            HORNS2_R[*t][*h][*ds][*ph][DsD0][*m]      = global_R;
            HORNS2_f[*t][*h][*ds][*ph][DsD0][*m]      = global_f;
            HORNS2_csi[*t][*h][*ds][*ph][DsD0][*m]    = global_csi;
            HORNS2_shift[*t][*h][*ds][*ph][DsD0][*m]  = global_shift;
            
            // DsD*0 gamma
            HILL2_a[*t][*h][*ds][*ph][DsD0][*m]       = DsDstar0_gamma_a;
            HILL2_b[*t][*h][*ds][*ph][DsD0][*m]       = DsDstar0_gamma_b;
            //HILL2_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HILL2_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma_simple[DsD0][*m];
            HILL2_R[*t][*h][*ds][*ph][DsD0][*m]       = global_R;
            HILL2_f[*t][*h][*ds][*ph][DsD0][*m]       = global_f;
            HILL2_csi[*t][*h][*ds][*ph][DsD0][*m]     = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsD0][*m]   = global_shift;            // PartReco shapes for DsD0 

            

            // Ds*D0* shape approximation 
            HILL_3_a[*t][*h][*ds][*ph][DsD0][*m]       = DsstarDstar0_gamma_pi0_a;
            HILL_3_b[*t][*h][*ds][*ph][DsD0][*m]       = DsstarDstar0_gamma_pi0_b;
            //HILL_3_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma[*t][*h][*ds][*ph][DsD0][*m][both][both];
            HILL_3_sigma[*t][*h][*ds][*ph][DsD0][*m]   = avg_sigma_simple[DsD0][*m];
            HILL_3_R[*t][*h][*ds][*ph][DsD0][*m]       = global_R;
            HILL_3_f[*t][*h][*ds][*ph][DsD0][*m]       = global_f;
            HILL_3_csi[*t][*h][*ds][*ph][DsD0][*m]     = global_csi_hill;
            HILL_3_shift[*t][*h][*ds][*ph][DsD0][*m]   = global_shift;            // PartReco shapes for DsD0 


            // ------------------------------------------
            // PartReco shapes for Ds*Phi
            // ------------------------------------------

            // missed pi0 components 
            // pi0 010
            HORNS_a[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_pi0_a;
            HORNS_b[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_pi0_b;
            //HORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both];
            HORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma_simple[DsPhi][*m];
            HORNS_R[*t][*h][*ds][*ph][DsPhi][*m]      = global_R;
            HORNS_f[*t][*h][*ds][*ph][DsPhi][*m]      = global_f;
            HORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;
            
            //pi0 101
            HILL2_a[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_pi0_a;
            HILL2_b[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_pi0_b;
            //HILL2_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both];
            HILL2_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma_simple[DsPhi][*m];
            HILL2_R[*t][*h][*ds][*ph][DsPhi][*m]      = global_R;
            HILL2_f[*t][*h][*ds][*ph][DsPhi][*m]      = global_f;
            HILL2_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;
            
            // Missed gamma components 
            //gamma 010
            HILL_a[*t][*h][*ds][*ph][DsPhi][*m]       = DsstarPhi_gamma_a;
            HILL_b[*t][*h][*ds][*ph][DsPhi][*m]       = DsstarPhi_gamma_b;
            //HILL_sigma[*t][*h][*ds][*ph][DsPhi][*m]   = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both]; 
            HILL_sigma[*t][*h][*ds][*ph][DsPhi][*m]   = avg_sigma_simple[DsPhi][*m]; 
            HILL_R[*t][*h][*ds][*ph][DsPhi][*m]       = global_R;
            HILL_f[*t][*h][*ds][*ph][DsPhi][*m]       = global_f;
            HILL_csi[*t][*h][*ds][*ph][DsPhi][*m]     = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsPhi][*m]   = global_shift;

            //gamma 101
            LITTLEHORNS_a[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_gamma_a;
            LITTLEHORNS_b[*t][*h][*ds][*ph][DsPhi][*m]      = DsstarPhi_gamma_b;
            //LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both]; 
            LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma_simple[DsPhi][*m]; 
            LITTLEHORNS_R[*t][*h][*ds][*ph][DsPhi][*m]      = global_R;
            LITTLEHORNS_f[*t][*h][*ds][*ph][DsPhi][*m]      = global_f;
            LITTLEHORNS_g[*t][*h][*ds][*ph][DsPhi][*m]      = global_G;
            LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi;
            LITTLEHORNS_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift;

            // ------------------------------------------
            // PartReco shapes for Ds*Phi Sideband
            // ------------------------------------------

            // missed pi0 components 
            // pi0 010
            HORNS_a[*t][*h][*ds][*ph][DsPhiSide][*m]            = DsstarPhi_pi0_a;
            HORNS_b[*t][*h][*ds][*ph][DsPhiSide][*m]            = DsstarPhi_pi0_b;
            //HORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]        = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both];
            HORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]        = avg_sigma_simple[DsPhiSide][*m];
            HORNS_R[*t][*h][*ds][*ph][DsPhiSide][*m]            = global_R;
            HORNS_f[*t][*h][*ds][*ph][DsPhiSide][*m]            = global_f;
            HORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]          = global_csi;
            HORNS_shift[*t][*h][*ds][*ph][DsPhiSide][*m]        = global_shift;
            
            //pi0 101
            HILL2_a[*t][*h][*ds][*ph][DsPhiSide][*m]            = DsstarPhi_pi0_a;
            HILL2_b[*t][*h][*ds][*ph][DsPhiSide][*m]            = DsstarPhi_pi0_b;
            //HILL2_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]        = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both];
            HILL2_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]        = avg_sigma_simple[DsPhiSide][*m];
            HILL2_R[*t][*h][*ds][*ph][DsPhiSide][*m]            = global_R;
            HILL2_f[*t][*h][*ds][*ph][DsPhiSide][*m]            = global_f;
            HILL2_csi[*t][*h][*ds][*ph][DsPhiSide][*m]          = global_csi_hill;
            HILL2_shift[*t][*h][*ds][*ph][DsPhiSide][*m]        = global_shift;
            
            // Missed gamma components 
            //gamma 010
            HILL_a[*t][*h][*ds][*ph][DsPhiSide][*m]             = DsstarPhi_gamma_a;
            HILL_b[*t][*h][*ds][*ph][DsPhiSide][*m]             = DsstarPhi_gamma_b;
            //HILL_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]         = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both]; 
            HILL_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]         = avg_sigma_simple[DsPhiSide][*m]; 
            HILL_R[*t][*h][*ds][*ph][DsPhiSide][*m]             = global_R;
            HILL_f[*t][*h][*ds][*ph][DsPhiSide][*m]             = global_f;
            HILL_csi[*t][*h][*ds][*ph][DsPhiSide][*m]           = global_csi_hill;
            HILL_shift[*t][*h][*ds][*ph][DsPhiSide][*m]         = global_shift;

            //gamma 101
            LITTLEHORNS_a[*t][*h][*ds][*ph][DsPhiSide][*m]      = DsstarPhi_gamma_a;
            LITTLEHORNS_b[*t][*h][*ds][*ph][DsPhiSide][*m]      = DsstarPhi_gamma_b;
            //LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both]; 
            LITTLEHORNS_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma_simple[DsPhiSide][*m]; 
            LITTLEHORNS_R[*t][*h][*ds][*ph][DsPhiSide][*m]      = global_R;
            LITTLEHORNS_f[*t][*h][*ds][*ph][DsPhiSide][*m]      = global_f;
            LITTLEHORNS_g[*t][*h][*ds][*ph][DsPhiSide][*m]      = global_G;
            LITTLEHORNS_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = global_csi;
            LITTLEHORNS_shift[*t][*h][*ds][*ph][DsPhiSide][*m]  = global_shift;


            // ===============

            HORNS_DsstKK_a[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsstKK_a_%s",     cat_name.c_str()), "HORNS_DsstKK a",  5051.478 );
            HORNS_DsstKK_b[*t][*h][*ds][*ph][DsPhi][*m]      = new RooRealVar(Form("HORNS_DsstKK_b_%s",     cat_name.c_str()), "HORNS_DsstKK b",  5128.577 );
            //HORNS_DsstKK_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhi][*m][both][both];
            HORNS_DsstKK_sigma[*t][*h][*ds][*ph][DsPhi][*m]  = avg_sigma_simple[DsPhi][*m];
            HORNS_DsstKK_R[*t][*h][*ds][*ph][DsPhi][*m]      = global_R;
            HORNS_DsstKK_f[*t][*h][*ds][*ph][DsPhi][*m]      = global_f;
            HORNS_DsstKK_csi[*t][*h][*ds][*ph][DsPhi][*m]    = global_csi;
            HORNS_DsstKK_shift[*t][*h][*ds][*ph][DsPhi][*m]  = global_shift_test;
            
            HORNS_DsstKK_a[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsstKK_a_%s",     cat_name.c_str()), "HORNS_DsstKK a",  5051.478 );
            HORNS_DsstKK_b[*t][*h][*ds][*ph][DsPhiSide][*m]      = new RooRealVar(Form("HORNS_DsstKK_b_%s",     cat_name.c_str()), "HORNS_DsstKK b",  5128.577 );
            //HORNS_DsstKK_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma[*t][*h][*ds][*ph][DsPhiSide][*m][both][both];
            HORNS_DsstKK_sigma[*t][*h][*ds][*ph][DsPhiSide][*m]  = avg_sigma_simple[DsPhiSide][*m];
            HORNS_DsstKK_R[*t][*h][*ds][*ph][DsPhiSide][*m]      = global_R;
            HORNS_DsstKK_f[*t][*h][*ds][*ph][DsPhiSide][*m]      = global_f;
            HORNS_DsstKK_csi[*t][*h][*ds][*ph][DsPhiSide][*m]    = global_csi;
            HORNS_DsstKK_shift[*t][*h][*ds][*ph][DsPhiSide][*m]  = global_shift_test;
            
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
   
            
            //HILL_a[*m]      = new RooRealVar(Form("HILL_a_%s",(*m).c_str()),     "HILL a",      4779.5,   4700,  4900);
            //HILL_b[*m]      = new RooRealVar(Form("HILL_b_%s",(*m).c_str()),     "HILL b",      5231.9,   5127,  5327);
            //HILL_csi[*m]    = new RooRealVar(Form("HILL_csi_%s",(*m).c_str()),   "HILL csi",     0.807,    -10,    10);
            //HILL_shift[*m]  = new RooRealVar(Form("HILL_shift_%s",(*m).c_str()), "HILL shift",   -14.2,   -100,   100);
            //HILL_sigma[*m]  = new RooRealVar(Form("HILL_sigma_%s",(*m).c_str()), "HILL sigma",   0.918,      0,    20);
            //HILL_R[*m]      = new RooRealVar(Form("HILL_R_%s",(*m).c_str()),     "HILL R",        11.3,      0,    20);
            //HILL_f[*m]      = new RooRealVar(Form("HILL_f_%s",(*m).c_str()),     "HILL f",       0.882,      0,     1);
            
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
              //comb_slope[*t][*h][*ds][*ph][*b][*m] = comb_slope[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()];

              for(std::vector<std::string>::iterator a=allmagnetList.begin();a!=allmagnetList.end();a++){
                for(std::vector<std::string>::iterator c=allchargeList.begin();c!=allchargeList.end();c++){
                  
                  // Use same mean value for all plots 
                  //mean_B[*t][*h][*ds][*ph][*b][*m][*c][*a] = mean_B[*typeList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()]; 
                  
                  // Allow sigma to vary between Ds decay modes
                  sigma[*t][*h][*ds][*ph][*b][*m][*c][*a]     = sigma[*t][both][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  sigma2[*t][*h][*ds][*ph][*b][*m][*c][*a]    = sigma2[*t][both][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  avg_sigma[*t][*h][*ds][*ph][*b][*m][*c][*a] = avg_sigma[*t][both][*ds][*ph][*b][*m][*chargeList.begin()][*a];

                  Low_Mass_total[*t][*h][*ds][*ph][*b][*m][*c][*a] = Low_Mass_total[*t][both][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  
                  sigma[*t][*h][*ds][*ph][*b][*m][*c][*a]     = sigma[*t][*h][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  sigma2[*t][*h][*ds][*ph][*b][*m][*c][*a]    = sigma2[*t][*h][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  avg_sigma[*t][*h][*ds][*ph][*b][*m][*c][*a] = avg_sigma[*t][*h][*ds][*ph][*b][*m][*chargeList.begin()][*a];

                  Low_Mass_total[*t][*h][*ds][*ph][*b][*m][*c][*a] = Low_Mass_total[*t][*h][*ds][*ph][*b][*m][*chargeList.begin()][*a];


                  yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a]     = yield_comb[*t][*h][*ds][*ph][*b][*m][*chargeList.begin()][*a];
                  B_yield[*t][*h][*ds][*ph][DsD0][*m][*c][*a]      = B_yield[*t][*h][*ds][*ph][DsD0][*m][*chargeList.begin()][*a];
                  
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
  RooAbsPdf *pdf_HILL_3 =0;
  RooAbsPdf *pdf_comb = 0;
  RooAbsPdf *pdf_Dsa1 = 0;
  RooAbsPdf *pdf_DsstD0 = 0;
  RooAbsPdf *pdf_DsDst0 = 0;
  RooAbsPdf *pdf_KKst = 0;
  RooAbsPdf *pdf_DsstKKst = 0;
  RooAbsPdf *pdf_DsstKK = 0;
  RooAbsPdf *pdf_DsD = 0;

  RooAbsPdf *pdf_Bs2DsDs   = 0; 
  RooAbsPdf *pdf_Bs2DsstDs = 0; 
  RooAbsPdf *pdf_B02DsD    = 0;
  RooAbsPdf *pdf_DD        = 0; 


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
                  //RooAbsReal* sig_CB  = sigma[*t][*h][*ds][*ph][*b][*m][*c][*a];
                  RooAbsReal* sig_CB  = sigma_simple[*b][*m];
                  RooAbsReal* n_CB    = ncb[*t][*h][*ds][*ph][*b][*m];
                  RooAbsReal* alp_CB  = alpha[*t][*h][*ds][*ph][*b][*m];
                  pdf_peak_cb1  = new RooCBShape( Form("pdf_peak_cb1_%s",tag.c_str()), "", *mB, *mub, *sig_CB,  *alp_CB, *n_CB );

                  //RooAbsReal* sig2_CB = sigma2[*t][*h][*ds][*ph][*b][*m][*c][*a];
                  RooAbsReal* sig2_CB = sigma2_simple[*b][*m];
                  pdf_peak_cb2  = new RooCBShape( Form("pdf_peak_cb2_%s",tag.c_str()), "", *mB, *mub, *sig2_CB, *alp_CB, *n_CB );
                  
                  //pdf_peak      = new RooAddPdf(  Form("pdf_peak_%s",    tag.c_str()), "", RooArgSet(*pdf_peak_cb1,*pdf_peak_cb2),*sig_frac[*t][*h][*ds][*ph][*b][*m]);
                  pdf_peak      = new RooAddPdf(  Form("pdf_peak_%s",    tag.c_str()), "", RooArgSet(*pdf_peak_cb1,*pdf_peak_cb2),*sig_frac_simple[*b][*m]);
                  
                  //comb bkg
                  RooRealVar* comb = comb_slope[*t][*h][*ds][*ph][*b][*m];
                  
                  if(varyAllowed&&par->variation[fixedBG_Comb_shape]){
                    RooPolynomial*  comb_flat  = new RooPolynomial(  Form("pdf_comb_flat_%s",tag.c_str()) , "", *mB );
                    RooExponential* comb_exp   = new RooExponential( Form("pdf_comb_exp_%s", tag.c_str()), "",*mB,*comb);
                    pdf_comb = new RooAddPdf(Form("pdf_comb_%s",tag.c_str()), "", RooArgSet(*comb_flat,*comb_exp),*comb_flat_fraction);
                  }else if(varyAllowed&&par->variation[fixedBG_Comb_shape_flat]){
                    pdf_comb = new RooPolynomial(  Form("pdf_comb_%s",tag.c_str()) , "", *mB );
                
                  }else if(varyAllowed&&par->variation[fixedBG_Comb_shape_line]){
                    pdf_comb  = new RooPolynomial(  Form("pdf_comb_%s",tag.c_str()) , "", *mB ,RooArgList(*comb_line_offset,*comb_line_slope),0);
                    //pdf_comb  = new RooChebychev(   Form("pdf_comb_%s",tag.c_str()) , "", *mB ,RooArgList(*comb_line_offset));
                  
                  } else{
                    pdf_comb = new RooExponential( Form("pdf_comb_%s",tag.c_str()), "",*mB,*comb); 
                  }
  
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

                    // Add Ds*D*0 shape

                    //Ds*D*0 gamma
                    RooAbsReal* HILL_3_var_a     = HILL_3_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_b     = HILL_3_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_csi   = HILL_3_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_shift = HILL_3_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_sigma = HILL_3_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_R     = HILL_3_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HILL_3_var_f     = HILL_3_f[*t][*h][*ds][*ph][*b][*m];
                    
                    if(varyAllowed&&par->variation[fixedBG_DsstDst0_shape]) {
                      RooCBShape* CB_1  = new RooCBShape( Form("pdf_DsstDst0_cb1_%s",tag.c_str()), "", *mB, *DsstDst0_mean, *DsstDst0_sigma1, *DsstDst0_alpha, *DsstDst0_n );
                      RooCBShape* CB_2  = new RooCBShape( Form("pdf_DsstDst0_cb2_%s",tag.c_str()), "", *mB, *DsstDst0_mean, *DsstDst0_sigma2, *DsstDst0_alpha, *DsstDst0_n );
                      
                      pdf_HILL_3        = new RooAddPdf(  Form("pdf_DsstDst0_%s",    tag.c_str()), "", RooArgSet(*CB_1,*CB_2),*DsstDst0_fraction);
                  
                    } else{
                      pdf_HILL_3 = new RooHILLdini(Form("pdf_DsstDst0_%s",tag.c_str()), "", *mB, *HILL_3_var_a, *HILL_3_var_b, *HILL_3_var_csi, *HILL_3_var_shift, *HILL_3_var_sigma, *HILL_3_var_R, *HILL_3_var_f);
                    
                    }



                    RooArgSet pdflist(  *pdf_comb
                                        ,*pdf_peak 
                                        ,*pdf_DsstD0 
                                        ,*pdf_DsDst0 
                                        ,*pdf_HILL_3
                                        );
                    RooArgSet nevents(  *yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*B_yield[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_dsstd0[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_dsdst0[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_DsstarDstar0[*t][*h][*ds][*ph][*b][*m][*c][*a] 
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
                    if(varyAllowed&&par->variation[fixedBG_Dsa1_smear]) {
                      pdf_Dsa1     = new RooFFTConvPdf(*(PartRecoPDF_Bs02Dsa1_Convolved[*m]),     Form("pdf_Dsa1_%s",tag.c_str()) );
                    } else { 
                      pdf_Dsa1     = new RooFFTConvPdf(*PartRecoPDF_Bs02Dsa1_Conv,     Form("pdf_Dsa1_%s",tag.c_str()) );
                      //pdf_Dsa1     = new RooFFTConvPdf(*PartRecoPDF_Bs02Dsa1_Conv_split[*b][*h],     Form("pdf_Dsa1_%s",tag.c_str()) );
                      //PartRecoPDF_Bs02Dsa1_Conv_split
                    }
                    
                    //Partially reconstructed Bs0 -> Ds* K Kst0
                    if(varyAllowed&&par->variation[fixedBG_DsstKKst_smear]) {
                      pdf_DsstKKst = new RooFFTConvPdf(*(PartRecoPDF_Bs02DsstKKst_Convolved[*m]), Form("pdf_DsstKKst_%s",tag.c_str()) );
                    }else{
                      pdf_DsstKKst = new RooFFTConvPdf(*PartRecoPDF_Bs02DsstKKst_Conv, Form("pdf_DsstKKst_%s",tag.c_str()) );
                    }                    
                    
                    pdf_KKst = new RooAddPdf(Form("pdf_KKst_%s", tag.c_str()),"",RooArgSet(*pdf_Dsa1,*pdf_DsstKKst), RooArgSet(*frac_Dsa1_DsstKKst) );

                    
                    // Test Horns shape 

                    RooAbsReal* HORNS_DsstKK_var_a     = HORNS_DsstKK_a[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_b     = HORNS_DsstKK_b[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_csi   = HORNS_DsstKK_csi[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_shift = HORNS_DsstKK_shift[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_sigma = HORNS_DsstKK_sigma[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_R     = HORNS_DsstKK_R[*t][*h][*ds][*ph][*b][*m];
                    RooAbsReal* HORNS_DsstKK_var_f     = HORNS_DsstKK_f[*t][*h][*ds][*ph][*b][*m];
                    pdf_DsstKK = new RooHORNSdini(Form("pdf_DsstKK_%s",tag.c_str()), "", *mB, *HORNS_DsstKK_var_a, *HORNS_DsstKK_var_b, *HORNS_DsstKK_var_csi, *HORNS_DsstKK_var_shift, *HORNS_DsstKK_var_sigma, *HORNS_DsstKK_var_R, *HORNS_DsstKK_var_f);
                
                    pdf_DsD = new RooFFTConvPdf(*(PartRecoPDF_B2DsD_Conv), Form("pdf_DsD_%s",tag.c_str()) );



                    // Make DD shapes 
                    pdf_Bs2DsDs   = new RooFFTConvPdf(*(PartRecoPDF_Bs2DsDs_Conv),   Form("pdf_Bs2DsDs_%s",  tag.c_str()) );
                    pdf_Bs2DsstDs = new RooFFTConvPdf(*(PartRecoPDF_Bs2DsstDs_Conv), Form("pdf_Bs2DsstDs_%s",tag.c_str()) );
                    pdf_B02DsD    = new RooFFTConvPdf(*(PartRecoPDF_B02DsD_Conv),    Form("pdf_B02DsD_%s",   tag.c_str()) );

                    pdf_DD        = new RooAddPdf(Form("pdf_DD_%s", tag.c_str()),"",RooArgSet(*pdf_Bs2DsDs,*pdf_Bs2DsstDs,*pdf_B02DsD), RooArgSet(*DD_DsDs_Fraction,*DD_DsstDs_Fraction) );

                    std::cout << "Made RooKeysPdf" << std::endl;

                    RooArgSet pdflist(   *pdf_peak 
                                        ,*pdf_comb
                                        ,*pdf_PartReco 
                                        //,*pdf_Dsa1, 
                                        ,*pdf_KKst
                                        //,*pdf_DsstKK
                                        //,*pdf_DsD
                                        ,*pdf_DD
                                       );
                                       // *pdf_DsstKKst );
                    std::cout << "Made pdf list" << std::endl;

                    RooArgSet nevents(   *B_yield[*t][*h][*ds][*ph][*b][*m][*c][*a]
                                        ,*yield_comb[*t][*h][*ds][*ph][*b][*m][*c][*a]
                                        ,*PR_total_yield[*t][*h][*ds][*ph][*b][*m][*c][*a]
                                        ,*yield_Dsa1[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        //,*yield_dsstKK[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        //,*yield_DsD[*t][*h][*ds][*ph][*b][*m][*c][*a] 
                                        ,*yield_DD[*t][*h][*ds][*ph][*b][*m][*c][*a]
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
                  //float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_dsa1 =  yield_Dsa1[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_partreco =   PR_total_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  
                  std::cout<< "\n---------------------------"<<std::endl;
                  std::cout<<(*b)<<", "<<(*m)<<", magnet: "<<(*a)<<", year: "<<(*t)<<"HelBin: " << (*h)<<std::endl;
                  //std::cout<< " B Yield: "<<y_peak<<std::endl;
                  //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                  std::cout<< " Comb:    "<<y_comb<<std::endl;
                  std::cout<< " Dsa1:    "<<y_dsa1<<std::endl;
                  std::cout<< " Ds*Phi:  "<<y_partreco<<std::endl;
                  std::cout<< "---------------------------\n"<<std::endl;
                  //B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->printValue();
                } else if((*b)==DsPhiSide){
                  if (needsBlinding) std::cout << "DsPhi Mode is BLIND: Printing blinded yield" << std::endl;
                  float y_peak =  B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  //float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_comb =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  float y_dsa1 =  yield_Dsa1[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                  
                  std::cout<< "\n---------------------------"<<std::endl;
                  std::cout<<(*b)<<", "<<(*m)<<", magnet: "<<(*a)<<", year: "<<(*t)<<", HelBin: " << (*h)<<std::endl;
                  //std::cout<< " B Yield: "<<y_peak<<std::endl;
                  //std::cout<< " Yield peak: "<<y_yp<<std::endl;
                  std::cout<< " Comb:    "<<y_comb<<std::endl;
                  std::cout<< " Dsa1:    "<<y_dsa1<<std::endl;

                  std::cout<< "---------------------------\n"<<std::endl;

                } else {
                  if(par->sumOverCharges){
                    float y_peak   =  B_yield[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    //float y_yp   =  yield_peak[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_comb   =  yield_comb[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_dsdst0 =  yield_dsdst0[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    float y_dsstd0 =  yield_dsstd0[*t][*h][*ds][*ph][*b][*m][both][*a]->getVal();
                    
                    std::cout<< "\n---------------------------"<<std::endl;
                    std::cout<<(*b)<<", "<<(*m)<<", magnet-"<<(*a)<<", year-"<<(*t)<<" :"<<", HelBin: " << (*h)<<std::endl;
                    //std::cout<< " B Yield: "<<y_peak<<std::endl;
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

  // Printing Branching fractions
  
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){ 
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            RooRealVar *DsPhi_temp = (RooRealVar*)Signal_total[*t][both][*ds][*ph][DsPhi][*m][both][*a];
            RooRealVar *DsD0_temp  = (RooRealVar*)yield_peak[*t][both][*ds][*ph][DsD0][*m][both][*a];
            RooRealVar *eff_temp   = (RooRealVar*)eff_ratio_rrv[*m];
            
            double dsphi_yield = DsPhi_temp->getVal();
            double dsd0_yield  = DsD0_temp->getVal();
            double eff_yield   = eff_temp->getVal();

            double dsphi_error_hi = DsPhi_temp->getAsymErrorHi();
            double dsd0_error_hi  = DsD0_temp->getAsymErrorHi();

            double dsphi_error_lo = DsPhi_temp->getAsymErrorLo();
            double dsd0_error_lo  = DsD0_temp->getAsymErrorLo();


            double dsphi_error = DsPhi_temp->getError();
            double dsd0_error  = DsD0_temp->getError();
            double eff_error   = eff_temp->getError();


            double branching_ratio = eff_yield * (dsphi_yield/dsd0_yield);
            double br_err    = branching_ratio * sqrt(pow(eff_error/eff_yield,2) +  pow(dsd0_error/dsd0_yield,2)    + pow(dsphi_error/dsphi_yield,2)    );
            double br_err_hi = branching_ratio * sqrt(pow(eff_error/eff_yield,2) +  pow(dsd0_error_hi/dsd0_yield,2) + pow(dsphi_error_hi/dsphi_yield,2) );
            double br_err_lo = branching_ratio * sqrt(pow(eff_error/eff_yield,2) +  pow(dsd0_error_lo/dsd0_yield,2) + pow(dsphi_error_lo/dsphi_yield,2) ); 

            std::cout<< "\n---------------------------"<<std::endl;
            std::cout<< "Ds mode: " << *m <<std::endl;
            std::cout<< "Branching Fraction: " << branching_ratio << " +/- " << br_err <<" (+ " << br_err_hi << ", -"<<br_err_lo << ")"<<std::endl;  
            std::cout<< "---------------------------"<<std::endl;
            
            if(par->fitBr){
              RooRealVar *Bf_temp = (RooRealVar*)Branching_fraction_all;
              
              std::cout<< "\n---------------------------"<<std::endl;
              std::cout<< "Fitted Branching Fraction: " << Bf_temp->getVal() << " +/- " << Bf_temp->getError() <<" (+ " << Bf_temp->getAsymErrorHi() << ", -"<<Bf_temp->getAsymErrorLo() << ")"<<std::endl;  
              std::cout<< "---------------------------"<<std::endl;
            }else if(par->fitFourBr){
              RooRealVar *Bf_PhiPi_temp  = (RooRealVar*)Branching_fraction[Ds2PhiPi];
              RooRealVar *Bf_KKPi_temp   = (RooRealVar*)Branching_fraction[Ds2KKPi];
              RooRealVar *Bf_KPiPi_temp  = (RooRealVar*)Branching_fraction[Ds2KPiPi];
              RooRealVar *Bf_PiPiPi_temp = (RooRealVar*)Branching_fraction[Ds2PiPiPi];
              
              std::cout<< "\n---------------------------"<<std::endl;
              std::cout<< "Fitted Branching Fraction PhiPi:  " << Bf_PhiPi_temp->getVal()  << " +/- " << Bf_PhiPi_temp->getError()  <<" (+ " << Bf_PhiPi_temp->getAsymErrorHi()  << ", -"<<Bf_PhiPi_temp->getAsymErrorLo()  << ")"<<std::endl;  
              std::cout<< "Fitted Branching Fraction KKPi:   " << Bf_KKPi_temp->getVal()   << " +/- " << Bf_KKPi_temp->getError()   <<" (+ " << Bf_KKPi_temp->getAsymErrorHi()   << ", -"<<Bf_KKPi_temp->getAsymErrorLo()   << ")"<<std::endl;  
              std::cout<< "Fitted Branching Fraction KPiPi:  " << Bf_KPiPi_temp->getVal()  << " +/- " << Bf_KPiPi_temp->getError()  <<" (+ " << Bf_KPiPi_temp->getAsymErrorHi()  << ", -"<<Bf_KPiPi_temp->getAsymErrorLo()  << ")"<<std::endl;  
              std::cout<< "Fitted Branching Fraction PiPiPi: " << Bf_PiPiPi_temp->getVal() << " +/- " << Bf_PiPiPi_temp->getError() <<" (+ " << Bf_PiPiPi_temp->getAsymErrorHi() << ", -"<<Bf_PiPiPi_temp->getAsymErrorLo() << ")"<<std::endl;  
              std::cout<< "---------------------------"<<std::endl;
            }

           
          }
        }
      }
    }
  }


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
