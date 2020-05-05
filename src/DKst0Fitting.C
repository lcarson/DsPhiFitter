#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <math.h>
#include <time.h>
#include <fstream>

#include "TF1.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TIterator.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TKey.h"
#include "TGraph2D.h"
#include "TLegend.h"

#include "RooHist.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "Roo1DTable.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooMCStudy.h"
#include "RooGenericPdf.h"

#include "DKst0Fitting.h"
#include "Parameters.h"
#include "CommonTools.h"
//#include "MultipleCandidates.h"

//**********************
//***----BLINDING----***
//**********************
//Bool_t isBlind2(kTRUE); //-> just using state==BLIND

DKst0Fitting::DKst0Fitting(Parameters* p, TApplication* app)
: Base()
, toys("toys")
, inputlist("contents of Final ntuple")
, fulllist("contents, including categories")
, storelist("contents of roodataset")
, reducedlist("reduced content")
  //, eventNumber("eventNumber","eventNumber",0,1e99)
, runNumber("runNumber","runNumber",0,1e99)
//, mB("Bu_D0constDconstPVconst_Bu_M","B mass with constraints", 5050 , 5900, "MeV/#font[12]{c}^{2}")
//, mB("Bu_DTF_M","B mass with constraints", 5050 , 5900, "MeV/#font[12]{c}^{2}")
, mB("Bu_DTF_M","B mass with constraints", 4900 , 5900, "MeV/#font[12]{c}^{2}")
//, mD0("D0_M","m(D)",1814.84, 1914.84,"MeV/c^{2}")
//, mDs("D_M","m(Ds)",1945, 1990,"MeV/c^{2}")
//, mKst0("Phi_M","m(Phi)",969.461, 1069.461,"MeV/c^{2}")
, mD0("D0_M","m(D)",1800, 1940,"MeV/c^{2}")
, mD("D_M","m(D)",1800, 1950,"MeV/c^{2}")
, mKst0("Kst0_M","m(Kst0)",650.0, 1100.0,"MeV/c^{2}")
  //, bach_dll("Bach_PIDK","#Delta(LL)",-200,200.)
//, bdt("D_BDT","MVA",0.4,0.9)
, BuIPCHI2("Bu_IPCHI2_OWNPV", "",-1000000 , 100000)
, DIPCHI2("D_IPCHI2_OWNPV", "",-1000000 , 100000)
, bdtD(  "D_BDT",        "MVA", -1.1, 1.1)
, bdtgD( "D_BDTG",       "MVA", -1.1, 1.1)
, bdtbD( "D_BDTB",       "MVA", -1.1, 1.1)
, bdtD0(  "D0_BDT",       "MVA", -1.1, 1.1)
, bdtgD0( "D0_BDTG",      "MVA", -1.1, 1.1)
, bdtbD0( "D0_BDTB",      "MVA", -1.1, 1.1)
, bdtKst0( "Kst0_BDT",     "MVA", -1.1, 1.1)
, bdtgKst0("Kst0_BDTG",    "MVA", -1.1, 1.1)
, bdtbKst0("Kst0_BDTB",    "MVA", -1.1, 1.1)
, bdtMC(  "Bu_MC_BDT",    "MVA", -1.1, 1.1)
, bdtgMC( "Bu_MC_BDTG",   "MVA", -1.1, 1.1)
, bdtbMC( "Bu_MC_BDTB",   "MVA", -1.1, 1.1)
// PID variables
, D0_K0_pidk(  "D0_K0_PIDK","PIDK", -100, 100)
, D0_K1_pidk(  "D0_K1_PIDK","PIDK", -100, 100)
, D_K_pidk(    "D_K_PIDK",  "PIDK", -100, 100)
, D_K0_pidk(   "D_K0_PIDK", "PIDK", -100, 100)
, D_K1_pidk(   "D_K1_PIDK", "PIDK", -100, 100)
, D_P_pidk(    "D_P_PIDK",  "PIDK", -100, 100)
, D_P0_pidk(   "D_P0_PIDK", "PIDK", -100, 100)
, D_P1_pidk(   "D_P1_PIDK", "PIDK", -100, 100)
, D_P2_pidk(   "D_P2_PIDK", "PIDK", -100, 100)
//, helicityD0(  "D0_LoKi_LV01",  "Angle",-1.1,1.1)
//, helicityKst0( "Kst0_LoKi_LV01","Angle",-1.1,1.1)
, helicityD0(  "helAngle",  "Angle",-1.1,1.1)
, helicityKst0( "helAngle","Angle",-1.1,1.1)
, mKK(         "D_KK_M",        "Ds_PhiPi_Mass",    0.0,2000.0)
, bid(         "Bu_ID",         "PDG_ID of B cand.",-pow(2,32),pow(2,32))
, L0Hadron_TOS("Bu_L0HadronDecision_TOS","",-0.5,1.5)
, L0_TIS(      "Bu_L0Global_TIS","",-0.5,1.5)
, Hlt2IncPhi_TOS(   "Bu_Hlt2IncPhiDecision_TOS","",-0.5,1.5)
, Hlt2PhiIncPhi_TOS("Bu_Hlt2PhiIncPhiDecision_TOS","",-0.5,1.5)
 // 2011/2012 Trigger Lines
, Hlt1TrackAllL0_TOS(   "Bu_Hlt1TrackAllL0Decision_TOS","",-0.5,1.5) //this is all == 1
, Hlt2Topo2BodyBBDT_TOS("Bu_Hlt2Topo2BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo3BodyBBDT_TOS("Bu_Hlt2Topo3BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo4BodyBBDT_TOS("Bu_Hlt2Topo4BodyBBDTDecision_TOS","",-0.5,1.5)
  // 2015 Trigger Lines
, Hlt1TrackMVA_TOS(   "Bu_Hlt1TrackMVADecision_TOS",   "",-0.5,1.5)
, Hlt1TwoTrackMVA_TOS("Bu_Hlt1TwoTrackMVADecision_TOS","",-0.5,1.5)
, Hlt2Topo2Body_TOS(  "Bu_Hlt2Topo2BodyDecision_TOS",  "",-0.5,1.5)
, Hlt2Topo3Body_TOS(  "Bu_Hlt2Topo3BodyDecision_TOS",  "",-0.5,1.5)
, Hlt2Topo4Body_TOS(  "Bu_Hlt2Topo4BodyDecision_TOS",  "",-0.5,1.5)
, KKPi_D_Veto(  "D_Veto",    "",-0.5,1.5)
, KKPi_Lc_Veto( "Lc_Veto",   "",-0.5,1.5)
, deltaMass(    "deltaMass", "",-0.5,1000000)
, deltaMass2(   "deltaMass2","",-0.5,100000)
, KPiPi_mPiKPi( "D_PiKPi_M", "", 0.0,100000)
, KPiPi_mPiPiPi("D_PiPiPi_M","", 0.0,100000)
, KPiPi_mKKPi(  "D_KKPi_M",  "", 0.0,100000)
, PiKPi_mKKPi(  "D_KKPi_M", "",  0.0,100000)
, PiKPi_mPiKK(  "D_PiKK_M", "",  0.0,100000)
, PiKPi_mKK1(   "D_KK1_M" , "",  0.0,100000)
, PiKPi_mKK2(   "D_KK2_M" , "",  0.0,100000)
, m1245("m1245","",0.0,100000)
, m145("m145","",0.0,100000)
, m245("m245","",0.0,100000)
, m345("m345","",0.0,100000)
, D_FDCHI2("D_FDCHI2_ORIVX","", 0.0,10000.0)
, D0_FDCHI2("D0_FDCHI2_ORIVX","", 0.0,10000.0)
, type("type","Type of data")
, mode("mode","D decay mode")
, Bmode("Bmode","B^{+} decay mode")
, magnet("magnet","dipole polarity")
, charge("charge","batchelor charge")
, bMassRegion("bMassRegion","m(B) band")
, helBin("helBin","Helicity Angle bin")
, m_mD0_Mu(1864.84)
, m_mDs_Mu(1968.30)
, m_mD_Mu(1869.58)
, m_mKst0_Mu(895.6)
, m_mBs0_Mu(5366.477)
, m_bdtCut(-0.9)
, state(1)
, typeList()
, yearList()
, modeList()
, BmodeList()
, chargeList()
, HelBinList()
{
  par=p;
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetNdivisions(205,"XYZ");
  gStyle->SetStatFont(132);
  gStyle->SetStatFontSize(0.08);
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTitleSize(0.077,"XYZ");
  gStyle->SetLabelSize(0.08,"XYZ");
  gStyle->SetTitleOffset(0.83,"X");
  gStyle->SetTitleOffset(1.05,"Y");
  gStyle->SetMarkerSize(0.5);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
 
  DefineModel();
  DefineRooCategories();
  if(par->readToys){
    DisplayToys();
  }else if(par->readSys){
    DisplaySys();
  }else if(par->doSensitivity){
    Sensitivity();
  }else{
    if(par->genToys){
     if(par->nToys==1){
        OrderToys(1);
        RunFullFit(true);
      }else{
        if(par->quickVersion){
          if(LoadDataSet()) return;
          RunFullFit(false);
        }
        RunManyToys();
        return;
      }
    }else{
      if(LoadDataSet()) return;
      //removed if(automatic) and else()
      if (par->manyFits) {
        RunManyFits();
      }else if (par->runEff){
        RunEfficiency();
      }else {
        RunFullFit((!par->batch)||par->variation[draw]); 
      }
    }
  }

  if(!par->batch){
    std::cout<<"Starting: app->Run()"<<std::endl;
    app->Run(true);
  }
}


void DKst0Fitting::DefineRooCategories()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::DefineRooCategories()"<<std::endl;
  
  if(par->splitYears) {
    if(par->dsetsReq[s21])   yearList.push_back(s21);
    if(par->dsetsReq[s21r1]) yearList.push_back(s21r1);
    if(par->dsetsReq[s24])   yearList.push_back(s24);
    if(par->dsetsReq[s26])   yearList.push_back(s26);
    if(0==yearList.size())   yearList.push_back("UNDEFTYPE");
  } else {
    yearList.push_back(both);
  }

  if(par->genToys||par->readToys)typeList.push_back(toy);
  if(par->dsetsReq[s21])   typeList.push_back(s21);
  if(par->dsetsReq[s21r1]) typeList.push_back(s21r1);
  if(par->dsetsReq[s24])   typeList.push_back(s24);
  if(par->dsetsReq[s26])   typeList.push_back(s26);

  if(0==typeList.size())   typeList.push_back("UNDEFTYPE"); 


  type.defineType(toy.c_str());
  type.defineType(s21.c_str());
  type.defineType(s21r1.c_str());
  type.defineType(s24.c_str());
  type.defineType(s26.c_str());
  type.defineType(both.c_str());

  mode.defineType(D2PiKPi.c_str());
  mode.defineType(D2KKPi.c_str());
  mode.defineType(D2PiPiPi.c_str());
  
  Bmode.defineType(DD0.c_str());
  Bmode.defineType(DKst0.c_str());
  Bmode.defineType(DKst0Side.c_str());

  if(par->modes[D2PiKPi])     modeList.push_back(D2PiKPi);
  if(par->modes[D2KKPi])      modeList.push_back(D2KKPi);
  if(par->modes[D2PiPiPi])    modeList.push_back(D2PiPiPi);

  if(par->Bmodes[DD0])        BmodeList.push_back(DD0);
  if(par->Bmodes[DKst0])      BmodeList.push_back(DKst0);
  if(par->Bmodes[DKst0Side])  BmodeList.push_back(DKst0Side);
 

  if(par->sumOverCharges){
    chargeList.push_back(both);
  }else{
    chargeList.push_back(minus);
    chargeList.push_back(plus);
  }
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

  if(par->splitHel) {
    HelBinList.push_back(Helbin1);
    HelBinList.push_back(Helbin2);
  } else {
    HelBinList.push_back(both);
  }
  
  helBin.defineType(Helbin1.c_str());
  helBin.defineType(Helbin2.c_str());
  helBin.defineType(both.c_str());

  bMassRegion.defineType("signal");
  bMassRegion.defineType("bckgrd");


  //Make titles
  std::map<std::string,std::string> chg;
  std::map<std::string,std::string> mag;
  std::map<std::string,std::string> mod;
  std::map<std::string,std::string> Bmod;
  std::map<std::string,std::string> hel;
  std::map<std::string,std::string> dsbdt;
  std::map<std::string,std::string> phibdt;
  std::map<std::string,std::string> year;

  Bmod[DD0]  ="D^{0}";            Bmod[DKst0]  ="K*^{0}";             Bmod[DKst0Side]  ="K*^{0}";
  mod[D2KKPi]="KK#pi";            mod[D2PiKPi] ="#piK#pi";            mod[D2PiPiPi]="#pi#pi#pi";
  chg[plus]  ="+";                chg[minus]   ="#font[122]{-}";      chg[both]    ="#font[122]{#pm}";
  mag[up]    ="UP";               mag[dn]      ="DN";                 mag[both]    ="";
  hel[both]  = "All cos(#theta)"; hel[Helbin1] = "|cos(#theta)|>0.4"; hel[Helbin2] = "|cos(#theta)|<0.4";
  year[s21]  = "2012";            year[s21r1]  = "2011";              year[s24]    = "2015";                
  year[s26]  = "2016";            year[both]   = "All years";
  
  std::vector<std::string> temp_list = yearList;
  if(par->genToys) temp_list = typeList; 
  for(std::vector<std::string>::iterator y=temp_list.begin();y!=temp_list.end();y++){ 
    for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
        for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
          for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
            for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
              title[*y][*h][*b][*m][*c][*a] = std::string(Form("B^{%s}#rightarrowD^{%s}#font[132]{(}#rightarrow%s#font[132]{)}%s#font[132]{(}#rightarrowK#pi#font[132]{)} %s", chg[*c].c_str(), chg[*c].c_str(),mod[*m].c_str(), Bmod[*b].c_str(), mag[*a].c_str()));
              bin_detail[*y][*h][*b][*m][*c][*a] = std::string(Form( "%s %s",year[*y].c_str(),hel[*h].c_str()));
              if(par->genToys) title[*y][*h][*b][*m][*c][*a]+=" TOY"; 
            }
          }
        }
      }
    }
  }
  //inputlist.add(eventNumber);
  //inputlist.add(bdt);
  inputlist.add(bid);
  inputlist.add(mD);
  inputlist.add(mD0);
  inputlist.add(mKst0);
  inputlist.add(BuIPCHI2);
  inputlist.add(DIPCHI2);
  inputlist.add(bdtD0);
  inputlist.add(bdtgD0);
  inputlist.add(bdtbD0);
  inputlist.add(bdtgD);
  inputlist.add(bdtbD);
  inputlist.add(bdtD);
  inputlist.add(bdtgKst0);
  inputlist.add(bdtbKst0);
  inputlist.add(bdtKst0);
  inputlist.add(bdtgMC);
  inputlist.add(bdtbMC);
  inputlist.add(bdtMC);
  inputlist.add(helicityD0);
  inputlist.add(helicityKst0);
  inputlist.add(mKK);

  //inputlist.add(bach_dll);
  inputlist.add(L0Hadron_TOS);
  inputlist.add(L0_TIS);
  inputlist.add(Hlt2IncPhi_TOS);
  inputlist.add(Hlt2PhiIncPhi_TOS);
  
  // 2011/2012 Trigger Lines
  inputlist.add(Hlt1TrackAllL0_TOS);
  inputlist.add(Hlt2Topo2Body_TOS);
  inputlist.add(Hlt2Topo3Body_TOS);
  inputlist.add(Hlt2Topo4Body_TOS);

  // 2015 Trigger Lines
  inputlist.add(Hlt1TrackMVA_TOS);
  inputlist.add(Hlt1TwoTrackMVA_TOS);
  inputlist.add(Hlt2Topo2BodyBBDT_TOS);
  inputlist.add(Hlt2Topo3BodyBBDT_TOS);
  inputlist.add(Hlt2Topo4BodyBBDT_TOS);

  // Veto variables
  inputlist.add(KKPi_D_Veto);
  inputlist.add(KKPi_Lc_Veto);
  inputlist.add(KPiPi_mPiKPi);
  inputlist.add(KPiPi_mPiPiPi);
  inputlist.add(KPiPi_mKKPi);
  inputlist.add(m1245);
  inputlist.add(m145);
  inputlist.add(m245);
  inputlist.add(m345);

  inputlist.add(PiKPi_mKKPi);
  inputlist.add(PiKPi_mPiKK);
  inputlist.add(PiKPi_mKK1);
  inputlist.add(PiKPi_mKK2);

  inputlist.add(deltaMass);
  inputlist.add(deltaMass2);
  // FDChi2

  inputlist.add(D_FDCHI2);
  inputlist.add(D0_FDCHI2);

  //PIDK variables
  inputlist.add(D0_K0_pidk);
  inputlist.add(D0_K1_pidk);
  inputlist.add(D_K_pidk);
  inputlist.add(D_K0_pidk);
  inputlist.add(D_K1_pidk);
  inputlist.add(D_P_pidk);
  inputlist.add(D_P0_pidk);
  inputlist.add(D_P1_pidk);
  inputlist.add(D_P2_pidk);

  storelist.add(mB);
  storelist.add(type);
  storelist.add(mode);
  storelist.add(Bmode);
  storelist.add(magnet);
  storelist.add(charge);
  storelist.add(helBin);

  reducedlist.add(type);
  reducedlist.add(mode);
  reducedlist.add(Bmode);
  reducedlist.add(magnet);
  reducedlist.add(charge);
  reducedlist.add(helBin);



  //reducedlist.add(mDs);
  //reducedlist.add(mD0);
  //reducedlist.add(mKst0);

  //storelist.add(mDs);
  //storelist.add(mD0);
  //storelist.add(mKst0);



  if(par->manyFits || par->runEff){
    reducedlist.add(bdtD0);
    reducedlist.add(bdtgD0);
    reducedlist.add(bdtbD0);
    reducedlist.add(bdtgD);
    reducedlist.add(bdtbD);
    reducedlist.add(bdtD);
    reducedlist.add(bdtgMC);
    reducedlist.add(bdtbMC);
    reducedlist.add(bdtMC);

    reducedlist.add(D0_K0_pidk);
    reducedlist.add(D0_K1_pidk);


    if(par->modes[D2KKPi]) {
      reducedlist.add(D_K0_pidk);
      reducedlist.add(D_K1_pidk);
      reducedlist.add(D_P_pidk);
    } else if(par->modes[D2PiKPi]){
      reducedlist.add(D_K_pidk);
      reducedlist.add(D_P0_pidk);
      reducedlist.add(D_P1_pidk);
    } else if(par->modes[D2PiPiPi]){
      reducedlist.add(D_P0_pidk);
      reducedlist.add(D_P1_pidk);
      reducedlist.add(D_P2_pidk);
    }

    storelist.add(bdtD0);
    storelist.add(bdtgD0);
    storelist.add(bdtbD0);
    storelist.add(bdtgD);
    storelist.add(bdtbD);
    storelist.add(bdtD);
    storelist.add(bdtgMC);
    storelist.add(bdtbMC);
    storelist.add(bdtMC);
    storelist.add(D0_K0_pidk);
    storelist.add(D0_K1_pidk);

    if(par->modes[D2KKPi]){
      storelist.add(D_K0_pidk);
      storelist.add(D_K1_pidk);
      storelist.add(D_P_pidk);
    } else if(par->modes[D2PiKPi]){
      storelist.add(D_K_pidk);
      storelist.add(D_P0_pidk);
      storelist.add(D_P1_pidk); 
    } else if(par->modes[D2PiPiPi]){
      storelist.add(D_P0_pidk);
      storelist.add(D_P1_pidk);
      storelist.add(D_P2_pidk);
    }
  }


  fulllist.add(storelist);
  fulllist.add(inputlist);
  fulllist.add(bMassRegion);

  reducedlist.add(mB);
  inputlist.add(mB);

  
  std::cout<< "inputList: " << std::endl;
  inputlist.Print();

  std::cout<< "reducedList: " << std::endl;
  reducedlist.Print();

  std::cout<< "fullList: " << std::endl;
  fulllist.Print();
  

} //end of funcn DefineRooCategories



std::string DKst0Fitting::reducedFileName(std::string t,std::string m,std::string Bm)
{
  gSystem->Exec("mkdir -p roodatasets");
  std::string reducedfile="roodatasets/";
  //reducedfile="Laurence_RooDataSets/";
  reducedfile+=Bm+underscore;
  reducedfile+=t+underscore;
  reducedfile+=m+underscore;
  reducedfile+=std::string( Form("RooDs_%i_%i.root",int(mB.getMin()),int(mB.getMax())) );//,int(m_bdtCut*10)));
  return reducedfile;
}



int DKst0Fitting::LoadDataSet()
{
  if(par->debug) std::cout<<"Running: LoadDataSet()"<<std::endl;
  data = new RooDataSet("data","Data data",storelist);
  for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){  
    for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++) {
      for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++) {
        std::string fn=reducedFileName(*t,*m,*b);
        RooDataSet* reducedData=0;
        if(par->quickVersion){
          TFile* in = TFile::Open(fn.c_str(),"read");
          if(in){
            std::cout<<"\nReading reduced file: '"<<fn<<"'\n"<<std::endl;
            reducedData = new RooDataSet(*(RooDataSet*)gDirectory->Get("reducedData"));
            in->Close();
          }
        }
        if(!reducedData){
          if(NULL==par->locationoption){
            std::cout<<"\nNo ntuple location specified! (Need the -L <directory> option)"<<std::endl;
            return 1;
          }
          reducedData = MakeDataSet(*t,*m,*b);
          if(NULL==reducedData) return 1;
        }
        data->append(*reducedData);
        delete reducedData;
      }
    }
  }
  AssignGlobalCategory();
  return 0;
}



RooDataSet* DKst0Fitting::MakeDataSet(std::string t,std::string m,std::string b)
{
  if(par->debug) std::cout<<"Running: MakeDataSet(): "<< t <<", "<< m <<", "<< b <<std::endl;
  RooDataSet* fulldata = new RooDataSet("fulldata","Data data",fulllist);
  
  std::vector<const char*> tupleName;
  //tupleName.push_back("TTT");
  tupleName.push_back("DecayTree");

  std::string s_mode(null);
  std::string s_Bmode(null);
  std::string s_type(null);

  for(std::vector<std::string>::iterator il=par->locations.begin();il!=par->locations.end();il++){

    std::vector<std::string> fileName;
    CommonTools::getdir((*il),fileName);
    
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
      std::string fullpathandname = (*il)+slash+(*it);
      if(par->debug) std::cout<<"The fullpathandname is: "<<fullpathandname<<"."<<std::endl; //for debugging
      bool accept=false;
      s_Bmode == "";
      if( (*it).find("D0D")    !=std::string::npos && b == DD0      ) { s_Bmode=DD0;      } 
      if( (*it).find("D0D")    !=std::string::npos && b == DKst0    ) { s_Bmode="";       } 
      if( (*it).find("D0D")    !=std::string::npos && b == DKst0Side) { s_Bmode="";       }    
      if( (*it).find("DKst0")  !=std::string::npos && b == DKst0    ) { s_Bmode=DKst0;    }    
      if( (*it).find("DKst0")  !=std::string::npos && b == DKst0Side) { s_Bmode=DKst0Side;}
      if( (*it).find("DKst0")  !=std::string::npos && b == DD0      ) { s_Bmode="";       }

      if(  (*it).find((m+underscore).c_str())!=std::string::npos && s_Bmode==b){ 
        accept=true;
        s_mode=m;
      }
      
      //Laurence
      if( (*it).find("2011")  !=std::string::npos ) { s_type=s21r1;}
      if( (*it).find("2012")  !=std::string::npos ) { s_type=s21; }
      if( (*it).find("2015")  !=std::string::npos ) { s_type=s24; }
      if( (*it).find("2016")  !=std::string::npos ) { s_type=s26; }



      if(!par->dsetsReq[s_type]) continue; 
      if(t!=s_type) continue;
      
      if(!par->Bmodes[s_Bmode]) continue; 
      if(b!=s_Bmode) continue;


      std::string s_magnet(both);
      if( (*it).find("MagUp")  !=std::string::npos ) { s_magnet=up; }
      if( (*it).find("MagDown")!=std::string::npos ) { s_magnet=dn; }    

      if(accept){
        TFile* tfile = TFile::Open(fullpathandname.c_str());
        if(tfile){
          //Laurence: cd() into folder inside the .root
          std::string tfileFolderName("");
          if (s_Bmode==DD0){
            if( s_mode==D2PiKPi )  { tfileFolderName = "B2D0D_D02KPi_D2PiKPi"; }
            if( s_mode==D2KKPi )   { tfileFolderName = "B2D0D_D02KPi_D2KKPi"; } 
            if( s_mode==D2PiPiPi ) { tfileFolderName = "B2D0D_D02KPi_D2PiPiPi"; }
          } else if(s_Bmode==DKst0||s_Bmode==DKst0Side){
            if( s_mode==D2PiKPi )  { tfileFolderName = "B2DKst0_Kst02KPi_D2PiKPi"; }
            if( s_mode==D2KKPi )   { tfileFolderName = "B2DKst0_Kst02KPi_D2KKPi"; } 
            if( s_mode==D2PiPiPi ) { tfileFolderName = "B2DKst0_Kst02KPi_D2PiPiPi"; }
          }
          const char *pathPlusFolderName = (fullpathandname+colon+slash+tfileFolderName).c_str();
          std::cout<<"The pathPlusFolderName is: "<<pathPlusFolderName<<"."<<std::endl; //for debugging
          tfile->cd(pathPlusFolderName);

          //not sure a loop is necessary, since the .root's only have one tuple in them, but we leave it as it is
          for(std::vector<const char*>::iterator itc=tupleName.begin(); itc!=tupleName.end(); itc++){
            if(gDirectory->FindObjectAny(*itc)){
              std::cout<<"   Reading tuple '"<<(*itc)<<"'  with " << ((TTree*)gDirectory->Get(*itc))->GetEntries() << " events" <<std::endl;
              RooDataSet *input = FinalDataSet(s_Bmode,s_type,s_mode,s_magnet,(TTree*)gDirectory->Get(*itc));
              fulldata->append(*input);
              delete input;
            }
          }
          tfile->Close();
          par->dsetsFound[s_type][s_mode]=true;
        }else{
          std::cout<<"No file at: "<<fullpathandname<<". Skipping"<<std::endl;
        }

      } //end of if(accept)

    } //end of loop over fileNames

  } //end of loop over locations

  if(!par->dsetsFound[t][m]){
    std::cout<<"Required dataset '"<<t<<"' for mode '"<<m<<"' not found."<<std::endl;
    std::cout<<"Abandoning. Check input locations and dataset/mode requirements."<<std::endl;
    return NULL;
  }

  std::string cutString("bMassRegion!=bMassRegion::bckgrd");
  RooDataSet* reducedData=(RooDataSet*)fulldata->reduce(cutString.c_str());
  //RooDataSet* reducedData=(RooDataSet*)reducedData2->reduce(cutString.c_str()); // Hack to easily change BDT cut
  delete fulldata;fulldata=reducedData;

  //IdMultipleCandidates(fulldata);
  //reducedData = (RooDataSet*)fulldata->reduce("clone!=clone::true");
  //std::cout<<"\n   Dataset built and checked for multiple candidates: "<<fulldata->sumEntries()-reducedData->sumEntries()<<" removed."<<std::endl;

  std::cout<<"--> Writing reduced file: '"<<reducedFileName(t,m,b)<<"' with "<<reducedData->sumEntries()<<" events."<<std::endl;
  reducedData = new RooDataSet("reducedData","",reducedData,storelist);
  TFile out(reducedFileName(t,m,b).c_str(),"recreate");
  reducedData->Write();
  out.Close();
  
  delete fulldata;
  return reducedData;
} //end of funcn MakeDataSet




RooDataSet* DKst0Fitting::FinalDataSet(const std::string s_Bmode, const std::string s_type, const std::string s_mode, std::string s_mag, TTree* tree)
{
  if(par->debug) std::cout<<"Running: FinalDataSet()"<<std::endl;
  if(!tree){ std::cout << "\n THE TREE IS A ZERO POINTER IN "<< s_Bmode <<" "<< s_type <<" "<<s_mode<<" ("<< s_mag <<")"<< std::endl; return 0; }
  RooDataSet *input = new RooDataSet(Form("%s_%s_%s",s_Bmode.c_str(),s_type.c_str(),s_mode.c_str()),Form("Final %s %s %s",s_Bmode.c_str(),s_type.c_str(),s_mode.c_str()),inputlist,RooFit::Import(*tree));
  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(type,helBin,Bmode,mode,magnet,charge,bMassRegion)); //removed "pid" from the RooArgSet
 

  int nTot=0;
  int nNotL0=0;
  int nNotL0TOS=0;
  int nNotTopoTOS=0;
  int nNotTrackTOS=0;
  int nTotIncPhi=0;
  int nTotOnlyIncPhi=0;
  for(int i=0;i<input->numEntries();i++){
    Bmode.setLabel(s_Bmode.c_str());
    mode.setLabel(s_mode.c_str());
    type.setLabel(s_type.c_str());
    magnet.setLabel(s_mag.c_str());

    
    const RooArgSet *row = input->get(i);
    //RooRealVar* D   = (RooRealVar*) row->find(bach_dll.GetName());
    
    // Split charges into categories
    RooRealVar* q   = (RooRealVar*) row->find(bid.GetName());
    if(q->getVal() > 0 && q->getVal() < pow(2,31)) charge.setLabel("plus");
    if(q->getVal() < 0 || q->getVal() > pow(2,31)) charge.setLabel("minus");

    // Split helicty angles into bins
    if(s_Bmode==DD0){
      RooRealVar* angle = (RooRealVar*) row->find(helicityD0.GetName());
      if(fabs(angle->getVal()) > 0.4 ) {
        helBin.setLabel(Helbin1.c_str());
      }else {
        helBin.setLabel(Helbin2.c_str());
      }
    }else{
      RooRealVar* angle = (RooRealVar*) row->find(helicityKst0.GetName());
      if(fabs(angle->getVal()) > 0.4 ) {
        helBin.setLabel(Helbin1.c_str());
      }else {
        helBin.setLabel(Helbin2.c_str());
      }
    }
  

    bMassRegion.setLabel("signal");


    if(s_Bmode==DD0){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* d   = (RooRealVar*) row->find(mD0.GetName());
        if(fabs(d->getVal()-m_mD0_Mu)>25) bMassRegion.setLabel("bckgrd");
        //if(fabs(d->getVal()-m_mD0_Mu)>50) bMassRegion.setLabel("bckgrd");
      }
    } else if(s_Bmode==DKst0){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* kst0m = (RooRealVar*) row->find(mKst0.GetName());
        if(fabs(kst0m->getVal()-m_mKst0_Mu)>50) bMassRegion.setLabel("bckgrd");
      }
    } else if(s_Bmode==DKst0Side){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* kst0m = (RooRealVar*) row->find(mKst0.GetName());
        if(fabs(kst0m->getVal()-m_mKst0_Mu)<50) bMassRegion.setLabel("bckgrd");
        if(fabs(kst0m->getVal()-m_mKst0_Mu)>150) bMassRegion.setLabel("bckgrd");
      }
    }

    if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* ds = (RooRealVar*) row->find(mD.GetName());
      if(fabs(ds->getVal()-m_mD_Mu)>25) bMassRegion.setLabel("bckgrd");
    }


    // Remove cross feed candidates
    
    if(s_mode==D2PiKPi){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* mmiss_kkpi = (RooRealVar*) row->find(PiKPi_mKKPi.GetName());
        RooRealVar* mmiss_pikk = (RooRealVar*) row->find(PiKPi_mPiKK.GetName());
        RooRealVar* mmiss_kk1  = (RooRealVar*) row->find(PiKPi_mKK1.GetName());
        RooRealVar* mmiss_kk2  = (RooRealVar*) row->find(PiKPi_mKK2.GetName());
        RooRealVar* p0_pidk    = (RooRealVar*) row->find(D_P0_pidk.GetName());
        RooRealVar* p1_pidk    = (RooRealVar*) row->find(D_P1_pidk.GetName());


        if(fabs(mmiss_kkpi->getVal()-m_mDs_Mu)<25 &&  (p0_pidk->getVal()> -10 || fabs(mmiss_kk1->getVal()-1020.0)<10) ) bMassRegion.setLabel("bckgrd");
        if(fabs(mmiss_pikk->getVal()-m_mDs_Mu)<25 &&  (p1_pidk->getVal()> -10 || fabs(mmiss_kk2->getVal()-1020.0)<10) ) bMassRegion.setLabel("bckgrd");
      }
    }


    /* 
    if(s_mode==D2KKPi){
      RooRealVar* dveto  = (RooRealVar*) row->find(KKPi_D_Veto.GetName());
      RooRealVar* lcveto = (RooRealVar*) row->find(KKPi_Lc_Veto.GetName());
      RooRealVar* deltam = (RooRealVar*) row->find(deltaMass.GetName());
      if (dveto->getVal()  == 1 ) {bMassRegion.setLabel("bckgrd");}
      if (lcveto->getVal() == 1 ) {bMassRegion.setLabel("bckgrd");}
      if (deltam->getVal() < 150 ) {bMassRegion.setLabel("bckgrd");}
    } 
    if(s_mode==D2PiKPi){
      RooRealVar* dveto  = (RooRealVar*) row->find(KKPi_D_Veto.GetName());
      RooRealVar* deltam = (RooRealVar*) row->find(deltaMass.GetName());
      if (dveto->getVal()  == 1 ) {bMassRegion.setLabel("bckgrd");}
      if (deltam->getVal() < 150 ) {bMassRegion.setLabel("bckgrd");}
    }
    if(s_mode==D2PiPiPi){
      RooRealVar* deltam  = (RooRealVar*) row->find(deltaMass.GetName());
      RooRealVar* deltam2 = (RooRealVar*) row->find(deltaMass2.GetName());
      if (deltam->getVal()  < 150) {bMassRegion.setLabel("bckgrd");}
      if (deltam2->getVal() < 150) {bMassRegion.setLabel("bckgrd");}
    }
  */
    /*
    //Veto Bs->phi phi events
    if((s_Bmode==DsPhi||s_Bmode==DsPhiSide)&&(s_mode==Ds2PhiPi||s_mode==Ds2KKPi)){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* mphiphi = (RooRealVar*) row->find(m1245.GetName());
        if(fabs(mphiphi->getVal()-m_mBs0_Mu)<50) bMassRegion.setLabel("bckgrd");
      }
    }

     //Veto Ds->KKPi in Ds->KKPi
    if((s_Bmode==DsPhi||s_Bmode==DsPhiSide)&&(s_mode==Ds2PhiPi||s_mode==Ds2KKPi)){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* mmisskkpi = (RooRealVar*) row->find(m345.GetName());
        if(fabs(mmisskkpi->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
      }
    }

    //Veto Ds->KKPi in Ds->PiPiPi
    if((s_Bmode==DsPhi||s_Bmode==DsPhiSide)&&s_mode==Ds2PiPiPi){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* mmisskkpi_1 = (RooRealVar*) row->find(m145.GetName());
        RooRealVar* mmisskkpi_2 = (RooRealVar*) row->find(m245.GetName());
        RooRealVar* mmisskkpi_3 = (RooRealVar*) row->find(m345.GetName());
        if(fabs(mmisskkpi_1->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_2->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_3->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_1->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_2->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_3->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
      }
    }    

    //Veto Ds->KKPi in Ds->KPiPi
    if((s_Bmode==DsPhi||s_Bmode==DsPhiSide)&&s_mode==Ds2KPiPi){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* mmisskkpi_2 = (RooRealVar*) row->find(m245.GetName());
        RooRealVar* mmisskkpi_3 = (RooRealVar*) row->find(m345.GetName());
        if(fabs(mmisskkpi_2->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_3->getVal()-m_mDs_Mu)<25) bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_2->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
        if(fabs(mmisskkpi_3->getVal()-m_mD_Mu)<25)  bMassRegion.setLabel("bckgrd");
      }
    }
    */


    // Test cut
    
    if(s_Bmode== DKst0 || s_Bmode == DKst0Side){
      RooRealVar* Bu_ipchi2    = (RooRealVar*) row->find(BuIPCHI2.GetName());
      RooRealVar* D_ipchi2     = (RooRealVar*) row->find(DIPCHI2.GetName());
      if(log10(Bu_ipchi2->getVal())>1.0)   bMassRegion.setLabel("bckgrd");
      if(log10(D_ipchi2->getVal() )<1.0)   bMassRegion.setLabel("bckgrd");
    }

    if(!par->manyFits){ // Only apply cuts if not doing optimisation
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* bdd     = (RooRealVar*) row->find(bdtgD.GetName());
        RooRealVar* bdkst0  = (RooRealVar*) row->find(bdtgKst0.GetName());;
        RooRealVar* bdd0    = (RooRealVar*) row->find(bdtgD0.GetName());
        if(s_Bmode==DD0      && bdd0->getVal()  < 0.7) {bMassRegion.setLabel("bckgrd");} //0.85
        if(s_Bmode==DKst0    && bdkst0->getVal()< 0.7) {bMassRegion.setLabel("bckgrd");}
        if(s_Bmode==DKst0Side&& bdkst0->getVal()< 0.7) {bMassRegion.setLabel("bckgrd");}

        // Main selection [BDT] cuts. Values chosen just by looking at D_BDT distribution, and removing a few events in the tail.
        // Same Phi/D0 Cut for all modes
        if(s_type==s21 || s_type==s21r1 ){
          if(s_mode==D2KKPi){if(bdd->getVal()< 0.5){bMassRegion.setLabel("bckgrd");}//+0.6
          } else if(s_mode==D2PiPiPi){         if(bdd->getVal()< 0.6){bMassRegion.setLabel("bckgrd");}  //+.30
          } else if(s_mode==D2PiKPi) {         if(bdd->getVal()< 0.6){bMassRegion.setLabel("bckgrd");}}  //+.54
        } else if(s_type==s24 || s_type==s26){
          if(s_mode==D2KKPi){if(bdd->getVal()< 0.4){bMassRegion.setLabel("bckgrd");}//+0.6
          } else if(s_mode==D2PiPiPi){        if(bdd->getVal()< 0.6){bMassRegion.setLabel("bckgrd");}  //+.30
          } else if(s_mode==D2PiKPi) {        if(bdd->getVal()< 0.6){bMassRegion.setLabel("bckgrd");}}
        }
      
      }
    }
    
    // Cuts to remove Charmless background 
    if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* D_fdchi2    = (RooRealVar*) row->find(D_FDCHI2.GetName());
      RooRealVar* D0_fdchi2   = (RooRealVar*) row->find(D0_FDCHI2.GetName());
      
      if((s_mode == D2KKPi) && s_Bmode == DD0){
        if(D_fdchi2->getVal() <0.0)  bMassRegion.setLabel("bckgrd");
        if(D0_fdchi2->getVal()<0.0)  bMassRegion.setLabel("bckgrd");
      
      } else if(s_mode == D2PiPiPi && s_Bmode == DD0){
        if(D_fdchi2->getVal()< 0.0)  bMassRegion.setLabel("bckgrd");
        if(D0_fdchi2->getVal()<0.0)   bMassRegion.setLabel("bckgrd");
      
      } else if(s_mode == D2PiKPi && s_Bmode == DD0){ 
        if(D_fdchi2->getVal() <0.0)  bMassRegion.setLabel("bckgrd");
        if(D0_fdchi2->getVal()<0.0)   bMassRegion.setLabel("bckgrd");
      
      } else if((s_mode == D2KKPi) && (s_Bmode == DKst0||s_Bmode == DKst0Side)){
        if(D_fdchi2->getVal()<0.0)   bMassRegion.setLabel("bckgrd");
      
      } else if(s_mode == D2PiPiPi && (s_Bmode == DKst0||s_Bmode == DKst0Side)){
        if(D_fdchi2->getVal()<0.0)   bMassRegion.setLabel("bckgrd");
      
      } else if(s_mode == D2PiKPi && (s_Bmode == DKst0||s_Bmode == DKst0Side)){
        if(D_fdchi2->getVal()<0.0)   bMassRegion.setLabel("bckgrd");
      }

    }

    if((s_type==s21||s_type==s21r1) && signal==std::string(bMassRegion.getLabel())){
      RooRealVar* _L0Global_TIS           = (RooRealVar*) row->find(L0_TIS.GetName());
      RooRealVar* _L0Hadron_TOS           = (RooRealVar*) row->find(L0Hadron_TOS.GetName());
      RooRealVar* _Hlt1TrackAllL0_TOS     = (RooRealVar*) row->find(Hlt1TrackAllL0_TOS.GetName());
      RooRealVar* _Hlt2Topo2BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo2BodyBBDT_TOS.GetName());
      RooRealVar* _Hlt2Topo3BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo3BodyBBDT_TOS.GetName());
      RooRealVar* _Hlt2Topo4BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo4BodyBBDT_TOS.GetName());
      RooRealVar* _Hlt2IncPhi_TOS         = (RooRealVar*) row->find(Hlt2IncPhi_TOS.GetName());
      if (i%1000 == 0){
        std::cout << std::endl;
        std::cout << " L0Global_TIS           "  << _L0Global_TIS->getVal()<< std::endl;
        std::cout << " L0Hadron_T0S           "  << _L0Hadron_TOS->getVal()<< std::endl;
        std::cout << " Hlt1TrackAllL0_TOS     "  << _Hlt1TrackAllL0_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2BodyBBDT_TOS  "  << _Hlt2Topo2BodyBBDT_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2BodyBBDT_TOS  "  << _Hlt2Topo2BodyBBDT_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2BodyBBDT_TOS  "  << _Hlt2Topo2BodyBBDT_TOS->getVal()<< std::endl;
        std::cout << " Hlt2IncPhi_TOS         "  << _Hlt2IncPhi_TOS->getVal()<< std::endl;    
        std::cout << std::endl;
        }
      nTot++;
      if(!(_L0Hadron_TOS->getVal() || _L0Global_TIS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotL0++;
      }else if(!(_L0Hadron_TOS->getVal())){
        nNotL0TOS++;
      }
      if(!_Hlt1TrackAllL0_TOS->getVal()){
        bMassRegion.setLabel("bckgrd");
        nNotTrackTOS++;
      }
      if(s_Bmode==DD0){
        if(!(_Hlt2Topo2BodyBBDT_TOS->getVal()||_Hlt2Topo3BodyBBDT_TOS->getVal()||_Hlt2Topo4BodyBBDT_TOS->getVal())){
          bMassRegion.setLabel("bckgrd");
          nNotTopoTOS++;
        }
      } else {
        if(!(_Hlt2Topo2BodyBBDT_TOS->getVal()||_Hlt2Topo3BodyBBDT_TOS->getVal()||_Hlt2Topo4BodyBBDT_TOS->getVal()||_Hlt2IncPhi_TOS->getVal())){
          bMassRegion.setLabel("bckgrd");
        }    
      }
      if(!(_Hlt2Topo2BodyBBDT_TOS->getVal()||_Hlt2Topo3BodyBBDT_TOS->getVal()||_Hlt2Topo4BodyBBDT_TOS->getVal())) nNotTopoTOS++;
      if(_Hlt2IncPhi_TOS->getVal()) nTotIncPhi++;
      if(!(_Hlt2Topo2BodyBBDT_TOS->getVal()||_Hlt2Topo3BodyBBDT_TOS->getVal()||_Hlt2Topo4BodyBBDT_TOS->getVal())&&_Hlt2IncPhi_TOS->getVal()) nTotOnlyIncPhi++;

    
    }

    if((s_type==s24|| s_type==s26) && signal==std::string(bMassRegion.getLabel())){

      RooRealVar* _L0Global_TIS        = (RooRealVar*) row->find(L0_TIS.GetName());
      RooRealVar* _L0Hadron_TOS        = (RooRealVar*) row->find(L0Hadron_TOS.GetName());
      RooRealVar* _Hlt1TrackMVA_TOS    = (RooRealVar*) row->find(Hlt1TrackMVA_TOS.GetName());
      RooRealVar* _Hlt1TwoTrackMVA_TOS = (RooRealVar*) row->find(Hlt1TwoTrackMVA_TOS.GetName());
      RooRealVar* _Hlt2Topo2Body_TOS   = (RooRealVar*) row->find(Hlt2Topo2Body_TOS.GetName());
      RooRealVar* _Hlt2Topo3Body_TOS   = (RooRealVar*) row->find(Hlt2Topo3Body_TOS.GetName());
      RooRealVar* _Hlt2Topo4Body_TOS   = (RooRealVar*) row->find(Hlt2Topo4Body_TOS.GetName());  
      RooRealVar* _Hlt2IncPhi_TOS;
      if (s_type==s24){ 
        _Hlt2IncPhi_TOS = (RooRealVar*) row->find(Hlt2IncPhi_TOS.GetName());
      }else { 
        _Hlt2IncPhi_TOS = (RooRealVar*) row->find(Hlt2PhiIncPhi_TOS.GetName());
      }
      if (i%1000 == 0){
        std::cout << std::endl;
        std::cout << " L0Global_TIS           "  << _L0Global_TIS->getVal()<< std::endl;
        std::cout << " L0Hadron_T0S           "  << _L0Hadron_TOS->getVal()<< std::endl;
        std::cout << " Hlt1TrackMVA_TOS       "  << _Hlt1TrackMVA_TOS->getVal()<< std::endl;
        std::cout << " Hlt1TwoTrackMVA_TOS    "  << _Hlt1TwoTrackMVA_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2Body_TOS      "  << _Hlt2Topo2Body_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2Body_TOS      "  << _Hlt2Topo2Body_TOS->getVal()<< std::endl;
        std::cout << " Hlt2Topo2Body_TOS      "  << _Hlt2Topo2Body_TOS->getVal()<< std::endl;
        std::cout << " Hlt2IncPhi_TOS         "  << _Hlt2IncPhi_TOS->getVal()<< std::endl;    
        std::cout << std::endl;
        }
      nTot++;

      if(!(_L0Hadron_TOS->getVal() || _L0Global_TIS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotL0++;
      }else if(!(_L0Hadron_TOS->getVal())){
        nNotL0TOS++;
      }


      if(!(_Hlt1TrackMVA_TOS->getVal()||_Hlt1TwoTrackMVA_TOS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotTrackTOS++;
      }

      if(s_Bmode==DD0){
        if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal())){
          bMassRegion.setLabel("bckgrd");
          nNotTopoTOS++;
        }
      } else {
        if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal()||_Hlt2IncPhi_TOS->getVal())){
          bMassRegion.setLabel("bckgrd");
          
        } 

      if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal())) nNotTopoTOS++;
      if(_Hlt2IncPhi_TOS->getVal()) nTotIncPhi++;
      if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal())&&_Hlt2IncPhi_TOS->getVal()) nTotOnlyIncPhi++;

      }

    }

    extra->add(RooArgSet(type,helBin,Bmode,mode,magnet,charge,bMassRegion));
  } //end of loop over entries
  

  std::cout <<"   Trigger: notL0TISTOS: "<<100*nNotL0/float(nTot)<<"% (TIS/L0Acc="<<100*nNotL0TOS/float(nTot-nNotL0)<<"%. Not TrackTOS: "<<100*nNotTrackTOS/float(nTot)<<"%. Not TopoTOS "<<100*nNotTopoTOS/float(nTot)<<"%. IncPhi "<<100*nTotIncPhi/float(nTot)<<"%. OnlyIncPhi "<<100*nTotOnlyIncPhi/float(nTot)<<"%."<<std::endl;
  input->merge(extra);
  delete extra;
  return input;
} //end of funcn FinalDataSet




void DKst0Fitting::AssignGlobalCategory()
{
  if(par->debug) std::cout<<"Running: AssignGlobalCategory()"<<std::endl;

  if(par->polarity==up){RooDataSet* reducedData = (RooDataSet*)data->reduce("magnet==magnet::up"); delete data; data=reducedData; }
  if(par->polarity==dn){RooDataSet* reducedData = (RooDataSet*)data->reduce("magnet==magnet::dn"); delete data; data=reducedData; }

  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(*model->Cat()));
  for(int i=0;i<data->numEntries();i++){
    const RooArgSet *row    = data->get(i);
    RooCategory* _Bmode     = (RooCategory*) row->find("Bmode");
    RooCategory* _mode      = (RooCategory*) row->find("mode");
    RooCategory* _charge    = (RooCategory*) row->find("charge");
    RooCategory* _magnet    = (RooCategory*) row->find("magnet");
    RooCategory* _helBin    = (RooCategory*) row->find("helBin");
    RooCategory* _type      = (RooCategory*) row->find("type");
    if(not par->splitYears || par->genToys) _type->setLabel(both.c_str());
    if(par->polarity==both) _magnet->setLabel(both.c_str());
    if(par->sumOverCharges) _charge->setLabel(both.c_str());
    if( not par->splitHel)  _helBin->setLabel(both.c_str());

    model->Cat()->setLabel(Form("%s_%s_%s_%s_%s_%s",_type->getLabel(),_helBin->getLabel(),_Bmode->getLabel(),_mode->getLabel(),_charge->getLabel(),_magnet->getLabel()));
    if(par->debug&&0==i%1000) std::cout << "Event " << i <<"==> Magnet: "<<_magnet->getLabel()<<", Year: "<<_type->getLabel()<<", Helicity: "<<_helBin->getLabel()<<", B Mode: "<<_Bmode->getLabel()<<", D Mode: "<<_mode->getLabel()<<", Charge: "<<_charge->getLabel()<<". Cat: "<<model->Cat()->getLabel()<<std::endl;
    extra->add(RooArgSet(*model->Cat()));
  }
  data->merge(extra);

  RooDataSet* reducedData = new RooDataSet("reducedData","",data,reducedlist);
  delete data; data=reducedData;
}



void DKst0Fitting::PrintDataSet()
{
  if(par->debug) std::cout<<"Running: PrintDataSet()"<<std::endl;
  bool verbose = 1 ;

  if(verbose){
    for(int i=0;i<data->numEntries();i++){

    if(par->debug&&0==i%1000) std::cout<<"Running: PrintDataSet() ==> Printing event " << i <<std::endl;
      const RooArgSet *row = data->get(i);
      RooRealVar*  v1 = (RooRealVar*)  row->find(mB.GetName());
      RooCategory* c0 = (RooCategory*) row->find("magnet");
      RooCategory* c1 = (RooCategory*) row->find("charge");
      RooCategory* c3 = (RooCategory*) row->find("type");
      RooCategory* c4 = (RooCategory*) row->find("mode");
      RooCategory* c5 = (RooCategory*) row->find("Bmode");
      RooCategory* c6 = (RooCategory*) row->find("helBin");
      if(!v1) std::cout << "Can't find v1 " << std::endl;
      if(!c0) std::cout << "Can't find c0 " << std::endl;
      if(!c1) std::cout << "Can't find c1 " << std::endl;
      if(!c3) std::cout << "Can't find c3 " << std::endl;
      if(!c4) std::cout << "Can't find c4 " << std::endl;
      if(!c5) std::cout << "Can't find c5 " << std::endl;
      if(!c6) std::cout << "Can't find c6 " << std::endl;
      if(par->debug&&0==i%1000) std::cout << "Event " << i <<": magnet "<<c0->getLabel()<<", "<<c3->getLabel()<<": "<<c1->getLabel()<<" "<<c4->getLabel()<<" "<<c5->getLabel()<< " "<<c6->getLabel()<<"\t mB="<<v1->getVal()<<"   "<<std::endl;
    }
  }

  data->table(*model->Cat())->Print("v");
  if(par->genToys){
    std::cout<<"   NUMBER OF GENERATED TOY EVENTS = "<<data->sumEntries()<<std::endl;
  }else{
    std::cout<<"                            TOTAL = "<<data->sumEntries()<<std::endl;
  }
}



void DKst0Fitting::DefineModel()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::DefineModel()"<<std::endl;
  /**********************************************/
  /**********************************************/
  /****/           state=BLIND;             /****/
  /**********************************************/
  /**********************************************/
  if(par->genToys) state=UNBLIND;
  model_gen = new DKst0Model(par,&mB,UNBLIND,1,1);
  model     = new DKst0Model(par,&mB,state,1,0);

  if(par->vary) model_novary  = new DKst0Model(par,&mB,UNBLIND,0,0);
  std::cout << "DKst0Model built; continuing..." << std::endl << std::endl;
  if(par->runEff){
    //model_phi = new PhiModel(par,&mKst0,&mD0,state);
    //std::cout << "PhiModel built; continuing..." << std::endl << std::endl;
    //model_ds  = new DsModel( par,&mD,state);
    //std::cout << "DsModel built; continuing..." << std::endl << std::endl;
  }
  reducedlist.add(*model->Cat());
}



void DKst0Fitting::RunFullFit(bool draw=true){
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunFullFit()"<<std::endl;
  PrintDataSet();
  double binWidth=10;//MeV/c2
  mB.setBins((mB.getMax()-mB.getMin())/binWidth);
  RooCategory* cat=model->Cat();
  RooCategory *cat_KKPi, *cat_PiKPi, *cat_PiPiPi;
  cat_KKPi  = new RooCategory( *(model->Cat()), "cat_KKPi"  );
  cat_PiKPi = new RooCategory( *(model->Cat()), "cat_PiKPi" );
  cat_PiPiPi= new RooCategory( *(model->Cat()), "cat_PiPiPi");

  RooArgSet tempcat = RooArgSet(*cat);
  std::cout << "=============> Printed data set " <<std::endl;
  tempcat.Print();

  RooSimultaneous*  sim=model->Pdf();
  RooDataHist* hist;
  if(par->binned) hist = data->binnedClone();

  if(par->doFit){
    RooAbsData* abs=(RooAbsData*) data;
    if(par->binned) abs=(RooAbsData*) hist;
    RooFitResult* result = 0;
    
    std::cout <<"Running: DKst0Fitting::RunFullFit() ---> Printing variables...\n"<<std::endl;
    
    std::cout << "============= Constant variables ==============\n" << std::endl ;
    RooRealVar* var=0;
    RooArgSet* vars = sim->getVariables();
    TIterator* it = vars->createIterator();
    int nv = 0; 
    while((var = (RooRealVar*)it->Next())) {
      if (var->isConstant()){
        nv++;
        std::cout << "Variable " << nv << "\t -> Initial Value: " << var->getVal();
        std::cout<< "\t" << var->GetName()<<std::endl;
      }
      Initial_value[var->GetName()] = var->getVal();
    }
    std::cout<<std::endl;
    std::cout << "============= Variable variables ==============\n" << std::endl ;
    var=0;
    vars = sim->getVariables();
    it = vars->createIterator();
    nv = 0; 
    while((var = (RooRealVar*)it->Next())) {
      if (!var->isConstant()){
        nv++;
        std::cout << "Variable " << nv << "\t -> Initial Value: " << var->getVal();
        std::cout<< "\t" << var->GetName()<<std::endl;
      }
      Initial_value[var->GetName()] = var->getVal();
    }
    std::cout<<std::endl;




    if(state==BLIND){
      // Silence RooFit to prevent it printing the fit result :
      RooMsgService::instance().setSilentMode( kTRUE );
      // And suppress INFO messages for plot range yields :
      //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    }
    
    std::cout <<"Running: DKst0Fitting::RunFullFit() ---> Starting fitTo"<<std::endl;
    result = sim->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
    TCanvas* corrCanv = new TCanvas("corrCanv","",1000,1000); corrCanv->cd();
    //result->correlationHist()->Draw("colz");
    TH2* corr = result->correlationHist();
    corr->GetXaxis()->SetLabelSize(0.01);
    corr->Draw("colz");
    corr->Draw("text same");
    //result->correlationHist()->DrawCopy("text same");
    std::cout << "===================================" << std::endl ;
    std::cout <<"Running: DKst0Fitting::RunFullFit() ---> Printing fit quality..."<<std::endl;  
    std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result->minNll()<<" covQual = " << result->covQual() << std::endl;

    if(state==BLIND){
      // If we are blind , we should print the fit status information
      // and the values of the non - blinded params :
      std::cout << std::endl ;
      std::cout << "===================================" << std::endl ;
      std::cout << " Fit complete " << std::endl ;
      std::cout << " covQual : " << result->covQual() << std::endl ;
      std::cout << " EDM : " << result->edm() << std::endl ;
      std::cout << " FCN at min : " << result->minNll() << std::endl ;
      std::cout << std::endl ;
    } else {
      result->Print("v");    
    }

    std::cout <<"Running: DKst0Fitting::RunFullFit() ---> Printing final variables...\n"<<std::endl;
    std::cout << "============= Constant variables ==============\n" << std::endl ;
    var=0;
    vars = sim->getVariables();
    it = vars->createIterator();
    nv = 0; 
    while((var = (RooRealVar*)it->Next())) {
      if (var->isConstant()){
        nv++;
        std::cout << "Variable " << nv << "\t -> Initial Value: " << var->getVal();
        std::cout<< "\t" << var->GetName()<<std::endl;
      }
    }
    std::cout<<std::endl;
    std::cout << "============= Variable variables ==============\n" << std::endl ;
    var=0;
    vars = sim->getVariables();
    it = vars->createIterator();
    nv = 0; 
    while((var = (RooRealVar*)it->Next())) {
      if (!var->isConstant()){
        nv++;
        std::cout << "Variable " << nv << "\t -> Initial Value: " << Initial_value[var->GetName()];
        if(0==std::string(var->GetName()).find("yield_peak_total_DKst0")){
          std::cout << "\t BLIND";
        } else{  
          std::cout << "\t Final Value: " << var->getVal();
        }
        std::cout<< "\t" << var->GetName()<<std::endl;
      }
    }
    std::cout<<std::endl;

    std::cout << "===================================" << std::endl ;
    model->PrintResult();

    gSystem->Exec("mkdir -p results");
    std::string lastfit_name = "lastRooFitResult";
    if (par->sumOverCharges) lastfit_name += "_summed";
    if (par->splitHel)       lastfit_name += "_splitHel";
    if(par->dsetsReq[s21])   lastfit_name += "_"+s21;
    if(par->dsetsReq[s21r1]) lastfit_name += "_"+s21r1;
    if(par->dsetsReq[s24])   lastfit_name += "_"+s24;
    if(par->dsetsReq[s26])   lastfit_name += "_"+s26;
    if(par->modes[D2PiKPi]) lastfit_name += "_"+D2PiKPi;
    if(par->modes[D2KKPi])  lastfit_name += "_"+D2KKPi;
    if(par->modes[D2PiPiPi])lastfit_name += "_"+D2PiPiPi;    if(!par->vary){
      TFile f(Form("results/%s.root",lastfit_name.c_str()),"RECREATE");
      result->Write(); 
      f.Close();
    } else {
      int seed = 20110602;
      if(par->useSeed) seed = par->seed;
      TFile f(Form("%s/vary_%s_%s_SEED_%i_%i.root",par->variationoption.c_str(),lastfit_name.c_str(),par->toVary.c_str(),seed,int(sim->expectedEvents(*cat))),"RECREATE");
      result->Write(); 
      f.Close();  
    }
  } // closes if(par->doFit)

  if(!draw) return;
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunFullFit() ----> Drawing Histograms..."<<std::endl;

	mB.setRange("lowersideband", 4900, 5220);
  mB.setRange("uppersideband", 5340, 5900);
  mB.setRange("wholesideband", 4900, 5900);
	//TString lowersideband_string = Form("%s>5050&&%s<5230",mB.GetName(),mB.GetName());//"mB<5230";
	//TString uppersideband_string = Form("%s>5330&&%s<5900",mB.GetName(),mB.GetName());//"mB>5330";
  //TString wholerange_string    = Form("%s>5050&&%s<5900",mB.GetName(),mB.GetName());//"mB>5330";

  std::vector<std::string> mergedList;
  mergedList.push_back("merged");
  if(par->doMerge) modeList = mergedList;

  // -------------------------------------------------------------------------
  std::vector<std::string> temp_list = yearList;
  if(par->genToys) temp_list = typeList; 
  for(std::vector<std::string>::iterator y=temp_list.begin();y!=temp_list.end();y++){ 
    for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
      for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
        //if(par->debug) std::cout<<"Mode: "<<(*b).c_str()<<", "<<(*m).c_str()<<", Year: "<<(*y).c_str()<<std::endl;
        std::string canvas_name = Form("canvas_%s_%s_%s",(*b).c_str(),(*m).c_str(),(*y).c_str());
        if (par->sumOverCharges) canvas_name += "_summed";
        if (par->splitHel)       canvas_name += "_splitHel";
        if(par->dsetsReq[s21])   canvas_name += "_"+s21;
        if(par->dsetsReq[s21r1]) canvas_name += "_"+s21r1;
        if(par->dsetsReq[s24])   canvas_name += "_"+s24;
        if(par->dsetsReq[s26])   canvas_name += "_"+s26;
        if( (*b).c_str()==DKst0||(*b).c_str()==DKst0Side){      
          binWidth=40;//MeV/c2
          mB.setBins((mB.getMax()-mB.getMin())/binWidth);
        }

        TCanvas* canvas=new TCanvas(canvas_name.c_str(),Form("%s_%s_%s",(*b).c_str(),(*m).c_str(),(*y).c_str()),40,0,(HelBinList.size()*chargeList.size()*600),600);
        TCanvas* canRes=0;
        RooHist* hresid=0;
        RooHist* hresid_2=0;
        //if( par->doFit )
        std::string canvasres_name = Form("residuals_%s_%s_%s",(*b).c_str(),(*m).c_str(),(*y).c_str());
        if (par->sumOverCharges) canvasres_name += "_summed";
        if (par->splitHel)       canvasres_name += "_splitHel";
        if(par->dsetsReq[s21])   canvasres_name += "_"+s21;
        if(par->dsetsReq[s21r1]) canvasres_name += "_"+s21r1;
        if(par->dsetsReq[s24])   canvasres_name += "_"+s24;
        if(par->dsetsReq[s26])   canvasres_name += "_"+s26;
        canRes = new TCanvas(canvasres_name.c_str(),Form("%s_%s",(*b).c_str(),(*m).c_str()),40,0,(HelBinList.size()*chargeList.size()*600),300);
        float maxH = 0.0;
        std::map<std::string,std::map<std::string,std::map<std::string,RooPlot*> > > plot;

        for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              std::string cat_str = Form("%s_%s_%s_%s_%s_%s",(*y).c_str(),(*h).c_str(),(*b).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());
              std::cout << "------------------------------------------" <<std::endl;
              std::cout << " Plotting Category: " << cat_str << std::endl; 
              std::cout << "------------------------------------------" <<std::endl;
              RooArgSet cat_argset;
              if(par->doMerge){
                std::cout << "Label: " << cat->getLabel() << std::endl;
              }else{
                cat->setLabel(cat_str.c_str());
                cat_argset.add(RooArgSet(*cat));
              }
             
              std::string cat_cut = Form("cat==cat::%s",cat_str.c_str());

              std::string cat_str_KKPi, cat_str_PiKPi, cat_str_PiPiPi;
              // PDF names
              std::string str_comb     = Form("pdf_comb_%s",     cat_str.c_str());    
              std::string str_peak     = Form("pdf_peak_%s",     cat_str.c_str());
              std::string str_PartReco = Form("pdf_PartReco_%s", cat_str.c_str());
              std::string str_DstD0    = Form("pdf_DstD0_%s",   cat_str.c_str());  
              std::string str_DDst0    = Form("pdf_DDst0_%s",   cat_str.c_str());  
              std::string str_DKPiPi   = Form("pdf_Bs2DKPiPi_%s", cat_str.c_str()); 
              std::string str_DstKPi   = Form("pdf_Bu2DstKPi_%s", cat_str.c_str()); 
              std::string str_Dst0Kst0 = Form("pdf_B02Dst0Kst0_%s", cat_str.c_str()); 
              std::string str_DKst0pi  = Form("pdf_Bs2DKst0pi_%s", cat_str.c_str());
              std::string str_DKst0pi_2= Form("pdf_Bs2DKst0pi_2_%s", cat_str.c_str());  

              if(par->doMerge) {
                cat_str_KKPi   = Form("%s_%s_%s_%s_%s_%s",(*y).c_str(),(*h).c_str(),(*b).c_str(),D2KKPi.c_str()  ,(*c).c_str(),(*a).c_str());
                cat_str_PiKPi  = Form("%s_%s_%s_%s_%s_%s",(*y).c_str(),(*h).c_str(),(*b).c_str(),D2PiKPi.c_str() ,(*c).c_str(),(*a).c_str());
                cat_str_PiPiPi = Form("%s_%s_%s_%s_%s_%s",(*y).c_str(),(*h).c_str(),(*b).c_str(),D2PiPiPi.c_str(),(*c).c_str(),(*a).c_str());

                cat_cut = Form("(cat==cat::%s||cat==cat::%s||cat==cat::%s)",cat_str_KKPi.c_str(),cat_str_PiKPi.c_str(),cat_str_PiPiPi.c_str());

                
                helBin.setLabel(   Form("%s",(*h).c_str()));
                Bmode.setLabel(    Form("%s",(*b).c_str()));

                std::cout << "Type: " << type.getLabel() << std::endl;
                std::cout << "Hel : " << helBin.getLabel()<< std::endl;
                std::cout << "Bmod: " << Bmode.getLabel()<< std::endl;
                std::cout << "mag : " << magnet.getLabel()<< std::endl;
                std::cout << "char: " << charge.getLabel()<< std::endl;

                RooArgSet test(Bmode,helBin);
                cat_argset.add(test);
                std::cout << "=============>> ARGSET" <<std::endl;
                cat_argset.Print();
                std::cout << "=============" <<std::endl;
                const RooArgSet *row = data->get(1);
                row->Print();
                RooCategory* c0 = (RooCategory*) row->find("magnet");
                RooCategory* c1 = (RooCategory*) row->find("charge");
                RooCategory* c2 = (RooCategory*) row->find("cat");
                RooCategory* c3 = (RooCategory*) row->find("type");
                RooCategory* c4 = (RooCategory*) row->find("mode");
                RooCategory* c5 = (RooCategory*) row->find("Bmode");
                RooCategory* c6 = (RooCategory*) row->find("helBin");
                std::cout << "Type: " << c3->getLabel() << std::endl;
                std::cout << "Hel : " << c6->getLabel()<< std::endl;
                std::cout << "Bmod: " << c5->getLabel()<< std::endl;
                std::cout << "mag : " << c0->getLabel()<< std::endl;
                std::cout << "char: " << c1->getLabel()<< std::endl;
                std::cout << "cat : " << c2->getLabel()<< std::endl;

                std::cout << "=============" <<std::endl;
                
                str_comb     = Form("pdf_comb_%s,pdf_comb_%s,pdf_comb_%s",                     cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );    
                str_peak     = Form("pdf_peak_%s,pdf_peak_%s,pdf_peak_%s",                     cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );    
                str_PartReco = Form("pdf_PartReco_%s,pdf_PartReco_%s,pdf_PartReco_%s",         cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );
                str_DstD0    = Form("pdf_DstD0_%s,pdf_DstD0_%s,pdf_DstD0_%s",                  cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );  
                str_DDst0    = Form("pdf_DDst0_%s,pdf_DDst0_%s,pdf_DDst0_%s",                  cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );  
                str_DKPiPi   = Form("pdf_B2sDKPiPi_%s,pdf_B2sDKPiPi_%s,pdf_B2sDKPiPi_%s",      cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() ); 
                str_Dst0Kst0 = Form("pdf_B02Dst0Kst0_%s,pdf_B02Dst0Kst0_%s,pdf_B02Dst0Kst0_%s",cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );
                str_DKst0pi  = Form("pdf_Bs2DKst0pi_%s,pdf_Bs2DKst0pi_%s,pdf_Bs2DKst0pi_%s",   cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );
                str_DKst0pi_2= Form("pdf_Bs2DKst0pi_2_%s,pdf_Bs2DKst0pi_2_%s,pdf_Bs2DKst0pi_2_%s",cat_str_KKPi.c_str(), cat_str_PiKPi.c_str(), cat_str_PiPiPi.c_str() );
              }

              double n_sideband_upper = data->sumEntries("1","uppersideband");
              double n_sideband_lower = data->sumEntries("1","lowersideband");          
              
              plot[*h][*c][*a] = mB.frame();
              
              if( ((*b).c_str()==DKst0||(*b).c_str()==DKst0Side ) and state==BLIND){     
                std::cout << "Keeping Bu region blind!" << std::endl;

                data->plotOn( plot[*h][*c][*a],
                              RooFit::Cut(cat_cut.c_str()), 
                              RooFit::DrawOption("PZ"),
                              RooFit::CutRange("lowersideband,uppersideband"),
                              RooFit::Name(Form("data_%s", cat_str.c_str() ))  );

                if(par->doDraw){                
                  if((*b).c_str()==DKst0||(*b).c_str()==DKst0Side ){
                    // ------------------------------------------
                    // Combinatoric + D*Kst0 + DKPiPi+ DstKPi + Dst0Kst0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s,%s,%s,%s,%s",str_DKst0pi_2.c_str(),str_DKst0pi.c_str(),str_Dst0Kst0.c_str(),str_DstKPi.c_str(),str_DKPiPi.c_str(),str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kYellow), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("DKst0pi_2") );
                    // ------------------------------------------
                    // Combinatoric + D*Kst0 + DKPiPi+ DstKPi + Dst0Kst0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s,%s,%s,%s,%s",str_DKst0pi.c_str(),str_Dst0Kst0.c_str(),str_DstKPi.c_str(),str_DKPiPi.c_str(),str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kOrange), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("DKst0pi") );

                    // ------------------------------------------
                    // Combinatoric + D*Kst0 + DKPiPi+ DstKPi + Dst0Kst0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s,%s,%s,%s",str_Dst0Kst0.c_str(),str_DstKPi.c_str(),str_DKPiPi.c_str(),str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kMagenta), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("Dst0Kst0") );
                    
                    // ------------------------------------------
                    // Combinatoric + D*Kst0 + DKPiPi+ DstKPi
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s,%s,%s",str_DstKPi.c_str(),str_DKPiPi.c_str(),str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kBlue), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("DstKPi") );
                    
                    // ------------------------------------------
                    // Combinatoric + D*Kst0 + DKPiPi
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s,%s",str_DKPiPi.c_str(),str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kRed), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("DKPiPi") );
                    // ------------------------------------------
                    // Combinatoric + D*Kst0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s,%s",str_PartReco.c_str(),str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ), 
                                RooFit::FillColor(kGreen), //kAzure+5
                                RooFit::Precision(1.E-7),
                                RooFit::DrawOption("B"), 
                                RooFit::LineWidth(0),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper),RooAbsReal::NumEvent),
                                RooFit::Name("DstKst0") );
                    // ------------------------------------------
                    // Combinatoric
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s",str_comb.c_str())),
                                RooFit::Range("lowersideband,uppersideband" ),         
                                RooFit::FillColor(kAzure-7),
                                RooFit::Precision(1.E-6),
                                RooFit::LineWidth(0), 
                                RooFit::DrawOption("B"),
                                RooFit::Normalization((n_sideband_lower+n_sideband_upper), RooAbsReal::NumEvent),
                                RooFit::Name("Combinatoric"));
                   
                  }

                  if(canRes) {
                    auto dataHist  = (RooHist*)  plot[*h][*c][*a]->getHist(Form("data_%s",cat_str.c_str())); 
                    auto curve1    = (RooCurve*) plot[*h][*c][*a]->getObject(1);  // 1 is index in the list of RooPlot items (see printout from plot[*h][*c][*a]->Print("V")  
                    auto curve2    = (RooCurve*) plot[*h][*c][*a]->getObject(2);
                    hresid         = dataHist->makePullHist(*curve1,true);
                    hresid_2       = dataHist->makePullHist(*curve2,true);
                  }
                } //closes if(par->doDraw)

                // Plot data again so it's on top
                data->plotOn( plot[*h][*c][*a], 
                              RooFit::Cut(cat_cut.c_str()), 
                              RooFit::DrawOption("PZ"),
                              RooFit::CutRange("lowersideband,uppersideband"));
              } else { // if Not Blind

                data->plotOn( plot[*h][*c][*a],
                              RooFit::Cut(cat_cut.c_str()), 
                              RooFit::DrawOption("PZ") );
                if(par->doDraw){
                  
                  if((*b).c_str() == DD0) { 

                    // ------------------------------------------
                    // Combinatoric + D*D0 + DD*0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                  RooFit::Slice(cat_argset), 
                                  RooFit::ProjWData(cat_argset,*data),
                                  RooFit::Components(Form("%s,%s,%s",str_DstD0.c_str(),str_DDst0.c_str(),str_comb.c_str())),
                                  RooFit::LineWidth(0),
                                  RooFit::FillStyle(1002),
                                  RooFit::FillColor(kAzure+5), 
                                  RooFit::DrawOption("F"), 
                                  RooFit::Name("DDst0"));
                
                    // ------------------------------------------
                    // Combinatoric + DD*0
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                  RooFit::Slice(cat_argset), 
                                  RooFit::ProjWData(cat_argset,*data) ,
                                  RooFit::Components(Form("%s,%s",str_DstD0.c_str(),str_comb.c_str())),
                                  RooFit::LineWidth(0),
                                  RooFit::FillStyle(1002),                                  
                                  RooFit::FillColor(kCyan-1), 
                                  RooFit::DrawOption("F"),
                                  RooFit::Name("DstD0") );

                  }else if((*b).c_str() == DKst0||(*b).c_str() == DKst0Side) {

                    // ------------------------------------------
                    // Combinatoric + Dsa1
                    // ------------------------------------------
                    sim->plotOn(  plot[*h][*c][*a],
                                  RooFit::Slice(cat_argset), 
                                  RooFit::ProjWData(cat_argset,*data) ,
                                  RooFit::Components(Form("%s,%s",str_PartReco.c_str(),str_comb.c_str())),
                                  RooFit::LineWidth(0),
                                  RooFit::LineStyle(kSolid) ,
                                  RooFit::FillColor(kGreen), 
                                  RooFit::DrawOption("F"),
                                  RooFit::Name("DstKst0") );                                            
                  }
                  
                  // ------------------------------------------
                  // Combinatoric
                  // ------------------------------------------
                  sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s",str_comb.c_str())),
                                RooFit::LineWidth(0),
                                RooFit::LineStyle(kDotted),
                                RooFit::FillColor(kAzure-7), 
                                RooFit::DrawOption("F"),
                                RooFit::Name("Combinatoric") );
                  
                  // ------------------------------------------
                  // Total
                  // ------------------------------------------
                  sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::LineWidth(2),
                                RooFit::LineColor(kBlack));

                  if(canRes) hresid = plot[*h][*c][*a]->pullHist();
                  

                  // ------------------------------------------
                  // Signal
                  // ------------------------------------------
                  sim->plotOn(  plot[*h][*c][*a],
                                RooFit::Slice(cat_argset), 
                                RooFit::ProjWData(cat_argset,*data) ,
                                RooFit::Components(Form("%s",str_peak.c_str())),
                                RooFit::LineWidth(3),
                                RooFit::LineStyle(kSolid) ,
                                RooFit::LineColor(kRed)  ,
                                RooFit::FillColor(kRed),
                                RooFit::Name("Signal") );
                } //closes if(par->doDraw)
                
                data->plotOn( plot[*h][*c][*a],
                              RooFit::Cut(cat_cut.c_str()), 
                              RooFit::DrawOption("PZ") );
      	      } //closes else
                
                


              float xgutter=0.010;
              float leftmargin=0.08;
              if (par->sumOverCharges) leftmargin=0.13;
              float bottommargin=0.09;
              //float nY=pidList.size();
              //float nY= 1.0;
              float nY= 1.0;
              float nX=chargeList.size()*magnetList.size()*HelBinList.size();
              float padwidth =(0.999-  leftmargin)/nX;
              float padheight=(0.999-bottommargin)/nY;
              float ygutter=xgutter*canvas->GetWindowWidth()/canvas->GetWindowHeight();
              int iX=HelBinList.size()*(chargeList.size()*(a-magnetList.begin()) + (c-chargeList.begin())) + (h-HelBinList.begin());
              //int iY=1-(p-pidList.begin()); 
              //int iY = 0; // we redefine iY as 0, otherwise (if iY is 1) then y1 is ~2 which doesn't work.
              int iY = 0;
              float x0=iX*padwidth+(iX?leftmargin:0);
              float x1=(iX+1)*padwidth+leftmargin;
              float y0=iY*padheight+(iY?bottommargin:0);
              float y1=(iY+1)*padheight+bottommargin;
              //std::cout<<"Laurence: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
              canpad[*h][*b][*m][*c][*a] = new TPad(Form("pad_%s",cat_str.c_str()),"",x0,y0,x1,y1);
              canpad[*h][*b][*m][*c][*a]->SetRightMargin( xgutter/(x1-x0));
              canpad[*h][*b][*m][*c][*a]->SetTopMargin(   ygutter/(y1-y0));
              canpad[*h][*b][*m][*c][*a]->SetLeftMargin(  xgutter+(iX?0:(leftmargin/(x1-x0))));
              canpad[*h][*b][*m][*c][*a]->SetBottomMargin(ygutter+(iY?0:(bottommargin/(y1-y0))));
              TPad* resPad = (TPad*)canpad[*h][*b][*m][*c][*a]->Clone();
              canvas->cd(); canpad[*h][*b][*m][*c][*a]->Draw();
              canpad[*h][*b][*m][*c][*a]->cd();
           
              plot[*h][*c][*a]->Draw();
              plot[*h][*c][*a]->SetTitle("");
           
              if( (iY or (nX-iX>1) ) ){ plot[*h][*c][*a]->GetXaxis()->SetTitle("");}else if((*b).c_str()==DKst0Side){plot[*h][*c][*a]->GetXaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*c][*a]->GetXaxis()->SetTitle("K*^{0} Sideband #font[12]{m}(#font[12]{DK*^{0}}) [MeV/#font[12]{c}^{2}] ");}
              if( (iY or (nX-iX>1) ) ){ plot[*h][*c][*a]->GetXaxis()->SetTitle("");}else if((*b).c_str()==DKst0){plot[*h][*c][*a]->GetXaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*c][*a]->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{DK*^{0}}) [MeV/#font[12]{c}^{2}] ");}
              if( (iY or (nX-iX>1) ) ){ plot[*h][*c][*a]->GetXaxis()->SetTitle("");}else if((*b).c_str()==DD0 ){plot[*h][*c][*a]->GetXaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*c][*a]->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{DD^{0}}) [MeV/#font[12]{c}^{2}] ");}
              if( iX or (nY-iY>1)){ plot[*h][*c][*a]->GetYaxis()->SetTitle("");}else{plot[*h][*c][*a]->GetYaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*c][*a]->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));}
              if(iY){ plot[*h][*c][*a]->GetXaxis()->SetLabelColor(0); }else{plot[*h][*c][*a]->GetXaxis()->SetLabelSize(0.045/(y1-y0));}
              if(iX){ plot[*h][*c][*a]->GetYaxis()->SetLabelColor(0); }else{plot[*h][*c][*a]->GetYaxis()->SetLabelSize(0.045/(y1-y0));}
              plot[*h][*c][*a]->GetXaxis()->SetTitleOffset(0.9);
        	    plot[*h][*c][*a]->GetYaxis()->SetTitleOffset(1.0);
              // Draw B decay title
              x0=((iX?0:leftmargin)+padwidth*0.79)/(padwidth+(iX?0:leftmargin));
              x1=((iX?0:leftmargin)+padwidth*0.99)/(padwidth+(iX?0:leftmargin));
              y0=((iY?0:bottommargin)+padheight*0.89)/(padheight+(iY?0:bottommargin));
              y1=((iY?0:bottommargin)+padheight*0.99)/(padheight+(iY?0:bottommargin));
              //std::cout<<"Title: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
              TPaveLabel *pav = new TPaveLabel(x0,y0,x1,y1,title[*y][*h][*b][*m][*c][*a].c_str(),"NDC");
              pav->SetBorderSize(0);
              pav->SetFillStyle(0);
              pav->SetTextFont(12);
              pav->SetTextSize(0.5);
              pav->SetTextAlign(31);
              pav->Draw();
              // Draw LHCb unofficial 
              x0=((iX?0:leftmargin)+padwidth*0.75)/(padwidth+(iX?0:leftmargin));
              x1=((iX?0:leftmargin)+padwidth*0.95)/(padwidth+(iX?0:leftmargin));
              y0=((iY?0:bottommargin)+padheight*0.79)/(padheight+(iY?0:bottommargin));
              y1=((iY?0:bottommargin)+padheight*0.89)/(padheight+(iY?0:bottommargin));
              //std::cout<<"LHCb: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging       
              TPaveLabel* lhcblabel = new TPaveLabel(x0,y0,x1,y1,"LHCb","NDC");
              lhcblabel->SetBorderSize(0);
              lhcblabel->SetFillStyle(0);
              lhcblabel->SetTextSize(0.5);
              lhcblabel->SetTextFont(62); 
              lhcblabel->SetTextAlign(31);  
              if(!par->genToys||state==UNBLIND) lhcblabel->Draw();
              // Draw 
              x0=((iX?0:leftmargin)+padwidth*0.79)/(padwidth+(iX?0:leftmargin));
              x1=((iX?0:leftmargin)+padwidth*0.99)/(padwidth+(iX?0:leftmargin));
              y0=((iY?0:bottommargin)+padheight*0.69)/(padheight+(iY?0:bottommargin));
              y1=((iY?0:bottommargin)+padheight*0.79)/(padheight+(iY?0:bottommargin));
              //std::cout<<"Details: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
              TPaveLabel *pav2 = new TPaveLabel(x0,y0,x1,y1,bin_detail[*y][*h][*b][*m][*c][*a].c_str(),"NDC");
              pav2->SetBorderSize(0);
              pav2->SetFillStyle(0);
              pav2->SetTextFont(12);
              pav2->SetTextSize(0.5);
              pav2->SetTextAlign(31);
              pav2->Draw();

              
              x0=((iX?0:leftmargin)+padwidth*0.49)/(padwidth+(iX?0:leftmargin));
              x1=((iX?0:leftmargin)+padwidth*0.89)/(padwidth+(iX?0:leftmargin));
              y0=((iY?0:bottommargin)+padheight*0.20)/(padheight+(iY?0:bottommargin));
              y1=((iY?0:bottommargin)+padheight*0.50)/(padheight+(iY?0:bottommargin));
              gStyle->SetLegendBorderSize(0);
              TLegend* leg = new TLegend(x0,y0,x1,y1);
              if((*b).c_str() == DD0){
                leg->AddEntry(plot[*h][*c][*a]->findObject("Signal"),    "B #rightarrow D D^{0}","l"); 
                leg->AddEntry(plot[*h][*c][*a]->findObject("DstD0"),     "B #rightarrow D* D^{0}","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("DDst0"),     "B #rightarrow D D*^{0}","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("Combinatoric"),"Combinatorics","f");
                leg->SetTextFont(12);
                if( !(iY or (nX-iX>1) ) )leg->Draw();
              } else {
                if (state == UNBLIND) leg->AddEntry(plot[*h][*c][*a]->findObject("Signal"),      "D K*^{0}","l"); 
                leg->AddEntry(plot[*h][*c][*a]->findObject("DstKPi"),      "B^{+} #rightarrow D* K #pi","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("DKPiPi"),      "B_{s}^{0} #rightarrow D^{-} K^{-} #pi^{+} #pi^{+}","f"); 
                leg->AddEntry(plot[*h][*c][*a]->findObject("DstKst0"),     "B^{+} #rightarrow D* K*^{0}","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("Dst0Kst0"),    "B^{0} #rightarrow D*^{0} K*^{0}","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("DKst0pi"),     "B_{s}^{0} #rightarrow D^{+} K*^{0} [#pi^{-}]","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("DKst0pi_2"),   "B_{s}^{0} #rightarrow D^{+} K*^{0} #pi^{-}","f");
                leg->AddEntry(plot[*h][*c][*a]->findObject("Combinatoric"),"Combinatorics","f");
                leg->SetTextFont(12);
                if( !(iY or (nX-iX>1) ) )leg->Draw();
              }

              if(plot[*h][*c][*a]->GetMaximum()>maxH) maxH=plot[*h][*c][*a]->GetMaximum();
              if(hresid){
                canRes->cd();
                resPad->Draw();
                resPad->cd();
                RooPlot* frame = mB.frame(RooFit::Title("Residual Distribution"));
                hresid->SetFillColor(4);
                hresid->SetLineColor(10);
                if ((*b)==DKst0Side)frame->GetXaxis()->SetTitle("K*^{0} Sideband #font[12]{m}(#font[12]{D K*^{0}}) [MeV/#font[12]{c}^{2}] ");
                if ((*b)==DKst0)    frame->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D K*^{0}}) [MeV/#font[12]{c}^{2}] ");
                if ((*b)==DD0)      frame->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D D^{0}}) [MeV/#font[12]{c}^{2}] ");
                //frame->GetYaxis()->SetTitle("Residual");
                frame->GetXaxis()->SetTitleOffset(0.9);
                frame->GetXaxis()->SetTitleSize(0.045/(y1-y0));
                frame->GetYaxis()->SetTitleOffset(0.9);
                frame->GetYaxis()->SetTitleSize(0.045/(y1-y0));

                frame->addPlotable(hresid,"B");
                if(hresid_2) hresid_2->SetFillColor(4);
                if(hresid_2) hresid_2->SetLineColor(10);
                if(hresid_2) frame->addPlotable(hresid_2,"B");
                frame->Draw();
                frame->SetTitle("");
                frame->GetYaxis()->SetRangeUser(-5,5); // range chosen such that 1/1000 chance of Gaussian error fluctuating to the edge
              } //end of loop over magnetList    
            } //end of loop over chargeList
          }
        }



        for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                if(0==plot[*h][*c][*a]) continue;
                plot[*h][*c][*a]->SetTitleFont(132,"X"); plot[*h][*c][*a]->SetLabelFont(132,"X");
                plot[*h][*c][*a]->SetTitleFont(132,"Y"); plot[*h][*c][*a]->SetLabelFont(132,"Y");
                plot[*h][*c][*a]->SetMaximum(maxH);
                plot[*h][*c][*a]->SetMinimum(0.01);
                
                if(par->genToys&&((*b)==DKst0||(*b)==DKst0Side )) plot[*h][*c][*a]->SetMaximum(100);

                if( ((*b)==DKst0 ||(*b)==DKst0Side)and state==BLIND ){  
                  canpad[*h][*b][*m][*c][*a]->cd();
                  maxH = 600.0;
                  plot[*h][*c][*a]->SetMaximum(maxH);
                  TPaveLabel *blindpav = new TPaveLabel(5220,0.000001,5340,maxH-0.5,"BLIND","");
                  blindpav->SetBorderSize(0); blindpav->SetTextSize(0.1); blindpav->SetTextAngle(30); 
                  blindpav->SetTextColor(2); 
                  blindpav->SetFillColor(10);
                  blindpav->SetFillStyle(0);
                  blindpav->Draw();
                }
                //std::cout<<(*m)<<" "<<(*c)<<" "<<(*a)<<"  chi2 = "<<plot[*h][*ds][*ph[*h]c][*a]->chiSquare()<<std::endl;
            }
          }
        }

        
        canvas->Print(Form("results/%s.eps",canvas->GetName()));
        canvas->Print(Form("results/%s.pdf",canvas->GetName()));
        //canvas->Print(Form("results/%s.png",canvas->GetName()));
        if(canRes){
         canRes->Print(Form("results/%s.eps",canRes->GetName()));
         canRes->Print(Form("results/%s.pdf",canRes->GetName()));
         //canRes->Print(Form("results/%s.png",canRes->GetName()));
        }

      } //end of loop over modeList 
    }//end of loop over BmodeList 
  }//end of loop over typeList 
}






void DKst0Fitting::RunManyFits()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunManyFits()"<<std::endl;
  
  std::string mva(par->MVAMethod);
  std::string mvatype(par->MVAType); 


  if(par->debug) std::cout<<"Running: DKst0Fitting::RunManyFits() --> Chosen setup: " << mvatype << ", " << mva<<std::endl;

  if(HelBinList.size()*BmodeList.size()*chargeList.size()*magnetList.size()*modeList.size()>1) {
    std::cout << "ERROR Can only run many fits over a single mode at a time" << std::endl;
    return;
  }

  std::string dsmode;
  if(par->modes[D2PiKPi]) dsmode = D2PiKPi;
  if(par->modes[D2KKPi])  dsmode = D2KKPi;
  if(par->modes[D2PiPiPi])dsmode = D2PiPiPi;

  std::string years="";
  if(par->dsetsReq[s26] )   years="_2016";
  if(par->dsetsReq[s24] )   years="_2015";
  if(par->dsetsReq[s21] )   years="_2012";
  if(par->dsetsReq[s21r1] ) years="_2011";
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1]) years = "_Run1"; 
  if(par->dsetsReq[s24] && par->dsetsReq[s26]) years = "_Run2";
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1] && par->dsetsReq[s24] && par->dsetsReq[s26]) years = "_All"; 

  if( mva!=BDT && mva!=BDTG && mva!= BDTB) {
    std::cout << " MVA method must be one of BDT, BDTG, or BDTB" << std::endl;
    return;
  }

  // Set up fit range
  int n_p = 5;
  if (mvatype==MC) n_p = 2;
  if (par->nBDTPoints>0) n_p = par->nBDTPoints;
  
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunManyFits() --> Number of points: " << n_p<<std::endl;

  std::map<std::string,std::map<std::string,double>> low;
  std::map<std::string,std::map<std::string,double>> high;
  

  std::map<std::string,std::map<std::string,std::string>> var_name;
  var_name[BDT][D]    = bdtD.GetName();
  var_name[BDTG][D]   = bdtgD.GetName();
  var_name[BDTB][D]   = bdtbD.GetName();
  var_name[BDT][Kst0] = bdtD0.GetName();
  var_name[BDTG][Kst0]= bdtgD0.GetName();
  var_name[BDTB][Kst0]= bdtbD0.GetName();
  var_name[BDT][MC]   = bdtMC.GetName();
  var_name[BDTG][MC]  = bdtgMC.GetName();
  var_name[BDTB][MC]  = bdtbMC.GetName();

  low[BDT][Kst0]  = -0.4;  high[BDT][Kst0] =  1.1;  
  low[BDT][D]     = -0.4;  high[BDT][D]    =  1.1;
  
  low[BDTG][Kst0] = -0.6;  high[BDTG][Kst0]=  1.1;
  low[BDTG][D]    = -1.0;  high[BDTG][D]   =  1.1;
  
  low[BDTB][Kst0] =  0.5;  high[BDTB][Kst0]=  1.05;
  low[BDTB][D]    = -0.5;  high[BDTB][D]   =  1.05;
  
  low[BDT][MC]   = -1.0;  high[BDT][MC]  =  0.45;
  low[BDTG][MC]  = -1.0;  high[BDTG][MC] =  1.05;
  low[BDTB][MC]  = -1.0;  high[BDTB][MC] =  1.05;
  double cut_value_MC[n_p];
  double cut_value_Kst0[n_p];
  double cut_value_D[n_p]; 

  double N_signal[n_p][n_p], N_background[n_p][n_p], N_background_red[n_p][n_p];
  //double N_sigerr[n_p][n_p];

  double N_signalMC[n_p], N_backgroundMC[n_p], N_background_redMC[n_p];

  double binWidth=10;//MeV/c2
  mB.setBins((mB.getMax()-mB.getMin())/binWidth);
  RooCategory* cat=model->Cat();
  RooSimultaneous*  sim=model->Pdf();

  RooFitResult* result = 0;
  if(mvatype==DATA){
    for (int i = 0; i < n_p ; i++){
      for (int j = 0; j < n_p ; j++){
        cut_value_D[i]    = low[mva][D]    + i*(high[mva][D]   -low[mva][D]   )/n_p;
        cut_value_Kst0[j] = low[mva][Kst0] + j*(high[mva][Kst0]-low[mva][Kst0])/n_p;
        
        
        std::string cut_BDT  = Form("%s>%f&&%s>%f", var_name[mva][D].c_str(), cut_value_D[i], var_name[mva][Kst0].c_str(), cut_value_Kst0[j]);
        
        if(par->debug) std::cout << "--- Running fit with: " << cut_BDT << std::endl;
        
        //cutdata = data;
        //RooAbsData* abs=(RooAbsData*) data;

        if(par->debug) std::cout << "--- Reducing Data: " << std::endl;
        //RooDataSet* BDTdataset=(RooDataSet*)abs->reduce(cut_BDT.c_str());
        //RooAbsData* BDTdataset=(RooAbsData*)abs->reduce(cut_BDT.c_str());
        //delete data; data=BDTdataset;
        //RooDataSet* data_plot = (RooDataSet*) BDTdataset;
        //RooDataHist* hist = BDTdataset->binnedClone();
        
        if(par->debug) std::cout << "--- Reducing RooDataSet: " << std::endl;
        if(par->debug) std::cout << "--- No. events before: " << data->numEntries() << std::endl;
        cutdata = (RooDataSet*) data->reduce(cut_BDT.c_str());     
        if(par->debug) std::cout << "--- No. events after: " << cutdata->numEntries() << std::endl;

        RooAbsData* abs=(RooAbsData*) cutdata;
        RooDataHist* hist=0;
        if(par->binned) hist = cutdata->binnedClone();
        if(par->binned) abs=(RooAbsData*) hist;

        if(par->debug) std::cout << "--- Doing Fit: " << std::endl;
        result = sim->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
        
        // Simple plots
        std::string name_string=  Form("%s_%s_Ds_%f_Phi_%f",dsmode.c_str(), mva.c_str(),cut_value_D[i], cut_value_Kst0[j] );
        TCanvas* canvas = new TCanvas(name_string.c_str(),"",40,0,600,600);
        RooPlot* simple_plot = mB.frame();
        std::string cat_str = Form("%s_%s_%s_%s_%s_%s",
                                   (*yearList.begin()).c_str(),
                                   (*HelBinList.begin()).c_str(),
                                   (*BmodeList.begin()).c_str(),
                                   (*modeList.begin()).c_str(),
                                   (*chargeList.begin()).c_str(),
                                   (*magnetList.begin()).c_str() );
        std::string decay_title = title[*yearList.begin()][*HelBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()];    
        cat->setLabel(cat_str.c_str());
        cutdata->plotOn( simple_plot,
                      RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                      RooFit::DrawOption("PZ") );
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::Components(Form("pdf_DstD0_%s,pdf_DDst0_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                      RooFit::FillColor(kAzure+5), 
                      RooFit::DrawOption("F") );
                      // Comb + Hills
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::Components(Form("pdf_DDst0_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
                      RooFit::FillColor(kCyan-1), 
                      RooFit::DrawOption("F") );
        // Comb
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())),
                      RooFit::FillColor(kAzure-7), 
                      RooFit::DrawOption("F") ); 
        // Total Lower
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::LineWidth(2),
                      RooFit::LineColor(kBlack) ); 
        cutdata->plotOn( simple_plot,
                      RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                      RooFit::DrawOption("PZ") );
        simple_plot->Draw();
        simple_plot->Draw();
        simple_plot->SetTitle("");

        simple_plot->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D D^{0}}) [MeV/#font[12]{c}^{2}] ");
        simple_plot->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));
        simple_plot->GetXaxis()->SetTitleOffset(0.9);
        simple_plot->GetYaxis()->SetTitleOffset(1.0);
        simple_plot->GetXaxis()->SetTitleSize(0.045);
        simple_plot->GetYaxis()->SetTitleSize(0.045);      
        simple_plot->GetXaxis()->SetLabelSize(0.045);
        simple_plot->GetYaxis()->SetLabelSize(0.045);

        fit_results = model->GetResult();
        std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result->minNll()<<" covQual = " << result->covQual() << std::endl;

        if(state==BLIND){

        } else {
          result->Print("v");
          model->PrintResult(); 
        }

        TFile f(Form("results/manyfits/%s/%s/%s%s_result.root",dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str()),"RECREATE");
        result->Write(); 
        f.Close();

        N_signal[i][j]     = fit_results[both][both][DD0][dsmode][both][both]["Peak"];
        //N_sigerr[i][j]     = fit_results[both][DD0][dsmode][both][both]["Peakerr"];
        N_background[i][j] = fit_results[both][both][DD0][dsmode][both][both]["Comb"];
        N_background_red[i][j] = fit_results[both][both][DD0][dsmode][both][both]["Comb_Reduced"];
        std::cout << "            N Signal: " << N_signal[i][j] << std::endl;
        std::cout << "        N Background: " << N_background[i][j] << std::endl;
        std::cout << "N Background reduced: " << N_background_red[i][j] << std::endl;
        std::string text_filename  = Form("results/manyfits/text/%s/%s/DATA/FitResult_%s%s.txt",dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str());
        std::ofstream output_file(text_filename.c_str());
        output_file << cut_value_D[i] <<":"<< cut_value_Kst0[j] <<":"<< N_signal[i][j] << ":" << N_background_red[i][j];
        output_file.close();

        TPaveLabel *pav = new TPaveLabel(0.85,0.85,0.95,0.95,decay_title.c_str(),"NDC");
        pav->SetBorderSize(0);
        pav->SetFillStyle(0);
        pav->SetTextFont(12);
        pav->SetTextSize(0.5);
        pav->SetTextAlign(31);
        pav->Draw();
        TPaveLabel *pav1 = new TPaveLabel(0.85,0.75,0.95,0.85,Form("Signal: %f",N_signal[i][j]),"NDC");
        pav1->SetBorderSize(0);
        pav1->SetFillStyle(0);
        pav1->SetTextFont(62);
        pav1->SetTextSize(0.3);
        pav1->SetTextAlign(31);
        pav1->Draw();
        TPaveLabel *pav2 = new TPaveLabel(0.85,0.65,0.95,0.75,Form("Comb BG: %f",N_background_red[i][j]),"NDC");
        pav2->SetBorderSize(0);
        pav2->SetFillStyle(0);
        pav2->SetTextFont(62);
        pav2->SetTextSize(0.3);
        pav2->SetTextAlign(31);
        pav2->Draw();
        canvas->Print(Form("results/manyfits/%s/%s/%s%s.eps",dsmode.c_str(), mva.c_str(),canvas->GetName(),years.c_str()));
        canvas->Print(Form("results/manyfits/%s/%s/%s%s.pdf",dsmode.c_str(), mva.c_str(),canvas->GetName(),years.c_str()));
        //canvas->Print(Form("results/manyfits/%s/%s/%s.png",dsmode.c_str(), mva.c_str(),canvas->GetName()));
        

        //delete BDTdataset;
        delete cutdata;
      }   
    }
    // Calculate Punzi and significance
    double Punzi_2D[n_p*n_p], Significance_2D[n_p*n_p];
    double cut_value_1D_Kst0[n_p*n_p], cut_value_1D_D[n_p*n_p]; 

    for (int i = 0; i < n_p ; i++){
      for (int j = 0; j < n_p ; j++){
        Punzi_2D[i*n_p +j]          = N_signal[i][j]/ (2.5 + sqrt(N_background_red[i][j]));
        Significance_2D[i*n_p+j]    = N_signal[i][j]/sqrt(N_signal[i][j]+N_background_red[i][j]);
        cut_value_1D_D[i*n_p+j]     =  cut_value_D[i];
        cut_value_1D_Kst0[i*n_p +j] =  cut_value_Kst0[j];
      }
    }

    std::cout << std::endl;
    for (int i=0; i<n_p; i++ ) {
      for (int j=0;j<n_p; j++){
        std::cout << "Cut Value Ds: " << cut_value_1D_D[i*n_p+j] <<"\t Cut Value Phi: " << cut_value_1D_Kst0[i*n_p+j] << " \t N Sig: " << N_signal[i][j] <<" \t N Bg: " << N_background_red[i][j] << " \t Punzi: " << Punzi_2D[i*n_p +j] <<" \t Significance: " << Significance_2D[i*n_p+j] << std::endl;
      }
    }

    std::cout << std::endl;

    std::cout << "double cut_value_1D_D []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << cut_value_1D_D[i];
      if(i!=(n_p*n_p-1))std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;
    
    std::cout << "double cut_value_1D_Kst0 []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << cut_value_1D_Kst0[i];
      if(i!=(n_p*n_p-1))std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;

    std::cout << "double Punzi_2D []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << Punzi_2D[i];
      if(i != (n_p*n_p-1)) std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;

    std::cout << "double Significance_2D []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << Significance_2D[i];
      if(i!=(n_p*n_p-1)) std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;
    std::cout << std::endl;
    // Plot Graphs of mva respose as function of cut 
    //gStyle->Reset(); 
    /*
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.03);
    //gStyle->SetNdivisions(205,"XYZ");
    //gStyle->SetStatFont(132);
    //gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(132,"XYZ");
    gStyle->SetLabelFont(132,"XYZ");
    gStyle->SetTitleSize(0.077,"XYZ");
    gStyle->SetLabelSize(0.08,"XYZ");
    gStyle->SetTitleOffset(0.83,"X");
    gStyle->SetTitleOffset(1.05,"Y");
    gStyle->SetMarkerSize(0.5);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    */

    TCanvas *cav_punzi = new TCanvas(Form("%s_%s_DATA_punzi",dsmode.c_str(), mva.c_str() ),Form("%s %s Punzi as function of cut value",dsmode.c_str(), mva.c_str()),200,10,700,500);
    TGraph2D *gr_punzi = new TGraph2D(n_p*n_p, cut_value_1D_D, cut_value_1D_Kst0, Punzi_2D);
    SetLHCbStyle(surf);
    gr_punzi->Draw("surf1z");
    gr_punzi->SetTitle(Form("Punzi    S/(5/2 + sqrt(B));D %s Cut Value;D0 %s Cut Value;Punzi",mva.c_str(),mva.c_str()));    
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_surface.eps",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_surface.pdf",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    //cav_punzi->Print(Form("results/manyfits/%s/%s/%s_surface.png",dsmode.c_str(), mva.c_str(),cav_punzi->GetName()));
    TFile f1(Form("results/manyfits/%s/%s/%s%s_surface_plots.root",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()),"RECREATE");
    cav_punzi->Write();
    gr_punzi->Write(); 
    f1.Close(); 

    SetLHCbStyle(cont);
    //gr_punzi->Draw("contz");
    gr_punzi->Draw("colz");  
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_contour.pdf",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_contour.eps",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    //cav_punzi->Print(Form("results/manyfits/%s/%s/%s_contour.png",dsmode.c_str(), mva.c_str(),cav_punzi->GetName()));
    TFile f2(Form("results/manyfits/%s/%s/%s%s_contour_plots.root",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()),"RECREATE");
    cav_punzi->Write();
    gr_punzi->Write(); 
    f2.Close(); 


    TCanvas *cav_sig = new TCanvas(Form("%s_%s_DATA_significance",dsmode.c_str(), mva.c_str() ),Form("%s %s Significance as function of cut value",dsmode.c_str(), mva.c_str()),200,10,700,500);
    TGraph2D *gr_sig = new TGraph2D(n_p*n_p, cut_value_1D_D, cut_value_1D_Kst0, Significance_2D);
    SetLHCbStyle(surf);
    gr_sig->Draw("surf1z");
    gr_sig->SetTitle(Form("Significance    S/sqrt(S+B);D %s Cut Value;D0 %s Cut Value;Significance",mva.c_str(),mva.c_str()));    
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_surface.eps",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()));
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_surface.pdf",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()));
    //cav_sig->Print(Form("results/manyfits/%s/%s/%s_surface.png",dsmode.c_str(), mva.c_str(),cav_sig->GetName()));
    TFile f3(Form("results/manyfits/%s/%s/%s%s_surface_plots.root",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()),"RECREATE");
    cav_sig->Write();
    gr_sig->Write(); 
    f3.Close(); 

    SetLHCbStyle(cont);
    gr_sig->Draw("colz");
    //gr_sig->Draw("contz");  
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_contour.eps",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()));
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_contour.pdf",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str())); 
    //cav_sig->Print(Form("results/manyfits/%s/%s/%s_contour.png",dsmode.c_str(), mva.c_str(),cav_sig->GetName()));
    TFile f4(Form("results/manyfits/%s/%s/%s%s_contour_plots.root",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()),"RECREATE");
    cav_sig->Write();
    gr_sig->Write(); 
    f4.Close(); 
  } 
  /*
  else if (mvatype==MC){

    for (int i = 0; i < n_p ; i++){
      

      cut_value_MC[i]  = low[mva][MC]  + i*(high[mva][MC] -low[mva][MC] )/n_p;
      
      std::string cut_BDT  = Form("%s>%f", var_name[mva][MC].c_str(), cut_value_MC[i]);
      if(dsmode==D2KKPi)   cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_K0_PIDK >1  && D_K1_PIDK >1  && D_P_PIDK  <-1 &&" + cut_BDT;
      if(dsmode==D2PiKPi)  cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_K_PIDK  >1  && D_P0_PIDK <-1 && D_P1_PIDK <-1 &&" + cut_BDT;
      if(dsmode==D2PiPiPi) cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_P0_PIDK <-1 && D_P1_PIDK >1  && D_P2_PIDK <-1 &&" + cut_BDT;

      if(par->debug) std::cout << "--- Reducing RooDataSet: "<< cut_BDT << std::endl;
      if(par->debug) std::cout << "--- No. events before: " << data->numEntries() << std::endl;
      cutdata = (RooDataSet*) data->reduce(cut_BDT.c_str());     
      if(par->debug) std::cout << "--- No. events after: " << cutdata->numEntries() << std::endl;

      RooAbsData* abs=(RooAbsData*) cutdata;
      RooDataHist* hist=0;
      if(par->binned) hist = cutdata->binnedClone();
      if(par->binned) abs=(RooAbsData*) hist;

      if(par->debug) std::cout << "--- Doing Fit: " << std::endl;
      result = sim->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
      
      // Simple plots
      std::string name_string=  Form("%s_%s_MC_%f",dsmode.c_str(), mva.c_str(),cut_value_MC[i] );
      TCanvas* canvas = new TCanvas(name_string.c_str(),"",40,0,600,600);
      RooPlot* simple_plot = mB.frame();
      std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s_%s",
                                 (*yearList.begin()).c_str(),
                                 (*HelBinList.begin()).c_str(),
                                 (*DsBDTBinList.begin()).c_str(),
                                 (*PhiBDTBinList.begin()).c_str(),
                                 (*BmodeList.begin()).c_str(),
                                 (*modeList.begin()).c_str(),
                                 (*chargeList.begin()).c_str(),
                                 (*magnetList.begin()).c_str() );
      std::string decay_title = title[*yearList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()];    
      cat->setLabel(cat_str.c_str());
      cutdata->plotOn( simple_plot,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::Components(Form("pdf_HORNS_%s,pdf_HILL_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                    RooFit::FillColor(kAzure+5), 
                    RooFit::DrawOption("F") );
                    // Comb + Hills
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::Components(Form("pdf_HILL_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
                    RooFit::FillColor(kCyan-1), 
                    RooFit::DrawOption("F") );
      // Comb
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())),
                    RooFit::FillColor(kAzure-7), 
                    RooFit::DrawOption("F") ); 
      // Total Lower
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::LineWidth(2),
                    RooFit::LineColor(kBlack) ); 
      cutdata->plotOn( simple_plot,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      simple_plot->Draw();
      simple_plot->Draw();
      simple_plot->SetTitle("");

      simple_plot->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D_{s}D^{0}}) [MeV/#font[12]{c}^{2}] ");
      simple_plot->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));
      simple_plot->GetXaxis()->SetTitleOffset(0.9);
      simple_plot->GetYaxis()->SetTitleOffset(1.0);
      simple_plot->GetXaxis()->SetTitleSize(0.045);
      simple_plot->GetYaxis()->SetTitleSize(0.045);      
      simple_plot->GetXaxis()->SetLabelSize(0.045);
      simple_plot->GetYaxis()->SetLabelSize(0.045);

      fit_results = model->GetResult();
      std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result->minNll()<<" covQual = " << result->covQual() << std::endl;

      result->Print("v");
      model->PrintResult();

      TFile f(Form("results/manyfits/%s/%s/%s%s_result.root",dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str()),"RECREATE");
      result->Write(); 
      f.Close();

      N_signalMC[i]         = fit_results[both][both][DD0][dsmode][both][both]["Peak"];
      N_backgroundMC[i]     = fit_results[both][both][DD0][dsmode][both][both]["Comb"];
      N_background_redMC[i] = fit_results[both][both][DD0][dsmode][both][both]["Comb_Reduced"];
      std::cout << "            N Signal: " << N_signalMC[i] << std::endl;
      std::cout << "        N Background: " << N_backgroundMC[i] << std::endl;
      std::cout << "N Background reduced: " << N_background_redMC[i] << std::endl;
      
      // Write fit result to file
      
      std::string text_filename  = Form("results/manyfits/text/%s/%s/MC/FitResult_%s%s.txt",dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str());
      std::ofstream output_file(text_filename.c_str());
      output_file << cut_value_MC[i] <<":"<< N_signalMC[i] <<":"<< N_background_redMC[i];
      output_file.close();

      TPaveLabel *pav = new TPaveLabel(0.85,0.85,0.95,0.95,decay_title.c_str(),"NDC");
      pav->SetBorderSize(0);
      pav->SetFillStyle(0);
      pav->SetTextFont(12);
      pav->SetTextSize(0.5);
      pav->SetTextAlign(31);
      pav->Draw();
      TPaveLabel *pav1 = new TPaveLabel(0.85,0.75,0.95,0.85,Form("Signal: %f",N_signalMC[i]),"NDC");
      pav1->SetBorderSize(0);
      pav1->SetFillStyle(0);
      pav1->SetTextFont(62);
      pav1->SetTextSize(0.3);
      pav1->SetTextAlign(31);
      pav1->Draw();
      TPaveLabel *pav2 = new TPaveLabel(0.85,0.65,0.95,0.75,Form("Comb BG: %f",N_background_redMC[i]),"NDC");
      pav2->SetBorderSize(0);
      pav2->SetFillStyle(0);
      pav2->SetTextFont(62);
      pav2->SetTextSize(0.3);
      pav2->SetTextAlign(31);
      pav2->Draw();
      canvas->Print(Form("results/manyfits/%s/%s/%s%s.eps",dsmode.c_str(), mva.c_str(),canvas->GetName(),years.c_str()));
      canvas->Print(Form("results/manyfits/%s/%s/%s%s.pdf",dsmode.c_str(), mva.c_str(),canvas->GetName(),years.c_str()));
      //canvas->Print(Form("results/manyfits/%s/%s/%s.png",dsmode.c_str(), mva.c_str(),canvas->GetName()));
      

      //delete BDTdataset;
      delete cutdata;
       
    }
    // Calculate Punzi and significance
    double Punzi_1D[n_p], Significance_1D[n_p];

    for (int i = 0; i < n_p ; i++){
      Punzi_1D[i]          = N_signalMC[i]/ (2.5 + sqrt(N_background_redMC[i]));
      Significance_1D[i]   = N_signalMC[i]/sqrt(N_signalMC[i]+N_background_redMC[i]);
      
    }
    // Setup plotting options
    SetLHCbStyle(norm);

    std::cout << std::endl;
    for (int i=0; i<n_p; i++ ) {
      std::cout << "Cut Value: " << cut_value_MC[i] << " \t N Sig: " << N_signalMC[i] <<" \t N Bg: " << N_background_redMC[i] << " \t Punzi: " << Punzi_1D[i] <<" \t Significance: " << Significance_1D[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "double cut_value_MC []={";
    for (int i=0; i<n_p; i++ ) {
      std::cout << cut_value_MC[i];
      if(i!=(n_p-1))std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;

    std::cout << "double Punzi_1D []={";
    for (int i=0; i<n_p; i++ ) {
      std::cout << Punzi_1D[i];
      if(i != (n_p-1)) std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;

    std::cout << "double Significance_1D []={";
    for (int i=0; i<n_p; i++ ) {
      std::cout << Significance_1D[i];
      if(i!=(n_p-1)) std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;
  

    TCanvas *cav_punzi = new TCanvas(Form("%s_%s_punzi",dsmode.c_str(), mva.c_str() ),Form("%s %s Punzi as function of cut value",dsmode.c_str(), mva.c_str()),200,10,700,500);
    TGraph  *gr_punzi  = new TGraph(n_p, cut_value_MC, Punzi_1D);
    gr_punzi->Draw("ALP");
    gr_punzi->SetTitle(Form("Punzi    S/(5/2 + sqrt(B));MC %s Cut Value;Punzi",mva.c_str()));    
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_MC.eps",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    cav_punzi->Print(Form("results/manyfits/%s/%s/%s%s_MC.pdf",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()));
    //cav_punzi->Print(Form("results/manyfits/%s/%s/%s_MC.png",dsmode.c_str(), mva.c_str(),cav_punzi->GetName()));
    TFile fa(Form("results/manyfits/%s/%s/%s%s_MC_plots.root",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()),"RECREATE");
    cav_punzi->Write();
    gr_punzi->Write(); 
    fa.Close();

    TCanvas *cav_sig = new TCanvas(Form("%s_%s_significance",dsmode.c_str(), mva.c_str() ),Form("%s %s Significance as function of cut value",dsmode.c_str(), mva.c_str()),200,10,700,500);
    TGraph  *gr_sig  = new TGraph(n_p, cut_value_MC, Significance_1D);
    gr_sig->Draw("ALP");
    gr_sig->SetTitle(Form("Significance    S/sqrt(S+B);MC %s Cut Value;Significance",mva.c_str()));    
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_MC.eps",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()));
    cav_sig->Print(Form("results/manyfits/%s/%s/%s%s_MC.pdf",dsmode.c_str(), mva.c_str(),cav_sig->GetName(),years.c_str()));
    //cav_sig->Print(Form("results/manyfits/%s/%s/%s_MC.png",dsmode.c_str(), mva.c_str(),cav_sig->GetName()));
    TFile fb(Form("results/manyfits/%s/%s/%s%s_MC_plots.root",dsmode.c_str(), mva.c_str(),cav_punzi->GetName(),years.c_str()),"RECREATE");
    cav_punzi->Write();
    gr_punzi->Write(); 
    fb.Close();
  }
*/

}

void DKst0Fitting::RunEfficiency()
{

  if(par->debug) std::cout<<"Running: DKst0Fitting::RunEfficiency() -- Not implemented"<<std::endl;
/*
  std::string mva(BDTG);
  std::string mvatype(DATA);
  
  if( mva!=BDT && mva!=BDTG && mva!= BDTB) {
    std::cout << " MVA method must be one of BDT, BDTG, or BDTB" << std::endl;
    return;
  }

  if(par->debug) std::cout<<"Running: DKst0Fitting::RunEfficiency() --> Chosen setup: " << mvatype << ", " << mva<<std::endl;

  if(yearList.size()*HelBinList.size()*DsBDTBinList.size()*PhiBDTBinList.size()*BmodeList.size()*chargeList.size()*magnetList.size()*modeList.size()>1) {
    std::cout << "ERROR Can only run efficiency fit over a single mode at a time" << std::endl;
    return;
  }

  std::string dsmode;
  if(par->modes[Ds2PhiPi]) dsmode = Ds2KKPi;  
  if(par->modes[Ds2KKPi])  dsmode = Ds2KKPi;
  if(par->modes[Ds2PiPiPi])dsmode = Ds2PiPiPi;
  if(par->modes[Ds2KPiPi]) dsmode = Ds2KPiPi;

  std::string bmode;
  if(par->Bmodes[DsD0])   bmode = DsD0;  
  if(par->Bmodes[DsPhi])  bmode = DsPhi;

  std::string years="";
  if(par->dsetsReq[s26] )   years="_2016";
  if(par->dsetsReq[s24] )   years="_2015";
  if(par->dsetsReq[s21] )   years="_2012";
  if(par->dsetsReq[s21r1] ) years="_2011";
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1]) years = "_Run1"; 
  if(par->dsetsReq[s24] && par->dsetsReq[s26]) years = "_Run2"; 
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1] && par->dsetsReq[s24]&& par->dsetsReq[s26]) years = "_All"; 

  if(par->debug) std::cout<<"Running: DKst0Fitting::RunEfficiency() --> Chosen setup: " << bmode << ", " << dsmode<<", " << years<<std::endl;


  // Set up fit range
  int n_p = 2;
  //if (par->nBDTPoints>0) n_p = par->nBDTPoints;
  
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunEfficiency() --> Number of points: " << n_p<<std::endl;

  std::map<std::string,std::map<std::string,double>> low;
  std::map<std::string,std::map<std::string,double>> high;
  

  std::map<std::string,std::map<std::string,std::string>> var_name;
  var_name[BDT][D]  = bdtD.GetName();
  var_name[BDTG][D] = bdtgD.GetName();
  var_name[BDTB][Ds] = bdtbD.GetName();
  if (bmode==DsD0){
    var_name[BDT][Kst0] = bdtD0.GetName();
    var_name[BDTG][Phi]= bdtgD0.GetName();
    var_name[BDTB][Phi]= bdtbD0.GetName();
  } else {
    var_name[BDT][Phi] = bdtKst0.GetName();
    var_name[BDTG][Phi]= bdtgKst0.GetName();
    var_name[BDTB][Phi]= bdtbKst0.GetName();
  }


  low[BDT][Phi]  = -1.0;  high[BDT][Phi] =  1.1;  
  low[BDT][Ds]   = -1.0;  high[BDT][Ds]  =  1.1;
  low[BDTG][Phi] = -1.0;  high[BDTG][Phi]=  1.1;
  low[BDTG][Ds]  = -1.0;  high[BDTG][Ds] =  1.1;
  low[BDTB][Phi] = -1.0;  high[BDTB][Phi]=  1.1;
  low[BDTB][Ds]  = -1.0;  high[BDTB][Ds] =  1.1;

  double cut_value_Kst0[n_p];
  double cut_value_D[n_p]; 

  double N_signal_Phi[n_p][n_p], N_background_Phi[n_p][n_p], N_background_red_Phi[n_p][n_p];
  double N_signal_Ds[n_p][n_p],  N_background_Ds[n_p][n_p],  N_background_red_Ds[n_p][n_p];

  double binWidth=2;//MeV/c2
  mD0.setBins( (mD0.getMax() -mD0.getMin() )/binWidth);
  mKst0.setBins((mKst0.getMax()-mKst0.getMin())/binWidth);
  mDs.setBins( (mDs.getMax() -mDs.getMin() )/binWidth);

  RooCategory*      cat_ds=model_ds->Cat();
  RooSimultaneous*  sim_ds=model_ds->Pdf();

  RooCategory*      cat_phi=model_phi->Cat();
  RooSimultaneous*  sim_phi=model_phi->Pdf();

  RooFitResult* result_phi = 0;
  RooFitResult* result_ds  = 0;
  for (int i = 0; i < n_p ; i++){
    for (int j = 0; j < n_p ; j++){
      cut_value_Ds[i]  = low[mva][Ds]  + i*(high[mva][Ds] -low[mva][Ds] )/n_p;
      cut_value_Phi[j] = low[mva][Phi] + j*(high[mva][Phi]-low[mva][Phi])/n_p;
           
      std::string cut_BDT  = Form("%s>%f&&%s>%f", var_name[mva][Ds].c_str(), cut_value_Ds[i], var_name[mva][Phi].c_str(), cut_value_Phi[j]);
      
      if(par->debug) std::cout << "--- Running fit with: " << cut_BDT << std::endl;
      if(par->debug) std::cout << "--- No. events before: " << data->numEntries() << std::endl;
      cutdata = (RooDataSet*) data->reduce(cut_BDT.c_str());     
      if(par->debug) std::cout << "--- No. events after: " << cutdata->numEntries() << std::endl;

      RooAbsData* abs=(RooAbsData*) cutdata;

      if(par->debug) std::cout << "--- Doing Fit to Phi mass: " << std::endl;
      result_phi = sim_phi->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
      
      if(par->debug) std::cout << "--- Doing Fit to Ds mass: " << std::endl;
      result_ds  = sim_ds->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
      
      // Plotting Phi mass
      std::string name_string=  Form("Phi_fit_%s_%s_Ds_%f_Phi_%f",dsmode.c_str(), mva.c_str(),cut_value_Ds[i], cut_value_Phi[j] );
      TCanvas* canvas_phi = new TCanvas(name_string.c_str(),"",40,0,600,600);
      RooPlot* simple_plot_phi;
      if(bmode==DsD0) { simple_plot_phi = mD0.frame();
      }else {simple_plot_phi = mKst0.frame();}

      std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s_%s",
                                 (*yearList.begin()).c_str(),
                                 (*HelBinList.begin()).c_str(),
                                 (*DsBDTBinList.begin()).c_str(),
                                 (*PhiBDTBinList.begin()).c_str(),
                                 (*BmodeList.begin()).c_str(),
                                 (*modeList.begin()).c_str(),
                                 (*chargeList.begin()).c_str(),
                                 (*magnetList.begin()).c_str() );
      std::string decay_title = title[*yearList.begin()][*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()];    
      cat_phi->setLabel(cat_str.c_str());
      cutdata->plotOn( simple_plot_phi,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      // Comb
      sim_phi->plotOn(  simple_plot_phi,
                    RooFit::Slice(RooArgSet(*cat_phi)), 
                    RooFit::ProjWData(RooArgSet(*cat_phi),*cutdata) ,
                    RooFit::Components(Form("pdf_Phi_comb_%s",cat_str.c_str())),
                    RooFit::FillColor(kAzure-7), 
                    RooFit::DrawOption("F") ); 
      // Total
      sim_phi->plotOn(  simple_plot_phi,
                    RooFit::Slice(RooArgSet(*cat_phi)), 
                    RooFit::ProjWData(RooArgSet(*cat_phi),*cutdata) ,
                    RooFit::LineWidth(2),
                    RooFit::LineColor(kBlack) ); 
      cutdata->plotOn( simple_plot_phi,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      simple_plot_phi->Draw();
      simple_plot_phi->Draw();
      simple_plot_phi->SetTitle("");

      if(par->Bmodes[DsD0])  simple_plot_phi->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D^{0}}) [MeV/#font[12]{c}^{2}] ");
      if(par->Bmodes[DsPhi]) simple_plot_phi->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{#phi) [MeV/#font[12]{c}^{2}] ");
      simple_plot_phi->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));
      simple_plot_phi->GetXaxis()->SetTitleOffset(0.9);
      simple_plot_phi->GetYaxis()->SetTitleOffset(1.0);
      simple_plot_phi->GetXaxis()->SetTitleSize(0.045);
      simple_plot_phi->GetYaxis()->SetTitleSize(0.045);      
      simple_plot_phi->GetXaxis()->SetLabelSize(0.045);
      simple_plot_phi->GetYaxis()->SetLabelSize(0.045);

      fit_results_phi = model_phi->GetResult();
      std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result_phi->minNll()<<" covQual = " << result_phi->covQual() << std::endl;

      result_phi->Print("v");
      model_phi->PrintResult();

      TFile f(Form("results/efficiency/%s/%s/%s/%s%s_Phi_result.root",bmode.c_str(),dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str()),"RECREATE");
      result_phi->Write(); 
      f.Close();

      N_signal_Phi[i][j]         = fit_results_phi[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Peak"];
      N_background_Phi[i][j]     = fit_results_phi[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb"];
      N_background_red_Phi[i][j] = fit_results_phi[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb_Reduced"];
      std::cout << "            N Signal Phi: " << N_signal_Phi[i][j] << std::endl;
      std::cout << "        N Background Phi: " << N_background_Phi[i][j] << std::endl;
      std::cout << "N Background reduced Phi: " << N_background_red_Phi[i][j] << std::endl;
      std::string text_filename  = Form("results/efficiency/text/%s/%s/%s/FitResult_Phi_%s%s.txt",bmode.c_str(),dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str());
      std::ofstream output_file(text_filename.c_str());
      output_file << cut_value_Ds[i] <<":"<< cut_value_Phi[j] <<":"<< N_signal_Phi[i][j] << ":" << N_background_red_Phi[i][j];
      output_file.close();

      TPaveLabel *pav = new TPaveLabel(0.85,0.85,0.95,0.95,decay_title.c_str(),"NDC");
      pav->SetBorderSize(0);
      pav->SetFillStyle(0);
      pav->SetTextFont(12);
      pav->SetTextSize(0.5);
      pav->SetTextAlign(31);
      pav->Draw();
      TPaveLabel *pav1 = new TPaveLabel(0.85,0.75,0.95,0.85,Form("Signal: %f",N_signal_Phi[i][j]),"NDC");
      pav1->SetBorderSize(0);
      pav1->SetFillStyle(0);
      pav1->SetTextFont(62);
      pav1->SetTextSize(0.3);
      pav1->SetTextAlign(31);
      pav1->Draw();
      TPaveLabel *pav2 = new TPaveLabel(0.85,0.65,0.95,0.75,Form("Comb BG: %f",N_background_Phi[i][j]),"NDC");
      pav2->SetBorderSize(0);
      pav2->SetFillStyle(0);
      pav2->SetTextFont(62);
      pav2->SetTextSize(0.3);
      pav2->SetTextAlign(31);
      pav2->Draw();
      canvas_phi->Print(Form("results/efficiency/%s/%s/%s/%s%s.eps",bmode.c_str(),dsmode.c_str(), mva.c_str(),canvas_phi->GetName(),years.c_str()));
      canvas_phi->Print(Form("results/efficiency/%s/%s/%s/%s%s.pdf",bmode.c_str(),dsmode.c_str(), mva.c_str(),canvas_phi->GetName(),years.c_str()));
      
      // ----- plot Ds fit --------
      name_string=  Form("Ds_fit_%s_%s_Ds_%f_Phi_%f",dsmode.c_str(), mva.c_str(),cut_value_Ds[i], cut_value_Phi[j] );
      TCanvas* canvas_ds = new TCanvas(name_string.c_str(),"",40,0,600,600);
      RooPlot* simple_plot_ds = mDs.frame();

      cat_ds->setLabel(cat_str.c_str());
      cutdata->plotOn( simple_plot_ds,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      // Comb
      sim_ds->plotOn(  simple_plot_ds,
                    RooFit::Slice(RooArgSet(*cat_ds)), 
                    RooFit::ProjWData(RooArgSet(*cat_ds),*cutdata) ,
                    RooFit::Components(Form("pdf_Ds_comb_%s",cat_str.c_str())),
                    RooFit::FillColor(kAzure-7), 
                    RooFit::DrawOption("F") ); 
      // Total
      sim_ds->plotOn(  simple_plot_ds,
                    RooFit::Slice(RooArgSet(*cat_ds)), 
                    RooFit::ProjWData(RooArgSet(*cat_ds),*cutdata) ,
                    RooFit::LineWidth(2),
                    RooFit::LineColor(kBlack) ); 
      cutdata->plotOn( simple_plot_ds,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      simple_plot_ds->Draw();
      simple_plot_ds->Draw();
      simple_plot_ds->SetTitle("");

      if(par->Bmodes[DsD0])  simple_plot_ds->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D_{s}}) [MeV/#font[12]{c}^{2}] ");
      simple_plot_ds->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));
      simple_plot_ds->GetXaxis()->SetTitleOffset(0.9);
      simple_plot_ds->GetYaxis()->SetTitleOffset(1.0);
      simple_plot_ds->GetXaxis()->SetTitleSize(0.045);
      simple_plot_ds->GetYaxis()->SetTitleSize(0.045);      
      simple_plot_ds->GetXaxis()->SetLabelSize(0.045);
      simple_plot_ds->GetYaxis()->SetLabelSize(0.045);

      fit_results_ds = model_ds->GetResult();
      std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result_ds->minNll()<<" covQual = " << result_ds->covQual() << std::endl;

      result_ds->Print("v");
      model_ds->PrintResult();

      TFile f_ds(Form("results/efficiency/%s/%s/%s/%s%s_Ds_result.root",bmode.c_str(),dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str()),"RECREATE");
      result_ds->Write(); 
      f_ds.Close();

      N_signal_Ds[i][j]         = fit_results_ds[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Peak"];
      N_background_Ds[i][j]     = fit_results_ds[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb"];
      N_background_red_Ds[i][j] = fit_results_ds[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb_Reduced"];
      std::cout << "            N Signal Ds: " << N_signal_Ds[i][j] << std::endl;
      std::cout << "        N Background Ds: " << N_background_Ds[i][j] << std::endl;
      std::cout << "N Background reduced Ds: " << N_background_red_Ds[i][j] << std::endl;
      text_filename  = Form("results/efficiency/text/%s/%s/%s/FitResult_Ds_%s%s.txt",bmode.c_str(),dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str());
      std::ofstream output_file_ds(text_filename.c_str());
      output_file_ds << cut_value_Ds[i] <<":"<< cut_value_Phi[j] <<":"<< N_signal_Ds[i][j] << ":" << N_background_red_Ds[i][j];
      output_file_ds.close();

      TPaveLabel *pav_ds = new TPaveLabel(0.85,0.85,0.95,0.95,decay_title.c_str(),"NDC");
      pav_ds->SetBorderSize(0);
      pav_ds->SetFillStyle(0);
      pav_ds->SetTextFont(12);
      pav_ds->SetTextSize(0.5);
      pav_ds->SetTextAlign(31);
      pav_ds->Draw();
      TPaveLabel *pav1_ds = new TPaveLabel(0.85,0.75,0.95,0.85,Form("Signal: %f",N_signal_Ds[i][j]),"NDC");
      pav1_ds->SetBorderSize(0);
      pav1_ds->SetFillStyle(0);
      pav1_ds->SetTextFont(62);
      pav1_ds->SetTextSize(0.3);
      pav1_ds->SetTextAlign(31);
      pav1_ds->Draw();
      TPaveLabel *pav2_ds = new TPaveLabel(0.85,0.65,0.95,0.75,Form("Comb BG: %f",N_background_Ds[i][j]),"NDC");
      pav2_ds->SetBorderSize(0);
      pav2_ds->SetFillStyle(0);
      pav2_ds->SetTextFont(62);
      pav2_ds->SetTextSize(0.3);
      pav2_ds->SetTextAlign(31);
      pav2_ds->Draw();
      canvas_ds->Print(Form("results/efficiency/%s/%s/%s/%s%s.eps",bmode.c_str(),dsmode.c_str(), mva.c_str(),canvas_ds->GetName(),years.c_str()));
      canvas_ds->Print(Form("results/efficiency/%s/%s/%s/%s%s.pdf",bmode.c_str(),dsmode.c_str(), mva.c_str(),canvas_ds->GetName(),years.c_str()));
      

      delete cutdata;
    }   
  } 
*/
}


void DKst0Fitting::RunManyToys()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::RunManyToys()"<<std::endl;
  RooArgSet* initialPars= model_gen->GetParameters();
  //RooArgSet* fittedPars = model->GetParameters();
  /*
  std::cout<<"Copying the last fitted parameters to the toy generator:"<<std::endl;
  TIterator* initialIter= initialPars->createIterator();
  TObject* initialObj;  TObject* fittedObj;
  while((initialObj=initialIter->Next())) {
    RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
    if(!initialRrv) continue;
    TIterator* fittedIter = fittedPars->createIterator();
    while((fittedObj=fittedIter->Next())) {
      RooRealVar* fittedRrv = dynamic_cast<RooRealVar*>(fittedObj);
      if(!fittedRrv) continue;
      if(std::string(initialRrv->GetName())==std::string(fittedRrv->GetName())){
        std::cout<<" "<<initialRrv->GetName()<<" <--- "<<fittedRrv->getVal()<<" (default: "<<initialRrv->getVal()<<")"<<std::endl;
        initialRrv->setVal(fittedRrv->getVal());
      }
    }
    delete fittedIter;
  }
  delete initialIter;
  */

  // Set Generator values to RooRealVar initial values and write them to a file

  std::string gen_location=par->toylocation+"/gen_vals";
  TIterator* initialIter= initialPars->createIterator();
  TObject* initialObj; 
  while((initialObj=initialIter->Next())) {
    RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
    if(!initialRrv) continue;
    std::string text_filename  = Form("%s/Init_Val_%s.txt", gen_location.c_str(),initialRrv->GetName());

    std::cout<<" "<<initialRrv->GetName()<<" <--- "<<initialRrv->getVal()<<" --> File: " << text_filename <<std::endl;
    
    std::ofstream output_file(text_filename.c_str());
    output_file << initialRrv->getVal();
    output_file.close();
  }
  delete initialIter;

  std::cout<<" "<<std::endl;
  OrderToys(par->nToys);
}



void DKst0Fitting::OrderToys(int n)
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::OrderToys()"<<std::endl;
  RooAbsPdf* sim = model_gen->Pdf();
  RooAbsCategory* cat = model_gen->Cat();
  if(par->debug) std::cout<<"Running: DKst0Fitting::OrderToys() ----> Created Model"<<std::endl;
  bool DB=false;
  int secs=time(NULL);
  if(1==n) secs=20110602;
  if(par->useSeed) secs = par->seed;
  std::cout <<  "====================== " << std::endl;
  std::cout <<  "Using Seed: " << secs << std::endl;
  std::cout <<  "====================== " << std::endl;

  RooRandom::randomGenerator()->SetSeed(secs);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  if(par->debug) std::cout<<"Running: DKst0Fitting::OrderToys() ----> Creating RooMCStudy:"<<std::endl;
  RooMCStudy* mcstudy = new RooMCStudy(*sim ,reducedlist,RooFit::Binned(par->binned),RooFit::FitOptions(RooFit::Save(true),RooFit::PrintLevel(DB?1:-1)));

  if(par->debug) std::cout<<"Running: DKst0Fitting::OrderToys() ----> Created RooMCStudy"<<std::endl;
  int nEvtsPerSample=sim->expectedEvents(*cat);
 
  if(par->debug) std::cout<<"Running: DKst0Fitting::OrderToys() ----> Expected Events: " << nEvtsPerSample <<std::endl; 
  if(1==n&&!DB){
    mcstudy->generate(n,nEvtsPerSample,true);
    data=(RooDataSet*)mcstudy->genData(n-1);
  }else{
    // Fit using fit model 
    if(true){
      std::cout<<"Generating "<<n<<" toys of "<<nEvtsPerSample<<" events."<<std::endl;
      mcstudy->generate(n,nEvtsPerSample,true);
      
      gSystem->Exec((std::string("mkdir -p ")+par->toylocation).c_str());
      std::string toyfile=par->toylocation+"/toy_";
      for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){toyfile+=(*m)+underscore;}
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){toyfile+=(*b)+underscore;}
      toyfile+=std::string(Form("toyFile.1f_%i.root",secs));  
      TFile f(toyfile.c_str(),"RECREATE");
      
      for(int i=0;i<n;i++){
        // Reset fit initial values to the generated ones...
        RooArgSet* initialPars= model_gen->GetParameters();
        std::string gen_location=par->toylocation+"/gen_vals";
        TIterator* initialIter= initialPars->createIterator();
        TObject* initialObj; 
        while((initialObj=initialIter->Next())) {
          RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
          if(!initialRrv) continue;
          std::string text_filename  = Form("%s/Init_Val_%s.txt", gen_location.c_str(),initialRrv->GetName());

          std::ifstream input;
          input.open(text_filename.c_str(),std::ifstream::in);
          if(!input){
            std::cout<<" Can't find file " << text_filename << std::endl; return;
          }
          std::string line = "";
          getline(input,line);
          initialRrv->setVal(atof(line.c_str()));
          input.close();
          }
        delete initialIter;


        std::cout<<"Fitting toy "<< i+1 <<" of "<<n<<"."<<std::endl;
        RooDataSet* mcdata=(RooDataSet*)mcstudy->genData(i);
        RooFitResult* mcresult = sim->fitTo(*mcdata,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(-1),RooFit::Optimize(false),RooFit::Timer(true));
        mcresult->Write(Form("toy%i",i));
      }
      gDirectory->ls();
      f.Close();
    }
    // Fit using MC model 
    if(false){
      std::cout<<"Generating "<<n<<" toys of "<<nEvtsPerSample<<" events."<<std::endl;
      mcstudy->generateAndFit(n,nEvtsPerSample);
      //if(DB){return;}
      
      gSystem->Exec((std::string("mkdir -p ")+par->toylocation).c_str());
      std::string toyfile=par->toylocation+"/toy_";
      for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){toyfile+=(*m)+underscore;}
      for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){toyfile+=(*b)+underscore;}
      toyfile+=std::string(Form("toyFile.1f_%i.root",secs));  
      TFile f(toyfile.c_str(),"RECREATE");
      for(int i=0;i<n;i++){
        mcstudy->fitResult(i)->Write(Form("toy%i",i));
      }
      mcstudy->fitParDataSet().Write("toyPars");
      gDirectory->ls();
      f.Close();
    }

  }
}



void DKst0Fitting::DisplayToys()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::DisplayToys()"<<std::endl;
  RooSimultaneous*  sim=model_gen->Pdf();
  std::vector<RooFitResult> Results; 
  RooMCStudy* mcstudy      = new RooMCStudy(*sim,reducedlist,RooFit::FitOptions("r"),RooFit::Extended());
  RooMCStudy* mcstudy_test = new RooMCStudy(*sim,reducedlist,RooFit::Extended());
  int nFits=0;
  int nBad=0;
  int nBadMINOS=0;
  int nGood=0;
  int nFPD=0;
  std::string toyfile("toy_");
  for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){toyfile+=(*m)+underscore;}
  for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){toyfile+=(*b)+underscore;}
  //toyfile+=std::string(Form("pid%3.1f_",m_pidCut));
  toyfile+=std::string("toyFile.1f_");  //need to check what the name should really be...
  std::cout <<" Looking for files like: " << par->toylocation+slash+toyfile <<"*"<< std::endl;
  std::vector<std::string> filename;
  CommonTools::getdir(par->toylocation.c_str(),filename);
  for(unsigned int i=0;i<filename.size();i++){
    if(filename[i].find(toyfile) == std::string::npos) continue;
    TFile f((par->toylocation+slash+filename[i]).c_str(),"READ");
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    while((key=(TKey*)nextkey())){
      if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
      RooFitResult *r = (RooFitResult*)key->ReadObj();
      nFits++;      
      if(r->covQual()<2){nBad++;continue;}
      if(r->covQual()<3){nFPD++;continue;}
      if(r->status()!=0){nBadMINOS++;}//continue;}
      nGood++;
      //std::cout << "minNLL: " <<  r->minNll() << std::endl;
      //std::cout << "status: " <<  r->status() << std::endl;
      //std::cout << "covQual: " <<  r->covQual() << std::endl;
      std::cout <<"\r   Collecting "<<nFits<<" toys, of which "<<nBad<<" don't converge and "<<100*float(nFPD)/nFits<<"% are forced positive definite and "<<100*float(nBadMINOS)/nFits<<"% have MINOS problems. "<< nGood << " are good."<<std::flush;
      //std::cout << "Added? "<< std::endl;
      mcstudy->addFitResult(*r);
      mcstudy_test->addFitResult(*r);
      Results.push_back(*r);

    }
  }

  std::cout<<"\n"<<std::endl;
  
  if(1==nFits){
    std::cout<<"SOMETHING WRONG. NO TOYS LOADED"<<std::endl;
    return;
  }

  // Get initial values from files...

  if(par->debug) std::cout<<"Running: DKst0Fitting::DisplayToys()--> Getting Generated values from files"<<std::endl;
  std::map<std::string,double> Generated_values;
  RooArgSet* initialPars= model_gen->GetParameters();
  std::string gen_location=par->toylocation+"/gen_vals";
  TIterator* initialIter= initialPars->createIterator();
  TObject* initialObj; 
  while((initialObj=initialIter->Next())) {
    RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
    if(!initialRrv) continue;
    std::string text_filename  = Form("%s/Init_Val_%s.txt", gen_location.c_str(),initialRrv->GetName());

    std::ifstream input;
    input.open(text_filename.c_str(),std::ifstream::in);
    if(!input){
      std::cout<<" Can't find file " << text_filename << std::endl; return;
    }
    std::string line = "";
    getline(input,line);
    Generated_values[initialRrv->GetName()] = atof(line.c_str()); 
    std::cout<<" "<<initialRrv->GetName()<<" <--- "<< line <<" --> From file: " << text_filename <<std::endl;
    
    input.close();
  }
  delete initialIter;


  if(true){
    if(par->debug) std::cout<<"Running: DKst0Fitting::DisplayToys()--> Calculating Asymmetric Error Pulls"<<std::endl;
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    std::vector<TCanvas*> v_canvas;
    std::map<int,int> Dodgy_pulls;
    RooArgSet* vars = sim->getVariables();
    std::vector<RooRealVar*> variables;
    RooRealVar* temp_var=0;
    int nv=0;
    int n_canvas=0;
    TIterator* it3 = vars->createIterator();
    while((temp_var = (RooRealVar*)it3->Next())) {
      if(temp_var->isConstant() or temp_var->InheritsFrom("RooAbsCategory") or 0==std::string(temp_var->GetName()).find(mB.GetName())) continue;
      nv++;
      variables.push_back(temp_var);      
      std::vector<double> v_values, v_errors, v_pulls;
      
      //RooRealVar* fit_ini = (RooRealVar*)initialPars->find(temp_var->GetName());
      std::cout << "Making pull plots for variable: " << nv << "  " << temp_var->GetName() <<std::endl;
      //double init_val = fit_ini->getVal();

      
      for(unsigned int i = 0; i< Results.size(); i++){
        RooFitResult rfr_temp = Results[i];
        RooRealVar* param_final = (RooRealVar*) (rfr_temp.floatParsFinal()).find(temp_var->GetName());
        //RooRealVar* param_init  = (RooRealVar*) (rfr_temp.floatParsInit()).find(temp_var->GetName());
        double final_value = param_final->getVal();
        //double init_value  = param_init->getVal();
        double init_value = Generated_values[temp_var->GetName()];

        double error  = param_final->getError();
        double err_hi = param_final->getAsymErrorHi();
        double err_lo = param_final->getAsymErrorLo();

        double pull;
        if(final_value < init_value) pull = (init_value - final_value)/err_hi;
        else                         pull = (final_value - init_value)/err_lo;

        if(fabs(pull)>5){
          std::cout << "Pull outside range---> Toy Number: "<< i  << "\t Pull: " << pull << "\t Value: " << final_value << "\t Error Low: " << err_lo<< "\t Error High: " << err_hi << std::endl;
          Dodgy_pulls[i]++;
        } else {
          v_values.push_back(final_value);

          //if(final_value < init_value) v_errors.push_back(err_hi);
          //else                         v_errors.push_back(err_lo);
          v_errors.push_back(error);
          v_pulls.push_back(pull);
        }        
      }

      double max_val = *max_element(v_values.begin(), v_values.end());
      double min_val = *min_element(v_values.begin(), v_values.end());
      
      double max_err = *max_element(v_errors.begin(), v_errors.end());
      double min_err = *min_element(v_errors.begin(), v_errors.end());
      
      double max_pul = *max_element(v_pulls.begin(), v_pulls.end());
      double min_pul = *min_element(v_pulls.begin(), v_pulls.end());
      max_pul =  5;
      min_pul = -5;
      SetLHCbStyle(norm);
      gStyle->SetOptStat(0);
      gStyle->SetMarkerSize(0.5);

      TCanvas* c  = new TCanvas(Form("canv_%s",temp_var->GetName()),temp_var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;
      TH1F* h_val = new TH1F(Form("Title;%s Value;Number of toys;", temp_var->GetTitle()), "Vars", 20, min_val-fabs(min_val)*0.1, max_val+fabs(max_val)*0.1);
      TH1F* h_err = new TH1F(Form("Title;%s Error; Number of toys;",temp_var->GetTitle()), "Errs", 20, min_err-fabs(min_err)*0.1, max_err+fabs(max_err)*0.1);
      TH1F* h_pul = new TH1F(Form("Title;%s Pulls; Number of toys;",temp_var->GetTitle()), "Pull", 20, min_pul, max_pul);

      for(unsigned int i = 0; i< v_values.size(); i++){
        h_val->Fill(v_values[i]);
        h_err->Fill(v_errors[i]);
        h_pul->Fill(v_pulls[i]);     
      } 

      c->cd(1) ; gPad->SetLeftMargin(0.15); h_val->GetXaxis()->SetTitle(Form("%s Value",temp_var->GetTitle())); h_val->GetXaxis()->SetNdivisions(5); h_val->GetYaxis()->SetTitleOffset(1.4); h_val->Draw("PE0") ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); h_err->GetXaxis()->SetTitle(Form("%s Error",temp_var->GetTitle())); h_err->GetXaxis()->SetNdivisions(5); h_err->GetYaxis()->SetTitleOffset(1.4); h_err->Draw("PE0") ;
      c->cd(3) ; gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->GetXaxis()->SetNdivisions(5); h_pul->GetYaxis()->SetTitleOffset(1.4); h_pul->Draw("PE0") ;
      double max=h_pul->GetMaximum();
      h_pul->SetMaximum(max*1.3);
      TF1 *gaussian = new TF1("Gaussian","gaus",min_pul,max_pul);
      h_pul->Fit("Gaussian","Q0R");
      double mean  = gaussian->GetParameter(1);
      double sigma = gaussian->GetParameter(2);
      const double* err=gaussian->GetParErrors();
      double meanerr=err[1];
      double sigmaerr=err[2];
      c->cd(3) ; gaussian->SetLineColor(kBlue); gaussian->Draw("same");

      TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
      TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
      pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
      pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
      pav1->SetFillColor(0);    pav2->SetFillColor(0); 
      pav1->SetTextFont(12);    pav2->SetTextFont(12);  
      pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
      pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
      pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
      pav1->Draw();             pav2->Draw();            

      gStyle->SetTitleFontSize(0.1);
      if ((nv-1)%20==0) {
        std::cout << "Creating Canvas.. " << n_canvas  << std::endl;
        v_canvas.push_back( new TCanvas(Form("canvas_all_%d_pull",n_canvas),"Pulls",20,20,800,800)  );
        v_canvas[n_canvas]->Divide(5,4);
        n_canvas++;
      }

      std::cout << "Moving to Canvas section  " << nv - (n_canvas-1)*20 << std::endl;
      v_canvas[n_canvas-1]->cd( nv - (n_canvas-1)*20 );
      gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->Draw("PE0"); gaussian->Draw("same");pav1->Draw();pav2->Draw(); 
      gPad->Update();
      c->Print(Form("%s/plots/Asym_%s.pdf",(par->toylocation).c_str(), c->GetName()));
    }
    for(unsigned int i = 0; i < v_canvas.size(); i++){
      v_canvas[i]->Print( Form("%s/plots/Asym_%s.pdf",(par->toylocation).c_str(),v_canvas[i]->GetName()) );
    }
    std::cout << "Dodgy Pull list: " <<std::endl;
    for(std::map<int,int>::iterator itF=Dodgy_pulls.begin(); itF!=Dodgy_pulls.end(); itF++) {
      std::cout << "Pull Number: " << itF->first << " count: " << itF->second <<std::endl;   
    }
  }

  // Symmetric errors 
  if(true){
    if(par->debug) std::cout<<"Running: DKst0Fitting::DisplayToys()--> Calculating Symmetric Error Pulls"<<std::endl;
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    std::vector<TCanvas*> v_canvas;
    std::map<int,int> Dodgy_pulls;
    RooArgSet* vars = sim->getVariables();
    std::vector<RooRealVar*> variables;
    RooRealVar* temp_var=0;
    int nv=0;
    int n_canvas=0;
    TIterator* it3 = vars->createIterator();
    while((temp_var = (RooRealVar*)it3->Next())) {
      if(temp_var->isConstant() or temp_var->InheritsFrom("RooAbsCategory") or 0==std::string(temp_var->GetName()).find(mB.GetName())) continue;
      nv++;
      variables.push_back(temp_var);      
      std::vector<double> v_values, v_errors, v_pulls;
      
      std::cout << "Making pull plots for variable: " << nv << "  " << temp_var->GetName() <<std::endl;
      
      for(unsigned int i = 0; i< Results.size(); i++){
        RooFitResult rfr_temp = Results[i];
        RooRealVar* param_final = (RooRealVar*) (rfr_temp.floatParsFinal()).find(temp_var->GetName());
        //RooRealVar* param_init  = (RooRealVar*) (rfr_temp.floatParsInit()).find(temp_var->GetName());
        double final_value = param_final->getVal();
        double init_value = Generated_values[temp_var->GetName()];

        double error  = param_final->getError();
        double err_hi = param_final->getAsymErrorHi();
        double err_lo = param_final->getAsymErrorLo();

        double pull = (init_value - final_value)/error;

        if(abs(pull)>5){
          std::cout << "Pull outside range---> Toy Number: "<< i  << "\t Pull: " << pull << "\t Value: " << final_value << "\t Error Low: " << err_lo<< "\t Error High: " << err_hi << std::endl;
          Dodgy_pulls[i]++;
        } else {
          v_values.push_back(final_value);
          v_errors.push_back(error);
          v_pulls.push_back(pull);

          if(par->debug&&i%100==0) std::cout << "DEBUG::  Value: " << final_value  << ", Gen Value: " << init_value<<", Error: " << error << ", Pull = " << pull << std::endl;
        }
        
      }

      double max_val = *max_element(v_values.begin(), v_values.end());
      double min_val = *min_element(v_values.begin(), v_values.end());
      
      double max_err = *max_element(v_errors.begin(), v_errors.end());
      double min_err = *min_element(v_errors.begin(), v_errors.end());
      
      double max_pul = *max_element(v_pulls.begin(), v_pulls.end());
      double min_pul = *min_element(v_pulls.begin(), v_pulls.end());
      max_pul =  5;
      min_pul = -5;
      SetLHCbStyle(norm);
      gStyle->SetOptStat(0);
      gStyle->SetMarkerSize(0.5);

      TCanvas* c  = new TCanvas(Form("canv_%s",temp_var->GetName()),temp_var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;
      TH1F* h_val = new TH1F(Form("Title;%s Value;Number of toys;", temp_var->GetTitle()), "Vars", 20, min_val-fabs(min_val)*0.1, max_val+fabs(max_val)*0.1);
      TH1F* h_err = new TH1F(Form("Title;%s Error; Number of toys;",temp_var->GetTitle()), "Errs", 20, min_err-fabs(min_err)*0.1, max_err+fabs(max_err)*0.1);
      TH1F* h_pul = new TH1F(Form("Title;%s Pulls; Number of toys;",temp_var->GetTitle()), "Pull", 20, min_pul, max_pul);

      for(unsigned int i = 0; i< v_values.size(); i++){
        h_val->Fill(v_values[i]);
        h_err->Fill(v_errors[i]);
        h_pul->Fill(v_pulls[i]);     
      } 

      c->cd(1) ; gPad->SetLeftMargin(0.15); h_val->GetXaxis()->SetTitle(Form("%s Value",temp_var->GetTitle())); h_val->GetXaxis()->SetNdivisions(5); h_val->GetYaxis()->SetTitleOffset(1.4); h_val->Draw("PE0") ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); h_err->GetXaxis()->SetTitle(Form("%s Error",temp_var->GetTitle())); h_err->GetXaxis()->SetNdivisions(5); h_err->GetYaxis()->SetTitleOffset(1.4); h_err->Draw("PE0") ;
      c->cd(3) ; gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->GetXaxis()->SetNdivisions(5); h_pul->GetYaxis()->SetTitleOffset(1.4); h_pul->Draw("PE0") ;
      double max=h_pul->GetMaximum();
      h_pul->SetMaximum(max*1.3);

      TF1 *gaussian = new TF1("Gaussian","gaus",min_pul,max_pul);
      h_pul->Fit("Gaussian","Q0R");
      double mean  = gaussian->GetParameter(1);
      double sigma = gaussian->GetParameter(2);
      const double* err=gaussian->GetParErrors();
      double meanerr=err[1];
      double sigmaerr=err[2];
      c->cd(3) ; gaussian->SetLineColor(kBlue); gaussian->Draw("same");
      TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
      TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
      pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
      pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
      pav1->SetFillColor(0);    pav2->SetFillColor(0); 
      pav1->SetTextFont(12);    pav2->SetTextFont(12);  
      pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
      pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
      pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
      pav1->Draw();             pav2->Draw();            

      gStyle->SetTitleFontSize(0.1);
      if ((nv-1)%20==0) {
        std::cout << "Creating Canvas.. " << n_canvas  << std::endl;
        v_canvas.push_back( new TCanvas(Form("canvas_all_%d_pull",n_canvas),"Pulls",20,20,800,800)  );
        v_canvas[n_canvas]->Divide(5,4);
        n_canvas++;
      }

      std::cout << "Moving to Canvas section  " << nv - (n_canvas-1)*20 << std::endl;
      v_canvas[n_canvas-1]->cd( nv - (n_canvas-1)*20 );
      gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->Draw("PE0"); gaussian->Draw("same");pav1->Draw();pav2->Draw(); 
      gPad->Update();
      c->Print(Form("toysDir/plots/Sym_%s.pdf",c->GetName()));
    }
    for(unsigned int i = 0; i < v_canvas.size(); i++){
      v_canvas[i]->Print( Form("toysDir/plots/Sym_%s.pdf",v_canvas[i]->GetName()) );
    }
    std::cout << "Dodgy Pull list: " <<std::endl;
    for(std::map<int,int>::iterator itF=Dodgy_pulls.begin(); itF!=Dodgy_pulls.end(); itF++) {
      std::cout << "Pull Number: " << itF->first << " count: " << itF->second <<std::endl;   
    }
  }  





  if(true){
    if(par->debug) std::cout<<"Running: DKst0Fitting::DisplayToys()--> Calculating Symmetric Error Pulls with RooMCStudy"<<std::endl;    
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    std::vector<TCanvas*> v_canvas;

    const RooDataSet mcs = mcstudy->fitParDataSet();
    std::map<int,int> Dodgy_pulls;

    RooArgSet* vars = sim->getVariables();
    std::vector<RooRealVar*> variables;
    RooRealVar* temp_var=0;
    int nv=0;
    int n_canvas=0;
    TIterator* it3 = vars->createIterator();
    while((temp_var = (RooRealVar*)it3->Next())) {
      if(temp_var->isConstant() or temp_var->InheritsFrom("RooAbsCategory") or 0==std::string(temp_var->GetName()).find(mB.GetName())) continue;
      nv++;
      variables.push_back(temp_var);
      
      std::vector<double> v_values, v_errors, v_pulls;
      
      double init_value = Generated_values[temp_var->GetName()];
      std::cout << "Making pull plots for variable: " << nv << "  " << temp_var->GetName() << "\t Generated value: " << init_value <<std::endl;
      
      for (int i = 0; i< mcs.numEntries(); i++){
        const RooArgSet* event = mcs.get(i);  
        RooRealVar* mcs_var = (RooRealVar*)event->find(temp_var->GetName());
        RooRealVar* mcs_err = (RooRealVar*)event->find(Form("%serr",temp_var->GetName()));
        //RooRealVar* mcs_gen = (RooRealVar*)event->find(Form("%sgen",temp_var->GetName()));
        double final_value = mcs_var->getVal();
        double error = mcs_err->getVal();
        double pull = (init_value-final_value)/error; 
        if(fabs(pull)>5){
          std::cout << "Pull outside range---> Toy Number: "<< i  << "\t Pull: " << pull << "\t Value: " << mcs_var->getVal() << "\t Error: " << mcs_err->getVal() << std::endl;
          Dodgy_pulls[i]++;
        } else {
          v_values.push_back(final_value);
          v_errors.push_back(error);
          v_pulls.push_back(pull);
          if(par->debug&&i%100==0) std::cout << "DEBUG::  Value: " << final_value << "Gen Value: " << init_value << ", Error: " << error << ", Pull = " << pull << std::endl;
        }
      }
        
      double max_val = *max_element(v_values.begin(), v_values.end());
      double min_val = *min_element(v_values.begin(), v_values.end());
      
      double max_err = *max_element(v_errors.begin(), v_errors.end());
      double min_err = *min_element(v_errors.begin(), v_errors.end());
      
      double max_pul = *max_element(v_pulls.begin(), v_pulls.end());
      double min_pul = *min_element(v_pulls.begin(), v_pulls.end());
      max_pul =  5;
      min_pul = -5;
      SetLHCbStyle(norm);
      gStyle->SetOptStat(0);
      gStyle->SetMarkerSize(0.5);

      TCanvas* c  = new TCanvas(Form("canv_%s",temp_var->GetName()),temp_var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;
      TH1F* h_val = new TH1F(Form("Title;%s Value;Number of toys;", temp_var->GetTitle()), "Vars", 20, min_val, max_val);
      TH1F* h_err = new TH1F(Form("Title;%s Error; Number of toys;",temp_var->GetTitle()), "Errs", 20, min_err, max_err);
      TH1F* h_pul = new TH1F(Form("Title;%s Pulls; Number of toys;",temp_var->GetTitle()), "Pull", 20, min_pul, max_pul);

      for(unsigned int i = 0; i< v_values.size(); i++){
        h_val->Fill(v_values[i]);
        h_err->Fill(v_errors[i]);
        h_pul->Fill(v_pulls[i]);     
      } 

      c->cd(1) ; gPad->SetLeftMargin(0.15); h_val->GetXaxis()->SetTitle(Form("%s Value",temp_var->GetTitle())); h_val->GetXaxis()->SetNdivisions(5); h_val->GetYaxis()->SetTitleOffset(1.4); h_val->Draw("PE0") ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); h_err->GetXaxis()->SetTitle(Form("%s Error",temp_var->GetTitle())); h_err->GetXaxis()->SetNdivisions(5); h_err->GetYaxis()->SetTitleOffset(1.4); h_err->Draw("PE0") ;
      c->cd(3) ; gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->GetXaxis()->SetNdivisions(5); h_pul->GetYaxis()->SetTitleOffset(1.4); h_pul->Draw("PE0") ;
      double max=h_pul->GetMaximum();
      h_pul->SetMaximum(max*1.3);

      TF1 *gaussian = new TF1("Gaussian","gaus",min_pul,max_pul);
      h_pul->Fit("Gaussian","Q0R");
      double mean  = gaussian->GetParameter(1);
      double sigma = gaussian->GetParameter(2);
      const double* err=gaussian->GetParErrors();
      double meanerr=err[1];
      double sigmaerr=err[2];
      c->cd(3) ; gaussian->SetLineColor(kBlue); gaussian->Draw("same");

      TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
      TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
      pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
      pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
      pav1->SetFillColor(0);    pav2->SetFillColor(0); 
      pav1->SetTextFont(12);    pav2->SetTextFont(12);  
      pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
      pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
      pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
      pav1->Draw();             pav2->Draw();            

      gStyle->SetTitleFontSize(0.1);
      if ((nv-1)%20==0) {
        std::cout << "Creating Canvas.. " << n_canvas  << std::endl;
        v_canvas.push_back( new TCanvas(Form("canvas_all_%d_pull",n_canvas),"Pulls",20,20,800,800)  );
        v_canvas[n_canvas]->Divide(5,4);
        n_canvas++;
      }

      std::cout << "Moving to Canvas section  " << nv - (n_canvas-1)*20 << std::endl;
      v_canvas[n_canvas-1]->cd( nv - (n_canvas-1)*20 );
      gPad->SetLeftMargin(0.15); h_pul->GetXaxis()->SetTitle(Form("%s Pulls",temp_var->GetTitle())); h_pul->Draw("PE0"); gaussian->Draw("same");pav1->Draw();pav2->Draw(); 
      gPad->Update();
      c->Print(Form("toysDir/plots/test_%s.pdf",c->GetName()));
    }
    for(unsigned int i = 0; i < v_canvas.size(); i++){
      v_canvas[i]->Print( Form("toysDir/plots/test_%s.pdf",v_canvas[i]->GetName()) );
    }
    std::cout << "Dodgy Pull list: " <<std::endl;
    for(std::map<int,int>::iterator itF=Dodgy_pulls.begin(); itF!=Dodgy_pulls.end(); itF++) {
      std::cout << "Pull Number: " << itF->first << " count: " << itF->second <<std::endl;   
    }

  }


  if(false){
    int nv=0;
    RooRealVar* var=0;
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    TCanvas* c_all = new TCanvas("canvas_all_pull","Pulls",0,0,800,800) ; 
    c_all->Divide(5,4) ;
    TCanvas* c_all2 = new TCanvas("canvas_all2_pull","Pulls",0,0,800,800) ; 
    c_all2->Divide(5,4) ;
    TCanvas* c_all3 = new TCanvas("canvas_all3_pull","Pulls",0,0,800,800) ; 
    c_all3->Divide(5,4) ;
    RooArgSet* vars = sim->getVariables();
    TIterator* it = vars->createIterator();
    while((var = (RooRealVar*)it->Next())) {
      if(var->isConstant() or var->InheritsFrom("RooAbsCategory") or 0==std::string(var->GetName()).find(mB.GetName())) {
        std::cout <<"Not Plotting variable: " << var->GetName()<<std::endl;
        continue;
      }
      nv++;
      std::cout<< std::endl<< nv <<". Processing variable: " << var->GetName()<<std::endl;
      //if(0==std::string(var->GetName()).find("R_d0kst_pass")) continue;
      //if((0==std::string(var->GetName()).find("n_")) and (0!=std::string(var->GetName()).find("n_b2dpi_"))) continue;
      //if((0!=std::string(var->GetName()).find("R_")) and 
      //   (0!=std::string(var->GetName()).find("A_")) and 
      //   (0!=std::string(var->GetName()).find("effPID")) and 
      //   (0!=std::string(var->GetName()).find("alphaL_dpi_")) and 
      //   (0!=std::string(var->GetName()).find("alphaR_dpi_")) and 
      //   (0!=std::string(var->GetName()).find("rb")) and 
      //   (0!=std::string(var->GetName()).find("rd")) and 
      //   (0!=std::string(var->GetName()).find("delta")) and 
      //   (0!=std::string(var->GetName()).find("gamma")) and 
      //   (0!=std::string(var->GetName()).find("MC_"))) continue;
      std::cout<<"  ------->  Getting PULLs " << std::endl;
      RooPlot* frame1 = mcstudy->plotParam(*var,RooFit::Bins(20));
      RooPlot* frame2 = mcstudy->plotError(*var,RooFit::Bins(20));
      //RooPlot* frame3 = mcstudy->plotPull( *var,RooFit::Bins(20),RooFit::Range(-5,5),RooFit::FitGauss(true));frame3->SetTitle(var->GetTitle());  
      //RooPlot* frame3 = mcstudy->plotPull( *var,RooFit::Bins(20),RooFit::Range(-5,5));frame3->SetTitle(var->GetTitle());
      RooPlot* frame3 = mcstudy->plotPull( *var, -5.0, 5.0, 20, kTRUE);
      //frame3.Draw()


      TCanvas* c = new TCanvas(Form("canv_%s",var->GetName()),var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;      
      std::cout<<"  ------->  Plotting PULLs " << std::endl;
      c->cd(1) ; gPad->SetLeftMargin(0.15); frame1->GetXaxis()->SetNdivisions(5); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw() ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); frame2->GetXaxis()->SetNdivisions(5); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw() ;
      c->cd(3) ; gPad->SetLeftMargin(0.15); frame3->GetXaxis()->SetNdivisions(5); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw() ;
      /*
      RooHist* hist = frame3->getHist();
      RooAbsData* hist_abs=(RooAbsData*) hist;
      RooRealVar pullMean("pullMean","Mean of pull",0,-5,5) ;
      RooRealVar pullSigma("pullSigma","Width of pull",1,0,5) ;
      RooGenericPdf pullGauss("pullGauss",
                              "Gaussian of pull",
                              "exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2))",
                              RooArgSet(*var,pullMean,pullSigma)  ) ;

      pullGauss.fitTo(*hist_abs) ;
      pullGauss.plotOn(frame3) ;
      //hist->Fit('gaus', 'LME');
      //RooGaussian* gaus = hist->GetFunction('gaus')    
      //gaus->Draw();
      */
      gStyle->SetTitleFontSize(0.1);
      if(nv<=20)        c_all->cd(nv);
      if(nv>20&&nv<=40) c_all2->cd(nv-20);
      if(nv>40)        c_all3->cd(nv-40);
      gPad->SetLeftMargin(0.15); frame3->Draw();
      ((TPaveText*)gPad->FindObject("pullGauss_paramBox"))->SetTextSize(0.04); //was 0.07
      gPad->Update();
      c->Print(Form("toysDir/plots/%s.pdf",c->GetName()));
    }

    TCanvas* c_nll = new TCanvas("canvas_nll","NLL",0,0,400,400) ;
    c_nll->cd(1)->SetLeftMargin(0.15); 
    mcstudy->plotNLL(RooFit::Bins(20))->Draw();  

    c_all->Print( Form("toysDir/plots/%s.pdf",c_all->GetName()) );
    c_all2->Print(Form("toysDir/plots/%s.pdf",c_all2->GetName()));
    c_all3->Print(Form("toysDir/plots/%s.pdf",c_all3->GetName()));
    c_nll->Print( Form("toysDir/plots/%s.pdf",c_nll->GetName()) );
    std::cout<< std::endl;
  }

  if(true){
    // Correlation Histograms
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetStatFont(132);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(132,"XYZ"); 
    gStyle->SetLabelFont(132,"XYZ");
    gStyle->SetTitleSize(0.08,"XYZ"); 
    gStyle->SetLabelSize(0.08,"XYZ");
    gStyle->SetLabelOffset(0.01,"XYZ");
    gStyle->SetPadTopMargin(0.01);
    gStyle->SetPadBottomMargin(0.2);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.01);
    
    int ipad=1;
    std::map<int,RooRealVar*> var;
    TCanvas* c = new TCanvas("canvas_correl_1","",0,0,1200,800);
    TCanvas* c2 = new TCanvas("canvas_correl_2","",0,0,1200,800);
    TCanvas* c3 = new TCanvas("canvas_correl_3","",0,0,1200,800);
    c->Divide(7,7); //first three vars => 45 plots
    c2->Divide(7,7); //next four vars => 46 plots
    c3->Divide(7,7); //remaining vars => 45 plots
    //Put variables here  
    /*
    var[0] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);          var[0]->setRange(1100, 1300);
    var[1] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);        var[1]->setRange(250, 400);
    var[2] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);         var[2]->setRange(70, 160);
    var[3] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);          var[3]->setRange(450, 950);
    var[4] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);        var[4]->setRange(250, 550);
    var[5] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);         var[5]->setRange(100, 350);
    var[6] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);      var[6]->setRange(1800, 2300);
    var[7] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);    var[7]->setRange(350, 700);
    var[8] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);     var[8]->setRange(100, 300);

    var[9] = new RooRealVar(*(RooRealVar*)model_gen->mean_B[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);              var[9]->setRange(5277, 5281);
    var[10] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);              var[10]->setRange(12, 16);
    var[11] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);            var[11]->setRange(11, 17);
    var[12] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);             var[12]->setRange(8, 20);
    var[13] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi]);                     var[13]->setRange(-0.0055, -0.003);
    var[14] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi]);                   var[14]->setRange(-0.006, -0.002);
    var[15] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi]);                    var[15]->setRange(-0.007, -0.0015);
    var[16] = new RooRealVar(*(RooRealVar*)model_gen->frac[toy][both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);               var[16]->setRange(0.4, 0.6);

    var[17] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KKPi][both][both]);        var[17]->setRange(1100, 1300);
    var[18] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2PiPiPi][both][both]);      var[18]->setRange(250, 400);
    var[19] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KPiPi][both][both]);       var[19]->setRange(70, 160);
    var[20] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KKPi][both][both]);        var[20]->setRange(100, 350);
    var[21] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KKPi][both][both]);    var[21]->setRange(1800, 2300);
    var[22] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2PiPiPi][both][both]);  var[22]->setRange(350, 700);
    var[23] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KPiPi][both][both]);   var[23]->setRange(100, 300);

    var[24] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KKPi][both][both]);             var[24]->setRange(12, 16);
    var[25] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2PiPiPi][both][both]);           var[25]->setRange(11, 17);
    var[26] = new RooRealVar(*(RooRealVar*)model_gen->sigma[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KPiPi][both][both]);            var[26]->setRange(8, 20);
    var[27] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[toy][both][DsBDTbin1][PhiBDTbin1][DsPhi][Ds2KKPi]);                    var[27]->setRange(-0.0055, -0.003);
    */
    RooArgSet* vars = sim->getVariables();
    std::vector<RooRealVar*> variables;
    RooRealVar* temp_var=0;
    TIterator* it2 = vars->createIterator();
    while((temp_var = (RooRealVar*)it2->Next())) {
      if(temp_var->isConstant() or temp_var->InheritsFrom("RooAbsCategory") or 0==std::string(temp_var->GetName()).find(mB.GetName())) {
        continue;
      }
      variables.push_back(temp_var);
    }

    std::map<int, std::map<int,TH1*> > hist;
    for(unsigned int i=0;i<variables.size()-1;i++){
      for(unsigned int j=i+1;j<variables.size();j++){
        hist[i][j] = mcstudy->fitParDataSet().createHistogram(Form("cor_%i_%i",i,j),*variables[i],RooFit::YVar(*variables[j])) ;
        //hist[i][j] = mcstudy->fitParDataSet().createHistogram(Form("cor_%i_%i",i,j),*var[i],RooFit::YVar(*var[j])) ;
        hist[i][j]->SetMarkerStyle(20);
        hist[i][j]->SetMarkerSize(0.02);
        hist[i][j]->SetTitle("");
        if(i < 3) { c->cd(ipad); 
        hist[i][j]->SetMarkerColor(kBlue-7+i);
        hist[i][j]->SetLineColor(kBlue-7+i); }
        else if(i < 7) { c2->cd(ipad-45); 
        hist[i][j]->SetMarkerColor(kBlue-7+i);
        hist[i][j]->SetLineColor(kBlue-7+i); }
        else { c3->cd(ipad-91); 
        hist[i][j]->SetMarkerColor(kBlue-11+i);
        hist[i][j]->SetLineColor(kBlue-11+i); }
        gPad->SetLeftMargin(0.25);
        gPad->SetBottomMargin(0.25);
        hist[i][j]->GetYaxis()->SetTitleOffset(1.4);
        hist[i][j]->GetXaxis()->SetTitleOffset(1.2);
        hist[i][j]->Draw("box");
        ipad++;
      }
    }
    
    c->Print(Form("toysDir/plots/%s.pdf",c->GetName()));
    c2->Print(Form("toysDir/plots/%s.pdf",c2->GetName()));
    c3->Print(Form("toysDir/plots/%s.pdf",c3->GetName()));

  } //end of if(true)


} //end of DisplayToys() funcn



void DKst0Fitting::DisplaySys()
{
  if(par->debug) std::cout<<"Running: DKst0Fitting::DisplaySys()"<<std::endl;
  //RooSimultaneous*  sim=model_gen->Pdf();
  std::vector<RooFitResult> Results; 
  int nFits=1;
  int nBad=0;
  int nGood=0;
  int nFPD=0;
  std::string lastfit_name = "lastRooFitResult";
  if (par->sumOverCharges) lastfit_name += "_summed";
  if (par->splitHel)       lastfit_name += "_splitHel";
  if(par->dsetsReq[s21])   lastfit_name += "_"+s21;
  if(par->dsetsReq[s21r1]) lastfit_name += "_"+s21r1;
  if(par->dsetsReq[s24])   lastfit_name += "_"+s24;
  if(par->dsetsReq[s26])   lastfit_name += "_"+s26;
  if(par->modes[D2PiKPi]) lastfit_name += "_"+D2PiKPi;
  if(par->modes[D2KKPi])  lastfit_name += "_"+D2KKPi;
  if(par->modes[D2PiPiPi])lastfit_name += "_"+D2PiPiPi;
  std::string sysfile("vary_");
  sysfile+=lastfit_name;  //need to check what the name should really be...
  //Form("%s/vary_%s_%s_SEED_%i_%i.root",par->variationoption.c_str(),lastfit_name.c_str(),par->toVary.c_str(),par->seed,int(sim->expectedEvents(*cat))

  std::cout <<" Looking for files like: " << par->variationoption+slash+sysfile <<"*"<< std::endl;
  std::vector<std::string> filename;
  CommonTools::getdir(par->variationoption.c_str(),filename);
  for(unsigned int i=0;i<filename.size();i++){
    if(filename[i].find(sysfile) == std::string::npos) continue;
    TFile f((par->variationoption+slash+filename[i]).c_str(),"READ");
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    while((key=(TKey*)nextkey())){
      if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
      RooFitResult *r = (RooFitResult*)key->ReadObj();
      if(r->covQual()<2){nBad++;continue;}
      if(r->covQual()<3){nFPD++;continue;}
      nGood++;
      std::cout <<"\r   Collecting "<<nFits<<" systematic fits, of which "<<nBad<<" don't converge and "<<100*float(nFPD)/nFits<<"% are forced positive definite. "<< nGood << " are good."<<std::flush;
      Results.push_back(*r);
      nFits++;
    }
  }
  std::cout<<"\n"<<std::endl;
  if(1==nFits){
    std::cout<<"SOMETHING WRONG. NO SYSTEMATICS LOADED"<<std::endl;
    return;
  }
  std::vector<std::string> sys_vars;

  sys_vars.push_back("yield_peak_DD0_D2PiKPi_both_both_both_both");
  sys_vars.push_back("yield_peak_DD0_D2KKPi_both_both_both_both");
  sys_vars.push_back("yield_peak_DD0_D2PiPiPi_both_both_both_both");

  
  if(par->fitBr){
    sys_vars.push_back("Branching_fraction");
  }else if(par->fitFourBr){
    sys_vars.push_back("Branching_fraction_D2PiKPi");
    sys_vars.push_back("Branching_fraction_D2KKPi");
    sys_vars.push_back("Branching_fraction_D2PiPiPi");
  } else{
    sys_vars.push_back("yield_peak_total_DKst0_D2PiKPi_both_both_both_both");
    sys_vars.push_back("yield_peak_total_DKst0_D2KKPi_both_both_both_both");
    sys_vars.push_back("yield_peak_total_DKst0_D2PiPiPi_both_both_both_both");

  }

  if(true){
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    std::vector<TCanvas*> v_canvas;
    //std::map<int,int> Dodgy_pulls;
    //RooArgSet* vars = sim->getVariables();
    std::vector<RooRealVar*> variables;
    RooRealVar* temp_var=0;
    int nv=0;
    for(unsigned int j = 0; j<sys_vars.size(); j++){
      nv++;
      std::vector<double> v_values, v_errors;
      std::cout << "Printing Variable: " << sys_vars[j] << std::endl;
      std::string var_title;
      std::string var_name;
      for(unsigned int i = 0; i< Results.size(); i++){
        RooFitResult rfr_temp = Results[i];
        temp_var = (RooRealVar*) (rfr_temp.floatParsFinal()).find(sys_vars[j].c_str());
        v_values.push_back(temp_var->getVal());
        v_errors.push_back(temp_var->getError());
        var_title = temp_var->GetTitle();
        var_name = temp_var->GetName();
      }

      double max_val = *max_element(v_values.begin(), v_values.end());
      double min_val = *min_element(v_values.begin(), v_values.end());
      double range_val = max_val - min_val;
      std::cout << "Min value: " << min_val<< std::endl;
      std::cout << "Min value: " << max_val<< std::endl;
      std::cout << "Range:     " << range_val<< std::endl;
      double max_err = *max_element(v_errors.begin(), v_errors.end());
      double min_err = *min_element(v_errors.begin(), v_errors.end());
      double range_error = max_err - min_err;
      std::cout << "Min error: " << min_err  << std::endl;
      std::cout << "Min error: " << max_err  << std::endl;
      std::cout << "Range:     " << range_error<< std::endl;

      SetLHCbStyle(norm);
      gStyle->SetOptStat(0);
      gStyle->SetMarkerSize(0.5);    
      TCanvas* c  = new TCanvas(Form("canv_%s",var_name.c_str()),var_name.c_str(),20*nv,20*nv,600,300); c->Divide(2) ;
      TH1F* h_val = new TH1F(Form("Title;%s Value;Number of fits;", var_title.c_str()), "Vars", 20, min_val-fabs(range_val)*0.25, max_val+fabs(range_val)*0.25);
      TH1F* h_err = new TH1F(Form("Title;%s Error; Number of fits;",var_title.c_str()), "Errs", 20, min_err-fabs(range_error)*0.25, max_err+fabs(range_error)*0.25);
      std::cout << "Min plot value: " << min_val-fabs(range_val)*0.25<< std::endl;
      std::cout << "Max plot value: " << max_val+fabs(range_val)*0.25<< std::endl;
      
      for(unsigned int i = 0; i< v_values.size(); i++){
        h_val->Fill(v_values[i]);
        h_err->Fill(v_errors[i]);     
      } 
      c->cd(1) ; gPad->SetLeftMargin(0.15); h_val->GetXaxis()->SetTitle(Form("%s Value",var_title.c_str())); h_val->GetXaxis()->SetNdivisions(5); h_val->GetYaxis()->SetTitleOffset(1.4); h_val->Draw("PE0") ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); h_err->GetXaxis()->SetTitle(Form("%s Error",var_title.c_str())); h_err->GetXaxis()->SetNdivisions(5); h_err->GetYaxis()->SetTitleOffset(1.4); h_err->Draw("PE0") ;
  
      TF1 *gaussian = new TF1("Gaussian","gaus",min_val-fabs(range_val)*0.25,max_val+fabs(range_val)*0.25);
      h_val->Fit("Gaussian","Q0R");
      double mean  = gaussian->GetParameter(1);
      double sigma = gaussian->GetParameter(2);
      const double* err=gaussian->GetParErrors();
      double meanerr=err[1];
      double sigmaerr=err[2];
      c->cd(1) ; gaussian->SetLineColor(kBlue); gaussian->Draw("same");

      TPaveLabel *pav1 = new TPaveLabel(0.2,0.89,0.93,0.94,Form("#bf{Mean: } %f #pm %f", mean,  meanerr),"NDC");
      TPaveLabel *pav2 = new TPaveLabel(0.2,0.81,0.93,0.86,Form("#bf{Sigma:} %f #pm %f",sigma, sigmaerr),"NDC");
      pav1->SetBorderSize(0);   pav2->SetBorderSize(0);
      pav1->SetFillStyle(1001); pav2->SetFillStyle(1001);
      pav1->SetFillColor(0);    pav2->SetFillColor(0); 
      pav1->SetTextFont(12);    pav2->SetTextFont(12);  
      pav1->SetTextSize(0.9);   pav2->SetTextSize(0.9);
      pav1->SetTextAlign(31);   pav2->SetTextAlign(31);
      pav1->SetTextColor(kRed); pav2->SetTextColor(kRed);
      pav1->Draw();             pav2->Draw();            
      c->Print(Form("%s/plots/%s.pdf",par->variationoption.c_str(),c->GetName()));
    }
  }
} //end of DisplaySys() funcn


void DKst0Fitting::Sensitivity(){

  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity()"<<std::endl;


  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Created Model"<<std::endl;
  RooAbsPdf* sim_withBr = model_gen->Pdf();
  RooAbsPdf* sim_noBr   = model->Pdf();
  RooAbsCategory* cat   = model_gen->Cat();
  RooArgSet* initialPars= model_gen->GetParameters();
  int secs=20110602;
  if(par->useSeed) secs = par->seed;
  std::cout <<  "====================== " << std::endl;
  std::cout <<  "Using Seed: " << secs << std::endl;
  std::cout <<  "====================== " << std::endl;

  RooRandom::randomGenerator()->SetSeed(secs);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Saving generated values to text files"<<std::endl;
  gSystem->Exec((std::string("mkdir -p ")+"sensitivityDir").c_str());
  gSystem->Exec((std::string("mkdir -p ")+"sensitivityDir/gen_vals").c_str());
  gSystem->Exec((std::string("mkdir -p ")+"sensitivityDir/text_vals").c_str());

  std::string gen_location="sensitivityDir/gen_vals";
  std::string text_location="sensitivityDir/text_vals";

  TIterator* initialIter= initialPars->createIterator();
  TObject* initialObj; 
  while((initialObj=initialIter->Next())) {
    RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
    if(!initialRrv) continue;
    std::string text_filename  = Form("%s/Init_Val_%s.txt", gen_location.c_str(),initialRrv->GetName());

    std::cout<<" "<<initialRrv->GetName()<<" <--- "<<initialRrv->getVal()<<" --> File: " << text_filename <<std::endl;
    
    std::ofstream output_file(text_filename.c_str());
    output_file << initialRrv->getVal();
    output_file.close();
  }
  delete initialIter;

  std::cout<<" "<<std::endl;


  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Creating RooMCStudy:"<<std::endl;
  RooMCStudy* mcstudy = new RooMCStudy(*sim_withBr ,reducedlist,RooFit::Binned(par->binned),RooFit::FitOptions(RooFit::Save(true),RooFit::PrintLevel(-1)));

  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Created RooMCStudy"<<std::endl;
  int nEvtsPerSample=sim_withBr->expectedEvents(*cat);
 
  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Expected Events: " << nEvtsPerSample <<std::endl; 
  // Fit using fit model 
  
  int nToys = par->sensitivityN;

  if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Generating " << nToys << " toys of "<<nEvtsPerSample<<" events."<<std::endl;
  mcstudy->generate(nToys,nEvtsPerSample,true);

  for(int i = 0; i<nToys;i++){

    std::string toyfile="sensitivityDir/";
    //for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){toyfile+=(*m)+underscore;}
    //for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){toyfile+=(*b)+underscore;}
    toyfile+=std::string(Form("sensitivity_%f_%i_%i_.root",par->sensitivityBR,secs,i));  
    TFile f(toyfile.c_str(),"RECREATE");
    
    
    // Reset fit initial values to the generated ones...
    //RooArgSet* initialPars= model_gen->GetParameters();
    //std::string gen_location=par->toylocation+"/gen_vals";
    TIterator* initialIter= initialPars->createIterator();
    TObject* initialObj; 
    while((initialObj=initialIter->Next())) {
      RooRealVar* initialRrv = dynamic_cast<RooRealVar*>(initialObj);
      if(!initialRrv) continue;
      std::string text_filename  = Form("%s/Init_Val_%s.txt", gen_location.c_str(),initialRrv->GetName());

      std::ifstream input;
      input.open(text_filename.c_str(),std::ifstream::in);
      if(!input){
        std::cout<<" Can't find file " << text_filename << std::endl; return;
      }
      std::string line = "";
      getline(input,line);
      initialRrv->setVal(atof(line.c_str()));
      input.close();
      }
    delete initialIter;
    

    if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Fitting toy "<< i <<" with signal..."<<std::endl;
    RooDataSet* mcdata=(RooDataSet*)mcstudy->genData(i);
    RooFitResult* mcresult_withSig = sim_withBr->fitTo(*mcdata,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(-1),RooFit::Optimize(false),RooFit::Timer(true));
    mcresult_withSig->Write("toy_withSig");
    if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> With Signal: " << mcresult_withSig->minNll()<<std::endl;
    
    if(false){
      if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Fitting toy "<< i <<" without signal..."<<std::endl;   
      RooFitResult* mcresult_noSig = sim_noBr->fitTo(*mcdata,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(-1),RooFit::Optimize(false),RooFit::Timer(true));
      mcresult_noSig->Write("toy_noSig");  
      if(par->debug) std::cout<<"Running: DKst0Fitting::Sensitivity() ----> Without Signal: " << mcresult_noSig->minNll()<<std::endl;

      std::string minLL_filename  = Form("%s/minLL_%f_%i_%i.txt",text_location.c_str(),par->sensitivityBR,secs,i);    
      std::ofstream minLL_file(minLL_filename.c_str());
      minLL_file << mcresult_withSig->minNll() << ":" << mcresult_noSig->minNll();
      minLL_file.close();
    }
    gDirectory->ls();
    f.Close();  

  }

}

void DKst0Fitting::SetLHCbStyle(std::string type = ""){

  if(type!=norm && type!= cont && type!= surf)
  {
    std::cout <<  "Type must be one of 1D, cont or surf. Defaulting to 1D" <<std::endl;
    type = norm;
  }

  // Copy and paste from lhcbstyle.C:

  // Use times new roman, precision 2 
  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
  // Line thickness
  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
  // Text size
  Double_t lhcbTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *lhcbStyle= new TStyle("lhcbStyle","LHCb plots style");
  
  //lhcbStyle->SetErrorX(0); //  don't suppress the error bar along X

  lhcbStyle->SetFillColor(1);
  lhcbStyle->SetFillStyle(1001);   // solid
  lhcbStyle->SetFrameFillColor(0);
  lhcbStyle->SetFrameBorderMode(0);
  lhcbStyle->SetPadBorderMode(0);
  lhcbStyle->SetPadColor(0);
  lhcbStyle->SetCanvasBorderMode(0);
  lhcbStyle->SetCanvasColor(0);
  lhcbStyle->SetStatColor(0);
  lhcbStyle->SetLegendBorderSize(0);

  // If you want the usual gradient palette (blue -> red)
  lhcbStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  lhcbStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  lhcbStyle->SetPaperSize(20,26);
  lhcbStyle->SetPadTopMargin(0.05);
  lhcbStyle->SetPadRightMargin(0.15); // increase for colz plots
  lhcbStyle->SetPadBottomMargin(0.16);
  lhcbStyle->SetPadLeftMargin(0.14);
  
  // use large fonts
  lhcbStyle->SetTextFont(lhcbFont);
  lhcbStyle->SetTextSize(lhcbTSize);
  lhcbStyle->SetLabelFont(lhcbFont,"x");
  lhcbStyle->SetLabelFont(lhcbFont,"y");
  lhcbStyle->SetLabelFont(lhcbFont,"z");
  lhcbStyle->SetLabelSize(lhcbTSize,"x");
  lhcbStyle->SetLabelSize(lhcbTSize,"y");
  lhcbStyle->SetLabelSize(lhcbTSize,"z");
  lhcbStyle->SetTitleFont(lhcbFont);
  lhcbStyle->SetTitleFont(lhcbFont,"x");
  lhcbStyle->SetTitleFont(lhcbFont,"y");
  lhcbStyle->SetTitleFont(lhcbFont,"z");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"x");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"y");
  lhcbStyle->SetTitleSize(1.2*lhcbTSize,"z");

  // use medium bold lines and thick markers
  lhcbStyle->SetLineWidth(lhcbWidth);
  lhcbStyle->SetFrameLineWidth(lhcbWidth);
  lhcbStyle->SetHistLineWidth(lhcbWidth);
  lhcbStyle->SetFuncWidth(lhcbWidth);
  lhcbStyle->SetGridWidth(lhcbWidth);
  lhcbStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  lhcbStyle->SetMarkerStyle(20);
  lhcbStyle->SetMarkerSize(1.0);

  // label offsets
  lhcbStyle->SetLabelOffset(0.010,"X");
  lhcbStyle->SetLabelOffset(0.010,"Y");

  // by default, do not display histogram decorations:
  //lhcbStyle->SetOptStat(0);  
  lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  lhcbStyle->SetStatFormat("6.3g"); // specified as c printf options
  lhcbStyle->SetOptTitle(0);
  lhcbStyle->SetOptFit(0);
  //lhcbStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  lhcbStyle->SetTitleOffset(0.95,"X");
  lhcbStyle->SetTitleOffset(0.95,"Y");
  lhcbStyle->SetTitleOffset(1.2,"Z");
  lhcbStyle->SetTitleFillColor(0);
  lhcbStyle->SetTitleStyle(0);
  lhcbStyle->SetTitleBorderSize(0);
  lhcbStyle->SetTitleFont(lhcbFont,"title");
  lhcbStyle->SetTitleX(0.0);
  lhcbStyle->SetTitleY(1.0); 
  lhcbStyle->SetTitleW(1.0);
  lhcbStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  lhcbStyle->SetStatBorderSize(0);
  lhcbStyle->SetStatFont(lhcbFont);
  lhcbStyle->SetStatFontSize(0.05);
  lhcbStyle->SetStatX(0.9);
  lhcbStyle->SetStatY(0.9);
  lhcbStyle->SetStatW(0.25);
  lhcbStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  lhcbStyle->SetPadTickX(1);
  lhcbStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  lhcbStyle->SetNdivisions(505,"x");
  lhcbStyle->SetNdivisions(510,"y");
 
   // Extra options
  if(type == cont||type==surf) lhcbStyle->SetPalette(112);
  if(type == surf) lhcbStyle->SetTitleOffset(1.5,"Y");
  if(type == surf) lhcbStyle->SetTitleOffset(1.5,"X");
  if(type == surf) lhcbStyle->SetTitleOffset(0.83,"Z");
  if(type == norm) lhcbStyle->SetPadRightMargin(0.05);
  

  gROOT->SetStyle("lhcbStyle");
  gROOT->ForceStyle();

  /*
  // add LHCb label
  lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                           0.87 - gStyle->GetPadTopMargin(),
                           gStyle->GetPadLeftMargin() + 0.20,
                           0.95 - gStyle->GetPadTopMargin(),
                           "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
  */
  TText *lhcbLabel = new TText();
  lhcbLabel->SetTextFont(lhcbFont);
  lhcbLabel->SetTextColor(1);
  lhcbLabel->SetTextSize(lhcbTSize);
  lhcbLabel->SetTextAlign(12);

  TLatex *lhcbLatex = new TLatex();
  lhcbLatex->SetTextFont(lhcbFont);
  lhcbLatex->SetTextColor(1);
  lhcbLatex->SetTextSize(lhcbTSize);
  lhcbLatex->SetTextAlign(12);

  std::cout << "-------------------------" << std::endl;  
  std::cout << "Set LHCb Style - Feb 2012" << std::endl;  
  std::cout << "Mode : " << type           << std::endl;
  std::cout << "-------------------------" << std::endl;  


}
