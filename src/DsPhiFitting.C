#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <math.h>
#include <time.h>
#include <fstream>

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

#include "DsPhiFitting.h"
#include "Parameters.h"
#include "CommonTools.h"
//#include "MultipleCandidates.h"
std::string cutStringBDT("D_BDTG>0.0 && D0_BDTG>0.0");

//**********************
//***----BLINDING----***
//**********************
//Bool_t isBlind2(kTRUE); //-> just using state==BLIND

DsPhiFitting::DsPhiFitting(Parameters* p, TApplication* app)
: Base()
, toys("toys")
, inputlist("contents of Final ntuple")
, fulllist("contents, including categories")
, storelist("contents of roodataset")
, reducedlist("reduced content")
  //, eventNumber("eventNumber","eventNumber",0,1e99)
, runNumber("runNumber","runNumber",0,1e99)
//, mB("Bu_D0constDconstPVconst_Bu_M","B mass with constraints", 5050 , 5900, "MeV/#font[12]{c}^{2}")
, mB("Bu_DTF_M","B mass with constraints", 5050 , 5900, "MeV/#font[12]{c}^{2}")
, mD0("D0_M","m(D)",1840, 1892,"MeV/c^{2}")
, mDs("D_M","m(Ds)",1944.3, 1994.3,"MeV/c^{2}")
, mPhi("Phi_M","m(Phi)",1000, 1040,"MeV/c^{2}")
  //, bach_dll("Bach_PIDK","#Delta(LL)",-200,200.)
//, bdt("D_BDT","MVA",0.4,0.9)
, bdtDs( "D_BDT","MVA",-1.1,1.1)
, bdtgDs("D_BDTG","MVA",-1.1,1.1)
, bdtbDs("D_BDTB","MVA",-1.1,1.1)
, bdtD0( "D0_BDT","MVA",-1.1,1.1)
, bdtgD0("D0_BDTG","MVA",-1.1,1.1)
, bdtbD0("D0_BDTB","MVA",-1.1,1.1)
, bdtPhi( "Phi_BDT","MVA",-1.1,1.1)
, bdtgPhi("Phi_BDTG","MVA",-1.1,1.1)
, bdtbPhi("Phi_BDTB","MVA",-1.1,1.1)
, bdtMC( "Bu_MC_BDT","MVA",-1.1,1.1)
, bdtgMC("Bu_MC_BDTG","MVA",-1.1,1.1)
, bdtbMC("Bu_MC_BDTB","MVA",-1.1,1.1)
// PID variables
, D0_K0_pidk("D0_K0_PIDK","PIDK",-100,100)
, D0_K1_pidk("D0_K1_PIDK","PIDK",-100,100)
, D_K_pidk("D_K_PIDK","PIDK",-100,100)
, D_K0_pidk("D_K0_PIDK","PIDK",-100,100)
, D_K1_pidk("D_K1_PIDK","PIDK",-100,100)
, D_P_pidk("D_P_PIDK","PIDK",-100,100)
, D_P0_pidk("D_P0_PIDK","PIDK",-100,100)
, D_P1_pidk("D_P1_PIDK","PIDK",-100,100)
, D_P2_pidk("D_P2_PIDK","PIDK",-100,100)
, helicityD0("D0_LoKi_LV01","Angle",-1.1,1.1)
, helicityPhi("Phi_LoKi_LV01","Angle",-1.1,1.1)
, mKK("D_KK_M","Ds_PhiPi_Mass",0.0,2000.0)
, bid("Bu_ID","PDG_ID of B cand.",-pow(2,32),pow(2,32))
, L0Hadron_TOS("Bu_L0HadronDecision_TOS","",-0.5,1.5)
, L0_TIS("Bu_L0Global_TIS","",-0.5,1.5)
, Hlt2IncPhi_TOS("Bu_Hlt2IncPhiDecision_TOS","",-0.5,1.5)
 // 2011/2012 Trigger Lines
, Hlt1TrackAllL0_TOS("Bu_Hlt1TrackAllL0Decision_TOS","",-0.5,1.5) //this is all == 1
, Hlt2Topo2BodyBBDT_TOS("Bu_Hlt2Topo2BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo3BodyBBDT_TOS("Bu_Hlt2Topo3BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo4BodyBBDT_TOS("Bu_Hlt2Topo4BodyBBDTDecision_TOS","",-0.5,1.5)
  // 2015 Trigger Lines
, Hlt1TrackMVA_TOS("Bu_Hlt1TrackMVADecision_TOS","",-0.5,1.5)
, Hlt1TwoTrackMVA_TOS("Bu_Hlt1TwoTrackMVADecision_TOS","",-0.5,1.5)
, Hlt2Topo2Body_TOS("Bu_Hlt2Topo2BodyDecision_TOS","",-0.5,1.5)
, Hlt2Topo3Body_TOS("Bu_Hlt2Topo3BodyDecision_TOS","",-0.5,1.5)
, Hlt2Topo4Body_TOS("Bu_Hlt2Topo4BodyDecision_TOS","",-0.5,1.5)
, KKPi_D_Veto("D_Veto","",-0.5,1.5)
, KKPi_Lc_Veto("Lc_Veto","",-0.5,1.5)
, KPiPi_mPiKPi("D_PiKPi_M","",0.0,100000)
, KPiPi_mPiPiPi("D_PiPiPi_M","",0.0,100000)
, KPiPi_mKKPi("D_KKPi_M","",0.0,100000)
, type("type","Type of data")
, mode("mode","D_{s} decay mode")
, Bmode("Bmode","B^{+} decay mode")
, magnet("magnet","dipole polarity")
, charge("charge","batchelor charge")
, bMassRegion("bMassRegion","m(B) band")
, helBin("helBin","Helicity Angle bin")
, DsBDTBin("DsBDTBin","Ds BDT bin")
, PhiBDTBin("PhiBDTBin","Phi BDT bin")
, m_mD0_Mu(1864.84)
, m_mDs_Mu(1968.30)
, m_mPhi_Mu(1019.461)
, m_bdtCut(-0.9)
, state(1)
, typeList()
, modeList()
, BmodeList()
, chargeList()
, HelBinList()
, DsBDTBinList()
, PhiBDTBinList()
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
      if (!par->manyFits) {
        RunFullFit(!par->batch); 
      }else {
        RunManyFits();
      }
    }
  }
  if(!par->batch){
    std::cout<<"Starting: app->Run()"<<std::endl;
    app->Run(true);
  }
}


void DsPhiFitting::DefineRooCategories()
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::DefineRooCategories()"<<std::endl;
  if(par->genToys)      typeList.push_back(toy);
  if(par->dsetsReq[s21])   typeList.push_back(s21);
  if(par->dsetsReq[s21r1]) typeList.push_back(s21r1);
  if(par->dsetsReq[s24]) typeList.push_back(s24);
  if(0==typeList.size()) typeList.push_back("UNDEFTYPE");

  type.defineType(toy.c_str());
  type.defineType(s21.c_str());
  type.defineType(s21r1.c_str());
  type.defineType(s24.c_str());

  mode.defineType(Ds2PhiPi.c_str());
  mode.defineType(Ds2KKPi.c_str());
  mode.defineType(Ds2PiPiPi.c_str());
  mode.defineType(Ds2KPiPi.c_str());
  
  Bmode.defineType(DsD0.c_str());
  Bmode.defineType(DsPhi.c_str());

  if(par->modes[Ds2PhiPi])     modeList.push_back(Ds2PhiPi);
  if(par->modes[Ds2KKPi])      modeList.push_back(Ds2KKPi);
  if(par->modes[Ds2PiPiPi])    modeList.push_back(Ds2PiPiPi);
  if(par->modes[Ds2KPiPi])     modeList.push_back(Ds2KPiPi);

  if(par->Bmodes[DsD0])        BmodeList.push_back(DsD0);
  if(par->Bmodes[DsPhi])       BmodeList.push_back(DsPhi);
 

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

  DsBDTBinList.push_back(DsBDTbin1);
  if(par->nDsBDTBins>1) DsBDTBinList.push_back(DsBDTbin2);

  PhiBDTBinList.push_back(PhiBDTbin1);
  if(par->nPhiBDTBins>1) PhiBDTBinList.push_back(PhiBDTbin2);
  
  DsBDTBin.defineType(DsBDTbin1.c_str());
  DsBDTBin.defineType(DsBDTbin2.c_str());

  PhiBDTBin.defineType(PhiBDTbin1.c_str());
  PhiBDTBin.defineType(PhiBDTbin2.c_str());

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

  Bmod[DsD0]="D^{0}";   Bmod[DsPhi]="#phi";
  mod[Ds2PhiPi]="#phi#pi"; mod[Ds2KKPi]="KK#pi"; mod[Ds2KPiPi]="K#pi#pi";    mod[Ds2PiPiPi]="#pi#pi#pi";
  chg[plus]="+";        chg[minus]="#font[122]{-}"; chg[both]="#font[122]{#pm}";
  mag[up]="UP";         mag[dn]="DN";               mag[both]="";
  hel[both] = "All cos(#theta)"; hel[Helbin1] = "|cos(#theta)|>0.4"; hel[Helbin2] = "|cos(#theta)|<0.4";
  dsbdt[DsBDTbin1]  = "BDT_{D_s} > 0.9";  dsbdt[DsBDTbin2] =  " 0.8 < BDT_{D_s} < 0.9"; 
  phibdt[DsBDTbin1] = "BDT_{#phi} > 0.9"; phibdt[DsBDTbin2] = " 0.8 < BDT_{#phi} < 0.9";
  if(par->nDsBDTBins ==1)  dsbdt[DsBDTbin1]   = "";
  if(par->nPhiBDTBins ==1) phibdt[PhiBDTbin1] = "";

  for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
    for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
      for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
        for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
          for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
            for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
              for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
                title[*h][*ds][*ph][*b][*m][*c][*a] = std::string(Form("B^{%s}#rightarrow D_{s}^{%s} #font[132]{(}#rightarrow %s#font[132]{)} %s #font[132]{(}#rightarrow KK#font[132]{)} %s", chg[*c].c_str(), chg[*c].c_str(),mod[*m].c_str(), Bmod[*b].c_str(), mag[*a].c_str()));
                bin_detail[*h][*ds][*ph][*b][*m][*c][*a] = std::string(Form( "%s %s %s",hel[*h].c_str(),dsbdt[*ds].c_str(),phibdt[*ph].c_str()));
                if(par->genToys) title[*h][*ds][*ph][*b][*m][*c][*a]+=" TOY"; 
              }
            }
          }
        }
      }
    }
  }
  //inputlist.add(eventNumber);
  //inputlist.add(bdt);
  inputlist.add(bid);
  inputlist.add(mDs);
  inputlist.add(mD0);
  inputlist.add(mPhi);
  inputlist.add(bdtD0);
  inputlist.add(bdtgD0);
  inputlist.add(bdtbD0);
  inputlist.add(bdtgDs);
  inputlist.add(bdtbDs);
  inputlist.add(bdtDs);
  inputlist.add(bdtgPhi);
  inputlist.add(bdtbPhi);
  inputlist.add(bdtPhi);
  inputlist.add(bdtgMC);
  inputlist.add(bdtbMC);
  inputlist.add(bdtMC);
  inputlist.add(helicityD0);
  inputlist.add(helicityPhi);
  inputlist.add(mKK);

  //inputlist.add(bach_dll);
  inputlist.add(L0Hadron_TOS);
  inputlist.add(L0_TIS);
  inputlist.add(Hlt2IncPhi_TOS);
  
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
  storelist.add(DsBDTBin);
  storelist.add(PhiBDTBin);

  fulllist.add(storelist);
  fulllist.add(inputlist);
  fulllist.add(bMassRegion);

  reducedlist.add(mB);
  inputlist.add(mB);

  if(par->manyFits){
    reducedlist.add(bdtD0);
    reducedlist.add(bdtgD0);
    reducedlist.add(bdtbD0);
    reducedlist.add(bdtgDs);
    reducedlist.add(bdtbDs);
    reducedlist.add(bdtDs);
    reducedlist.add(bdtgMC);
    reducedlist.add(bdtbMC);
    reducedlist.add(bdtMC);


    reducedlist.add(D0_K0_pidk);
    reducedlist.add(D0_K1_pidk);

    if(par->modes[Ds2KKPi]) {
      reducedlist.add(D_K0_pidk);
      reducedlist.add(D_K1_pidk);
      reducedlist.add(D_P_pidk);
    } else if(par->modes[Ds2KPiPi]){
      reducedlist.add(D_K_pidk);
      reducedlist.add(D_P0_pidk);
      reducedlist.add(D_P1_pidk);
    } else if(par->modes[Ds2KPiPi]){
      reducedlist.add(D_P0_pidk);
      reducedlist.add(D_P1_pidk);
      reducedlist.add(D_P2_pidk);
    }

    storelist.add(bdtD0);
    storelist.add(bdtgD0);
    storelist.add(bdtbD0);
    storelist.add(bdtgDs);
    storelist.add(bdtbDs);
    storelist.add(bdtDs);
    storelist.add(bdtgMC);
    storelist.add(bdtbMC);
    storelist.add(bdtMC);
    storelist.add(D0_K0_pidk);
    storelist.add(D0_K1_pidk);

    if(par->modes[Ds2KKPi]){
      storelist.add(D_K0_pidk);
      storelist.add(D_K1_pidk);
      storelist.add(D_P_pidk);
    } else if(par->modes[Ds2KPiPi]){
      storelist.add(D_K_pidk);
      storelist.add(D_P0_pidk);
      storelist.add(D_P1_pidk); 
    } else if(par->modes[Ds2KPiPi]){
      storelist.add(D_P0_pidk);
      storelist.add(D_P1_pidk);
      storelist.add(D_P2_pidk);
    }

  }
} //end of funcn DefineRooCategories



std::string DsPhiFitting::reducedFileName(std::string t,std::string m,std::string Bm)
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



int DsPhiFitting::LoadDataSet()
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



RooDataSet* DsPhiFitting::MakeDataSet(std::string t,std::string m,std::string b)
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
      //if(par->debug) std::cout<<"The fullpathandname is: "<<fullpathandname<<"."<<std::endl; //for debugging
      bool accept=false;
      if( (*it).find("D0D")   !=std::string::npos ) { s_Bmode=DsD0; }    
      if( (*it).find("PhiD")  !=std::string::npos ) { s_Bmode=DsPhi;}

      if( s_Bmode==b && ((m == Ds2PhiPi && (*it).find((Ds2KKPi+underscore).c_str())!=std::string::npos ) 
                          || (*it).find((m+underscore).c_str())!=std::string::npos) ) {
      //if(  (*it).find((m+underscore).c_str())!=std::string::npos && s_Bmode==b){ 
        accept=true;
        s_mode=m;
      }
      
      //Laurence
      if( (*it).find("2011")  !=std::string::npos ) { s_type=s21r1;}
      if( (*it).find("2012")  !=std::string::npos ) { s_type=s21; }
      if( (*it).find("2015")  !=std::string::npos ) { s_type=s24; }



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
          if (s_Bmode==DsD0){
            if( s_mode==Ds2PhiPi )  { tfileFolderName = "B2D0D_D02KK_D2KKPi"; } //NB it's D2KKPi, not Ds2KKPi
            if( s_mode==Ds2KKPi )   { tfileFolderName = "B2D0D_D02KK_D2KKPi"; } //NB it's D2KKPi, not Ds2KKPi
            if( s_mode==Ds2PiPiPi ) { tfileFolderName = "B2D0D_D02KK_D2PiPiPi"; }
            if( s_mode==Ds2KPiPi )  { tfileFolderName = "B2D0D_D02KK_D2KPiPi"; }
          } else if(s_Bmode==DsPhi){
            if( s_mode==Ds2PhiPi )  { tfileFolderName = "B2PhiD_Phi2KK_D2KKPi"; } //NB it's D2KKPi, not Ds2KKPi
            if( s_mode==Ds2KKPi )   { tfileFolderName = "B2PhiD_Phi2KK_D2KKPi"; } //NB it's D2KKPi, not Ds2KKPi
            if( s_mode==Ds2PiPiPi ) { tfileFolderName = "B2PhiD_Phi2KK_D2PiPiPi"; }
            if( s_mode==Ds2KPiPi )  { tfileFolderName = "B2PhiD_Phi2KK_D2KPiPi"; }
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




RooDataSet* DsPhiFitting::FinalDataSet(const std::string s_Bmode, const std::string s_type, const std::string s_mode, std::string s_mag, TTree* tree)
{
  if(par->debug) std::cout<<"Running: FinalDataSet()"<<std::endl;
  if(!tree){ std::cout << "\n THE TREE IS A ZERO POINTER IN "<< s_Bmode <<" "<< s_type <<" "<<s_mode<<" ("<< s_mag <<")"<< std::endl; return 0; }
  RooDataSet *input = new RooDataSet(Form("%s_%s_%s",s_Bmode.c_str(),s_type.c_str(),s_mode.c_str()),Form("Final %s %s %s",s_Bmode.c_str(),s_type.c_str(),s_mode.c_str()),inputlist,RooFit::Import(*tree));
  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(type,helBin,DsBDTBin,PhiBDTBin,Bmode,mode,magnet,charge,bMassRegion)); //removed "pid" from the RooArgSet
 

  int nTot=0;
  int nNotL0=0;
  int nNotL0TOS=0;
  int nNotTopoTOS=0;
  int nNotTrackTOS=0;
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
    if(s_Bmode==DsD0){
      RooRealVar* angle = (RooRealVar*) row->find(helicityD0.GetName());
      if(fabs(angle->getVal()) > 0.4 ) {
        helBin.setLabel(Helbin1.c_str());
      }else {
        helBin.setLabel(Helbin2.c_str());
      }
    }else{
      RooRealVar* angle = (RooRealVar*) row->find(helicityPhi.GetName());
      if(fabs(angle->getVal()) > 0.4 ) {
        helBin.setLabel(Helbin1.c_str());
      }else {
        helBin.setLabel(Helbin2.c_str());
      }
    }
  
  
    // Split into BDT bins
    if(par->nDsBDTBins==1){
      DsBDTBin.setLabel(DsBDTbin1.c_str());
    } else {
      RooRealVar* bd  = (RooRealVar*) row->find(bdtgDs.GetName());
      if(bd->getVal()> 0.65) DsBDTBin.setLabel(DsBDTbin1.c_str());
      if(bd->getVal()< 0.65) DsBDTBin.setLabel(DsBDTbin2.c_str());
    }
    if(par->nPhiBDTBins==1){
      PhiBDTBin.setLabel(PhiBDTbin1.c_str());
    }else{
      if(s_Bmode==DsD0){
        RooRealVar* bdd0 = (RooRealVar*) row->find(bdtgD0.GetName());
        if(bdd0->getVal()> 0.90) PhiBDTBin.setLabel(PhiBDTbin1.c_str());
        if(bdd0->getVal()< 0.90) PhiBDTBin.setLabel(PhiBDTbin2.c_str());
      } else {       
        RooRealVar* bdphi = (RooRealVar*) row->find(bdtgPhi.GetName());
        if(bdphi->getVal()> 0.90) PhiBDTBin.setLabel(PhiBDTbin1.c_str());
        if(bdphi->getVal()< 0.90) PhiBDTBin.setLabel(PhiBDTbin2.c_str());
      }
    }


    bMassRegion.setLabel("signal");

    if(par->modes[Ds2PhiPi]){
      if(s_mode==Ds2PhiPi && signal==std::string(bMassRegion.getLabel())){
          RooRealVar* d   = (RooRealVar*) row->find(mKK.GetName());
          if(fabs(d->getVal()-m_mPhi_Mu)>10) bMassRegion.setLabel("bckgrd");
      } else if(s_mode==Ds2KKPi && signal==std::string(bMassRegion.getLabel()) ) {
          RooRealVar* d   = (RooRealVar*) row->find(mKK.GetName());
          if(fabs(d->getVal()-m_mPhi_Mu)<=10) bMassRegion.setLabel("bckgrd");
      }
    } 

    if(s_Bmode==DsD0){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* d   = (RooRealVar*) row->find(mD0.GetName());
        if(fabs(d->getVal()-m_mD0_Mu)>25) bMassRegion.setLabel("bckgrd");
      }
    } else if(s_Bmode==DsPhi){
      if(signal==std::string(bMassRegion.getLabel())){
        RooRealVar* phim = (RooRealVar*) row->find(mPhi.GetName());
        if(fabs(phim->getVal()-m_mPhi_Mu)>25) bMassRegion.setLabel("bckgrd");
      }

    }


    // Remove cross feed candidtes 
    if(s_mode==Ds2KKPi){
      RooRealVar* dveto  = (RooRealVar*) row->find(KKPi_D_Veto.GetName());
      RooRealVar* lcveto = (RooRealVar*) row->find(KKPi_Lc_Veto.GetName());
      if (dveto->getVal() == 1 ) {bMassRegion.setLabel("bckgrd");}
      if (lcveto->getVal() == 1 ) {bMassRegion.setLabel("bckgrd");}
    }


    if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* bd  = (RooRealVar*) row->find(bdtgDs.GetName());
      // Main selection [BDT] cuts. Values chosen just by looking at D_BDT distribution, and removing a few events in the tail.
      if(s_mode==Ds2KKPi){
        if(s_type==s21){   if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }//+.40
        if(s_type==s21r1){ if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }
        if(s_type==s24){  if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }
      }
      if(s_mode==Ds2PiPiPi){
        if(s_type==s21){   if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} } //+.30
        if(s_type==s21r1){ if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} } 
        if(s_type==s24){   if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }
      }
      if(s_mode==Ds2KPiPi){
        if(s_type==s21){   if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} } //+.54
        if(s_type==s21r1){ if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }
        if(s_type==s24){   if(bd->getVal()<-11.0){bMassRegion.setLabel("bckgrd");} }
      }
    }
   
     if((s_type==s21||s_type==s21r1) && signal==std::string(bMassRegion.getLabel())){
      RooRealVar* _L0Global_TIS           = (RooRealVar*) row->find(L0_TIS.GetName());
      RooRealVar* _L0Hadron_TOS           = (RooRealVar*) row->find(L0Hadron_TOS.GetName());
      RooRealVar* _Hlt1TrackAllL0_TOS     = (RooRealVar*) row->find(Hlt1TrackAllL0_TOS.GetName());
      RooRealVar* _Hlt2Topo2BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo2BodyBBDT_TOS.GetName());
      RooRealVar* _Hlt2Topo3BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo3BodyBBDT_TOS.GetName());
      RooRealVar* _Hlt2Topo4BodyBBDT_TOS  = (RooRealVar*) row->find(Hlt2Topo4BodyBBDT_TOS.GetName());
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
      if(!(_Hlt2Topo2BodyBBDT_TOS->getVal()||_Hlt2Topo3BodyBBDT_TOS->getVal()||_Hlt2Topo4BodyBBDT_TOS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotTopoTOS++;
      }
    }

    if(s_type==s24 && signal==std::string(bMassRegion.getLabel())){
      RooRealVar* _L0Global_TIS        = (RooRealVar*) row->find(L0_TIS.GetName());
      RooRealVar* _L0Hadron_TOS        = (RooRealVar*) row->find(L0Hadron_TOS.GetName());
      RooRealVar* _Hlt1TrackMVA_TOS    = (RooRealVar*) row->find(Hlt1TrackMVA_TOS.GetName());
      RooRealVar* _Hlt1TwoTrackMVA_TOS = (RooRealVar*) row->find(Hlt1TwoTrackMVA_TOS.GetName());
      RooRealVar* _Hlt2Topo2Body_TOS   = (RooRealVar*) row->find(Hlt2Topo2Body_TOS.GetName());
      RooRealVar* _Hlt2Topo3Body_TOS   = (RooRealVar*) row->find(Hlt2Topo3Body_TOS.GetName());
      RooRealVar* _Hlt2Topo4Body_TOS   = (RooRealVar*) row->find(Hlt2Topo4Body_TOS.GetName());
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
      if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotTopoTOS++;
      }
    }

    extra->add(RooArgSet(type,helBin,DsBDTBin,PhiBDTBin,Bmode,mode,magnet,charge,bMassRegion));
  } //end of loop over entries
 
  std::cout <<"   Trigger: notL0TISTOS: "<<100*nNotL0/float(nTot)<<"% (TIS/L0Acc="<<100*nNotL0TOS/float(nTot-nNotL0)<<"%. Not TrackTOS: "<<100*nNotTrackTOS/float(nTot)<<"%. Not TopoTOS "<<100*nNotTopoTOS/float(nTot)<<"%."<<std::endl;
  input->merge(extra);
  delete extra;
  return input;
} //end of funcn FinalDataSet




void DsPhiFitting::AssignGlobalCategory()
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
    RooCategory* _DsBDTBin  = (RooCategory*) row->find("DsBDTBin");
    RooCategory* _PhiBDTBin = (RooCategory*) row->find("PhiBDTBin");
    if(par->polarity==both) _magnet->setLabel(both.c_str());
    if(par->sumOverCharges) _charge->setLabel(both.c_str());
    if( not par->splitHel)  _helBin->setLabel(both.c_str());

    model->Cat()->setLabel(Form("%s_%s_%s_%s_%s_%s_%s",_helBin->getLabel(),_DsBDTBin->getLabel(),_PhiBDTBin->getLabel(),_Bmode->getLabel(),_mode->getLabel(),_charge->getLabel(),_magnet->getLabel()));
    if(par->debug&&0==i%1000) std::cout << "Event " << i <<": magnet "<<_magnet->getLabel()<<", "<<_helBin->getLabel()<<", "<<_DsBDTBin->getLabel()<<", "<<_PhiBDTBin->getLabel()<<", "<<_Bmode->getLabel()<<", "<<_mode->getLabel()<<": "<<_charge->getLabel()<<" "<<model->Cat()->getLabel()<<std::endl;
    extra->add(RooArgSet(*model->Cat()));
  }
  data->merge(extra);

  RooDataSet* reducedData = new RooDataSet("reducedData","",data,reducedlist);
  delete data; data=reducedData;

}



void DsPhiFitting::PrintDataSet()
{
  if(par->debug) std::cout<<"Running: PrintDataSet()"<<std::endl;
  bool verbose = 0 ;

  if(verbose){
    for(int i=0;i<data->numEntries();i++){
      const RooArgSet *row = data->get(i);
      RooRealVar*  v1 = (RooRealVar*)  row->find(mB.GetName());
      RooCategory* c0 = (RooCategory*) row->find("magnet");
      RooCategory* c1 = (RooCategory*) row->find("charge");
      RooCategory* c3 = (RooCategory*) row->find("type");
      RooCategory* c4 = (RooCategory*) row->find("mode");
      RooCategory* c5 = (RooCategory*) row->find("Bmode");
      RooCategory* c6 = (RooCategory*) row->find("helBin");
      RooCategory* c7 = (RooCategory*) row->find("DsBDTBin");
      RooCategory* c8 = (RooCategory*) row->find("PhiBDTBin");
      if(par->debug&&0==i%1000) std::cout << "Event " << i <<": magnet "<<c0->getLabel()<<", "<<c3->getLabel()<<": "<<c1->getLabel()<<" "<<c4->getLabel()<<" "<<c5->getLabel()<< " "<<c6->getLabel()<<" "<<c7->getLabel()<<" "<<c8->getLabel()<<"\t mB="<<v1->getVal()<<"   "<<std::endl;
    }
  }

  data->table(*model->Cat())->Print("v");
  if(par->genToys){
    std::cout<<"   NUMBER OF GENERATED TOY EVENTS = "<<data->sumEntries()<<std::endl;
  }else{
    std::cout<<"                            TOTAL = "<<data->sumEntries()<<std::endl;
  }
}



void DsPhiFitting::DefineModel()
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::DefineModel()"<<std::endl;
  /**********************************************/
  /**********************************************/
  /****/           state=BLIND;             /****/
  /**********************************************/
  /**********************************************/
  if(par->genToys) state=UNBLIND;
  model_gen = new DsPhiModel(par,&mB,UNBLIND);
  model = new DsPhiModel(par,&mB,state);
  std::cout << "Model built; continuing..." << std::endl << std::endl;
  reducedlist.add(*model->Cat());
}



void DsPhiFitting::RunFullFit(bool draw=true)
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::RunFullFit()"<<std::endl;
  PrintDataSet();
  double binWidth=10;//MeV/c2
  mB.setBins((mB.getMax()-mB.getMin())/binWidth);
  RooCategory* cat=model->Cat();
  RooSimultaneous*  sim=model->Pdf();
  RooDataHist* hist = data->binnedClone();

  if(par->doFit){
    RooAbsData* abs=(RooAbsData*) data;
    if(par->binned) abs=(RooAbsData*) hist;
    RooFitResult* result = 0;
    std::cout <<"Starting fitTo"<<std::endl;
    //if(par->fundamental){
      //if(par->automatic){
        //result = sim->fitTo(*abs,RooFit::Save(),RooFit::ExternalConstraints(model->Constraints()),RooFit::Hesse(false),RooFit::PrintLevel(-1),RooFit::Optimize(false));
	//}else{
	//   result = sim->fitTo(*abs,RooFit::Save(),RooFit::ExternalConstraints(model->Constraints()),RooFit::Optimize(false));
	// }
      //}else{
      result = sim->fitTo(*abs,RooFit::Save(),RooFit::Minos(par->minos),RooFit::PrintLevel(par->batch?-1:1),RooFit::Optimize(false),RooFit::Timer(true));
      //}
    TCanvas* corrCanv = new TCanvas("corrCanv","",1000,1000); corrCanv->cd();
    result->correlationHist()->Draw("colz");
    result->correlationHist()->DrawCopy("text same");
    std::cout <<"\n-------- minNLL = "; std::cout.precision(10); std::cout << result->minNll()<<" covQual = " << result->covQual() << std::endl;


    result->Print("v");
    model->PrintResult();

    gSystem->Exec("mkdir -p results");
    TFile f("results/lastRooFitResult.root","RECREATE");
    result->Write(); 
    f.Close(); 

  } // closes if(par->doFit)

  if(!draw) return;
  if(par->debug) std::cout<<"Drawing Histograms..."<<std::endl;

  //---------- Start of drawing --------------

// -------------------------------------------------------------------------

	mB.setRange("lowersideband", 5050, 5230);
	mB.setRange("uppersideband", 5330, 5900);
//	mB.setRange("fullfitrange", 5050, 5900);

//	TCut lowersideband = TCut(TString(convertToString(mB) + "<" + convertToString(5230)));
//	TCut uppersideband = TCut(TString(convertToString(mB) + ">" + convertToString(5330)));
	TString lowersideband_string = Form("%s>5050&&%s<5230",mB.GetName(),mB.GetName());//"mB<5230";
	TString uppersideband_string = Form("%s>5330&&%s<5900",mB.GetName(),mB.GetName());//"mB>5330";
  TString wholerange_string    = Form("%s>5050&&%s<5900",mB.GetName(),mB.GetName());//"mB>5330";
  int n_sideband_low, n_sideband_high, n_wholerange;
	//Double_t Nentries(data->sumEntries());

  // -------------------------------------------------------------------------
  for(std::vector<std::string>::iterator b=BmodeList.begin();b!=BmodeList.end();b++){
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
      if(par->debug) std::cout<<"Mode: "<<(*b).c_str()<<", "<<(*m).c_str()<<std::endl;
      std::string canvas_name = Form("canvas_%s_%s",(*b).c_str(),(*m).c_str());
      if (par->sumOverCharges) canvas_name += "_summed";
      if (par->splitHel) canvas_name += "_splitHel";
      TCanvas* canvas=new TCanvas(canvas_name.c_str(),Form("%s_%s",(*b).c_str(),(*m).c_str()),40,0,(HelBinList.size()*chargeList.size()*600),600);
      TCanvas* canRes=0;
      RooHist* hresid=0;
      RooHist* hresid_2=0;
      //if( par->doFit )
        std::string canvasres_name = Form("canvas_%s_%s",(*b).c_str(),(*m).c_str());
        if (par->sumOverCharges) canvas_name += "_summed";
        if (par->splitHel) canvas_name += "_splitHel";
        canRes = new TCanvas(canvasres_name.c_str(),Form("%s_%s",(*b).c_str(),(*m).c_str()),40,0,(HelBinList.size()*chargeList.size()*600),600);
      //std::map<std::string,float> maxH;
      float maxH = 0.0;
      std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,RooPlot*> > > > > plot;
      //if((*p)==pass) continue;

      for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
        for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
          for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s",(*h).c_str(),(*ds).c_str(),(*ph).c_str(),(*b).c_str(),(*m).c_str(),(*c).c_str(),(*a).c_str());
                plot[*h][*ds][*ph][*c][*a] = mB.frame();
                if(par->debug) std::cout<<"Conditions: "<<(*a).c_str()<<", "<<(*c).c_str()<<std::endl;
                
                if( (*b).c_str()==DsPhi and state==BLIND){      
                  std::cout << "Keeping Bu region blind!" << std::endl;
                  TString cat_string = Form("&&cat==cat::%s",cat_str.c_str());
                  n_sideband_low =(data->reduce(lowersideband_string+cat_string))->numEntries();
                  n_sideband_high=(data->reduce(uppersideband_string+cat_string))->numEntries();
                  n_wholerange   =(data->reduce(wholerange_string+cat_string))->numEntries();

                  if(par->debug) std::cout<<"No. in Lower Sideband: "<< n_sideband_low << std::endl;
                  if(par->debug) std::cout<<"No. in Upper Sideband: "<< n_sideband_high << std::endl;
                  if(par->debug) std::cout<<"No. in Whole Range:    "<< n_wholerange << std::endl;

                  if(par->debug) std::cout<<"Frac. in Lower Sideband: "<< (double)n_sideband_low/n_wholerange << std::endl;
                  if(par->debug) std::cout<<"Frac. in Upper Sideband: "<< (double)n_sideband_high/n_wholerange << std::endl;

                  cat->setLabel(cat_str.c_str());
                  data->plotOn( plot[*h][*ds][*ph][*c][*a],
                                RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                                RooFit::DrawOption("PZ") );
                                //RooFit::CutRange("lowersideband,uppersideband")

                  if(par->doDraw){
                    // Comb + Hills + Horns
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_XXst_%s,pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                                  //RooFit::NormRange("lowersideband,uppersideband"),
                                  RooFit::FillColor(kAzure+5), 
                                  RooFit::DrawOption("F") );
                                  //RooFit::Normalization(n_sideband_low+n_sideband_high,RooAbsReal::NumEvent));

                    // Comb + Hills
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
                                  //RooFit::NormRange("lowersideband,uppersideband"),
                                  RooFit::FillColor(kCyan-1), 
                                  RooFit::DrawOption("F") );
                                  //RooFit::Normalization(n_sideband_low+n_sideband_high, RooAbsReal::NumEvent));
                    // Comb
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())),
                                  //RooFit::NormRange("lowersideband,uppersideband"),
                                  RooFit::FillColor(kAzure-7), 
                                  RooFit::DrawOption("F") ); 
                                  //RooFit::Normalization(n_sideband_low+n_sideband_high, RooAbsReal::NumEvent));
                    // Total Lower
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::LineWidth(2),
                                  RooFit::LineColor(kBlack) );    

                    if(canRes) hresid = plot[*h][*ds][*ph][*c][*a]->pullHist();

                    /*
                    // Total Lower
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::LineWidth(2),
                                  RooFit::LineColor(kBlack),
                                  RooFit::Range("lowersideband"),
                                  //RooFit::NormRange("lowersideband"),
                                  RooFit::Normalization(1,RooAbsReal::Relative));  
                                  //RooFit::Normalization(n_sideband_low,RooAbsReal::NumEvent) );     
                    if(canRes) hresid = plot[*h][*ds][*ph][*c][*a]->pullHist();
                    
                    // Total Upper
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::LineWidth(2),
                                  RooFit::LineColor(kBlack),
                                  RooFit::Range("uppersideband"),
                                  //RooFit::NormRange("uppersideband"),
                                  RooFit::Normalization(1,RooAbsReal::Relative)); 
                                  //RooFit::Normalization(n_sideband_high,RooAbsReal::NumEvent));    
                    if(canRes) hresid_2 = plot[*h][*ds][*ph][*c][*a]->pullHist();
                    */

                  } //closes if(par->doDraw)

                  data->plotOn( plot[*h][*ds][*ph][*c][*a],
                                RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                                RooFit::DrawOption("PZ") );
                                //RooFit::CutRange("lowersideband,uppersideband")
                } else {
                  if(par->debug) std::cout<<"Cat Label: "<< cat_str<<std::endl;
                  cat->setLabel(cat_str.c_str());
                  data->plotOn( plot[*h][*ds][*ph][*c][*a],
                                RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                                RooFit::DrawOption("PZ") );
                  if(par->doDraw){
                    // Comb + Hills + Horns
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_XXst_%s,pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                                  RooFit::LineWidth(3),
                                  RooFit::LineStyle(kSolid) ,
                                  RooFit::FillColor(kAzure+5), 
                                  RooFit::DrawOption("F") );
                    
                    // Comb + Hills
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
                                  RooFit::LineWidth(1),
                                  RooFit::LineStyle(kDotted),
                                  RooFit::FillStyle(1002),
                                  RooFit::FillColor(kCyan-1), 
                                  RooFit::DrawOption("F") );
                    // Comb
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_comb_%s",cat_str.c_str())),
                                  RooFit::LineWidth(2),
                                  RooFit::LineStyle(kDotted),
                                  RooFit::FillColor(kAzure-7), 
                                  RooFit::DrawOption("F") );
                    // Total
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],
                                  RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::LineWidth(2),
                                  RooFit::LineColor(kBlack));

                    if(canRes) hresid = plot[*h][*ds][*ph][*c][*a]->pullHist();
                    
                    sim->plotOn(  plot[*h][*ds][*ph][*c][*a],RooFit::Slice(RooArgSet(*cat)), 
                                  RooFit::ProjWData(RooArgSet(*cat),*data) ,
                                  RooFit::Components(Form("pdf_peak_%s",cat_str.c_str())),
                                  RooFit::LineWidth(3),
                                  RooFit::LineStyle(kSolid) ,
                                  RooFit::LineColor(kRed)  ,
                                  RooFit::FillColor(kRed) );
                  } //closes if(par->doDraw)
                  data->plotOn( plot[*h][*ds][*ph][*c][*a],
                                RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                                RooFit::DrawOption("PZ") );
        	      } //closes else
                  

                  float xgutter=0.010;
                  float leftmargin=0.08;
                  if (par->sumOverCharges) leftmargin=0.13;
                  float bottommargin=0.09;
                  //float nY=pidList.size();
                  //float nY= 1.0;
                  float nY= DsBDTBinList.size()*PhiBDTBinList.size();
                  float nX=chargeList.size()*magnetList.size()*HelBinList.size();
                  float padwidth =(0.999-  leftmargin)/nX;
                  float padheight=(0.999-bottommargin)/nY;
                  float ygutter=xgutter*canvas->GetWindowWidth()/canvas->GetWindowHeight();
                  int iX=HelBinList.size()*(chargeList.size()*(a-magnetList.begin()) + (c-chargeList.begin())) + (h-HelBinList.begin());
                  //int iY=1-(p-pidList.begin()); 
                  //int iY = 0; // we redefine iY as 0, otherwise (if iY is 1) then y1 is ~2 which doesn't work.
                  int iY = PhiBDTBinList.size()*(ds-DsBDTBinList.begin()) + (ph-PhiBDTBinList.begin());
                  float x0=iX*padwidth+(iX?leftmargin:0);
                  float x1=(iX+1)*padwidth+leftmargin;
                  float y0=iY*padheight+(iY?bottommargin:0);
                  float y1=(iY+1)*padheight+bottommargin;
                  //std::cout<<"Laurence: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
                  canpad[*h][*ds][*ph][*b][*m][*c][*a] = new TPad(Form("pad_%s",cat_str.c_str()),"",x0,y0,x1,y1);
                  canpad[*h][*ds][*ph][*b][*m][*c][*a]->SetRightMargin( xgutter/(x1-x0));
                  canpad[*h][*ds][*ph][*b][*m][*c][*a]->SetTopMargin(   ygutter/(y1-y0));
                  canpad[*h][*ds][*ph][*b][*m][*c][*a]->SetLeftMargin(  xgutter+(iX?0:(leftmargin/(x1-x0))));
                  canpad[*h][*ds][*ph][*b][*m][*c][*a]->SetBottomMargin(ygutter+(iY?0:(bottommargin/(y1-y0))));
                  TPad* resPad = (TPad*)canpad[*h][*ds][*ph][*b][*m][*c][*a]->Clone();
                  canvas->cd(); canpad[*h][*ds][*ph][*b][*m][*c][*a]->Draw();
                  canpad[*h][*ds][*ph][*b][*m][*c][*a]->cd();
               
                  plot[*h][*ds][*ph][*c][*a]->Draw();
                  plot[*h][*ds][*ph][*c][*a]->SetTitle("");
               
                  if( iY or (nX-iX>1)){ plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetTitle("");}else{plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D_{s}D^{0}}) [MeV/#font[12]{c}^{2}] ");}
                  if( iX or (nY-iY>1)){ plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetTitle("");}else{plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetTitleSize(0.045/(y1-y0));plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));}
                  if(iY){ plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetLabelColor(0); }else{plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetLabelSize(0.045/(y1-y0));}
                  if(iX){ plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetLabelColor(0); }else{plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetLabelSize(0.045/(y1-y0));}
                  plot[*h][*ds][*ph][*c][*a]->GetXaxis()->SetTitleOffset(0.9);
            	    plot[*h][*ds][*ph][*c][*a]->GetYaxis()->SetTitleOffset(1.0);
                  // Draw B decay title
                  x0=((iX?0:leftmargin)+padwidth*0.79)/(padwidth+(iX?0:leftmargin));
                  x1=((iX?0:leftmargin)+padwidth*0.99)/(padwidth+(iX?0:leftmargin));
                  y0=((iY?0:bottommargin)+padheight*0.89)/(padheight+(iY?0:bottommargin));
                  y1=((iY?0:bottommargin)+padheight*0.99)/(padheight+(iY?0:bottommargin));
                  std::cout<<"Title: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
                  TPaveLabel *pav = new TPaveLabel(x0,y0,x1,y1,title[*h][*ds][*ph][*b][*m][*c][*a].c_str(),"NDC");
                  pav->SetBorderSize(0);
                  pav->SetFillStyle(0);
                  pav->SetTextFont(12);
                  pav->SetTextSize(0.5);
                  pav->SetTextAlign(31);
                  pav->Draw();
                  // Draw LHCb unofficial 
                  x0=((iX?0:leftmargin)+padwidth*0.79)/(padwidth+(iX?0:leftmargin));
                  x1=((iX?0:leftmargin)+padwidth*0.99)/(padwidth+(iX?0:leftmargin));
                  y0=((iY?0:bottommargin)+padheight*0.79)/(padheight+(iY?0:bottommargin));
                  y1=((iY?0:bottommargin)+padheight*0.89)/(padheight+(iY?0:bottommargin));
                  std::cout<<"LHCb: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging       
                  TPaveLabel* lhcblabel = new TPaveLabel(x0,y0,x1,y1,"LHCb Unofficial","NDC");
                  lhcblabel->SetBorderSize(0);
                  lhcblabel->SetFillStyle(0);
                  lhcblabel->SetTextSize(0.5);
                  lhcblabel->SetTextFont(62); 
                  lhcblabel->SetTextAlign(31);  
                  if(!par->genToys) lhcblabel->Draw();
                  // Draw 
                  x0=((iX?0:leftmargin)+padwidth*0.79)/(padwidth+(iX?0:leftmargin));
                  x1=((iX?0:leftmargin)+padwidth*0.99)/(padwidth+(iX?0:leftmargin));
                  y0=((iY?0:bottommargin)+padheight*0.69)/(padheight+(iY?0:bottommargin));
                  y1=((iY?0:bottommargin)+padheight*0.79)/(padheight+(iY?0:bottommargin));
                  std::cout<<"Details: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
                  TPaveLabel *pav2 = new TPaveLabel(x0,y0,x1,y1,bin_detail[*h][*ds][*ph][*b][*m][*c][*a].c_str(),"NDC");
                  pav2->SetBorderSize(0);
                  pav2->SetFillStyle(0);
                  pav2->SetTextFont(12);
                  pav2->SetTextSize(0.5);
                  pav2->SetTextAlign(31);
                  pav2->Draw();


                  if(plot[*h][*ds][*ph][*c][*a]->GetMaximum()>maxH) maxH=plot[*h][*ds][*ph][*c][*a]->GetMaximum();
                  if(hresid){
                    canRes->cd();
                    resPad->Draw();
                    resPad->cd();
                    RooPlot* frame = mB.frame(RooFit::Title("Residual Distribution"));
                    hresid->SetFillColor(4);
                    hresid->SetLineColor(10);
                    frame->addPlotable(hresid,"B");
                    if(hresid_2) hresid_2->SetFillColor(4);
                    if(hresid_2) hresid_2->SetLineColor(10);
                    if(hresid_2) frame->addPlotable(hresid_2,"B");
                    frame->Draw();
                    frame->SetTitle("");
                    frame->GetYaxis()->SetRangeUser(-5,5); // range chosen such that 1/1000 chance of Gaussian error fluctuating to the edge
                  }

              } //end of loop over magnetList
             
            } //end of loop over chargeList
          }
        }
      }

      for(std::vector<std::string>::iterator h=HelBinList.begin();h!=HelBinList.end();h++){ 
        for(std::vector<std::string>::iterator ds=DsBDTBinList.begin();ds!=DsBDTBinList.end();ds++){ 
          for(std::vector<std::string>::iterator ph=PhiBDTBinList.begin();ph!=PhiBDTBinList.end();ph++){
            for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
              for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
                  if(0==plot[*h][*ds][*ph][*c][*a]) continue;
                  plot[*h][*ds][*ph][*c][*a]->SetTitleFont(132,"X"); plot[*h][*ds][*ph][*c][*a]->SetLabelFont(132,"X");
                  plot[*h][*ds][*ph][*c][*a]->SetTitleFont(132,"Y"); plot[*h][*ds][*ph][*c][*a]->SetLabelFont(132,"Y");
                  plot[*h][*ds][*ph][*c][*a]->SetMaximum(maxH);
                  plot[*h][*ds][*ph][*c][*a]->SetMinimum(0.01);
                  if( (*b)==DsPhi and state==BLIND ){  
                    canpad[*h][*ds][*ph][*b][*m][*c][*a]->cd();
                    TPaveLabel *blindpav = new TPaveLabel(5229,0.000001,5329,maxH-5,"BLIND","");
                    blindpav->SetBorderSize(0); blindpav->SetTextSize(0.1); blindpav->SetTextAngle(30); 
                    blindpav->SetTextColor(2); blindpav->SetFillColor(10);
                    blindpav->Draw();
                  }
                  //std::cout<<(*m)<<" "<<(*c)<<" "<<(*a)<<"  chi2 = "<<plot[*h][*ds][*ph[*h][*ds][*ph]c][*a]->chiSquare()<<std::endl;
              }
            }
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

} //end of RunFullFit() funcn


void DsPhiFitting::RunManyFits()
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::RunManyFits()"<<std::endl;
  
  std::string mva(par->MVAMethod);
  std::string mvatype(par->MVAType); 


  if(par->debug) std::cout<<"Running: DsPhiFitting::RunManyFits() --> Chosen setup: " << mvatype << ", " << mva<<std::endl;

  if(HelBinList.size()*DsBDTBinList.size()*PhiBDTBinList.size()*BmodeList.size()*chargeList.size()*magnetList.size()*modeList.size()>1) {
    std::cout << "ERROR Can only run many fits over a single mode at a time" << std::endl;
    return;
  }

  std::string dsmode;
  if(par->modes[Ds2PhiPi]) dsmode = Ds2KKPi;  
  if(par->modes[Ds2KKPi])  dsmode = Ds2KKPi;
  if(par->modes[Ds2PiPiPi])dsmode = Ds2PiPiPi;
  if(par->modes[Ds2KPiPi]) dsmode = Ds2KPiPi;

  std::string years="";
  if(par->dsetsReq[s24] )   years="_2015";
  if(par->dsetsReq[s21] )   years="_2012";
  if(par->dsetsReq[s21r1] ) years="_2011";
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1]) years = "_Run1"; 
  if(par->dsetsReq[s21] && par->dsetsReq[s21r1] && par->dsetsReq[s24]) years = "_All"; 

  if( mva!=BDT && mva!=BDTG && mva!= BDTB) {
    std::cout << " MVA method must be one of BDT, BDTG, or BDTB" << std::endl;
    return;
  }

  // Set up fit range
  int n_p = 5;
  if (mvatype==MC) n_p = 2;
  if (par->nBDTPoints>0) n_p = par->nBDTPoints;
  
  if(par->debug) std::cout<<"Running: DsPhiFitting::RunManyFits() --> Number of points: " << n_p<<std::endl;

  std::map<std::string,std::map<std::string,double>> low;
  std::map<std::string,std::map<std::string,double>> high;
  

  std::map<std::string,std::map<std::string,std::string>> var_name;
  var_name[BDT][Ds]  = bdtDs.GetName();
  var_name[BDTG][Ds] = bdtgDs.GetName();
  var_name[BDTB][Ds] = bdtbDs.GetName();
  var_name[BDT][Phi] = bdtD0.GetName();
  var_name[BDTG][Phi]= bdtgD0.GetName();
  var_name[BDTB][Phi]= bdtbD0.GetName();
  var_name[BDT][MC]  = bdtMC.GetName();
  var_name[BDTG][MC] = bdtgMC.GetName();
  var_name[BDTB][MC] = bdtbMC.GetName();

  low[BDT][Phi]  = -0.4;  high[BDT][Phi] =  1.1;  
  low[BDT][Ds]   = -0.4;  high[BDT][Ds]  =  1.1;
  low[BDTG][Phi] =  0.5;  high[BDTG][Phi]=  1.1;
  low[BDTG][Ds]  =  0.0;  high[BDTG][Ds] =  1.1;
  low[BDTB][Phi] =  0.5;  high[BDTB][Phi]=  1.05;
  low[BDTB][Ds]  = -0.5;  high[BDTB][Ds] =  1.05;
  
  low[BDT][MC]   = -1.0;  high[BDT][MC]  =  0.45;
  low[BDTG][MC]  = -1.0;  high[BDTG][MC] =  1.05;
  low[BDTB][MC]  = -1.0;  high[BDTB][MC] =  1.05;
  double cut_value_MC[n_p];
  double cut_value_Phi[n_p];
  double cut_value_Ds[n_p]; 

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
        cut_value_Ds[i]  = low[mva][Ds]  + i*(high[mva][Ds] -low[mva][Ds] )/n_p;
        cut_value_Phi[j] = low[mva][Phi] + j*(high[mva][Phi]-low[mva][Phi])/n_p;
        
        
        std::string cut_BDT  = Form("%s>%f&&%s>%f", var_name[mva][Ds].c_str(), cut_value_Ds[i], var_name[mva][Phi].c_str(), cut_value_Phi[j]);
        
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
        std::string name_string=  Form("%s_%s_Ds_%f_Phi_%f",dsmode.c_str(), mva.c_str(),cut_value_Ds[i], cut_value_Phi[j] );
        TCanvas* canvas = new TCanvas(name_string.c_str(),"",40,0,600,600);
        RooPlot* simple_plot = mB.frame();
        std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s",
                                   (*HelBinList.begin()).c_str(),
                                   (*DsBDTBinList.begin()).c_str(),
                                   (*PhiBDTBinList.begin()).c_str(),
                                   (*BmodeList.begin()).c_str(),
                                   (*modeList.begin()).c_str(),
                                   (*chargeList.begin()).c_str(),
                                   (*magnetList.begin()).c_str() );
        std::string decay_title = title[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()];    
        cat->setLabel(cat_str.c_str());
        cutdata->plotOn( simple_plot,
                      RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                      RooFit::DrawOption("PZ") );
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::Components(Form("pdf_XXst_%s,pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                      RooFit::FillColor(kAzure+5), 
                      RooFit::DrawOption("F") );
                      // Comb + Hills
        sim->plotOn(  simple_plot,
                      RooFit::Slice(RooArgSet(*cat)), 
                      RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                      RooFit::Components(Form("pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
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

        N_signal[i][j]     = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Peak"];
        //N_sigerr[i][j]     = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Peakerr"];
        N_background[i][j] = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb"];
        N_background_red[i][j] = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb_Reduced"];
        std::cout << "            N Signal: " << N_signal[i][j] << std::endl;
        std::cout << "        N Background: " << N_background[i][j] << std::endl;
        std::cout << "N Background reduced: " << N_background_red[i][j] << std::endl;
        std::string text_filename  = Form("results/manyfits/text/%s/%s/DATA/FitResult_%s%s.txt",dsmode.c_str(), mva.c_str(),name_string.c_str(),years.c_str());
        std::ofstream output_file(text_filename.c_str());
        output_file << cut_value_Ds[i] <<":"<< cut_value_Phi[j] <<":"<< N_signal[i][j] << ":" << N_background_red[i][j];
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
    double cut_value_1D_Phi[n_p*n_p], cut_value_1D_Ds[n_p*n_p]; 

    for (int i = 0; i < n_p ; i++){
      for (int j = 0; j < n_p ; j++){
        Punzi_2D[i*n_p +j]         = N_signal[i][j]/ (2.5 + sqrt(N_background_red[i][j]));
        Significance_2D[i*n_p+j]   = N_signal[i][j]/sqrt(N_signal[i][j]+N_background_red[i][j]);
        cut_value_1D_Ds[i*n_p+j]   =  cut_value_Ds[i];
        cut_value_1D_Phi[i*n_p +j] =  cut_value_Phi[j];
      }
    }

    std::cout << std::endl;
    for (int i=0; i<n_p; i++ ) {
      for (int j=0;j<n_p; j++){
        std::cout << "Cut Value Ds: " << cut_value_1D_Ds[i*n_p+j] <<"\t Cut Value Phi: " << cut_value_1D_Phi[i*n_p+j] << " \t N Sig: " << N_signal[i][j] <<" \t N Bg: " << N_background_red[i][j] << " \t Punzi: " << Punzi_2D[i*n_p +j] <<" \t Significance: " << Significance_2D[i*n_p+j] << std::endl;
      }
    }

    std::cout << std::endl;

    std::cout << "double cut_value_1D_Ds []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << cut_value_1D_Ds[i];
      if(i!=(n_p*n_p-1))std::cout<< ", ";
    }
    std::cout <<"};"<< std::endl;
    
    std::cout << "double cut_value_1D_Phi []={";
    for (int i=0; i<n_p*n_p; i++ ) {
      std::cout << cut_value_1D_Phi[i];
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
    TGraph2D *gr_punzi = new TGraph2D(n_p*n_p, cut_value_1D_Ds, cut_value_1D_Phi, Punzi_2D);
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
    TGraph2D *gr_sig = new TGraph2D(n_p*n_p, cut_value_1D_Ds, cut_value_1D_Phi, Significance_2D);
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
  } else if (mvatype==MC){

    for (int i = 0; i < n_p ; i++){
      

      cut_value_MC[i]  = low[mva][MC]  + i*(high[mva][MC] -low[mva][MC] )/n_p;
      
      std::string cut_BDT  = Form("%s>%f", var_name[mva][MC].c_str(), cut_value_MC[i]);
      if(dsmode==Ds2KKPi)   cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_K0_PIDK >1  && D_K1_PIDK >1  && D_P_PIDK  <-1 &&" + cut_BDT;
      if(dsmode==Ds2KPiPi)  cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_K_PIDK  >1  && D_P0_PIDK <-1 && D_P1_PIDK <-1 &&" + cut_BDT;
      if(dsmode==Ds2PiPiPi) cut_BDT = "D0_K0_PIDK >1 && D0_K1_PIDK >1 && D_P0_PIDK <-1 && D_P1_PIDK >1  && D_P2_PIDK <-1 &&" + cut_BDT;

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
      std::string cat_str = Form("%s_%s_%s_%s_%s_%s_%s",
                                 (*HelBinList.begin()).c_str(),
                                 (*DsBDTBinList.begin()).c_str(),
                                 (*PhiBDTBinList.begin()).c_str(),
                                 (*BmodeList.begin()).c_str(),
                                 (*modeList.begin()).c_str(),
                                 (*chargeList.begin()).c_str(),
                                 (*magnetList.begin()).c_str() );
      std::string decay_title = title[*HelBinList.begin()][*DsBDTBinList.begin()][*PhiBDTBinList.begin()][*BmodeList.begin()][*modeList.begin()][*chargeList.begin()][*magnetList.begin()];    
      cat->setLabel(cat_str.c_str());
      cutdata->plotOn( simple_plot,
                    RooFit::Cut(Form("cat==cat::%s",cat_str.c_str())), 
                    RooFit::DrawOption("PZ") );
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::Components(Form("pdf_XXst_%s,pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str(),cat_str.c_str())),
                    RooFit::FillColor(kAzure+5), 
                    RooFit::DrawOption("F") );
                    // Comb + Hills
      sim->plotOn(  simple_plot,
                    RooFit::Slice(RooArgSet(*cat)), 
                    RooFit::ProjWData(RooArgSet(*cat),*cutdata) ,
                    RooFit::Components(Form("pdf_XstX_%s,pdf_comb_%s",cat_str.c_str(),cat_str.c_str())),
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

      N_signalMC[i]         = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Peak"];
      N_backgroundMC[i]     = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb"];
      N_background_redMC[i] = fit_results[both][DsBDTbin1][PhiBDTbin1][DsD0][dsmode][both][both]["Comb_Reduced"];
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
}


void DsPhiFitting::RunManyToys()
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::RunManyToys()"<<std::endl;
  RooArgSet* initialPars= model_gen->GetParameters();
  RooArgSet* fittedPars = model->GetParameters();
  
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
  std::cout<<" "<<std::endl;
  OrderToys(par->nToys);
}



void DsPhiFitting::OrderToys(int n)
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::OrderToys()"<<std::endl;
  RooAbsPdf* sim = model_gen->Pdf();
  RooAbsCategory* cat = model_gen->Cat();
  
  bool DB=false;
  int secs=time(NULL);
  if(1==n) secs=20110602;
  RooRandom::randomGenerator()->SetSeed(secs);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooMCStudy* mcstudy = new RooMCStudy(*sim ,reducedlist,RooFit::Binned(par->binned),RooFit::FitOptions(RooFit::Save(true),RooFit::PrintLevel(DB?1:-1)));

  int nEvtsPerSample=sim->expectedEvents(*cat);
  if(1==n&&!DB){
    mcstudy->generate(n,nEvtsPerSample,true);
    data=(RooDataSet*)mcstudy->genData(n-1);
  }else{
    std::cout<<"Generating "<<n<<"toys of "<<nEvtsPerSample<<" events."<<std::endl;
    mcstudy->generateAndFit(n,nEvtsPerSample);
    if(DB){return;}
    
    gSystem->Exec((std::string("mkdir -p ")+par->toylocation).c_str());
    std::string toyfile=par->toylocation+"/toy_";
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){toyfile+=(*m)+underscore;}
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



void DsPhiFitting::DisplayToys()
{
  if(par->debug) std::cout<<"Running: DsPhiFitting::DisplayToys()"<<std::endl;
  RooSimultaneous*  sim=model_gen->Pdf();
  
  RooMCStudy* mcstudy = new RooMCStudy(*sim,reducedlist,RooFit::FitOptions("r"),RooFit::Extended());
  int nFits=1;
  int nBad=0;
  int nFPD=0;
  std::string toyfile("toy_");
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){toyfile+=(*m)+underscore;}
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
      if(r->covQual()<2){nBad++;continue;}
      if(r->covQual()<3){nFPD++;continue;}
      std::cout <<"\r   Collecting "<<nFits<<" toys, of which "<<nBad<<" don't converge and "<<100*float(nFPD)/nFits<<"% are forced positive definite."<<std::flush;
      mcstudy->addFitResult(*r);
      nFits++;
    }
  }
  std::cout<<"\n"<<std::endl;
  if(1==nFits){
    std::cout<<"SOMETHING WRONG. NO TOYS LOADED"<<std::endl;
    return;
  }
  
  if(true){
    int nv=0;
    RooRealVar* var=0;
    gStyle->SetPalette(1);//gStyle->SetOptStat(0) ;
    gStyle->SetOptTitle(0) ; 
    TCanvas* c_all = new TCanvas("canvas_all_pull","Pulls",0,0,800,800) ; 
    c_all->Divide(5,4) ;
    RooArgSet* vars = sim->getVariables();
    TIterator* it = vars->createIterator();
    while((var = (RooRealVar*)it->Next())) {
      std::cout<< std::endl<< var->GetName()<<std::flush;
      if(var->isConstant()) continue;
      if(var->InheritsFrom("RooAbsCategory")) continue;
      if(0==std::string(var->GetName()).find(mB.GetName())) continue;
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
      std::cout<<"  ------->  Plotting PULL " << std::flush;
      RooPlot* frame1 = mcstudy->plotParam(*var,RooFit::Bins(20));
      RooPlot* frame2 = mcstudy->plotError(*var,RooFit::Bins(20));
      RooPlot* frame3 = mcstudy->plotPull( *var,RooFit::Bins(20),RooFit::Range(-5,5),RooFit::FitGauss(true));frame3->SetTitle(var->GetTitle());  
      nv++;
      TCanvas* c = new TCanvas(Form("canv_%s",var->GetName()),var->GetName(),20*nv,20*nv,900,300); c->Divide(3) ;
      c->cd(1) ; gPad->SetLeftMargin(0.15); frame1->GetXaxis()->SetNdivisions(5); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw() ;
      c->cd(2) ; gPad->SetLeftMargin(0.15); frame2->GetXaxis()->SetNdivisions(5); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw() ;
      c->cd(3) ; gPad->SetLeftMargin(0.15); frame3->GetXaxis()->SetNdivisions(5); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw() ;
      gStyle->SetTitleFontSize(0.1);
      c_all->cd(nv); gPad->SetLeftMargin(0.15); frame3->Draw();
      ((TPaveText*)gPad->FindObject("pullGauss_paramBox"))->SetTextSize(0.04); //was 0.07
      gPad->Update();
      c->Print(Form("toysDir/plots/%s.pdf",c->GetName()));
    }
    //c_all->cd(nv+1)->SetLeftMargin(0.15); 
    //mcstudy->plotNLL(RooFit::Bins(20))->Draw();  
    c_all->Print(Form("toysDir/plots/%s.pdf",c_all->GetName()));
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
    var[0] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);          var[0]->setRange(1100, 1300);
    var[1] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);        var[1]->setRange(250, 400);
    var[2] = new RooRealVar(*(RooRealVar*)model_gen->yield_peak[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);         var[2]->setRange(70, 160);
    var[3] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);          var[3]->setRange(450, 950);
    var[4] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);        var[4]->setRange(250, 550);
    var[5] = new RooRealVar(*(RooRealVar*)model_gen->yield_comb[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);         var[5]->setRange(100, 350);
    var[6] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);      var[6]->setRange(1800, 2300);
    var[7] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);    var[7]->setRange(350, 700);
    var[8] = new RooRealVar(*(RooRealVar*)model_gen->PR_total_yield[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);     var[8]->setRange(100, 300);

    var[9] = new RooRealVar(*(RooRealVar*)model_gen->mean_B[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);              var[9]->setRange(5277, 5281);
    var[10] = new RooRealVar(*(RooRealVar*)model_gen->sigma[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);              var[10]->setRange(12, 16);
    var[11] = new RooRealVar(*(RooRealVar*)model_gen->sigma[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi][both][both]);            var[11]->setRange(11, 17);
    var[12] = new RooRealVar(*(RooRealVar*)model_gen->sigma[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi][both][both]);             var[12]->setRange(8, 20);
    var[13] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi]);                     var[13]->setRange(-0.0055, -0.003);
    var[14] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2PiPiPi]);                   var[14]->setRange(-0.006, -0.002);
    var[15] = new RooRealVar(*(RooRealVar*)model_gen->comb_slope[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KPiPi]);                    var[15]->setRange(-0.007, -0.0015);
    var[16] = new RooRealVar(*(RooRealVar*)model_gen->frac[both][DsBDTbin1][PhiBDTbin1][DsD0][Ds2KKPi][both][both]);               var[16]->setRange(0.4, 0.6);

    std::map<int, std::map<int,TH1*> > hist;
    for(unsigned int i=0;i<var.size()-1;i++){
      for(unsigned int j=i+1;j<var.size();j++){
        hist[i][j] = mcstudy->fitParDataSet().createHistogram(Form("cor_%i_%i",i,j),*var[i],RooFit::YVar(*var[j])) ;
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

void DsPhiFitting::SetLHCbStyle(std::string type = ""){

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
