#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <math.h>
#include <time.h>

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

//**********************
//***----BLINDING----***
//**********************
Bool_t isBlind2(kFALSE);

DsPhiFitting::DsPhiFitting(Parameters* p, TApplication* app)
: Base()
, toys("toys")
, inputlist("contents of Final ntuple")
, fulllist("contents, including categories")
, storelist("contents of roodataset")
, reducedlist("reduced content")
  //, eventNumber("eventNumber","eventNumber",0,1e99)
, runNumber("runNumber","runNumber",0,1e99)
, mB("Bu_D0constDconstPVconst_Bu_M","B mass with constraints", 5050 , 5900, "MeV/#font[12]{c}^{2}")
, mD("D0_M","m(D)",1840, 1892,"MeV/c^{2}")
  //, bach_dll("Bach_PIDK","#Delta(LL)",-200,200.)
, bdt("D_BDT","MVA",0.4,0.9)
, bid("Bu_ID","PDG_ID of B cand.",-pow(2,32),pow(2,32))
, L0Hadron_TOS("Bu_L0HadronDecision_TOS","",-0.5,1.5)
, L0_TIS("Bu_L0Global_TIS","",-0.5,1.5)
, Hlt1TrackAllL0_TOS("Bu_Hlt1TrackAllL0Decision_TOS","",-0.5,1.5) //this is all == 1
, Hlt2Topo2Body_TOS("Bu_Hlt2Topo2BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo3Body_TOS("Bu_Hlt2Topo3BodyBBDTDecision_TOS","",-0.5,1.5)
, Hlt2Topo4Body_TOS("Bu_Hlt2Topo4BodyBBDTDecision_TOS","",-0.5,1.5)
, type("type","Type of data")
, mode("mode","D^{0} decay mode")
, magnet("magnet","dipole polarity")
, charge("charge","batchelor charge")
, bMassRegion("bMassRegion","m(B) band")
, m_mD_Mu(1866.6)
, m_bdtCut(-0.9)
, state(1)
, typeList()
, modeList()
, chargeList()
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
      RunFullFit(!par->batch); 
    }
  }
  if(!par->batch){
    std::cout<<"Starting: app->Run()"<<std::endl;
    app->Run(true);
  }
}


void DsPhiFitting::DefineRooCategories()
{
  if(par->genToys)      typeList.push_back(toy);
  if(par->dsetsReq[s21])   typeList.push_back(s21);
  if(par->dsetsReq[s21r1]) typeList.push_back(s21r1);
  if(0==typeList.size()) typeList.push_back("UNDEFTYPE");

  type.defineType(toy.c_str());
  type.defineType(s21.c_str());
  type.defineType(s21r1.c_str());

  mode.defineType(Ds2KKPi.c_str());
  mode.defineType(Ds2PiPiPi.c_str());
  mode.defineType(Ds2KPiPi.c_str());
  
  if(par->modes[Ds2KKPi])      modeList.push_back(Ds2KKPi);
  if(par->modes[Ds2PiPiPi])    modeList.push_back(Ds2PiPiPi);
  if(par->modes[Ds2KPiPi])     modeList.push_back(Ds2KPiPi);

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

  bMassRegion.defineType("signal");
  bMassRegion.defineType("bckgrd");


  //Make titles
  std::map<std::string,std::string> chg;
  //std::map<std::string,std::string> opp;
  std::map<std::string,std::string> mag;
  chg[plus]="+"; chg[minus]="#font[122]{-}"; chg[both]="#font[122]{#pm}";
  //opp[plus]="#font[122]{-}"; opp[minus]="+"; opp[both]="#font[122]{#mp}";
  mag[up]="UP";mag[dn]="DN";mag[both]="";
  for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
    for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
      title[Ds2KKPi][*c][*a]    = std::string(Form("B^{%s}#rightarrow D_{s}^{%s} #font[132]{(}#rightarrow KK#pi#font[132]{)} D^{0} #font[132]{(}#rightarrow KK#font[132]{)} %s", chg[*c].c_str(), chg[*c].c_str(), mag[*a].c_str()));
      title[Ds2PiPiPi][*c][*a]  = std::string(Form("B^{%s}#rightarrow D_{s}^{%s} #font[132]{(}#rightarrow #pi#pi#pi#font[132]{)} D^{0} #font[132]{(}#rightarrow KK#font[132]{)} %s", chg[*c].c_str(), chg[*c].c_str(), mag[*a].c_str()));
      title[Ds2KPiPi][*c][*a]  = std::string(Form("B^{%s}#rightarrow D_{s}^{%s} #font[132]{(}#rightarrow K#pi#pi#font[132]{)} D^{0} #font[132]{(}#rightarrow KK#font[132]{)} %s", chg[*c].c_str(), chg[*c].c_str(), mag[*a].c_str()));

      if(par->genToys){
        for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){  title[*m][*c][*a]+=" TOY"; }
        }
    }
  }

  //inputlist.add(eventNumber);
  inputlist.add(bdt);
  inputlist.add(bid);
  inputlist.add(mD);
  //inputlist.add(bach_dll);
  inputlist.add(L0Hadron_TOS);
  inputlist.add(L0_TIS);
  inputlist.add(Hlt1TrackAllL0_TOS);
  inputlist.add(Hlt2Topo2Body_TOS);
  inputlist.add(Hlt2Topo3Body_TOS);
  inputlist.add(Hlt2Topo4Body_TOS);

  storelist.add(mB);
  storelist.add(type);
  storelist.add(mode);
  storelist.add(magnet);
  storelist.add(charge);

  fulllist.add(storelist);
  fulllist.add(inputlist);
  fulllist.add(bMassRegion);

  reducedlist.add(mB);
  inputlist.add(mB);
} //end of funcn DefineRooCategories



std::string DsPhiFitting::reducedFileName(std::string t,std::string m)
{
  gSystem->Exec("mkdir -p roodatasets");
  std::string reducedfile="roodatasets/";
  reducedfile="Laurence_RooDataSets/";
  reducedfile+=t+underscore;
  reducedfile+=m+underscore;
  reducedfile+=std::string( Form("RooDs_%i_%i.root",int(mB.getMin()),int(mB.getMax())) );//,int(m_bdtCut*10)));
  return reducedfile;
}



int DsPhiFitting::LoadDataSet()
{
  data = new RooDataSet("data","Data data",storelist);
  for(std::vector<std::string>::iterator t=typeList.begin();t!=typeList.end();t++){
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
      std::string fn=reducedFileName(*t,*m);
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
        reducedData = MakeDataSet(*t,*m);
        if(NULL==reducedData) return 1;
      }
      data->append(*reducedData);
      delete reducedData;
    }
  }
  AssignGlobalCategory();
  return 0;
}



RooDataSet* DsPhiFitting::MakeDataSet(std::string t,std::string m)
{
  RooDataSet* fulldata = new RooDataSet("fulldata","Data data",fulllist);
  
  std::vector<const char*> tupleName;
  //tupleName.push_back("TTT");
  tupleName.push_back("DecayTree");

  std::string s_mode(null);
  std::string s_type(null);

  for(std::vector<std::string>::iterator il=par->locations.begin();il!=par->locations.end();il++){

    std::vector<std::string> fileName;
    CommonTools::getdir((*il),fileName);
    
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
      std::string fullpathandname = (*il)+slash+(*it);
      std::cout<<"The fullpathandname is: "<<fullpathandname<<"."<<std::endl; //for debugging
      bool accept=false;
      if(  (*it).find((m+underscore).c_str())!=std::string::npos ){
        accept=true;
        s_mode=m;
      }
      
      //Laurence
      if( (*it).find("2011")  !=std::string::npos ) { s_type=s21r1; }
      if( (*it).find("2012")  !=std::string::npos ) { s_type=s21; }
      if(!par->dsetsReq[s_type]) continue; 
      if(t!=s_type) continue;

      std::string s_magnet(both);
      if( (*it).find("MagUp")  !=std::string::npos ) { s_magnet=up; }
      if( (*it).find("MagDown")!=std::string::npos ) { s_magnet=dn; }    

      if(accept){
        TFile* tfile = TFile::Open(fullpathandname.c_str());
        if(tfile){
          //Laurence: cd() into folder inside the .root
          std::string tfileFolderName("");
          if( s_mode==Ds2KKPi )   { tfileFolderName = "B2D0D_D02KK_D2KKPi"; } //NB it's D2KKPi, not Ds2KKPi
          if( s_mode==Ds2PiPiPi ) { tfileFolderName = "B2D0D_D02KK_D2PiPiPi"; }
          if( s_mode==Ds2KPiPi )  { tfileFolderName = "B2D0D_D02KK_D2KPiPi"; }
        
          const char *pathPlusFolderName = (fullpathandname+colon+slash+tfileFolderName).c_str();
          std::cout<<"The pathPlusFolderName is: "<<pathPlusFolderName<<"."<<std::endl; //for debugging
          tfile->cd(pathPlusFolderName);

          //not sure a loop is necessary, since the .root's only have one tuple in them, but we leave it as it is
          for(std::vector<const char*>::iterator itc=tupleName.begin(); itc!=tupleName.end(); itc++){
            if(gDirectory->FindObjectAny(*itc)){
              std::cout<<"Reading tuple '"<<(*itc)<<"' from '"<<fullpathandname<<"'"<<std::endl;
              RooDataSet *input = FinalDataSet(s_type,s_mode,s_magnet,(TTree*)gDirectory->Get(*itc));
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
  delete fulldata;fulldata=reducedData;

  //IdMultipleCandidates(fulldata);
  //reducedData = (RooDataSet*)fulldata->reduce("clone!=clone::true");
  //std::cout<<"\n   Dataset built and checked for multiple candidates: "<<fulldata->sumEntries()-reducedData->sumEntries()<<" removed."<<std::endl;

  std::cout<<"\n   Writing reduced file: '"<<reducedFileName(t,m)<<"' with "<<reducedData->sumEntries()<<" events.\n"<<std::endl;
  reducedData = new RooDataSet("reducedData","",reducedData,storelist);
  TFile out(reducedFileName(t,m).c_str(),"recreate");
  reducedData->Write();
  out.Close();
  
  delete fulldata;
  return reducedData;
} //end of funcn MakeDataSet




RooDataSet* DsPhiFitting::FinalDataSet(const std::string s_type, const std::string s_mode, std::string s_mag, TTree* tree)
{
  if(!tree){ std::cout << "\n THE TREE IS A ZERO POINTER IN "<< s_type <<" "<<s_mode<<" ("<< s_mag <<")"<< std::endl; return 0; }
  RooDataSet *input = new RooDataSet(Form("%s_%s",s_type.c_str(),s_mode.c_str()),Form("Final %s %s",s_type.c_str(),s_mode.c_str()),inputlist,RooFit::Import(*tree));
  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(type,mode,magnet,charge,bMassRegion)); //removed "pid" from the RooArgSet
 

  int nTot=0;
  int nNotL0=0;
  int nNotL0TOS=0;
  int nNotTopoTOS=0;
  int nNotTrackTOS=0;
  for(int i=0;i<input->numEntries();i++){

    mode.setLabel(s_mode.c_str());
    type.setLabel(s_type.c_str());
    magnet.setLabel(s_mag.c_str());

    const RooArgSet *row = input->get(i);
    //RooRealVar* D   = (RooRealVar*) row->find(bach_dll.GetName());
    RooRealVar* q   = (RooRealVar*) row->find(bid.GetName());

    if(q->getVal() > 0 && q->getVal() < pow(2,31)) charge.setLabel("plus");
    if(q->getVal() < 0 || q->getVal() > pow(2,31)) charge.setLabel("minus");

  
    bMassRegion.setLabel("signal");

    if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* d   = (RooRealVar*) row->find(mD.GetName());
      if(fabs(d->getVal()-m_mD_Mu)>25) bMassRegion.setLabel("bckgrd");
    }


    if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* bd  = (RooRealVar*) row->find(bdt.GetName());
      // Main selection [BDT] cuts. Values chosen just by looking at D_BDT distribution, and removing a few events in the tail.
      if(s_mode==Ds2KKPi){
        if(s_type==s21){   if(bd->getVal()<+.40){bMassRegion.setLabel("bckgrd");} }
        if(s_type==s21r1){ if(bd->getVal()<+.40){bMassRegion.setLabel("bckgrd");} }
      }
      if(s_mode==Ds2PiPiPi){
        if(s_type==s21){   if(bd->getVal()<+.30){bMassRegion.setLabel("bckgrd");} }
        if(s_type==s21r1){ if(bd->getVal()<+.30){bMassRegion.setLabel("bckgrd");} }
      }
      if(s_mode==Ds2KPiPi){
        if(s_type==s21){   if(bd->getVal()<+.54){bMassRegion.setLabel("bckgrd");} }
        if(s_type==s21r1){ if(bd->getVal()<+.54){bMassRegion.setLabel("bckgrd");} }
      }
    }
   
     if(signal==std::string(bMassRegion.getLabel())){
      RooRealVar* _L0Global_TIS       = (RooRealVar*) row->find(L0_TIS.GetName());
      RooRealVar* _L0Hadron_TOS       = (RooRealVar*) row->find(L0Hadron_TOS.GetName());
      RooRealVar* _Hlt1TrackAllL0_TOS = (RooRealVar*) row->find(Hlt1TrackAllL0_TOS.GetName());
      RooRealVar* _Hlt2Topo2Body_TOS  = (RooRealVar*) row->find(Hlt2Topo2Body_TOS.GetName());
      RooRealVar* _Hlt2Topo3Body_TOS  = (RooRealVar*) row->find(Hlt2Topo3Body_TOS.GetName());
      RooRealVar* _Hlt2Topo4Body_TOS  = (RooRealVar*) row->find(Hlt2Topo4Body_TOS.GetName());
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
      if(!(_Hlt2Topo2Body_TOS->getVal()||_Hlt2Topo3Body_TOS->getVal()||_Hlt2Topo4Body_TOS->getVal())){
        bMassRegion.setLabel("bckgrd");
        nNotTopoTOS++;
      }
    }

    extra->add(RooArgSet(type,mode,magnet,charge,bMassRegion));
  } //end of loop over entries
 
  std::cout <<"   Trigger: notL0TISTOS: "<<100*nNotL0/float(nTot)<<"% (TIS/L0Acc="<<100*nNotL0TOS/float(nTot-nNotL0)<<"%. Not TrackTOS: "<<100*nNotTrackTOS/float(nTot)<<"%. Not TopoTOS "<<100*nNotTopoTOS/float(nTot)<<"%."<<std::endl;
  input->merge(extra);
  delete extra;
  return input;
} //end of funcn FinalDataSet




void DsPhiFitting::AssignGlobalCategory()
{
  if(par->polarity==up){RooDataSet* reducedData = (RooDataSet*)data->reduce("magnet==magnet::up"); delete data; data=reducedData; }
  if(par->polarity==dn){RooDataSet* reducedData = (RooDataSet*)data->reduce("magnet==magnet::dn"); delete data; data=reducedData; }

  RooDataSet *extra = new RooDataSet("extra","extra stuff",RooArgSet(*model->Cat()));
  for(int i=0;i<data->numEntries();i++){
    const RooArgSet *row = data->get(i);
    RooCategory* _mode   = (RooCategory*) row->find("mode");
    RooCategory* _charge = (RooCategory*) row->find("charge");
    RooCategory* _magnet = (RooCategory*) row->find("magnet");
    if(par->polarity==both) _magnet->setLabel(both.c_str());
    if(par->sumOverCharges) _charge->setLabel(both.c_str());
    model->Cat()->setLabel(Form("%s_%s_%s",_mode->getLabel(),_charge->getLabel(),_magnet->getLabel()));
    //std::cout <<"magnet "<<_magnet->getLabel()<<", "<<_mode->getLabel()<<": "<<_charge->getLabel()<<" "<<model->Cat()->getLabel()<<std::endl;
    extra->add(RooArgSet(*model->Cat()));
  }
  data->merge(extra);

  RooDataSet* reducedData = new RooDataSet("reducedData","",data,reducedlist);
  delete data; data=reducedData;

}



void DsPhiFitting::PrintDataSet()
{
  bool verbose = 0 ;

  if(verbose){
    for(int i=0;i<data->numEntries();i++){
      const RooArgSet *row = data->get(i);
      RooRealVar*  v1 = (RooRealVar*)  row->find(mB.GetName());
      RooCategory* c0 = (RooCategory*) row->find("magnet");
      RooCategory* c1 = (RooCategory*) row->find("charge");
      RooCategory* c3 = (RooCategory*) row->find("type");
      RooCategory* c4 = (RooCategory*) row->find("mode");
      std::cout <<"magnet "<<c0->getLabel()<<", "<<c3->getLabel()<<": "<<c1->getLabel()<<" "<<c4->getLabel()<<"\t mB="<<v1->getVal()<<"   "<<std::endl;
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
  

  //---------- Start of drawing --------------

// -------------------------------------------------------------------------

	mB.setRange("lowersideband", 5050, 5230);
	mB.setRange("uppersideband", 5330, 5900);
//	mB.setRange("fullfitrange", 5050, 5900);

//	TCut lowersideband = TCut(TString(convertToString(mB) + "<" + convertToString(5230)));
//	TCut uppersideband = TCut(TString(convertToString(mB) + ">" + convertToString(5330)));
	TString lowersideband = "mB<5230";
	TString uppersideband = "mB>5330";

	Double_t nd1(0.), nd2(0.);
	Double_t Nentries(data->sumEntries());

// -------------------------------------------------------------------------

  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
    TCanvas* canvas=new TCanvas(Form("canvas_%s%s",(*m).c_str(),(par->sumOverCharges?"_summed":"")),Form("%s",(*m).c_str()),40,0,(chargeList.size()*600),600);
    TCanvas* canRes=0;
    RooHist* hresid=0;
    if( par->doFit )
      canRes = new TCanvas(Form("canres_%s%s",(*m).c_str(),(par->sumOverCharges?"_summed":"")),Form("%s",(*m).c_str()),40,0,(chargeList.size()*600),600);
    //std::map<std::string,float> maxH;
    float maxH = 0.0;
    std::map<std::string,std::map<std::string,RooPlot*> > plot;
    //if((*p)==pass) continue;
    for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
        plot[*c][*a] = mB.frame();

// -------------------------------------------------------------------------
	if(isBlind2)
	{
// -------------------------------------------------------------------------


	std::cout << "Keeping Bu region blind!" << std::endl;
            nd1=(static_cast<Double_t>((data->reduce(lowersideband)->numEntries())))/(static_cast<Double_t>(Nentries)); 
            nd2=(static_cast<Double_t>((data->reduce(uppersideband)->numEntries())))/(static_cast<Double_t>(Nentries)); 

        cat->setLabel(Form("%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str()));
        data->plotOn(plot[*c][*a],RooFit::Cut(Form("cat==cat::%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())), RooFit::DrawOption("PZ"), RooFit::CutRange("lowersideband"));
        data->plotOn(plot[*c][*a],RooFit::Cut(Form("cat==cat::%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())), RooFit::DrawOption("PZ"), RooFit::CutRange("uppersideband"));

        if(par->doDraw){
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::LineWidth(1), RooFit::Range("lowersideband"), RooFit::Normalization(nd1,RooAbsReal::Relative));
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::LineWidth(1), RooFit::Range("uppersideband"),RooFit:: Normalization(nd2,RooAbsReal::Relative));

            if(canRes) hresid = plot[*c][*a]->pullHist();

            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsd0st_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(2),RooFit::LineStyle(kSolid) ,RooFit::LineColor(kGreen), RooFit::Range("lowersideband"), RooFit::Normalization(nd1,RooAbsReal::Relative));
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsd0st_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(2),RooFit::LineStyle(kSolid) ,RooFit::LineColor(kGreen), RooFit::Range("uppersideband"), RooFit::Normalization(nd2,RooAbsReal::Relative));

            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsstd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted),RooFit::FillStyle(3002), RooFit::DrawOption("F"), RooFit::Range("lowersideband"), RooFit::Normalization(nd1,RooAbsReal::Relative) );
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsstd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted),RooFit::FillStyle(3002), RooFit::DrawOption("F"), RooFit::Range("uppersideband"), RooFit::Normalization(nd2,RooAbsReal::Relative) );

            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_comb_dsd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted), RooFit::Range("lowersideband"),RooFit::Normalization(nd1,RooAbsReal::Relative) );
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_comb_dsd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted), RooFit::Range("uppersideband"), RooFit::Normalization(nd2,RooAbsReal::Relative) );

        } //closes if(par->doDraw)
	

// -------------------------------------------------------------------------
	} //closes if(isBlind)
	else
	{
        cat->setLabel(Form("%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str()));
        data->plotOn(plot[*c][*a],RooFit::Cut(Form("cat==cat::%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())), RooFit::DrawOption("PZ") );
        if(par->doDraw){
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::LineWidth(1));
            if(canRes) hresid = plot[*c][*a]->pullHist();
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(2),RooFit::LineStyle(kSolid) ,RooFit::LineColor(kRed)  ,RooFit::FillColor(kRed) );
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsd0st_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(2),RooFit::LineStyle(kSolid) ,RooFit::LineColor(kGreen));
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_dsstd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted),RooFit::FillStyle(3002), RooFit::DrawOption("F") );
            sim->plotOn( plot[*c][*a],RooFit::Slice(RooArgSet(*cat)), RooFit::ProjWData(RooArgSet(*cat),*data) ,RooFit::Components(Form("pdf_comb_dsd0_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str())),RooFit::LineWidth(1),RooFit::LineStyle(kDotted) );

        } //closes if(par->doDraw)
	} //closes else
// -------------------------------------------------------------------------


          float xgutter=0.010;
          float leftmargin=0.08;
          if (par->sumOverCharges) leftmargin=0.13;
          float bottommargin=0.09;
          //float nY=pidList.size();
          float nY= 1.0;
          float nX=chargeList.size()*magnetList.size();
          float padwidth =(0.999-  leftmargin)/nX;
          float padheight=(0.999-bottommargin)/nY;
          float ygutter=xgutter*canvas->GetWindowWidth()/canvas->GetWindowHeight();
          int iX=chargeList.size()*(a-magnetList.begin()) + (c-chargeList.begin());
          //int iY=1-(p-pidList.begin()); 
          int iY = 0; // we redefine iY as 0, otherwise (if iY is 1) then y1 is ~2 which doesn't work.
          float x0=iX*padwidth+(iX?leftmargin:0);
          float x1=(iX+1)*padwidth+leftmargin;
          float y0=iY*padheight+(iY?bottommargin:0);
          float y1=(iY+1)*padheight+bottommargin;
          //std::cout<<"Laurence: x0 = "<<x0<<", x1 = "<<x1<<", y0 = "<<y0<<", y1 = "<<y1<<std::endl; //for debugging
          canpad[*m][*c][*a] = new TPad(Form("pad_%s_%s_%s",(*m).c_str(),(*c).c_str(),(*a).c_str()),"",x0,y0,x1,y1);
          canpad[*m][*c][*a]->SetRightMargin( xgutter/(x1-x0));
          canpad[*m][*c][*a]->SetTopMargin(   ygutter/(y1-y0));
          canpad[*m][*c][*a]->SetLeftMargin(  xgutter+(iX?0:(leftmargin/(x1-x0))));
          canpad[*m][*c][*a]->SetBottomMargin(ygutter+(iY?0:(bottommargin/(y1-y0))));
          TPad* resPad = (TPad*)canpad[*m][*c][*a]->Clone();
          canvas->cd(); canpad[*m][*c][*a]->Draw();
          canpad[*m][*c][*a]->cd();
       
          plot[*c][*a]->Draw();
          plot[*c][*a]->SetTitle("");
       
          if( iY or (nX-iX>1)){ plot[*c][*a]->GetXaxis()->SetTitle("");}else{plot[*c][*a]->GetXaxis()->SetTitleSize(0.045/(y1-y0));plot[*c][*a]->GetXaxis()->SetTitle("#font[12]{m}(#font[12]{D_{s}D^{0}}) [MeV/#font[12]{c}^{2}] ");}
          if( iX or (nY-iY>1)){ plot[*c][*a]->GetYaxis()->SetTitle("");}else{plot[*c][*a]->GetYaxis()->SetTitleSize(0.045/(y1-y0));plot[*c][*a]->GetYaxis()->SetTitle(Form("Events / ( %i MeV/#font[12]{c}^{2} ) ",int(binWidth)));}
          if(iY){ plot[*c][*a]->GetXaxis()->SetLabelColor(0); }else{plot[*c][*a]->GetXaxis()->SetLabelSize(0.045/(y1-y0));}
          if(iX){ plot[*c][*a]->GetYaxis()->SetLabelColor(0); }else{plot[*c][*a]->GetYaxis()->SetLabelSize(0.045/(y1-y0));}
          plot[*c][*a]->GetXaxis()->SetTitleOffset(0.9);
    	  plot[*c][*a]->GetYaxis()->SetTitleOffset(0.7);

          x0=((iX?0:leftmargin)+padwidth*0.50)/(padwidth+(iX?0:leftmargin));
          x1=((iX?0:leftmargin)+padwidth*0.90)/(padwidth+(iX?0:leftmargin));
          y0=((iY?0:bottommargin)+padheight*0.40)/(padheight+(iY?0:bottommargin));
          y1=((iY?0:bottommargin)+padheight*0.50)/(padheight+(iY?0:bottommargin));
          TPaveLabel *pav = new TPaveLabel(x0,y0,x1,y1,title[*m][*c][*a].c_str(),"NDC");
          pav->SetBorderSize(0);
          pav->SetFillStyle(0);
          pav->SetTextFont(12);
          pav->SetTextSize(1);
          pav->SetTextAlign(31);
          pav->Draw();
          x0=((iX?0:leftmargin)+padwidth*0.60)/(padwidth+(iX?0:leftmargin));
          x1=((iX?0:leftmargin)+padwidth*0.90)/(padwidth+(iX?0:leftmargin));
          y0=((iY?0:bottommargin)+padheight*0.75)/(padheight+(iY?0:bottommargin));
          y1=((iY?0:bottommargin)+padheight*0.85)/(padheight+(iY?0:bottommargin));
          TPaveLabel* lhcblabel = new TPaveLabel(x0,y0,x1,y1,"LHCb","NDC");
          lhcblabel->SetBorderSize(0);
          lhcblabel->SetFillStyle(0);
          lhcblabel->SetTextSize(1);
          lhcblabel->SetTextFont(62); 
          lhcblabel->SetTextAlign(31);  
          if(!par->genToys) lhcblabel->Draw();
          if(plot[*c][*a]->GetMaximum()>maxH) maxH=plot[*c][*a]->GetMaximum();
          if(hresid){
            canRes->cd();
            resPad->Draw();
            resPad->cd();
            RooPlot* frame = mB.frame(RooFit::Title("Residual Distribution"));
            hresid->SetFillColor(4);
            hresid->SetLineColor(10);
            frame->addPlotable(hresid,"B");
            frame->Draw();
            frame->SetTitle("");
            frame->GetYaxis()->SetRangeUser(-5,5); // range chosen such that 1/1000 chance of Gaussian error fluctuating to the edge
          }

      } //end of loop over magnetList
     
    } //end of loop over chargeList


    for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
          if(0==plot[*c][*a]) continue;
          plot[*c][*a]->SetTitleFont(132,"X"); plot[*c][*a]->SetLabelFont(132,"X");
          plot[*c][*a]->SetTitleFont(132,"Y"); plot[*c][*a]->SetLabelFont(132,"Y");
          plot[*c][*a]->SetMaximum(maxH);
          plot[*c][*a]->SetMinimum(0.01);
          //if((((*m)==Ds2KKPi||(*m)==Ds2PiPiPi||(*m)==Ds2KPiPi)) and state==BLIND){
          //  canpad[*m][*c][*a]->cd();
          //  TPaveLabel *blindpav = new TPaveLabel(5229,0.000001,5329,maxH-5,"BLIND","");
          //  blindpav->SetBorderSize(0); blindpav->SetTextSize(0.1); blindpav->SetTextAngle(30); 
          //  blindpav->SetTextColor(2); blindpav->SetFillColor(10);
          //  blindpav->Draw();
          //}
          //std::cout<<(*m)<<" "<<(*c)<<" "<<(*a)<<"  chi2 = "<<plot[*c][*a]->chiSquare()<<std::endl;
      }
    }
    
    canvas->Print(Form("results/%s.eps",canvas->GetName()));
    canvas->Print(Form("results/%s.pdf",canvas->GetName()));
    canvas->Print(Form("results/%s.png",canvas->GetName()));
    if(canRes){
     canRes->Print(Form("results/%s.eps",canRes->GetName()));
     canRes->Print(Form("results/%s.pdf",canRes->GetName()));
     canRes->Print(Form("results/%s.png",canRes->GetName()));
    }

  } //end of loop over modeList 

} //end of RunFullFit() funcn



void DsPhiFitting::RunManyToys()
{
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
  RooSimultaneous*  sim=model_gen->Pdf();
  
  RooMCStudy* mcstudy = new RooMCStudy(*sim,reducedlist,RooFit::FitOptions("r"),RooFit::Extended());
  int nFits=1;
  int nBad=0;
  int nFPD=0;
  std::string toyfile("toy_");
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){toyfile+=(*m)+underscore;}
  //toyfile+=std::string(Form("pid%3.1f_",m_pidCut));
  toyfile+=std::string("toyName.1f_");  //need to check what the name should really be...
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
    TCanvas* c_all = new TCanvas("canvas_all_pull","Pulls",0,0,800,800) ; c_all->Divide(4,4) ;
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
      ((TPaveText*)gPad->FindObject("pullGauss_paramBox"))->SetTextSize(0.07);
      gPad->Update();
    }
    //c_all->cd(nv+1)->SetLeftMargin(0.15); 
    //mcstudy->plotNLL(RooFit::Bins(20))->Draw();  
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
    TCanvas* c = new TCanvas("canvas_corr","",0,0,1200,800);

    c->Divide(6,4);
    //Put variables here
    //var[0] = new RooRealVar(*(RooRealVar*)model_gen->R_cabibbo);var[0]->setRange( 0.07,0.08);
    //var[1] = new RooRealVar(*(RooRealVar*)model_gen->gamma);    var[1]->setRange( 0.0,3.1);
  
    std::map<int, std::map<int,TH1*> > hist;
    for(unsigned int i=0;i<var.size()-1;i++){
      for(unsigned int j=i+1;j<var.size();j++){
        hist[i][j] = mcstudy->fitParDataSet().createHistogram(Form("cor_%i_%i",i,j),*var[i],RooFit::YVar(*var[j])) ;
        hist[i][j]->SetMarkerStyle(20);
        hist[i][j]->SetMarkerSize(0.02);
        hist[i][j]->SetMarkerColor(kBlue-7+i);
        hist[i][j]->SetLineColor(kBlue-7+i);
        hist[i][j]->SetTitle("");
        c->cd(ipad);
        gPad->SetLeftMargin(0.25);
        gPad->SetBottomMargin(0.25);
        hist[i][j]->GetYaxis()->SetTitleOffset(1.4);
        hist[i][j]->Draw("box");
        ipad++;
      }
    }

  }


} //end of DisplayToys() funcn
