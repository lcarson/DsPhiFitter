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
//#include "RooHORNSdini.h"
//#include "RooHILLdini.h"

DsPhiModel::DsPhiModel(Parameters* p, RooRealVar* pmB, bool nB)
: Base()
  //, misID(p)
  //, calibration(p)
  //, pid("pid","PID boolean")
, type("type","Type of data") 
, mode("mode","D_{s} decay mode")
, charge("charge","batchelor charge")
, magnet("magnet","magnet polarity")
, blindString("BlindLeadingTheBlind")
, needsBlinding(nB)
  //, pidList()
, typeList()
, modeList()
, chargeList()
, magnetList()
  //, m_pidCut(-1000)
{
  par=p;
  mB=pmB;
  rand = new TRandom3(time(NULL));
  if(needsBlinding){
    RooRandom::randomGenerator()->SetSeed(20110602);
    std::cout<<"blind model envoked"<<std::endl;
  }
  //  m_pidCut=int(par->pidCut);
  DefineRooCategories();
  DefineModel();
}


void DsPhiModel::DefineRooCategories() //like D0h with "PIDcut" removed
{
  if(par->modes[Ds2KKPi])      modeList.push_back(Ds2KKPi);
  if(par->modes[Ds2PiPiPi])    modeList.push_back(Ds2PiPiPi);
  if(par->modes[Ds2KPiPi])     modeList.push_back(Ds2KPiPi);
 
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){       
      mode.defineType( (*m).c_str()); }
    
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

  cat = new RooCategory("cat","mode/charge/polarity");
  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
    for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
        std::stringstream str;
        str<<(*m)<<underscore<<(*c)<<underscore<<(*a);
        cat->defineType(str.str().c_str());
      }  
    }
  }
}


void DsPhiModel::DefineModel()
{

  // --------- Yields for DsD0, in an array (as line 447 of Model.C) ---------------
  for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
    for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){

      yield_dsd0[Ds2KKPi][*c][*a] =   new RooRealVar(Form("yield_dsd0_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "Yield D_{s}D0 (KK#pi)", 1660, 0, 8300);
      yield_dsd0[Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_dsd0_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "Yield D_{s}D0 (#pi#pi#pi)", 350, 0, 1750);
      yield_dsd0[Ds2KPiPi][*c][*a] =  new RooRealVar(Form("yield_dsd0_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "Yield D_{s}D0 (K#pi#pi)", 200, 0, 1000);

      yield_comb[Ds2KKPi][*c][*a] =   new RooRealVar(Form("yield_comb_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "Yield Comb. (KK#pi)", 3000, 0, 6000);
      yield_comb[Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_comb_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "Yield Comb. (#pi#pi#pi)", 3000, 0, 6000);
      yield_comb[Ds2KPiPi][*c][*a] =  new RooRealVar(Form("yield_comb_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "Yield Comb. (K#pi#pi)", 3000, 0, 6000);

      //yield_dsd0st[Ds2KKPi][*c][*a] =   new RooRealVar(Form("yield_dsd0st_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "Yield D_{s}D0^{*} (KK#pi)", 3000, 0, 6000);
      //yield_dsd0st[Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_dsd0st_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "Yield D_{s}D0^{*} (#pi#pi#pi)", 3000, 0, 6000);
      //yield_dsd0st[Ds2KPiPi][*c][*a] =  new RooRealVar(Form("yield_dsd0st_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "Yield D_{s}D0^{*} (K#pi#pi)", 3000, 0, 6000);

      //yield_dsstd0[Ds2KKPi][*c][*a] =   new RooRealVar(Form("yield_dsstd0_Ds2KKPi_%s_%s",(*c).c_str(),(*a).c_str()),   "Yield D_{s}^{*}D0 (KK#pi)", 3000, 0, 6000);
      //yield_dsstd0[Ds2PiPiPi][*c][*a] = new RooRealVar(Form("yield_dsstd0_Ds2PiPiPi_%s_%s",(*c).c_str(),(*a).c_str()), "Yield D_{s}^{*}D0 (#pi#pi#pi)", 3000, 0, 6000);
      //yield_dsstd0[Ds2KPiPi][*c][*a] =  new RooRealVar(Form("yield_dsstd0_Ds2KPiPi_%s_%s",(*c).c_str(),(*a).c_str()),  "Yield D_{s}^{*}D0 (K#pi#pi)", 3000, 0, 6000);

    } //end of loop over chargeList
  }  //end of loop over magnetList

  // --------- Mean B mass ---------------
  double default_mB=5282.5;
  std::map< std::string, std::map< std::string, std::map< std::string, RooRealVar* > > > mean_B; //could put this defn in the .h
  for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
      mean_B[*m][both][both] = new RooRealVar(Form("mean_B_%s"           ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][plus][up]   = new RooRealVar(Form("mean_B_plus_up_%s"   ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][minus][dn]  = new RooRealVar(Form("mean_B_minus_dn_%s"  ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][plus][dn]   = new RooRealVar(Form("mean_B_plus_dn_%s"   ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][minus][up]  = new RooRealVar(Form("mean_B_minus_up_%s"  ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][plus][both] = new RooRealVar(Form("mean_B_plus_both_%s" ,(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      mean_B[*m][minus][both]= new RooRealVar(Form("mean_B_minus_both_%s",(*m).c_str()), "Mean B mass", default_mB, 5270, 5290);
      //if(*m==d2kpi||*m==d2pik||*m==d2kk||*m==d2pipi||*m==d2kpipi0||*m==d2pikpi0){
      //  mean_B[*m][plus][both] = mean_B[*m][both][both];
      //  mean_B[*m][minus][both]= mean_B[*m][both][both];
      //}
      mean_B[*m][both][up]   = mean_B[*m][both][both];
      mean_B[*m][both][dn]   = mean_B[*m][both][both];
  }

  // --------- Signal width ---------------
  std::map<std::string,std::map<std::string,std::map<std::string,RooRealVar*> > > sigma_dsd0;
  for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
    sigma_dsd0[*m][both][both] = new RooRealVar(Form("sigma_dsd0_%s",(*m).c_str())          ,"D_{s}D0 sigma",14,  5,  20);
    sigma_dsd0[*m][plus][both] = sigma_dsd0[*m][both][both];
    sigma_dsd0[*m][minus][both] = sigma_dsd0[*m][both][both];
    sigma_dsd0[*m][plus][up]   = new RooRealVar(Form("sigma_dsd0_%s_plus_up",(*m).c_str())  ,"D_{s}D0 sigma",14,  5,  20);
    sigma_dsd0[*m][plus][dn]   = new RooRealVar(Form("sigma_dsd0_%s_plus_dn",(*m).c_str())  ,"D_{s}D0 sigma",14,  5,  20);
    sigma_dsd0[*m][minus][up]  = new RooRealVar(Form("sigma_dsd0_%s_minus_up",(*m).c_str()) ,"D_{s}D0 sigma",14,  5,  20);//sigma_dsd0[*m][minus][up] = sigma_dsd0[*m][plus][dn];
    sigma_dsd0[*m][minus][dn]  = new RooRealVar(Form("sigma_dsd0_%s_minus_dn",(*m).c_str()) ,"D_{s}D0 sigma",14,  5,  20);//sigma_dsd0[*m][minus][dn] = sigma_dsd0[*m][plus][up];
    sigma_dsd0[*m][both][up]   = sigma_dsd0[*m][both][both]; 
    sigma_dsd0[*m][both][dn]   = sigma_dsd0[*m][both][both];
  }

  // --------- Signal tail params ---------------
  std::map<std::string,RooRealVar*> ncb_dsd0;
  std::map<std::string,RooRealVar*> alpha_dsd0;
  for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
    ncb_dsd0[*m]     = new RooRealVar(Form("ncb_dsd0_%s",(*m).c_str())     , "D_{s}D0 nCB",   2.000, 0,  10);
    alpha_dsd0[*m]   = new RooRealVar(Form("alphaL_dsd0_%s",(*m).c_str())  , "D_{s}D0 alpha", 5.0, 0, 10.0); //in Model.C it was 0.113, 0, 0.2  
    alpha_dsd0[*m]->setConstant();
  }


  // --------- Comb. background slope ---------------
  std::map<std::string,RooRealVar*> comb_slope_dsd0;
  for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
    comb_slope_dsd0[*m] = new RooRealVar(Form("comb_slope_dsd0_%s",(*m).c_str())   ,"D_{s}D0 comb. slope", -0.1, -1.0, -0.0000000001);
  }


  // --------- PartReco bkg: RooHORNSDini for DsD0*  ---------------
  //std::map<std::string,RooRealVar*> HORNS_a;       //first horn peak mass
  //std::map<std::string,RooRealVar*> HORNS_b;       //second horn peak mass
  //std::map<std::string,RooRealVar*> HORNS_csi;     //slope parameter, <1 if peaks of unequal height
  //std::map<std::string,RooRealVar*> HORNS_shift;   //shift along B mass - rigid shift of full shape
  //std::map<std::string,RooRealVar*> HORNS_sigma;   //width of core resolution Gaussian
  //std::map<std::string,RooRealVar*> HORNS_R;       //ratio of widths of two Gaussians
  //std::map<std::string,RooRealVar*> HORNS_f;       //between 0 and 1, fraction of total area in core Gauss

  //for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
  //HORNS_a[*m]      = new RooRealVar(Form("HORNS_a_%s",(*m).c_str()), "HORNS a", 5052, 5000,  5100);
  //HORNS_b[*m]      = new RooRealVar(Form("HORNS_b_%s",(*m).c_str()), "HORNS b", 5124, 5100,  5180);
  //HORNS_csi[*m]    = new RooRealVar(Form("HORNS_csi_%s",(*m).c_str()), "HORNS csi", 0.8, 0,  1);
  //HORNS_shift[*m]  = new RooRealVar(Form("HORNS_shift_%s",(*m).c_str()), "HORNS shift", 4.8, 0,  20);
  //HORNS_sigma[*m]  = new RooRealVar(Form("HORNS_sigma_%s",(*m).c_str()), "HORNS sigma", 14, 0,  30);
  //HORNS_R[*m]      = new RooRealVar(Form("HORNS_R_%s",(*m).c_str()), "HORNS R", 4.9, 0,  50);
  //HORNS_f[*m]      = new RooRealVar(Form("HORNS_f_%s",(*m).c_str()), "HORNS f", 0.97, 0,  1);
    //}

  // --------- PartReco bkg: RooHILLDini for Ds*D0 ---------------
  //std::map<std::string,RooRealVar*> HILL_a;       
  //std::map<std::string,RooRealVar*> HILL_b;      
  //std::map<std::string,RooRealVar*> HILL_csi;    
  //std::map<std::string,RooRealVar*> HILL_shift;   
  //std::map<std::string,RooRealVar*> HILL_sigma; 
  //std::map<std::string,RooRealVar*> HILL_R;       
  //std::map<std::string,RooRealVar*> HILL_f;      

  //for(std::vector<std::string>::iterator m=allmodeList.begin();m!=allmodeList.end();m++){
  //HILL_a[*m]      = new RooRealVar(Form("HILL_a_%s",(*m).c_str()), "HILL a", 4800, 4700,  4900);
  //HILL_b[*m]      = new RooRealVar(Form("HILL_b_%s",(*m).c_str()), "HILL b", 5229, 5127,  5327);
  //HILL_csi[*m]    = new RooRealVar(Form("HILL_csi_%s",(*m).c_str()), "HILL csi", 0.8, -10,  10);
  //HILL_shift[*m]  = new RooRealVar(Form("HILL_shift_%s",(*m).c_str()), "HILL shift", -8.5, -100,  100);
  //HILL_sigma[*m]  = new RooRealVar(Form("HILL_sigma_%s",(*m).c_str()), "HILL sigma", 2.1, 0,  200);
  //HILL_R[*m]      = new RooRealVar(Form("HILL_R_%s",(*m).c_str()), "HILL R", 0.13, 0,  20);
  //HILL_f[*m]      = new RooRealVar(Form("HILL_f_%s",(*m).c_str()), "HILL f", 0.01, 0,  1);
  //}





  // --------- Define the simultaneous PDFs ---------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  sim = new RooSimultaneous("model","Simultaneous model",*cat);

  RooAbsPdf *pdf_dsd0 = 0;
  //RooAbsPdf *pdf_dsd0st = 0;
  //RooAbsPdf *pdf_dsstd0 = 0;
  RooAbsPdf *pdf_comb_dsd0 = 0;

  for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
    for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
      for(std::vector<std::string>::iterator c=chargeList.begin();c!=chargeList.end();c++){

        std::string tag=(*m)+underscore+(*c)+underscore+(*a);
        //DsD0
        RooAbsReal* mub = mean_B[*m][*c][*a];
        RooAbsReal* sig_dsd0 = sigma_dsd0[*m][*c][*a];
        RooAbsReal* n_dsd0 = ncb_dsd0[*m];
        RooAbsReal* alp_dsd0 = alpha_dsd0[*m];
        pdf_dsd0  = new RooCBShape( Form("pdf_dsd0_%s",tag.c_str()), "", *mB, *mub, *sig_dsd0, *alp_dsd0, *n_dsd0 );

        //comb bkg
        RooRealVar* comb_dsd0 = comb_slope_dsd0[*m];
        pdf_comb_dsd0 = new RooExponential( Form("pdf_comb_dsd0_%s",tag.c_str()), "",*mB,*comb_dsd0);
    
        //DsD0*
        //RooRealVar* dsd0st_a = HORNS_a[*m];
        //RooRealVar* dsd0st_b = HORNS_b[*m];
        //RooRealVar* dsd0st_csi = HORNS_csi[*m];
        //RooRealVar* dsd0st_shift = HORNS_shift[*m];
        //RooRealVar* dsd0st_sigma = HORNS_sigma[*m];
        //RooRealVar* dsd0st_R = HORNS_R[*m];
        //RooRealVar* dsd0st_f = HORNS_f[*m];
        //pdf_dsd0st = new RooHORNSdini(Form("pdf_dsd0st_%s",tag.c_str()), "", *mB, *dsd0st_a, *dsd0st_b, *dsd0st_csi, *dsd0st_shift, *dsd0st_sigma, *dsd0st_R, *dsd0st_f);

        //Ds*D0
        //RooRealVar* dsstd0_a = HILL_a[*m];
        //RooRealVar* dsstd0_b = HILL_b[*m];
        //RooRealVar* dsstd0_csi = HILL_csi[*m];
        //RooRealVar* dsstd0_shift = HILL_shift[*m];
        //RooRealVar* dsstd0_sigma = HILL_sigma[*m];
        //RooRealVar* dsstd0_R = HILL_R[*m];
        //RooRealVar* dsstd0_f = HILL_f[*m];
        //pdf_dsstd0 = new RooHILLdini(Form("pdf_dsstd0_%s",tag.c_str()), "", *mB, *dsstd0_a, *dsstd0_b, *dsstd0_csi, *dsstd0_shift, *dsstd0_sigma, *dsstd0_R, *dsstd0_f);




        //List of PDFs and yields
        //With PartReco
        //RooArgSet pdflist( *pdf_dsd0, *pdf_comb_dsd0, *pdf_dsd0st, *pdf_dsstd0 );
        //RooArgSet nevents( *yield_dsd0[*m][*c][*a], *yield_comb[*m][*c][*a], *yield_dsd0st[*m][*c][*a], *yield_dsstd0[*m][*c][*a] );
        //Without PartReco
        RooArgSet pdflist( *pdf_dsd0, *pdf_comb_dsd0 );
        RooArgSet nevents( *yield_dsd0[*m][*c][*a], *yield_comb[*m][*c][*a] );

        //Add to master PDF
        RooAddPdf* masterPdf = new RooAddPdf(Form("masterPdf_%s",tag.c_str()) ,"",pdflist, nevents);
        std::stringstream str;
        str<<(*m)<<underscore<<(*c)<<underscore<<(*a);
        sim->addPdf(*masterPdf,str.str().c_str());

      }
    }
  } //closing of for loops over m,a and c

} //end of funcn DefineModel()



RooArgSet* DsPhiModel::GetParameters()
{
  RooArgSet* pars = sim->getParameters(RooArgSet(*mB,mode,charge,magnet));
  return pars;
}



void DsPhiModel::PrintResult()
{
  if(not needsBlinding){
 
    for(std::vector<std::string>::iterator m=modeList.begin();m!=modeList.end();m++){
      for(std::vector<std::string>::iterator a=magnetList.begin();a!=magnetList.end();a++){
        if(par->sumOverCharges){
          float y_dsd0 =   yield_dsd0[*m][both][*a]->getVal();
          float y_comb =   yield_comb[*m][both][*a]->getVal();
          //float y_dsd0st = yield_dsd0st[*m][both][*a]->getVal();
          //float y_dsstd0 = yield_dsstd0[*m][both][*a]->getVal();
          std::cout<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
          std::cout<<" D_{s}D0: "<<y_dsd0<<"  comb: "<<y_comb<<std::endl;
          //std::cout<<" D_{s}D0^{*}: "<<y_dsd0st<<" D_{s}^{*}D0: "<<y_dsstd0<<std::endl;
          std::cout<<"---------------------------"<<std::endl;
        }else{
          float y_dsd0_minus =   yield_dsd0[*m][minus][*a]->getVal();
          float y_comb_minus =   yield_comb[*m][minus][*a]->getVal();
          //float y_dsd0st_minus = yield_dsd0st[*m][minus][*a]->getVal();
          //float y_dsstd0_minus = yield_dsstd0[*m][minus][*a]->getVal();
          float y_dsd0_plus =    yield_dsd0[*m][plus][*a]->getVal();
          float y_comb_plus =    yield_comb[*m][plus][*a]->getVal();
          //float y_dsd0st_plus =  yield_dsd0st[*m][plus][*a]->getVal();
          //float y_dsstd0_plus =  yield_dsstd0[*m][plus][*a]->getVal();

       
          std::cout<<(*m)<<", magnet-"<<(*a)<<" :"<<std::endl;
          //std::cout<<" D_{s}D0 minus: "<<y_dsd0_minus<<" comb. minus: "<<y_comb_minus<<" D_{s}D0^{*} minus: "<<y_dsd0st_minus<<" D_{s}^{*}D0: minus: "<<y_dsstd0_minus<<std::endl;
          //std::cout<<" D_{s}D0 plus: "<<y_dsd0_plus<<" comb. plus: "<<y_comb_plus<<" D_{s}D0^{*} plus: "<<y_dsd0st_plus<<" D_{s}^{*}D0: plus: "<<y_dsstd0_plus<<std::endl;
          //std::cout<<"--------------------------------------------------"<<std::endl;
          //std::cout<<" D_{s}D0 total: "<<y_dsd0_minus+y_dsd0_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<" D_{s}D0^{*} total: "<<y_dsd0st_minus+y_dsd0st_plus<<" D_{s}^{*}D0: total: "<<y_dsstd0_minus+y_dsstd0_plus<<std::endl;

          std::cout<<" D_{s}D0 minus: "<<y_dsd0_minus<<" comb. minus: "<<y_comb_minus<<std::endl;
          std::cout<<" D_{s}D0 plus: "<<y_dsd0_plus<<" comb. plus: "<<y_comb_plus<<std::endl;
          std::cout<<"--------------------------------------------------"<<std::endl;
          std::cout<<" D_{s}D0 total: "<<y_dsd0_minus+y_dsd0_plus<<" comb. total: "<<y_comb_minus+y_comb_plus<<std::endl;
        }
        std::cout<<std::endl;
      }
    }
  } //closes if(not needsBlinding)

  std::cout<<std::endl<<std::endl;
}

