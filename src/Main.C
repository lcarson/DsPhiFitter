#include <iostream>
#include <string>

#include "TApplication.h"
#include "TROOT.h"
#include "RooMsgService.h"

#include "CommonTools.h"
#include "Parameters.h"
#include "DsPhiFitting.h"

int main(int argc, char *argv[])
{
	gROOT->SetStyle("Plain");
	
  TApplication app("app",0,0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  Parameters* par = new Parameters();
  int sc=par->readCommandLine(argc,argv);
  if(sc) return sc;
  if(par->batch){
    gROOT->SetBatch(kTRUE);
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }
  if(par->debug) std::cout<<"Running: DsPhiFitting"<<std::endl;
  DsPhiFitting fit(par,&app);
  std::cout<<std::endl;
}
