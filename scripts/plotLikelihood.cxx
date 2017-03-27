void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotLikelihood(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "B #rightarrow D_{s} #phi";
    if (mode == "DKst0") name = "B #rightarrow D K^{*0}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string systematicsDir = "systematicsDir/";
    std::vector<std::string> fileName;
    getdir(systematicsDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> gen_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("vary")  ==std::string::npos ) continue;
        n_points++;
    }

    std::map<double,double> values;
    double* n_br    = new double[n_points];
    double* n_lik   = new double[n_points];

    int count = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("vary")  ==std::string::npos ) continue;
        std::cout << "Filename: " << *it << std::endl;

        TFile newfile((systematicsDir+*it).c_str(),"READ");
        if( !(newfile.GetListOfKeys()->Contains("fitresult_model_reducedData"))  ) continue; 
        RooFitResult* result = (RooFitResult*)newfile.Get("fitresult_model_reducedData");
        RooRealVar* param_final = (RooRealVar*) (result->constPars()).find("Branching_fraction");
        double br = param_final->getVal();
        double minNll = result->minNll();
        std::cout << "Br: " << br  << "\tMin Nll: " << minNll<<  std::endl;
        if(minNll>0){
          values[br] = minNll;
          count++ ;
        }
        

    }

    double minx = 24.0; 
    double miny = 38540.7;      // value and NLL of the approved result

    n_points = values.size();
    double* n_br_sorted    = new double[n_points];
    double* n_lik_sorted   = new double[n_points];
    count = 0;
    for (std::map<double,double>::iterator it=values.begin(); it!=values.end(); ++it){
      n_br_sorted[count]  = it->first;
      n_lik_sorted[count] = it->second - miny;
      count++;
    }



    SetLHCbStyle("norm");

    //TCanvas *cav = new TCanvas("Canvas","Likelihood",200,10,700,500);
  
    TGraph* gr_likelihood = new TGraph(count, n_br_sorted, n_lik_sorted);

    gr_likelihood->SetTitle(";Branching Fraction (#times 10^{-7}) ;#Delta(log(L))"); 
    gr_likelihood->Draw("APC");

    double MIN  = 18.0;         // minimum value of the scanned parameter
    double MAX  = 30.0;         // maximum value
    double Syst =  1.0;         // total systematic error



    double likelMAX     = -99999;
    double likelMAXCONV = -99999;
    double LInterpolation[5000];
    double xInterpolation[5000];
    int    nInterpolation = 100;
    double step = (MAX-MIN)/nInterpolation;
    double tmp = MIN;
    
    for(int i=0; i<=nInterpolation; i++){
      xInterpolation[i]=tmp;
      LInterpolation[i]=TMath::Exp(-gr_likelihood->Eval(tmp,0,"S"));
      //LInterpolation[i]=gr_likelihood->Eval(tmp,0,"S");
      if(LInterpolation[i]>likelMAX) likelMAX = LInterpolation[i];
      std::cout << "test: " << xInterpolation[i] << "\t " << LInterpolation[i] <<"\t Thing: " << gr_likelihood->Eval(tmp,0,"S")<<  std::endl;
      
      tmp = tmp+step;

    }

    TGraph *gL = new TGraph(nInterpolation,xInterpolation,LInterpolation);
    gL->SetLineColor(2);
    gL->SetLineWidth(4);
    gL->SetMarkerColor(2);
    gL->SetMarkerStyle(0);
    gL->SetTitle(" ");
    gL->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gL->GetYaxis()->SetTitle("Likelihood");
    gL->Draw("ACsame");


    double zConvMIN = 999999999;
    double LConvoluted[5000];
    double zConv[5000]; 
    double zConvs[5000];
    double t=MIN;
    for(int i=0; i<=nInterpolation; i++){
      xInterpolation[i]=t;
      double sum=0;
      for(double TAU=MIN;TAU<=MAX; TAU=TAU+step){
        double left = gL->Eval( TAU      ,0,"S")*TMath::Gaus(t- TAU      ,0,Syst,kTRUE);
        double right= gL->Eval((TAU+step),0,"S")*TMath::Gaus(t-(TAU+step),0,Syst,kTRUE);
        sum = sum + ((left+right)/2*step);
      }

      LConvoluted[i]=(sum);
      if(LConvoluted[i]>likelMAXCONV) likelMAXCONV=LConvoluted[i];
      zConv[i]=-TMath::Log(LConvoluted[i]);
      if(zConv[i]<zConvMIN) zConvMIN=zConv[i];
      t=t+step;
    }

    for(int i=0; i<=nInterpolation; i++) zConv[i]=zConv[i]-zConvMIN;
    for(int i=0; i<=nInterpolation; i++) zConvs[i]=sqrt(2*zConv[i]);
    
    TGraph *gConv = new TGraph(nInterpolation,xInterpolation,LConvoluted);
    gConv->SetLineColor(4);
    gConv->SetLineWidth(2);
    gConv->SetLineStyle(kDashed);
    gConv->SetMarkerColor(4);
    gConv->SetMarkerStyle(0);
    gConv->SetTitle("Likelihood Convoluted");
    gConv->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gConv->GetYaxis()->SetTitle("Likelihood Convoluted");
    TGraph *gzConv = new TGraph(nInterpolation,xInterpolation,zConv);
    TGraph *gzConvsv = new TGraph(nInterpolation,xInterpolation,zConvs);
    
    
    gzConv->SetLineColor(4);
    gzConv->SetLineWidth(2);
    gzConv->SetMarkerColor(4);
    gzConv->SetMarkerStyle(0);
    gzConv->SetTitle("Likelihood Convoluted");
    gzConv->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gzConv->GetYaxis()->SetTitle("Likelihood Convoluted");
    //////////////////////////////////////////////
    gzConvsv->SetLineColor(4);
    gzConvsv->SetLineWidth(2);
    gzConvsv->SetMarkerColor(4);
    gzConvsv->SetMarkerStyle(0);
    gzConvsv->SetTitle("Likelihood Convoluted");
    gzConvsv->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gzConvsv->GetYaxis()->SetTitle("Likelihood Convoluted");
  
    TCanvas *c1 = new TCanvas("likelihood","",800,600);
    c1->Divide(2);
    c1->cd(1)->SetGrid();  
    gr_likelihood->Draw("ACP");  gzConv->Draw("CPsame");
    //lhcbpreliminary->Draw();
    c1->cd(2)->SetGrid();
    //gs->Draw("ACP");gzConvsv->Draw("CPsame");
    gL->Draw("AC");
    gConv->Draw("Csame");

    c1->Print("results/Likelihood_plot.eps");
    c1->Print("results/Likelihood_plot.pdf");
}

void SetLHCbStyle(std::string type = ""){

  if(type!="norm" && type!= "cont" && type!= "surf")
  {
    std::cout <<  "Type must be one of 1D, cont or surf. Defaulting to 1D" <<std::endl;
    type = "norm";
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
  lhcbStyle->SetOptStat(0);  
  //lhcbStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
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
  if(type == "cont"||type=="surf") lhcbStyle->SetPalette(112);
  if(type == "surf") lhcbStyle->SetTitleOffset(1.5,"Y");
  if(type == "surf") lhcbStyle->SetTitleOffset(1.5,"X");
  if(type == "surf") lhcbStyle->SetTitleOffset(0.83,"Z");
  if(type == "norm") lhcbStyle->SetPadRightMargin(0.05);
  

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
