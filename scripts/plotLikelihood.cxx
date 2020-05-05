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

    double minx = 1.167895512; 
    double miny = 37002.25333;      // value and NLL of the approved result

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

    TCanvas *cav = new TCanvas("Canvas","Likelihood",200,10,700,500);
  
    TGraph* gr_loglikelihood = new TGraph(count, n_br_sorted, n_lik_sorted);

    gr_loglikelihood->SetTitle(";#it{B(B^{+}#rightarrow D_{s}^{+}#phi)} #times 10^{-7} ;#Delta(log(Likelihood))"); 
    //gr_loglikelihood->Draw("APC");


    // Quick checks
    
    double* likelihood_temp    = new double[n_points];
    for(int i=0; i<count; i++){
      
      likelihood_temp[i]=exp(-n_lik_sorted[i]);
      std::cout << n_br_sorted[i]<< "\t"<< likelihood_temp[i] << std::endl;
    }

    TGraph* gr_likelihood_temp = new TGraph(count, n_br_sorted, likelihood_temp);
    //gr_likelihood_temp->Draw("APC");




    double MIN  = -10.0;         // minimum value of the scanned parameter
    double MAX  =  20.0;         // maximum value
    double Syst = 0.788;         // total systematic error

    double* raw_significance    = new double[n_points];
    for(int i=0; i<count; i++) raw_significance[i]=sqrt(2*n_lik_sorted[i]);
  
    TGraph* gr_raw_significance = new TGraph(count, n_br_sorted, raw_significance);

    gr_raw_significance->SetTitle(";#it{B(B^{+}#rightarrow D_{s}^{+}#phi)} #times 10^{-7} ;#sqrt{2#Delta log(L)}"); 
    //gr_raw_significance->Draw("APC"); //REMOVE




    double likelMAX     = -99999;
    double likelMAXCONV = -99999;
    double LInterpolation[5000];
    double xInterpolation[5000];
    int    nInterpolation = 500;
    double step = (MAX-MIN)/nInterpolation;
    double tmp = MIN;
    double deltaLogLikelihood[5000]; 
    double significance[5000];
    
    //double LLMIN = 999999999;


    for(int i=0; i<=nInterpolation; i++){
      xInterpolation[i]=tmp;
      LInterpolation[i]=TMath::Exp(-gr_loglikelihood->Eval(tmp,0,"S"));
      //LInterpolation[i]=gr_loglikelihood->Eval(tmp,0,"S");
      if(LInterpolation[i]>likelMAX) likelMAX = LInterpolation[i];
      //if(LInterpolation[i]<LLMIN) LLMIN = LInterpolation[i];
      deltaLogLikelihood[i]=-TMath::Log(LInterpolation[i]);
      //std::cout << "test: " << xInterpolation[i] << "\t " << LInterpolation[i] <<"\t Thing: " << gr_loglikelihood->Eval(tmp,0,"S")<< "\t Thing test : " << TMath::Exp(-gr_loglikelihood->Eval(tmp,0,"S")) << std::endl; 
      std::cout << "delta LL: " << deltaLogLikelihood[i] <<std::endl;
      tmp = tmp+step;

    }

    
    for(int i=0; i<=nInterpolation; i++) significance[i]=sqrt(2*deltaLogLikelihood[i]);
    //for(int i=0; i<=nInterpolation; i++) LInterpolation[i] = LInterpolation[i]-LLMIN; 



    TGraph *gr_likelihood = new TGraph(nInterpolation,xInterpolation,LInterpolation);
    gr_likelihood->SetLineColor(kBlack);
    gr_likelihood->SetLineWidth(2);
    gr_likelihood->SetMarkerColor(2);
    gr_likelihood->SetMarkerStyle(0);
    gr_likelihood->SetTitle(" ");
    gr_likelihood->GetXaxis()->SetTitle("#it{B(B^{+}#rightarrow D_{s}^{+}#phi)} #times 10^{-7}");
    gr_likelihood->GetYaxis()->SetTitle("Likelihood");
    gr_likelihood->Draw("ACsame");


    double zConvMIN = 999999999;
    double LConvoluted[5000];
    double zConv[5000]; 
    double zConvs[5000];
    double t=MIN;
    for(int i=0; i<=nInterpolation; i++){
      xInterpolation[i]=t;
      double sum=0;
      for(double TAU=MIN;TAU<=MAX; TAU=TAU+step){
        double left = gr_likelihood->Eval( TAU      ,0,"S")*TMath::Gaus(t- TAU      ,0,Syst,kTRUE);
        double right= gr_likelihood->Eval((TAU+step),0,"S")*TMath::Gaus(t-(TAU+step),0,Syst,kTRUE);
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
    
    TGraph *gr_likelihood_Conv = new TGraph(nInterpolation,xInterpolation,LConvoluted);
    gr_likelihood_Conv->SetLineColor(kBlue);
    gr_likelihood_Conv->SetLineWidth(2);
    gr_likelihood_Conv->SetMarkerColor(4);
    gr_likelihood_Conv->SetMarkerStyle(0);
    gr_likelihood_Conv->SetTitle("Likelihood Convoluted");
    gr_likelihood_Conv->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gr_likelihood_Conv->GetYaxis()->SetTitle("Likelihood Convoluted");
    gr_likelihood_Conv->Draw("Csame");

    TGraph *gr_loglikelihood_Conv   = new TGraph(nInterpolation,xInterpolation,zConv);
    TGraph *gr_Conv_significance    = new TGraph(nInterpolation,xInterpolation,zConvs);
    TGraph *gr_significance         = new TGraph(nInterpolation,xInterpolation,significance);
    
    
    gr_loglikelihood_Conv->SetLineColor(4);
    gr_loglikelihood_Conv->SetLineWidth(2);
    gr_loglikelihood_Conv->SetMarkerColor(4);
    gr_loglikelihood_Conv->SetMarkerStyle(0);
    gr_loglikelihood_Conv->SetTitle("Likelihood Convoluted");
    gr_loglikelihood_Conv->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gr_loglikelihood_Conv->GetYaxis()->SetTitle("Likelihood Convoluted");

    //////////////////////////////////////////////
    gr_Conv_significance->SetLineColor(4);
    gr_Conv_significance->SetLineWidth(2);
    gr_Conv_significance->SetMarkerColor(4);
    gr_Conv_significance->SetMarkerStyle(0);
    gr_Conv_significance->SetTitle("Significance Convoluted");
    gr_Conv_significance->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gr_Conv_significance->GetYaxis()->SetTitle("Significance Convoluted");
    //////////////////////////////////////////////
    gr_significance->SetLineColor(4);
    gr_significance->SetLineWidth(2);
    gr_significance->SetMarkerColor(4);
    gr_significance->SetMarkerStyle(0);
    gr_significance->SetTitle("#sqrt{2#times#Delta log(L)}");
    gr_significance->GetXaxis()->SetTitle("Branching Fraction (#times 10^{-7})");
    gr_significance->GetYaxis()->SetTitle("#sqrt{2#times#Delta log(L)}");
  
    TCanvas *c1 = new TCanvas("likelihood","",900,500);
    c1->Divide(2);
    c1->cd(1)->SetGrid();  
    gr_loglikelihood->Draw("ACP");  gr_loglikelihood_Conv->Draw("CPsame");
    gr_loglikelihood->GetXaxis()->SetRangeUser(MIN,MAX);
    
    //lhcbpreliminary->Draw();
    c1->cd(2)->SetGrid();
    //gs
    //gr_loglikelihood->Draw("ACP");gr_Conv_significance->Draw("CPsame");
    gr_raw_significance->Draw("ACP"); gr_Conv_significance->Draw("CPsame");
    gr_raw_significance->GetXaxis()->SetRangeUser(MIN,MAX);
    //gr_likelihood->Draw("AC"); gr_likelihood_Conv->Draw("Csame");

    c1->Print("results/Likelihood_plot.eps");
    c1->Print("results/Likelihood_plot.pdf");


    TCanvas *c_integrals = new TCanvas("likelihood_integrals","",900,500);

    // Integrate TGraphs
    TF1 likelihood_raw_function( "likelihood_raw_function", [&](double *xInterpolation, double *){ return gr_likelihood->Eval(     xInterpolation[0]); },MIN,MAX,0);
    TF1 likelihood_Conv_function("likelihood_Conv_function",[&](double *xInterpolation, double *){ return gr_likelihood_Conv->Eval(xInterpolation[0]); },MIN,MAX,0);
    
    TGraph *g_integral_raw  = (TGraph*)likelihood_raw_function.DrawIntegral();
    TGraph *g_integral_Conv = (TGraph*)likelihood_Conv_function.DrawIntegral();


    double integral_raw_range = likelihood_raw_function.Integral(0.0,4.8);
    double integral_raw_total = likelihood_raw_function.Integral(0.0,20.0);
    std::cout << "Fraction raw:  " << integral_raw_range/integral_raw_total <<std::endl;

    double integral_Conv_range = likelihood_Conv_function.Integral(0.0,4.8);
    double integral_Conv_total = likelihood_Conv_function.Integral(0.0,20.0);
    std::cout << "Fraction Conv: " << integral_Conv_range/integral_Conv_total <<std::endl;

    std::cout << "Value at Zero: " << g_integral_raw->Eval(0.0) <<std::endl;
    std::cout << "Value at Zero: " << g_integral_Conv->Eval(0.0) <<std::endl;
    
    // Subtract offset 
    for (int i=0;i<g_integral_raw->GetN(); i++) g_integral_raw->GetY()[i]  = g_integral_raw->GetY()[i]  - g_integral_raw->Eval(0.0);
    for (int i=0;i<g_integral_Conv->GetN();i++) g_integral_Conv->GetY()[i] = g_integral_Conv->GetY()[i] - g_integral_Conv->Eval(0.0);
    

    for (int i=0;i<g_integral_raw->GetN(); i++) g_integral_raw->GetY()[i]  *= 1/g_integral_raw->Eval(18.0);
    for (int i=0;i<g_integral_Conv->GetN();i++) g_integral_Conv->GetY()[i] *= 1/g_integral_Conv->Eval(18.0);
    

    g_integral_raw->SetLineColor(kRed);
    g_integral_Conv->SetLineColor(kBlue);
    g_integral_Conv->SetLineWidth(1);
    g_integral_raw->SetLineWidth(1);

    g_integral_raw->GetXaxis()->SetRangeUser(0.0,5.0);
    g_integral_Conv->GetXaxis()->SetRangeUser(0.0,5.0);
    g_integral_raw->GetYaxis()->SetRangeUser(0.0,1.0);
    g_integral_Conv->GetYaxis()->SetRangeUser(0.0,1.0);

    g_integral_raw->Draw("ACP");
    g_integral_Conv->Draw("sameCP");

    g_integral_raw->GetXaxis()->SetRangeUser(MIN,MAX);
    g_integral_Conv->GetXaxis()->SetRangeUser(MIN,MAX);


    TLine *line_95 = new TLine(MIN,0.95,MAX,0.95);
    line_95->Draw();

    double before_95 = 0.0;
    double after_95  = 0.0;
    bool found_95 = false;
    for (int i=0;i<g_integral_raw->GetN(); i++){
      //std::cout << "("<< g_integral_raw->GetX()[i]<< ","<<g_integral_raw->GetY()[i] << ")"<<std::endl;

      if(g_integral_raw->GetY()[i]>0.95 && !found_95){
        before_95 = g_integral_raw->GetX()[i-1];
        after_95 =  g_integral_raw->GetX()[i];
        found_95 = true;
      }
    }

    std::cout << "RAW: " << before_95 << "\t" << after_95 << std::endl;

    double range = after_95- before_95;

    double before_interp_95 = 0.0;
    double after_interp_95  = 0.0;
    bool found_interp_95 = false;
    int No_points = 100;

    for(int i=0; i<No_points; i++ ){
      double x = before_95 + i*(range / No_points);

      if(g_integral_raw->Eval(x)>0.95 && !found_interp_95){
        before_interp_95 = x - range / No_points;
        after_interp_95 =  x;
        found_interp_95 = true;
      }

    }
    std::cout << "RAW: " << before_interp_95 << "\t" << after_interp_95 << std::endl;

    //---------------------
    before_95 = 0.0;
    after_95  = 0.0;
    found_95 = false;
    for (int i=0;i<g_integral_Conv->GetN(); i++){
      //std::cout << "("<< g_integral_Conv->GetX()[i]<< ","<<g_integral_Conv->GetY()[i] << ")"<<std::endl;

      if(g_integral_Conv->GetY()[i]>0.95 && !found_95){
        before_95 = g_integral_Conv->GetX()[i-1];
        after_95 =  g_integral_Conv->GetX()[i];
        found_95 = true;
      }
    }

    std::cout << "CONV: " << before_95 << "\t" << after_95 << std::endl;

    range = after_95- before_95;

    before_interp_95 = 0.0;
    after_interp_95  = 0.0;
    found_interp_95 = false;
    No_points = 100;

    for(int i=0; i<No_points; i++ ){
      double x = before_95 + i*(range / No_points);

      if(g_integral_Conv->Eval(x)>0.95 && !found_interp_95){
        before_interp_95 = x - range / No_points;
        after_interp_95 =  x;
        found_interp_95 = true;
      }

    }
    std::cout << "CONV: " << before_interp_95 << "\t" << after_interp_95 << std::endl;


    TCanvas *c_limits = new TCanvas("limits","Plots I want",900,500);
    c_limits->Divide(2);
    c_limits->cd(1);
    gr_loglikelihood->SetLineColor(kBlack);
    gr_loglikelihood_Conv->SetLineColor(kBlue);

    gr_loglikelihood->Draw("AC");  
    gr_loglikelihood_Conv->Draw("Csame");
    gr_loglikelihood->GetXaxis()->SetRangeUser(-5.0,10.0);
    gr_loglikelihood->GetYaxis()->SetRangeUser(0.0,20.0);
    c_limits->cd(2);

    gr_likelihood->GetXaxis()->SetRangeUser(-5.0,10.0);
    gr_likelihood_Conv->GetXaxis()->SetRangeUser(-5.0,10.0);

    gr_likelihood->SetLineColor(kBlack);
    gr_likelihood_Conv->SetLineColor(kBlue);

    gr_likelihood->Draw("AC");
    //gr_likelihood_Conv->Draw("Csame");

    int shaded_points = 41+3;
    TGraph *grshade = new TGraph(shaded_points);
    grshade->SetPoint(0,0.0,0.0);
    
    for (int i=0;i<(shaded_points-1);i++) {
      double x = i*0.1;
      grshade->SetPoint(i+1,x,gr_likelihood->Eval(x));
    }

    grshade->SetPoint(shaded_points-1,4.1,0.0);

    grshade->SetFillColor(kBlack);
    grshade->SetFillStyle(3005);
    grshade->SetFillStyle(3002);
    //grshade->SetFillColor(kGray);
    //grshade->SetFillStyle(1001);

    //grshade->Draw("Fsame");




    int shaded_points_Conv = 44+3;
    TGraph *grshade_Conv = new TGraph(shaded_points_Conv);
    grshade_Conv->SetPoint(0,0.0,0.0);
    
    for (int i=0;i<(shaded_points_Conv-1);i++) {
      double x = i*0.1;
      grshade_Conv->SetPoint(i+1,x,gr_likelihood_Conv->Eval(x));
    }

    grshade_Conv->SetPoint(shaded_points_Conv-1,4.4,0.0);

    //grshade_Conv->SetFillColor(kBlue);
    //grshade_Conv->SetFillStyle(3004);
    grshade_Conv->SetFillColor(kAzure-4);
    grshade_Conv->SetFillStyle(1001);
    grshade_Conv->Draw("Fsame");
    grshade->Draw("Fsame");

    gr_likelihood->Draw("Csame");
    gr_likelihood_Conv->Draw("Csame");

    gPad->RedrawAxis();

    c_limits->Print("results/Likelihood_limits.eps");
    c_limits->Print("results/Likelihood_limits.pdf");
    c_limits->Print("results/Likelihood_limits.png");
    c_limits->Print("results/Likelihood_limits.root");
    c_limits->Print("results/Likelihood_limits.C");

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
