void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotSensitivity(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "B #rightarrow D_{s} #phi";
    if (mode == "DKst0") name = "B #rightarrow D K^{*0}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string sensitivityDir = "sensitivityDir/";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> gen_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("sensitivity")  ==std::string::npos ) continue;
        n_points++;
    }

    double* n_gen_br   = new double[n_points];
    double* n_fit_br   = new double[n_points];

    int count = 0;
    int nFits=0;
    int nBad=0;
    int nFPD = 0;
    int nBadMINOS = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("sensitivity")  ==std::string::npos ) continue;
        std::cout << "Filename: " << *it << std::endl;
    
        std::vector<std::string> tokens;
        split((*it),tokens,"_");

        double br = atof(tokens[1].c_str());
        int seed  = atoi(tokens[2].c_str());
        int num   = atoi(tokens[3].c_str());


        TFile newfile((sensitivityDir+*it).c_str(),"READ");
        if( !(newfile.GetListOfKeys()->Contains("toy_withSig"))  ) continue; 
        RooFitResult* result = (RooFitResult*)newfile.Get("toy_withSig");
        RooRealVar* param_final = (RooRealVar*) (result->floatParsFinal()).find("Branching_fraction");
        double fitted_br = param_final->getVal();
        std::cout << "Br: " << br << " Fitted Br: " << fitted_br << " seed: " << seed << " Number: " << num << std::endl;
        n_gen_br[count] = br;
        n_fit_br[count] = fitted_br;
        count++ ;
        gen_count[br]++;
        gen_br[br].push_back(fitted_br);


        /*

        TIter nextkey(newfile.GetListOfKeys());
        TKey *key;
        while((key=(TKey*)nextkey())){
          if(std::string(key->GetClassName())!=std::string("RooFitResult")) continue;
          RooFitResult *r = (RooFitResult*)key->ReadObj();
          //if(std::string(r->GetName()) != std::string("toy_withSig")) continue;
          //if( !newfile->GetListOfKeys()->Contains("toy_withSig"))  ) continue; 
          std::cout << "RooFitResult Name: " << std::endl;
          r->Print();
          nFits++;      
          if(r->covQual()<2){nBad++;continue;}
          if(r->covQual()<3){nFPD++;continue;}
          if(r->status()!=0){nBadMINOS++;}//continue;}
          
          RooRealVar* param_final = (RooRealVar*) (r->floatParsFinal()).find("Branching_fraction");
          double fitted_br = param_final->getVal();

        

        
        TFile *newfile = TFile::Open((sensitivityDir+*it).c_str());
        if(!newfile){
          std::cout << "File not found: " << sensitivityDir+*it << std::endl;
          continue;
        }
        RooFitResult* result = (RooFitResult*)newfile->Get("toy_withSig");
        RooRealVar* param_final = (RooRealVar*) (result->floatParsFinal()).find("Branching_fraction");
        double fitted_br = param_final->getVal();
        

          std::cout << "Br: " << br << " Fitted Br: " << fitted_br << " seed: " << seed << " Number: " << num << std::endl;
          n_gen_br[count] = br;
          n_fit_br[count] = fitted_br;
          count++ ;
          gen_count[br]++;
          gen_br[br].push_back(fitted_br);
          //newfile->Close();
        }
      */

    }

    int n_gen_points = gen_br.size();
    std::cout << "Gen points: " << n_gen_points <<std::endl;
    double* n_gen_line_br   = new double[n_gen_points];
    double* n_fit_line_br   = new double[n_gen_points];
    
    double* n_gen_line2_br   = new double[n_gen_points];
    double* n_fit_line2_br   = new double[n_gen_points];
    
    double* n_gen_line_br_lower   = new double[n_gen_points];
    double* n_fit_line_br_lower   = new double[n_gen_points];
    
    double* n_gen_line2_br_lower   = new double[n_gen_points];
    double* n_fit_line2_br_lower   = new double[n_gen_points];
    count = 0;
    for(std::map<double,std::vector<double>>::iterator it=gen_br.begin();it!=gen_br.end();it++){ 
        std::vector<double> temp = it->second;
        double gen_value = it->first;
        std::sort (temp.begin(), temp.end());

        std::cout << "Gen BR: " << gen_value << " Length: " << temp.size()<< " Limit point: " << temp.size()*0.95 << " Rounded: " << floor(temp.size()*0.95 + 0.5) << " Value: " << temp[floor(temp.size()*0.95 + 0.5)-1] << std::endl;
        
        //std::cout << "Values: [" ; for(std::vector<double>::iterator itv=temp.begin();itv!=temp.end();itv++){std::cout << *itv <<", ";} std::cout << "]"<< std::endl;
        
        n_gen_line_br[count] = gen_value;
        n_fit_line_br[count] = temp[floor(temp.size()*0.95 + 0.5)-1];

        n_gen_line2_br[count] = gen_value;
        n_fit_line2_br[count] = temp[floor(temp.size()*0.90 + 0.5)-1];
        
        n_gen_line_br_lower[count] = gen_value;
        n_fit_line_br_lower[count] = temp[floor(temp.size()*0.05 + 0.5)-1];

        n_gen_line2_br_lower[count] = gen_value;
        n_fit_line2_br_lower[count] = temp[floor(temp.size()*0.10 + 0.5)-1]; 

        //std::cout << "n_gen_line_br[" << count << "] = " << n_gen_line_br[count] <<std::endl;
        //std::cout << "n_fit_line_br[" << count << "] = " << n_fit_line_br[count] <<std::endl;
        count++;
    
    }

    SetLHCbStyle("norm");

    TCanvas *cav = new TCanvas("Canvas","Sensitivity",200,10,700,500);
  
    TGraph* gr_sensitivity = new TGraph(n_points, n_fit_br, n_gen_br);

    gr_sensitivity->SetTitle(Form("Sensitivity %s;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})",name.c_str())); 
    //gr_sensitivity->Draw("AP*");
    
    //TCanvas *cav2 = new TCanvas("Canvas2","Sensitivity2",200,10,700,500);
  
    TGraph* gr_sensitivity_line = new TGraph(n_gen_points, n_fit_line_br, n_gen_line_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line2 = new TGraph(n_gen_points, n_fit_line2_br, n_gen_line2_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
  
    TGraph* gr_sensitivity_line_lower = new TGraph(n_gen_points, n_fit_line_br_lower, n_gen_line_br_lower);
    gr_sensitivity_line_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line2_lower = new TGraph(n_gen_points, n_fit_line2_br_lower, n_gen_line2_br_lower);
    gr_sensitivity_line_lower->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    

    TGraph *grshade = new TGraph(2*n_gen_points);
    for (i=0;i<n_gen_points;i++) {
      grshade->SetPoint(i,n_fit_line_br[i],n_gen_line_br[i]);
      grshade->SetPoint(n_gen_points+i,n_fit_line_br_lower[n_gen_points-i-1],n_gen_line_br_lower[n_gen_points-i-1]);
    }
    grshade->SetFillColor(kYellow);
    TGraph *grshade2 = new TGraph(2*n_gen_points);
    for (i=0;i<n_gen_points;i++) {
      grshade2->SetPoint(i,n_fit_line2_br[i],n_gen_line2_br[i]);
      grshade2->SetPoint(n_gen_points+i,n_fit_line2_br_lower[n_gen_points-i-1],n_gen_line2_br_lower[n_gen_points-i-1]);
    }

    grshade2->SetFillColor(kGreen);
/*
   //grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("l");
   gr_sensitivity_line->Draw("l");
   gr_sensitivity_line_lower->Draw("l");
*/
    //gr_sensitivity_line->Draw("AL");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Sensitivity %s;Fitted Br(%s) (#times10^{-7});Generated Br(%s) (#times10^{-7})",name.c_str(),name.c_str(),name.c_str()));

    mg->Add(grshade,"F");
    mg->Add(grshade2,"F");
    mg->Add(gr_sensitivity_line,"PL");
    mg->Add(gr_sensitivity_line2,"PL");
    mg->Add(gr_sensitivity_line_lower,"PL");
    mg->Add(gr_sensitivity_line2_lower,"PL");
    mg->Add(gr_sensitivity,"p");
    mg->Draw("A");
    

    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kYellow);
    gr_sensitivity_line_lower->SetLineColor(kYellow);

    gr_sensitivity_line2->SetLineColor(kGreen);
    gr_sensitivity_line2_lower->SetLineColor(kGreen);

    /*
    
    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kYellow);
    gr_sensitivity_line->SetFillColor(kYellow);
    gr_sensitivity_line_lower->SetFillColor(0);
    gr_sensitivity_line->Draw("AF");
    gr_sensitivity_line_lower->Draw("F");

    gr_sensitivity_line2->SetLineColor(kGreen);
    gr_sensitivity_line2_lower->SetLineColor(kGreen);

*/

    cav->Print("results/Sensitivity_plot.eps");
    cav->Print("results/Sensitivity_plot.pdf");
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
