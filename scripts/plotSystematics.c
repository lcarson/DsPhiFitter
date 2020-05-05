void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotSystematics(std::string mode = "DsPhi",double low = 0.0, double high = 30.0){


    std::string name;
    if (mode == "DsPhi") name = "B #rightarrow D_{s} #phi";
    if (mode == "DKst0") name = "B #rightarrow D K^{*0}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string sensitivityDir = "systematicsDir/";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> gen_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("vary")  ==std::string::npos ) continue;
        n_points++;
        std::cout << "Found file: " << *it <<std::endl;
    }

    double* n_gen_br   = new double[n_points];
    double* n_fit_br   = new double[n_points];

    TH1D* hist = new TH1D("Branching_fraction", "", 100 , low ,high );

    int count = 0;
    int nBad = 0;
    int nFPD = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("vary")  ==std::string::npos ) continue;
        std::cout << "Filename: " << *it << std::endl;
    
        //std::vector<std::string> tokens;
        //split((*it),tokens,"_");

        //double br = atof(tokens[1].c_str());
        //int seed  = atoi(tokens[2].c_str());
        //int num   = atoi(tokens[3].c_str());


        TFile *newfile = TFile::Open((sensitivityDir+*it).c_str());
        RooFitResult* result = (RooFitResult*)newfile->Get("fitresult_model_reducedData");

        if(result->covQual()<2){nBad++;continue;}
        if(result->covQual()<3){nFPD++;continue;}
        RooRealVar* param_final = (RooRealVar*) (result->floatParsFinal()).find("Branching_fraction");
        double fitted_br = param_final->getVal();


        std::cout << " Fitted Br: " << fitted_br << std::endl;
        //n_gen_br[count] = br;
        n_fit_br[count] = fitted_br;
        count++ ;
        //gen_count[br]++;
        //gen_br[br].push_back(fitted_br);
        hist->Fill(fitted_br);
        newfile->Close();
    }
    /*
    int n_gen_points = gen_br.size();
    std::cout << "Gen points: " << n_gen_points <<std::endl;
    double* n_gen_line_br   = new double[n_gen_points];
    double* n_fit_line_br   = new double[n_gen_points];
    
    double* n_gen_line2_br   = new double[n_gen_points];
    double* n_fit_line2_br   = new double[n_gen_points];
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

        //std::cout << "n_gen_line_br[" << count << "] = " << n_gen_line_br[count] <<std::endl;
        //std::cout << "n_fit_line_br[" << count << "] = " << n_fit_line_br[count] <<std::endl;
        count++;
    
    }
    */

    SetLHCbStyle("norm");

    TCanvas *cav = new TCanvas("Canvas","Systematics",200,10,700,500);
    hist->Draw();

    TF1 *gaussian = new TF1("Gaussian","gaus",0,30);
    hist->Fit("Gaussian","Q0R");
    double mean  = gaussian->GetParameter(1);
    double sigma = gaussian->GetParameter(2);
    const double* err=gaussian->GetParErrors();
    double meanerr=err[1];
    double sigmaerr=err[2];
    gaussian->SetLineColor(kBlue); gaussian->Draw("same");

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
    
    /*
    TGraph* gr_sensitivity = new TGraph(n_points, n_fit_br, n_gen_br);

    gr_sensitivity->SetTitle(Form("Sensitivity %s;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})",name.c_str())); 
    //gr_sensitivity->Draw("AP*");
    
    //TCanvas *cav2 = new TCanvas("Canvas2","Sensitivity2",200,10,700,500);
  
    TGraph* gr_sensitivity_line = new TGraph(n_gen_points, n_fit_line_br, n_gen_line_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");  
    
    TGraph* gr_sensitivity_line2 = new TGraph(n_gen_points, n_fit_line2_br, n_gen_line2_br);
    gr_sensitivity_line->SetTitle("Sensitivity;Fitted Branching Fraction (#times10^{-7});Generated Branching Fraction (#times10^{-7})"); 
    //gr_sensitivity_line->Draw("AL");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Sensitivity %s;Fitted Br(%s) (#times10^{-7});Generated Br(%s) (#times10^{-7})",name.c_str(),name.c_str(),name.c_str()));
    mg->Add(gr_sensitivity,"p");
    mg->Add(gr_sensitivity_line,"PL");
    mg->Add(gr_sensitivity_line2,"PL");
    mg->Draw("A");
    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kRed);
    gr_sensitivity_line2->SetLineColor(kGreen);

    cav->Print("results/Systematics_plot.eps");
    cav->Print("results/Systematics_plot.pdf");
    */
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
