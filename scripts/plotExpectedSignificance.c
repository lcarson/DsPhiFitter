void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotExpectedSignificance(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "B #rightarrow D_{s} #phi";
    if (mode == "DKst0") name = "B #rightarrow D K^{*0}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string sensitivityDir = "sensitivityDir/text_vals/";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> dLL_count;
    std::map<double,std::vector<double>> dLL_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("minLL")  ==std::string::npos ) continue;
        n_points++;
    }

    double* n_gen_br   = new double[n_points];
    double* n_fit_br   = new double[n_points];

    int count = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("minLL")  ==std::string::npos ) continue;
        std::cout << "Filename: " << *it << std::endl;
    
        std::vector<std::string> tokens;
        split((*it),tokens,"_");

        double br = atof(tokens[1].c_str());
        int seed  = atoi(tokens[2].c_str());
        int num   = atoi(tokens[3].c_str());


        std::ifstream input;
        input.open((sensitivityDir+*it).c_str(),std::ifstream::in);
        if(!input){
          std::cout<<" Can't find file " << sensitivityDir+*it << std::endl; return;
        }
        std::string line = "";
        getline(input,line);
        
        std::vector<std::string> tokens_2;
        split(line,tokens_2,":");
        double minLL1 = atof(tokens_2[0].c_str());
        double minLL2 = atof(tokens_2[1].c_str());

        double sigma = sqrt(-2 * (minLL1 - minLL2));
        std::cout << "MinLL 1: " << minLL1 << "\t MinLL 2: " << minLL2 << "\t sigma: " << sigma << std::endl;
        dLL_br[br].push_back(sigma);

        input.close();

        count++ ;
    }
    int n_br_vals = dLL_br.size();
    
    double* n_dLL_br   = new double[n_br_vals];
    double* n_dLL_3num = new double[n_br_vals];
    double* n_dLL_5num = new double[n_br_vals];
    double* n_dLL_3frac = new double[n_br_vals];
    double* n_dLL_5frac = new double[n_br_vals];
    int no_count = 0;
    for(std::map<double,std::vector<double>>::iterator it=dLL_br.begin();it!=dLL_br.end();it++){
      std::cout << "Br: " << it->first << "\t No points: " << (it->second).size() << std::endl;
      int above5 = 0;
      int above3 = 0;
      for(std::vector<double>::iterator itv=(it->second).begin();itv!=(it->second).end();itv++){
        if(*itv > 5.0 ) above5++;
        if(*itv > 3.0 ) above3++;
      }
      std::cout << "     Frac above 5 " << 100*(double)above5/(it->second).size() << "%"<< std::endl;
      std::cout << "     Frac above 3 " << 100*(double)above3/(it->second).size() << "%"<< std::endl;
      
      n_dLL_br[no_count]=(it->first);
      n_dLL_3frac[no_count]=((double)above3/(it->second).size());
      n_dLL_5frac[no_count]=((double)above5/(it->second).size());

      no_count++;
    }

    SetLHCbStyle("norm");

    TCanvas *cav = new TCanvas("Canvas","Sensitivity",200,10,700,500);
  
    TGraph* gr_sensitivity_line = new TGraph(n_br_vals, n_dLL_br, n_dLL_5frac );
    
    TGraph* gr_sensitivity_line2 = new TGraph(n_br_vals, n_dLL_br, n_dLL_3frac);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Significance %s;Generated Br(%s) (#times10^{-7});Fraction of toys with given #sigma",name.c_str(),name.c_str()));
    mg->Add(gr_sensitivity_line,"PL");
    mg->Add(gr_sensitivity_line2,"PL");
    mg->Draw("A");
    gr_sensitivity_line->SetLineColor(kRed);
    gr_sensitivity_line2->SetLineColor(kGreen);

    TLegend* leg = new TLegend(0.65,0.2,0.8,0.4);
    leg->AddEntry(gr_sensitivity_line,"5 #sigma","l");
    leg->AddEntry(gr_sensitivity_line2,"3 #sigma","l");
    leg->Draw();

    TLine *line = new TLine(18.7,0, 18.7,1.05);
    line->SetLineColor(kBlue+2);
    line->SetLineWidth(2.0);
    line->Draw();

    cav->Print("results/ExpectedSignificance_plot.eps");
    cav->Print("results/ExpectedSignificance_plot.pdf");

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

    cav->Print("results/Sensitivity_plot.eps");
    cav->Print("results/Sensitivity_plot.pdf");
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
