void SetLHCbStyle(std::string);
#include "tools.h"
#include <stdlib.h> 
// ===========================================================================

void plotLimitFromToy(std::string mode = "DsPhi"){


    std::string name;
    if (mode == "DsPhi") name = "B #rightarrow D_{s} #phi";
    if (mode == "DKst0") name = "B #rightarrow D K^{*0}";

    std::cout << "Using mode name: " << mode << " Name: " << name << std::endl;

    std::string sensitivityDir = "sensitivityDir/likelihood";
    std::vector<std::string> fileName;
    getdir(sensitivityDir,fileName);
    std::map<double,int> gen_count;
    std::map<double,std::vector<double>> gen_br;

    int n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("text_")  ==std::string::npos ) continue;
        n_points++;
    }

    double* n_gen_br    = new double[n_points];
    double* fraction_br = new double[n_points];

    int count = 0;
    int nFits=0;
    int nBad=0;
    int nFPD = 0;
    int nBadMINOS = 0;
    n_points = 0;
    for(std::vector<std::string>::iterator it=fileName.begin();it!=fileName.end();it++){
        if( (*it).find("text_")  ==std::string::npos ) continue;
    
        std::vector<std::string> tokens;
        //split((*it),tokens,"_");
        std::ifstream input;
        std::string fullname = sensitivityDir +"/" + (*it);
        input.open(fullname.c_str(),std::ifstream::in);
        std::string line = "";
        getline(input,line);
        split(line,tokens,":");
        if(tokens[3] == "nan" || tokens[3] == "-nan"){
          std::cout <<"Found bad value: " << (*it) << " with line: " << line <<std::endl; 
          continue;
        }

        double br        = atof(tokens[0].c_str());
        double fraction  = atof(tokens[3].c_str());

        // Check if toy was good
        tokens.clear();
        split((*it),tokens,"_");
        // text_likelihood_9.500000_70199_96_.txt
        //sensitivity_9.500000_70199_79_.root
        std::string result_filename = Form("sensitivityDir/sensitivity_%s_%s_%s_.root",
                                           tokens[2].c_str(),
                                           tokens[3].c_str(),
                                           tokens[4].c_str()
                                           );
        TFile newfile(result_filename.c_str(),"READ");
        if( !(newfile.GetListOfKeys()->Contains("toy_withSig"))  ) continue; 
        RooFitResult* r = (RooFitResult*)newfile.Get("toy_withSig");
        
        if(count%1000 == 0) std::cout << "Filename: " << *it << std::endl;
        if(count%1000 == 0) std::cout << "Finding results file: " << result_filename <<std::endl;
        if(count%1000 == 0) std::cout << "Result covQual: " << r->covQual() <<std::endl;
        if(r->covQual()<2){continue;}
        if(r->covQual()<3){continue;}
        if(fraction < 0 || fraction >1 || br<0 || br > 10){
          std::cout <<"Found bad value: " << fraction << " at br:" << br <<std::endl; 
          continue;
        } 

        n_gen_br[count]    = br;
        fraction_br[count] = fraction;
        count++ ;
        gen_count[br]++;
        gen_br[br].push_back(fraction);
        n_points++;

    }

    int n_gen_points = gen_br.size();
    std::cout << "Gen points: " << n_gen_points <<std::endl;
   
    double* n_gen_br_mean   = new double[n_gen_points];
    double* fraction_br_mean   = new double[n_gen_points];

    count = 0;
    for(std::map<double,std::vector<double>>::iterator it=gen_br.begin();it!=gen_br.end();it++){ 
        std::vector<double> temp = it->second;
        double gen_value = it->first;
        std::sort (temp.begin(), temp.end());

        double average = accumulate( temp.begin(), temp.end(), 0.0 )/ temp.size();


        TCanvas *cav = new TCanvas(Form("Canvas_%d",count),"Sensitivity",200,10,700,500);
        TH1D *h_temp = new TH1D(Form("h_%d",count), "", 100,0.0,1.0);
        for(int i = 0; i < temp.size(); i++){
          h_temp->Fill(temp[i]);
        }
        h_temp->Draw();


        n_gen_br_mean[count]    = gen_value;
        //fraction_br_mean[count] = average;
        fraction_br_mean[count] = h_temp->GetMean();

        std::cout << "Gen BR: " << gen_value;
        std::cout << " Length: " << temp.size();
        std::cout << " Middle?: " << (int)(temp.size()/2.0);
        std::cout << " Mean: "<< fraction_br_mean[count];
        std::cout << " Median: "<< temp[(int)(temp.size()/2.0)];
        std::cout<<std::endl;

        /*
        n_gen_line_br[count] = gen_value;
        n_fit_line_br[count] = temp[floor(temp.size()*0.95 + 0.5)-1];

        n_gen_line2_br[count] = gen_value;
        n_fit_line2_br[count] = temp[floor(temp.size()*0.683 + 0.5)-1];
        
        n_gen_line_br_lower[count] = gen_value;
        n_fit_line_br_lower[count] = temp[floor(temp.size()*0.05 + 0.5)-1];

        n_gen_line2_br_lower[count] = gen_value;
        n_fit_line2_br_lower[count] = temp[floor(temp.size()*0.317 + 0.5)-1]; 
*/
        //std::cout << "n_gen_line_br[" << count << "] = " << n_gen_line_br[count] <<std::endl;
        //std::cout << "n_fit_line_br[" << count << "] = " << n_fit_line_br[count] <<std::endl;
        count++;
    
    }

    SetLHCbStyle("norm");

    TCanvas *cav = new TCanvas("Canvas","Sensitivity",200,10,700,500);
  
    TGraph* gr_sensitivity = new TGraph(n_points,  n_gen_br,fraction_br);

    gr_sensitivity->SetTitle(Form("Sensitivity %s;Generated Branching Fraction (#times10^{-7});Fraction below fit value",name.c_str())); 
    gr_sensitivity->Draw("AP*");
    
    TGraph* gr_mean = new TGraph(n_gen_points,  n_gen_br_mean,fraction_br_mean);

    gr_mean->SetTitle(Form("Sensitivity %s;Generated Branching Fraction (#times10^{-7});Fraction below fit value",name.c_str())); 
    gr_mean->SetLineColor(kBlue);
    gr_mean->SetMarkerColor(kBlue);
    gr_mean->Draw("LPsame");
    //gr_mean->Draw("APC");

    TLine *line = new TLine(0.0,0.05,10.0,0.05);
    line->SetLineColor(kRed);
    line->Draw();

    cav->Print("results/limit_from_toys.eps");
    cav->Print("results/limit_from_toys.pdf");

    /*
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
    */
/*
   //grshade->SetFillStyle(3013);
   grshade->SetFillColor(16);
   grshade->Draw("l");
   gr_sensitivity_line->Draw("l");
   gr_sensitivity_line_lower->Draw("l");
*/
    //gr_sensitivity_line->Draw("AL");
/*
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(Form("Sensitivity %s;Fitted Br(%s) (#times10^{-7});Generated Br(%s) (#times10^{-7})",name.c_str(),name.c_str(),name.c_str()));

    mg->Add(grshade,"F");
    mg->Add(grshade2,"F");
    mg->Add(gr_sensitivity_line,"PL");
    mg->Add(gr_sensitivity_line2,"PL");
    mg->Add(gr_sensitivity_line_lower,"PL");
    mg->Add(gr_sensitivity_line2_lower,"PL");
    //mg->Add(gr_sensitivity,"p");
    mg->Draw("A");
    

    gr_sensitivity->SetMarkerStyle(2);
    gr_sensitivity_line->SetLineColor(kYellow);
    gr_sensitivity_line_lower->SetLineColor(kYellow);

    gr_sensitivity_line2->SetLineColor(kGreen);
    gr_sensitivity_line2_lower->SetLineColor(kGreen);
*/
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
/*
    double limit = 4.215;

    TLine *line = new TLine(1.16,0.0,1.16,limit);
    line->SetLineWidth(1);
    line->SetLineColor(kRed);
    line->Draw();    

    TLine *line2 = new TLine(-4.0,limit,1.16,limit);
    line2->SetLineStyle(kDashed);
    line2->SetLineWidth(1);
    line2->SetLineColor(kRed);
    line2->Draw();

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
