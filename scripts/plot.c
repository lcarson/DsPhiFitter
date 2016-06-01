#include <fstream>
#include <dirent.h>
#include <errno.h>
#include <iostream>
void SetLHCbStyle(std::string);
int getdir (std::string, std::vector<std::string> &);
void split(const std::string& ,std::vector<std::string>& ,const std::string&);

std::string oneD;
std::string cont;
std::string surf;

void plot_cont(int,double[],double[],double[],std::string);
void plot_1D(int,double[],double[],std::string);
void plot_surf(int,double[],double[],double[],std::string);
void plot(std::string, std::string, std::string, std::string);

/*void plot(int whichone = 0){

  if(whichone == 0) plot_cont();
  if(whichone == 1) plot_surf();
  if(whichone == 2) plot_1D();
  if(whichone == 3) plot_file();
}
*/

void plot(std::string datatype="MC",std::string method="BDT",std::string mode = "Ds2KKPi",std::string years= ""){
  oneD = "oneD";
  cont = "cont";
  surf = "surf";
  if (datatype!= "MC" && datatype!= "DATA"){
    std::cout << "Type must be MC or DATA, you typed: " << datatype << std::endl;
    return;
  }
 if (method!= "BDT" && method!= "BDTG" && method!= "BDTB"){
    std::cout << "Method must be BDT, BDTG, or BDTB, you typed: " << method << std::endl;
    return;
  }
  if (mode!= "Ds2KKPi" && mode!= "Ds2KPiPi" && mode!= "Ds2PiPiPi"){
    std::cout << "Mode must be Ds2KKPi, Ds2KPiPi, or Ds2PiPiPi, you typed: " << mode << std::endl;
    return;
  }

  std::string dir = Form("results/manyfits/text/%s/%s/%s/",mode.c_str(),method.c_str(),datatype.c_str());
  std::vector<std::string> files;
  getdir(dir,files);
  std::cout << "Number of files: " << files.size() << std::endl; 
  if (files.size()==0) return;

  std::vector<double> n_sig;
  std::vector<double> n_bkg;
  std::vector<double> cut_ds;
  std::vector<double> cut_phi;
  std::vector<double> cut_MC;
  double ratio = 0.023;
  int count =0;
  for(std::vector<std::string>::iterator it=files.begin();it!=files.end();it++){ 
    std::cout << "Processing file : " << (*it) << std::endl;
    std::string filename = (*it);
    std::vector<std::string> params;
    std::ifstream input;
    input.open((dir+filename).c_str(),std::ifstream::in);
    std::string line = "";
    getline(input,line);
    std::cout<< "Line read from file: "<< line << std::endl;
    split(line,params,":");
    if((*it).find("txt") !=std::string::npos && ((*it).find(years.c_str())!= std::string::npos || years=="") ){
      //cut_ds[count]  = atof(params[0].c_str());
      //cut_phi[count] = atof(params[1].c_str());
      //n_sig[count]   = atof(params[2].c_str());
      //n_bkg[count]   = atof(params[3].c_str());
      if (datatype=="DATA"){
        cut_ds.push_back( atof(params[0].c_str()));
        cut_phi.push_back(atof(params[1].c_str()));
        n_sig.push_back(  atof(params[2].c_str()));
        n_bkg.push_back(  atof(params[3].c_str()));
      } else {
        cut_MC.push_back(  atof(params[0].c_str()));  
        n_sig.push_back(   atof(params[1].c_str()));
        n_bkg.push_back(   atof(params[2].c_str()));
      }
      
      count++;
    }
  }
  double* n_sig_ar   = new double[count];
  double* n_bkg_ar   = new double[count];
  double* cut_ds_ar  = new double[count];
  double* cut_phi_ar = new double[count];
  double* cut_MC_ar = new double[count];

  double* n_punzi_ar   = new double[count];
  double* n_signi_ar   = new double[count];
  double* n_sigaj_ar   = new double[count];
  double* n_diff_ar    = new double[count];
  double* n_purity_ar  = new double[count];
  double* n_pursig_ar  = new double[count];
 
  for(int i=0; i<count; i++){
    n_sig_ar[i]   = n_sig[i];
    n_bkg_ar[i]   = n_bkg[i];
    if (datatype=="DATA"){
      cut_ds_ar[i]  = cut_ds[i];
      cut_phi_ar[i] = cut_phi[i];
    } else {
      cut_MC_ar[i]  = cut_MC[i];
    }

    n_punzi_ar[i] = n_sig[i] / (2.5 + sqrt(n_bkg[i]));
    if (n_bkg[i]+n_sig[i]>0){ 
      n_signi_ar[i] = n_sig[i] / sqrt(n_bkg[i]+n_sig[i]);
      n_sigaj_ar[i] = (ratio*n_sig[i]) / sqrt(n_bkg[i]+ratio*n_sig[i]);
    }else {
      n_signi_ar[i] = 0;
      n_sigaj_ar[i] = 0;
    }
    n_diff_ar[i] = n_punzi_ar[i]/n_signi_ar[i];
    n_purity_ar[i] = n_sig[i]/(n_sig[i]+n_bkg[i]);
    n_pursig_ar[i] = n_purity_ar[i]*n_punzi_ar[i];

    if (datatype=="DATA") std::cout << "Cut Ds:"<< cut_ds_ar[i] << ", \t Cut Phi: " << cut_phi_ar[i]<<", \t NSig: "<<n_sig_ar[i]<<", \t NBKG: " << n_bkg_ar[i] << ", \t Punzi: "<<n_punzi_ar[i]<<", \t Significance: " << n_signi_ar[i] <<", \t Significance adjusted: " << n_sigaj_ar[i] << std::endl;
    if (datatype=="MC")   std::cout << "Cut MC:"<< cut_MC_ar[i] << ", \t NSig: "<<n_sig_ar[i]<<", \t NBKG: " << n_bkg_ar[i] << ", \t Punzi: "<<n_punzi_ar[i]<<", \t Significance: " << n_signi_ar[i] <<", \t Significance adjusted: " << n_sigaj_ar[i] << std::endl;
  }
  if (count==0){
    std::cout<< "Found no files :(" <<std::endl;
    return;
  } else {
    std::cout << "Found " << count << " files..." <<std::endl;
  }
  if(datatype=="MC"){
    //plot_1D(count,cut_MC_ar,n_punzi_ar,dir);
    //plot_1D(count,cut_MC_ar,n_signi_ar,dir);
    plot_1D(count,cut_MC_ar,n_pursig_ar,dir);
  } else {
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_punzi_ar,dir);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_signi_ar,dir);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_sigaj_ar,dir);
    plot_cont(count,cut_ds_ar,cut_phi_ar,n_pursig_ar,dir);
  }

std::cout<<"Done" <<std::endl;
}

void plot_1D(int n_p,double cut_value_MC [],double Punzi_1D [], std::string dir){
  SetLHCbStyle(oneD);
  TCanvas *cav_punzi = new TCanvas("name","Title",200,10,700,500);
  TGraph  *gr_punzi  = new TGraph(n_p, cut_value_MC, Punzi_1D);
  gr_punzi->SetTitle("Punzi    S/(5/2 + sqrt(B));MC BDT Cut Value;Punzi");  
  gr_punzi->Draw("AP");
}

void plot_surf(int n_p,double cut_value_1D_Ds[],double cut_value_1D_Phi[],double Punzi_2D [], std::string dir){
  SetLHCbStyle(surf);
  std::cout << "Set Style" << std::endl;
  TCanvas *cav_punzi = new TCanvas("Ds2KKPi_BDT_punzi","Ds2KKPi BDT Punzi as function of cut value",200,10,700,500);
  std::cout << "Made Canvas" << std::endl;
  TGraph2D *gr_punzi = new TGraph2D(n_p, cut_value_1D_Ds, cut_value_1D_Phi, Punzi_2D); 
  std::cout << "Made TGraph2D" << std::endl;  
  gr_punzi->SetTitle("Punzi    S/(5/2 + sqrt(B));D BDT Cut Value;D0 BDT Cut Value;Punzi");  
  std::cout << "Set Title" << std::endl;
  gr_punzi->Draw("surf1z");
  //cav_punzi->Draw();

  TFile tfile(Form("%sPlot.root",dir.c_str()),"RECREATE");
  cav_punzi->Write();
  gr_punzi->Write(); 
  tfile.Close();

}

void plot_cont(int n_p,double cut_value_1D_Ds[],double cut_value_1D_Phi[],double Punzi_2D [], std::string dir){
  SetLHCbStyle(cont);
  TCanvas *cav_punzi3 = new TCanvas("Ds2KKPi_BDT_punzi_cont","Ds2KKPi BDT Punzi as function of cut value",200,10,700,500);
  TGraph2D *gr_punzi2 = new TGraph2D(n_p, cut_value_1D_Ds, cut_value_1D_Phi, Punzi_2D);   
  gr_punzi2->SetTitle("Punzi    S/(5/2 + sqrt(B));D BDTG Cut Value;D0 BDTG Cut Value;Punzi"); 
  gr_punzi2->Draw("colz");
}

void SetLHCbStyle(std::string type = ""){
    oneD = "oneD";
    cont = "cont";
    surf = "surf";
  if(type!=oneD && type!= cont && type!= surf)
  {
    std::cout <<  "Type must be one of 1D, cont or surf. Defaulting to 1D" <<std::endl;
    type = oneD;
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
  if ( type == cont){
    lhcbStyle->SetPalette(112);
    gROOT->SetStyle("lhcbStyle");
    gROOT->ForceStyle();

  } else if (type == surf){  
    lhcbStyle->SetTitleOffset(1.5,"Y");
    lhcbStyle->SetTitleOffset(1.5,"X");
    lhcbStyle->SetTitleOffset(0.83,"Z");
    lhcbStyle->SetPalette(112);
    gROOT->SetStyle("lhcbStyle");
    gROOT->ForceStyle();
  } else {
    lhcbStyle->SetPadRightMargin(0.05);
    gROOT->SetStyle("lhcbStyle");
    gROOT->ForceStyle();
  }


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
void split(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = " ")
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

int getdir (std::string dir, std::vector<std::string> &files)
{
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL) {
    std::cout << "Error(" << errno << ") opening " << dir << std::endl;
    return errno;
    }
  while ((dirp = readdir(dp)) != NULL) {
    std::string file(dirp->d_name);
    files.push_back(file);
  }
  closedir(dp);
  std::cout <<"\n" << std::endl;
  return 0;
}
