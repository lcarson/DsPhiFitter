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

void plot_cont(int,double[],double[],double[],std::map<std::string,std::string>);
void plot_1D(int,double[],double[],std::map<std::string,std::string>);
void plot_surf(int,double[],double[],double[],std::map<std::string,std::string>);
void plot(std::string, std::string, std::string, std::string);

/*void plot(int whichone = 0){

  if(whichone == 0) plot_cont();
  if(whichone == 1) plot_surf();
  if(whichone == 2) plot_1D();
  if(whichone == 3) plot_file();
}
*/

void plot_punzi(std::string datatype="MC",std::string method="BDT",std::string mode = "Ds2KKPi",std::string years= ""){
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
  if (mode!= "Ds2KKPi" && mode!= "Ds2KPiPi" && mode!= "Ds2PiPiPi" && mode!= "D2KKPi" && mode!= "D2PiKPi" && mode!= "D2PiPiPi"){
    std::cout << "Mode must be Ds2KKPi, Ds2KPiPi, or Ds2PiPiPi, you typed: " << mode << std::endl;
    return;
  }

  std::string dir = Form("results/manyfits/text/%s/%s/%s/",mode.c_str(),method.c_str(),datatype.c_str());

  std::map<std::string,std::string> config;
  config["datatype"] = datatype;
  config["method"]   = method;
  config["mode"]     = mode;
  config["years"]    = years;
  config["dir"]      = dir;

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
  
  // Get Signal effiencies 
  std::string dir_eff;
  std::string decay    = "Phi2KK";
  std::string decay_D0 = "D02KK";
  dir_eff = Form("/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Efficiency_Study/DV_v41r1/MC_Method/text/%s/",mode.c_str());
  
  if(mode == "D2PiKPi" || mode == "D2KKPi"||mode == "D2PiPiPi" ) {
    dir_eff = Form("/home/hadavizadeh/Bc_Analysis/DataStripping/B_KstD/Efficiency_Study/DV_v41r1/MC_Method/text/%s/",mode.c_str());
    decay = "Kst02KPi";
    decay_D0 = "D02KPi";
  }
  std::vector<std::string> files_eff;
  getdir(dir_eff,files_eff);

  std::cout << "Number of efficiency files: " << files_eff.size() << std::endl; 
  if (files_eff.size()==0) return;

  std::vector<double> sig_efficiency;
  std::vector<double> sig_efficiency_DsD0;
  std::ifstream input;
  for(int i=0; i<cut_ds.size();i++){
    std::cout << "Looking for signal efficiency of cuts: Ds" << cut_ds[i] << " and Phi: " << cut_phi[i] <<std::endl;
    int found = 0;
    if (years=="Run1"){
      double eff_2011_dn;
      double eff_2011_up;
      double eff_2012_dn;
      double eff_2012_up;
      int found_11_dn=0;
      int found_11_up=0;
      int found_12_dn=0;
      int found_12_up=0;
      // ========== Find 2011 MagDown Eff. ==================
      std::string text_filename_11_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2011",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_11_dn){
          found_11_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_11_dn << " \t" << atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2011_dn = atof(params[4].c_str());
        }
      }
      if(found_11_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_11_dn <<std::endl;
      }
      // ========== Find 2011 MagUp Eff. ==================
      std::string text_filename_11_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2011",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_11_up){
          found_11_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_11_up << " \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2011_up = atof(params[4].c_str());
        }
      }
      if(found_11_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_11_up <<std::endl;
      }
      // ========== Find 2012 MagDown Eff. ==================
      std::string text_filename_12_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2012",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_12_dn){
          found_12_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_12_dn << " \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2012_dn = atof(params[4].c_str());
        }
      }
      if(found_12_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_12_dn <<std::endl;
      }
      // ========== Find 2012 MagUp Eff. ==================
      std::string text_filename_12_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2012",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_12_up){
          found_12_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_12_up <<" \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2012_up = atof(params[4].c_str());
        }
      }
      if(found_12_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_12_up <<std::endl;
      }
      double average = (eff_2011_dn +eff_2011_up+eff_2012_dn+eff_2012_up)/4;
      sig_efficiency.push_back(average);
    } else if(years=="Run2"){
      double eff_2015_dn;
      double eff_2015_up;
      double eff_2016_dn;
      double eff_2016_up;
      int found_15_dn;
      int found_15_up;
      int found_16_dn;
      int found_16_up;
      
      // ========== Find 2015 MagDown Eff. ==================
      std::string text_filename_15_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2015",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_15_dn){
          found_15_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_15_dn << " \t" << atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2015_dn = atof(params[4].c_str());
        }
      }
      if(found_15_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_15_dn <<std::endl;
      }
      // ========== Find 2015 MagUp Eff. ==================
      std::string text_filename_15_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2015",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_15_up){
          found_15_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_15_up << " \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2015_up = atof(params[4].c_str());
        }
      }
      if(found_15_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_15_up <<std::endl;
      }
            // ========== Find 2016 MagDown Eff. ==================
      std::string text_filename_16_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2016",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_16_dn){
          found_16_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_16_dn << " \t" << atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2016_dn = atof(params[4].c_str());
        }
      }
      if(found_16_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_16_dn <<std::endl;
      }
      // ========== Find 2016 MagUp Eff. ==================
      std::string text_filename_16_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2016",  decay.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_16_up){
          found_16_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_16_up << " \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2016_up = atof(params[4].c_str());
        }
      }
      if(found_16_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_16_up <<std::endl;
      }

      double average = (eff_2015_dn+eff_2015_up+eff_2016_dn +eff_2016_up)/4;
      sig_efficiency.push_back(average);
    }
  }
  // ===============================// 
  //         Find eff for DsD0
  // ===============================//


  for(int i=0; i<cut_ds.size();i++){
    std::cout << "Looking for signal efficiency of cuts: Ds" << cut_ds[i] << " and Phi: " << cut_phi[i] <<std::endl;
    int found = 0;
    if (years=="Run1"){
      double eff_2011_dn;
      double eff_2011_up;
      double eff_2012_dn;
      double eff_2012_up;
      int found_11_dn=0;
      int found_11_up=0;
      int found_12_dn=0;
      int found_12_up=0;
      // ========== Find 2011 MagDown Eff. ==================
      std::string text_filename_11_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2011",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_11_dn){
          found_11_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_11_dn <<" \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2011_dn = atof(params[4].c_str());
        }
      }
      if(found_11_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_11_dn <<std::endl;
      }
      // ========== Find 2011 MagUp Eff. ==================
      std::string text_filename_11_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2011",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_11_up){
          found_11_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_11_up <<" \t"<<  atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2011_up = atof(params[4].c_str());
        }
      }
      if(found_11_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_11_up <<std::endl;
      }
      // ========== Find 2012 MagDown Eff. ==================
      std::string text_filename_12_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2012",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_12_dn){
          found_12_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_12_dn  <<" \t"<<atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2012_dn = atof(params[4].c_str());
        }
      }
      if(found_12_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_12_dn <<std::endl;
      }
      // ========== Find 2012 MagUp Eff. ==================
      std::string text_filename_12_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2012",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_12_up){
          found_12_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_12_up <<" \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2012_up = atof(params[4].c_str());
        }
      }
      if(found_12_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_12_up <<std::endl;
      }
      double average = (eff_2011_dn +eff_2011_up+eff_2012_dn+eff_2012_up)/4;
      sig_efficiency_DsD0.push_back(average);
    
    } else if(years=="Run2"){
      double eff_2015_dn;
      double eff_2015_up;
      double eff_2016_dn;
      double eff_2016_up;
      int found_15_dn;
      int found_15_up;
      int found_16_dn;
      int found_16_up;
      // ========== Find 2015 MagDown Eff. ==================
      std::string text_filename_15_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2015",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_15_dn){
          found_15_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_15_dn <<" \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2015_dn = atof(params[4].c_str());
        }
      }
      if(found_15_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_15_dn <<std::endl;
      }
      // ========== Find 2015 MagUp Eff. ==================
      std::string text_filename_15_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2015",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_15_up){
          found_15_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_15_up <<" \t"<<  atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2015_up = atof(params[4].c_str());
        }
      }
      if(found_15_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_15_up <<std::endl;
      }
      // ========== Find 2016 MagDown Eff. ==================
      std::string text_filename_16_dn  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Down", "2016",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_16_dn){
          found_16_dn++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file: " << text_filename_16_dn <<" \t"<< atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2016_dn = atof(params[4].c_str());
        }
      }
      if(found_16_dn==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_16_dn <<std::endl;
      }
      // ========== Find 2016 MagUp Eff. ==================
      std::string text_filename_16_up  = Form("Eff_%s_%s_%s_%s_%f_%f.txt", "Up", "2016",  decay_D0.c_str() ,mode.c_str(),cut_phi[i],cut_ds[i]);
      for(std::vector<std::string>::iterator it=files_eff.begin();it!=files_eff.end();it++){ 
        if((*it)==text_filename_16_up){
          found_16_up++;
          std::string filename = (*it);
          std::vector<std::string> params;
          std::ifstream input;
          input.open((dir_eff+filename).c_str(),std::ifstream::in);
          std::string line = "";
          getline(input,line);
          split(line,params,":");
          std::cout<< "=== Found file:   " << text_filename_16_up <<" \t"<<  atof(params[4].c_str()) * 100 << "%" <<std::endl;
          eff_2016_up = atof(params[4].c_str());
        }
      }
      if(found_16_up==0) {
        std::cout<< "WARNING =====> Missing file:" << text_filename_16_up <<std::endl;
      }


      double average = (eff_2015_dn +eff_2015_up+eff_2016_dn +eff_2016_up)/4;
      sig_efficiency_DsD0.push_back(average);

    }
  }


  std::cout<<std::endl;
  std::cout << "Ds cuts:" << std::endl;
  for(int i=0; i<cut_ds.size();i++){
    std::cout << cut_ds[i] << std::endl;
  }
  std::cout<<std::endl;

  std::cout << "D0 cuts:" << std::endl;
  for(int i=0; i<cut_phi.size();i++){
    std::cout << cut_phi[i] << std::endl;
  }
  std::cout<<std::endl;

  double* n_sig_eff_ar   = new double[count];
  double* n_DsD0_eff_ar   = new double[count];

  double* n_sig_ar   = new double[count];
  double* n_bkg_ar   = new double[count];
  double* n_sig_corrected_ar   = new double[count];
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
    n_sig_eff_ar[i]  = sig_efficiency[i];
    n_DsD0_eff_ar[i] = sig_efficiency_DsD0[i];
    n_sig_ar[i]   = n_sig[i];

    n_sig_corrected_ar[i] = (sig_efficiency[i]/sig_efficiency_DsD0[i]) * 0.002533749 * n_sig[i];

    n_bkg_ar[i]   = n_bkg[i];
    if (datatype=="DATA"){
      cut_ds_ar[i]  = cut_ds[i];
      cut_phi_ar[i] = cut_phi[i];
    } else {
      cut_MC_ar[i]  = cut_MC[i];
    }

    // -------- Definitions of FOMs ------------------
    //n_punzi_ar[i] = n_sig[i] / (2.5 + sqrt(n_bkg[i]));
    n_punzi_ar[i] = sig_efficiency[i] / (2.5 + sqrt(n_bkg[i]));
    if (n_bkg[i]+n_sig[i]>0){ 
      n_signi_ar[i] = n_sig[i] / sqrt(n_bkg[i]+n_sig[i]);
      n_sigaj_ar[i] = (n_sig_corrected_ar[i]) / sqrt(n_bkg[i]+n_sig_corrected_ar[i]);
    }else {
      n_signi_ar[i] = 0;
      n_sigaj_ar[i] = 0;
    }
    n_diff_ar[i] = n_punzi_ar[i]/n_signi_ar[i];
    n_purity_ar[i] = n_sig[i]/(n_sig[i]+n_bkg[i]);
    if((n_sig[i]+n_bkg[i])<0.01) n_purity_ar[i] = 0; 
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
    plot_1D(count,cut_MC_ar,n_sigaj_ar,config);
    //plot_1D(count,cut_MC_ar,n_signi_ar,dir);
    //plot_1D(count,cut_MC_ar,n_pursig_ar,dir);
  } else {
    plot_cont(count,cut_ds_ar,cut_phi_ar,n_punzi_ar,config);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_sigaj_ar,config);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_signi_ar,dir);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_sigaj_ar,dir);
    //plot_cont(count,cut_ds_ar,cut_phi_ar,n_pursig_ar,dir);
  }

std::cout<<"Done" <<std::endl;
}

void plot_1D(int n_p,double cut_value_MC [],double Punzi_1D [], std::map<std::string,std::string> config){
  SetLHCbStyle(oneD);
  TCanvas *cav_punzi = new TCanvas(Form("%s_%s_punzi_%s",config["mode"].c_str(),config["method"].c_str(),config["years"].c_str()),Form("%s %s Punzi as function of cut value",config["mode"].c_str(),config["method"].c_str()),200,10,700,500);
  TGraph  *gr_punzi  = new TGraph(n_p, cut_value_MC, Punzi_1D);
  gr_punzi->SetTitle(Form("Punzi    S/(5/2 + sqrt(B));MC %s Cut Value;S/sqrt(S+B)",config["method"].c_str()));  
  gr_punzi->Draw("AP");

  cav_punzi->Print(Form("%s%s_MC.eps",config["dir"].c_str(),cav_punzi->GetName()));
  cav_punzi->Print(Form("%s%s_MC.pdf",config["dir"].c_str(),cav_punzi->GetName()));
}

void plot_surf(int n_p,double cut_value_1D_Ds[],double cut_value_1D_Phi[],double Punzi_2D [], std::map<std::string,std::string> config){
  SetLHCbStyle(surf);
  std::cout << "Set Style" << std::endl;
  TCanvas *cav_punzi = new TCanvas(Form("%s_%s_punzi_%s",config["mode"].c_str(),config["method"].c_str(),config["years"].c_str()),Form("%s %s Punzi as function of cut value",config["mode"].c_str(),config["method"].c_str()),200,10,700,500);
  std::cout << "Made Canvas" << std::endl;
  TGraph2D *gr_punzi = new TGraph2D(n_p, cut_value_1D_Ds, cut_value_1D_Phi, Punzi_2D); 
  std::cout << "Made TGraph2D" << std::endl;  
  gr_punzi->SetTitle(Form("Punzi    S/(5/2 + sqrt(B));D %s Cut Value;D0 %s Cut Value;S/sqrt(S+B)",config["method"].c_str(),config["method"].c_str()));  
  std::cout << "Set Title" << std::endl;
  gr_punzi->Draw("surf1z");

  cav_punzi->Print(Form("%s%s_surf.eps",config["dir"].c_str(),cav_punzi->GetName()));
  cav_punzi->Print(Form("%s%s_surf.pdf",config["dir"].c_str(),cav_punzi->GetName()));
    
  //TFile tfile(Form("%sPlot.root",dir.c_str()),"RECREATE");
  //cav_punzi->Write();
  //gr_punzi->Write(); 
  //tfile.Close();  
}

void plot_cont(int n_p,double cut_value_1D_Ds[],double cut_value_1D_Phi[],double Punzi_2D [], std::map<std::string,std::string> config){
  SetLHCbStyle(cont);
  TCanvas *cav_punzi3 = new TCanvas(Form("%s_%s_punzi_%s",config["mode"].c_str(),config["method"].c_str(),config["years"].c_str()),Form("%s %s Punzi as function of cut value",config["mode"].c_str(),config["method"].c_str()),200,10,700,500);
  TGraph2D *gr_punzi2 = new TGraph2D(n_p, cut_value_1D_Ds, cut_value_1D_Phi, Punzi_2D);   
  gr_punzi2->SetTitle(Form("Punzi    S/(5/2 + sqrt(B));D %s Cut Value;D0 %s Cut Value;S/sqrt(S+B)",config["method"].c_str(),config["method"].c_str())); 
  gr_punzi2->Draw("colz");

  cav_punzi3->Print(Form("%s%s_cont.eps",config["dir"].c_str(),cav_punzi3->GetName()));
  cav_punzi3->Print(Form("%s%s_cont.pdf",config["dir"].c_str(),cav_punzi3->GetName()));
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
