#ifndef PARAMETERS_H 
#define PARAMETERS_H 1

// Include files
#include<iostream>
#include<string>
#include<vector>
#include<map>

#include "Base.h"

class Parameters : public Base {
 public: 
  /// Standard constructor
  Parameters( ); 
  int readCommandLine(unsigned int,char**);
  bool exists(std::string);
  void defineTypeColors();
  void help();

  virtual ~Parameters(){}

  std::string toylocation;
  //std::string autolocation;
  std::string locationoption;
  std::string variationoption;
  std::vector<std::string> locations;
  std::string moption;
  std::string doption;
  std::string coption;
  std::string soption;
  std::string BDToption;
  std::string polarity;
  std::string normcat;
  std::string toVary;
  //std::string toFix;
  bool readToys;
  bool genToys;
  int  nToys;
  int  nDsBDTBins;
  int  nPhiBDTBins;
  bool doDraw;
  bool doFit;
  bool sumOverCharges;
  bool batch;
  bool debug;
  bool splitHel;
  bool splitYears;
  bool manyFits;
  bool runEff;
  int  seed;
  bool useSeed;
  bool vary;
  bool readSys;
  bool binned;
  bool fitBr;
  bool fitFourBr;
  bool doMerge;
  bool doSensitivity;
  double sensitivityBR;
  bool doLikelihood;
  double likelihoodBR;
  int sensitivityN;
  std::string MVAMethod; 
  std::string MVAType;
  int nBDTPoints;
  std::map<std::string,bool> modes;
  std::map<std::string,bool> Bmodes;
  std::map<std::string,bool> dsetsReq;
  std::map<std::string,std::map<std::string,bool> > dsetsFound;
  bool quickVersion;
  std::map<std::string,bool> variation;
  bool minos;
};
#endif // PARAMETERS_H
