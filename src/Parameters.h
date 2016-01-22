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
  //std::string variationoption;
  std::vector<std::string> locations;
  std::string moption;
  std::string doption;
  std::string polarity;
  //std::string toVary;
  //std::string toFix;
  bool readToys;
  bool genToys;
  int  nToys;
  bool doDraw;
  bool doFit;
  bool sumOverCharges;
  bool batch;
  //bool vary;
  bool binned;
  std::map<std::string,bool> modes;
  std::map<std::string,bool> dsetsReq;
  std::map<std::string,std::map<std::string,bool> > dsetsFound;
  bool quickVersion;
  //std::map<std::string,bool> variation;
  bool minos;
};
#endif // PARAMETERS_H
