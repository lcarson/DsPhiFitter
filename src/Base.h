#ifndef BASE_H 
#define BASE_H 1

// Include files
#include<iostream>
#include<string>
#include<vector>
#include<map>

#include "TString.h"

class Ntuple;

class Base {
public: 
  /// Standard constructor
  Base(); 
  void isTypeUsed(std::string type,bool u){used[type]=u;}
  bool isTypeUsed(std::string type){return used[type];  }
  void defineTypeColors();
  virtual ~Base( ){}

  std::vector<std::string> sources;
  std::vector<std::string> types;
  std::map<std::string,bool> used;
  std::map<std::string,TString> Types;
  std::map<std::string,int> typeColor;
  
	std::string toy;
	std::string s21;
	std::string s21r1;
	std::string Ds2KKPi;
	std::string Ds2PiPiPi;
	std::string Ds2KPiPi;
	std::string signal;
	//std::string buffer;
	//std::string sidebd;
	std::string bckgrd;
	std::string minus;
	std::string plus;
	std::string up;
	std::string dn;
	std::string both;
	std::string sensitive;
	std::string dot;
	std::string slash;
        std::string colon;
	std::string underscore;
        std::string ws;
        std::string coda;
	std::string blind;
	std::string unblind;
	std::string null;

	double pi;
	double twopi;
	double pibytwo;
  
  bool BLIND;
  bool UNBLIND;

	std::vector<std::string> allmodeList;
  std::vector<std::string> allchargeList;
  std::vector<std::string> allmagnetList;
};



#endif // BASE_H
