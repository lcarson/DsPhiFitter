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
    std::string s24;
    std::string s26;
    std::string Ds2PhiPi;
    std::string Ds2KKPi;
    std::string Ds2PiPiPi;
    std::string Ds2KPiPi;
    std::string D2KKPi;
    std::string D2PiPiPi;
    std::string D2PiKPi;
    std::string DsD0;
    std::string DsPhi;
    std::string DD0;
    std::string DKst0;
    std::string DKst0Side;
    std::string DsPhiSide;
    std::string DsPhiSideWide;
    std::string Helbin1;
    std::string Helbin2;
    std::string DsBDTbin1;
    std::string DsBDTbin2;
    std::string PhiBDTbin1;
    std::string PhiBDTbin2;
    std::string signal;
    std::string Ds;
    std::string D;
    std::string D0;
    std::string Phi;
    std::string Kst0;
    std::string BDT;
    std::string BDTG;
    std::string BDTB;
    std::string DATA;
    std::string MC;
    std::string cont;
    std::string surf;
    std::string norm;
    std::string merged;
    std::string mergedDModes;
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
    std::vector<std::string> allBmodeList;
    std::vector<std::string> allchargeList;
    std::vector<std::string> allmagnetList;    
    std::vector<std::string> allHelbinList;
    std::vector<std::string> allDsBDTbinList;
    std::vector<std::string> allPhiBDTbinList;

    std::string fixedSig_n;
    std::string fixedSig_alpha;
    std::string fixedSig_Sigmaratio;
    std::string fixedSig_Sigmafrac;
    std::string fixedSig_NormSigma;
    std::string fixedSig_BinRatios;
    std::string fixedBG_DsD0;
    std::string fixedBG_DsPhi;
    std::string fixedBG_Dsa1;
    std::string fixedBG_hel;
    std::string fixedBG_noDsstPhi;
    std::string fixedBG_noDsa1;
    std::string fixedBG_slope;
    std::string fixedBG_Dsa1_smear;
    std::string fixedBG_DsstKKst_smear;
    std::string fixedBG_DD_Ratios;
    std::string fixedBG_DD_Fractions;
    std::string fixedBG_DsstDst0_endpoints;
    std::string fixedBG_DsstDst0_shape;
    std::string fixedBG_Comb_shape;
    std::string fixedBG_Comb_shape_flat;
    std::string fixedBG_Comb_shape_line;
    std::string fixedBG_DsKK_fractions;
    std::string draw;
    std::string fixedSig_double;
    std::string fixedBG_DD0;
    std::string fixedBG_DKst0;
    std::string doNothing;
};



#endif // BASE_H
