// Include files 

#include "Base.h"
#include "CommonTools.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Base
//
// 2009-06-19 : Malcolm John
//-----------------------------------------------------------------------------

Base::Base()
	: toy("toy")
  , s21("s21")
	, s21r1("s21r1")
  , s24("s24")
  , s26("s26")
  , Ds2PhiPi("Ds2PhiPi")
  , Ds2KKPi("Ds2KKPi")
	, Ds2PiPiPi("Ds2PiPiPi")
  , Ds2KPiPi("Ds2KPiPi")
  , DsD0("DsD0")
  , DsPhi("DsPhi")
  , DsPhiSide("DsPhiSide")
  , Helbin1("Helbin1")
  , Helbin2("Helbin2")
  , DsBDTbin1("DsBDTbin1")
  , DsBDTbin2("DsBDTbin2")
  , PhiBDTbin1("PhiBDTbin1")
  , PhiBDTbin2("PhiBDTbin2")
	, signal("signal")
  , Ds("Ds")
  , D0("D0")
  , Phi("Phi")
  , BDT("BDT")
  , BDTG("BDTG")
  , BDTB("BDTB")
  , DATA("DATA")
  , MC("MC")
  , cont("cont")
  , surf("surf")
  , norm("norm")
	  //, buffer("buffer")
	  //, sidebd("sidebd")
	, bckgrd("bckgrd")
	, minus("minus")
	, plus("plus")
	, up("up")
	, dn("dn")
	, both("both")
	, sensitive("sensitive")
	, dot(".")
	, slash("/")
  , colon(":")
	, underscore("_")
	, ws("WS")
	, coda(".root")
	, blind("blind")
	, unblind("unblind")
	, null("NULL")
	, pi(3.141592653589793)
	, twopi(6.283185307179586)
	, pibytwo(1.570796326794897)
  , BLIND(true)
	, UNBLIND(false)
  , allmodeList()
  , allBmodeList()
  , allchargeList()
  , allmagnetList()
  , allHelbinList()
  , allDsBDTbinList()
  , allPhiBDTbinList()
{
  types.push_back(signal);
  defineTypeColors();

  std::vector<std::string>::iterator iterType=types.begin();
  for(;iterType!=types.end();iterType++){
    std::string type = (*iterType);
    used.insert(std::make_pair(type,false));
  }

  allmodeList.push_back(Ds2PhiPi);
  allmodeList.push_back(Ds2KKPi);
  allmodeList.push_back(Ds2PiPiPi);
  allmodeList.push_back(Ds2KPiPi);
  allBmodeList.push_back(DsD0);
  allBmodeList.push_back(DsPhi);
  allBmodeList.push_back(DsPhiSide);
  allchargeList.push_back(both);
  allchargeList.push_back(minus);
  allchargeList.push_back(plus);
  allmagnetList.push_back(up);
  allmagnetList.push_back(dn);
  allmagnetList.push_back(both);
  allHelbinList.push_back(both);
  allHelbinList.push_back(Helbin1);
  allHelbinList.push_back(Helbin2);
  allDsBDTbinList.push_back(DsBDTbin1);
  allDsBDTbinList.push_back(DsBDTbin2);
  allPhiBDTbinList.push_back(PhiBDTbin1);
  allPhiBDTbinList.push_back(PhiBDTbin2);
}

void Base::defineTypeColors()
{
  Base::typeColor.insert(std::make_pair(signal,3));
}
