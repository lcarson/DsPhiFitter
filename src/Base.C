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
	, Ds2KKPi("Ds2KKPi")
	, Ds2PiPiPi("Ds2PiPiPi")
	, Ds2KPiPi("Ds2KPiPi")
	, signal("signal")
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
  , allchargeList()
  , allmagnetList()
{
  types.push_back(signal);
  defineTypeColors();

  std::vector<std::string>::iterator iterType=types.begin();
  for(;iterType!=types.end();iterType++){
    std::string type = (*iterType);
    used.insert(std::make_pair(type,false));
  }

  allmodeList.push_back(Ds2KKPi);
  allmodeList.push_back(Ds2PiPiPi);
  allmodeList.push_back(Ds2KPiPi);
  allchargeList.push_back(both);
  allchargeList.push_back(minus);
  allchargeList.push_back(plus);
  allmagnetList.push_back(up);
  allmagnetList.push_back(dn);
  allmagnetList.push_back(both);
}

void Base::defineTypeColors()
{
  Base::typeColor.insert(std::make_pair(signal,3));
}
