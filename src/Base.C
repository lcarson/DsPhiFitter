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
  , D2KKPi("D2KKPi")
  , D2PiPiPi("D2PiPiPi")
  , D2PiKPi("D2PiKPi")
  , DsD0("DsD0")
  , DsPhi("DsPhi")
  , DD0("DD0")
  , DKst0("DKst0")
  , DKst0Side("DKst0Side")
  , DsPhiSide("DsPhiSide")
  , DsPhiSideWide("DsPhiSideWide")
  , Helbin1("Helbin1")
  , Helbin2("Helbin2")
  , DsBDTbin1("DsBDTbin1")
  , DsBDTbin2("DsBDTbin2")
  , PhiBDTbin1("PhiBDTbin1")
  , PhiBDTbin2("PhiBDTbin2")
	, signal("signal")
  , Ds("Ds")
  , D("D")
  , D0("D0")
  , Phi("Phi")
  , Kst0("Kst0")
  , BDT("BDT")
  , BDTG("BDTG")
  , BDTB("BDTB")
  , DATA("DATA")
  , MC("MC")
  , cont("cont")
  , surf("surf")
  , norm("norm")
  , merged("merged")
  , mergedDModes("mergedDModes")
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
  //systematic variations
  ,fixedSig_n("fixedSig_n")
  ,fixedSig_alpha("fixedSig_alpha")
  ,fixedSig_Sigmaratio("fixedSig_Sigmaratio")
  ,fixedSig_Sigmafrac("fixedSig_Sigmafrac")
  ,fixedSig_NormSigma("fixedSig_NormSigma")
  ,fixedSig_BinRatios("fixedSig_BinRatios")
  ,fixedBG_DsD0("fixedBG_DsD0")
  ,fixedBG_DsPhi("fixedBG_DsPhi")
  ,fixedBG_Dsa1("fixedBG_Dsa1")
  ,fixedBG_hel("fixedBG_hel")
  ,fixedBG_noDsstPhi("fixedBG_noDsstPhi")
  ,fixedBG_noDsa1("fixedBG_noDsa1")
  ,fixedBG_slope("fixedBG_slope")
  ,fixedBG_Dsa1_smear("fixedBG_Dsa1_smear")
  ,fixedBG_DsstKKst_smear("fixedBG_DsstKKst_smear")
  ,fixedBG_DD_Ratios("fixedBG_DD_Ratios")
  ,fixedBG_DD_Fractions("fixedBG_DD_Fractions")
  ,fixedBG_DsstDst0_endpoints("fixedBG_DsstDst0_endpoints")
  ,fixedBG_DsstDst0_shape("fixedBG_DsstDst0_shape")
  ,fixedBG_Comb_shape("fixedBG_Comb_shape")
  ,fixedBG_Comb_shape_flat("fixedBG_Comb_shape_flat")
  ,fixedBG_Comb_shape_line("fixedBG_Comb_shape_line")
  ,fixedBG_DsKK_fractions("fixedBG_DsKK_fractions")
  ,draw("draw")
  ,fixedSig_double("fixedSig_double")   
  ,fixedBG_DD0("fixedBG_DD0")
  ,fixedBG_DKst0("fixedBG_DKst0")
  ,doNothing("doNothing")
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
  allBmodeList.push_back(DsPhiSideWide);
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
