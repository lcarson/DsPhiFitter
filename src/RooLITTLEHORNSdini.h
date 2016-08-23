
/*-----------------------------------------------------------------------------------
Author: Faye Cheung
Email:  faye.cheung@cern.ch

Double horn distribution + constant convoluted with a Gaussian resolution function
Based on class RooHORNSdini

-- v3 --
SHIFT included
EFFICIENCY CORRECTION included

-- v8 --
EFFICIENCY is a double GAUSSIAN: two extra parameters: ratio_sigma, fraction_sigma
-----------------------------------------------------------------------------------*/

#ifndef ROO_ROOLITTLEHORNSDINI
#define ROO_ROOLITTLEHORNSDINI

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "Rtypes.h"


class RooRealVar;

class  RooLITTLEHORNSdini : public RooAbsPdf {


public:
  RooLITTLEHORNSdini(const char *name, const char *title,
		  RooAbsReal& _m, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _csi, RooAbsReal& _shift, RooAbsReal& _sigma, RooAbsReal& _ratio_sigma, RooAbsReal& _fraction_sigma, RooAbsReal& _shiftg);
  RooLITTLEHORNSdini(const  RooLITTLEHORNSdini& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new  RooLITTLEHORNSdini(*this,newname); }
  inline virtual ~ RooLITTLEHORNSdini() {}


protected:
  RooRealProxy m;
  RooRealProxy a;
  RooRealProxy b;
  RooRealProxy csi;
  RooRealProxy shift;
  RooRealProxy sigma;
  RooRealProxy ratio_sigma;
  RooRealProxy fraction_sigma;
  RooRealProxy shiftg;
  Double_t evaluate() const;


private:
  ClassDef( RooLITTLEHORNSdini,2) //  RooLITTLEHORNSdini function PDF
};

#endif
