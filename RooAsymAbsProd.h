/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOASYMABSPROD
#define ROOASYMABSPROD

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooAsymAbsProd : public RooAbsReal {
public:
  RooAsymAbsProd() {} ; 
  RooAsymAbsProd(const char *name, const char *title,
	      RooAbsReal& _neg_factor,
	      RooAbsReal& _pos_factor,
	      RooAbsReal& _base_par,
              bool _same_initial_and_final=false );
  RooAsymAbsProd(const RooAsymAbsProd& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooAsymAbsProd(*this,newname); }
  inline virtual ~RooAsymAbsProd() { }

protected:

  RooRealProxy neg_factor ;
  RooRealProxy pos_factor ;
  RooRealProxy base_par ;
  bool same_initial_and_final ;

  Double_t evaluate() const ;

private:

  ClassDef(RooAsymAbsProd,1) // Your description goes here...
};
 
#endif