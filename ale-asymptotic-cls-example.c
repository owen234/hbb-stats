#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "RooFitResult.h"

#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace RooStats;


void RunAsyCLsSingle(TString fileName, double testVal, TString outFile) {

  gROOT->LoadMacro("RooBetaPdf.cxx+") ;
  gROOT->LoadMacro("RooRatio.cxx+") ;
  gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;

  // get relevant objects out of the "ws" file

  TFile *file = TFile::Open(fileName);
  if(!file){
    cout <<"file not found" << endl;
    return;
  } 

  RooWorkspace* w = (RooWorkspace*) file->Get("ws");
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }

  ModelConfig* mc = (ModelConfig*) w->obj("SbModel");
  RooAbsData* data = w->data("ra2b_observed_rds");

  if( !data || !mc ){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }

  RooRealVar* myPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  myPOI->setRange(0, 1000.);

  ModelConfig* bModel = (ModelConfig*) w->obj("BModel");
  ModelConfig* sbModel = (ModelConfig*) w->obj("SbModel");

  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  profll.SetOneSided(1);

  AsymptoticCalculator asymptotic_calculator(*data, *bModel, *sbModel);
  HypoTestInverter calc_asymptotic_calculator(asymptotic_calculator);
  
  calc_asymptotic_calculator.SetConfidenceLevel(0.95);
  calc_asymptotic_calculator.SetVerbose(false);
  calc_asymptotic_calculator.UseCLs(true);

  HypoTestInverterResult* res_asymptotic_calculator ;
  calc_asymptotic_calculator.SetFixedScan( 1 , testVal , testVal) ;
  res_asymptotic_calculator = calc_asymptotic_calculator.GetInterval();


  // get expected limits:

  SamplingDistribution * s = res_asymptotic_calculator->GetExpectedPValueDist(0);
  if ( !s)  return ; 
  const std::vector<double> & values = s->GetSamplingDistribution();

  double maxSigma = res_asymptotic_calculator->fgAsymptoticMaxSigma;
  double dsig = 2* maxSigma/ (values.size() -1) ;         

  int  im = (int) TMath::Floor ( ( -1 + maxSigma )/dsig + 0.5);
  int  i0 = (int) TMath::Floor ( ( maxSigma )/dsig + 0.5);
  int  ip = (int) TMath::Floor ( ( +1 + maxSigma )/dsig + 0.5);

  double CLs_expM = values[im];
  double CLs_exp0 = values[i0];
  double CLs_expP = values[ip];

  // dump results string to output file
  ofstream outStream ;
  outStream.open(outFile,ios::app) ;
  
  outStream << "testVal = " << testVal 
	    << "   CLs = " << res_asymptotic_calculator->CLs(0) 
	    << "   CLb = " << res_asymptotic_calculator->CLb(0) 
	    << "   CLsplusb = " << res_asymptotic_calculator->CLsplusb(0) 
	    << "   CLs_exp = " << CLs_exp0 
	    << "   CLs_exp(-1s) = " << CLs_expM 
	    << "   CLs_exp(+1s) = " << CLs_expP << endl ; 
  
  outStream.close() ;

  return ;

}
