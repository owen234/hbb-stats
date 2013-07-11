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
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace RooStats;


void RunToyScan5(TString fileName, double startVal, double stopVal, TString outFile) {

  int nToys = 1000 ;

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
  TestStatistic * testStat = &profll;

  HypoTestCalculatorGeneric *  hc = 0;
  hc = new FrequentistCalculator(*data, *bModel, *sbModel);
  
  ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
  toymcs->SetNEventsPerToy(1);
  toymcs->SetTestStatistic(testStat);
  
  ((FrequentistCalculator *)hc)->SetToys(nToys,nToys);
  
  HypoTestInverter calc(*hc);
  calc.SetConfidenceLevel(0.95);
  calc.UseCLs(true);
  calc.SetVerbose(true);

  calc.SetFixedScan(5,startVal,stopVal);
  HypoTestInverterResult * res_toysCLs_calculator = calc.GetInterval();


  // dump results string to output file
  ofstream outStream ;
  outStream.open(outFile,ios::app) ;
  
  outStream << "CLs = " << res_toysCLs_calculator->UpperLimit() 
	    << "   CLs_exp = " << res_toysCLs_calculator->GetExpectedUpperLimit(0) 
	    << "   CLs_exp(-1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(-1) 
	    << "   CLs_exp(+1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(1) << endl ;
  
  outStream.close() ;

  return ;

}
