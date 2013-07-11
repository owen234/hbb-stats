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

void ws_run_cls_asymptotic( const char* wsfilename, double testVal, TString outFile ) {

   TFile* wsfile = TFile::Open( wsfilename ) ;
   if ( wsfile==0x0 ) {
      printf("\n\n\n *** problem opening workspace file %s\n\n", wsfilename ) ;
      return ;
   } else {
      printf("\n\n Opened workspace file %s\n", wsfilename ) ;
   }

   RooWorkspace* w = (RooWorkspace*) wsfile -> Get("ws") ;
   if ( w==0x0 ) {
      printf("\n\n\n *** Did not find workspace in %s\n\n\n", wsfilename ) ;
      return ;
   } else {
      printf(" Found workspace.\n" ) ;
   }

   ModelConfig* mc = (ModelConfig*) w -> obj( "SbModel" ) ;
   if ( mc==0x0 ) {
      printf("\n\n\n *** Did not find SbModel ModelConfig.\n\n") ;
      return ;
   } else {
      printf(" Found SbModel ModelConfig.\n") ;
   }

   RooAbsData* data = w -> data( "hbb_observed_rds" ) ;
   if ( data==0x0 ) {
      printf("\n\n\n *** Did not find RooDataset hbb_observed_rds.\n\n") ;
      return ;
   } else {
      printf(" Found RooDataset hbb_observed_rds.\n") ;
   }

   RooRealVar* myPOI = (RooRealVar*) mc -> GetParametersOfInterest() -> first();
   if ( myPOI==0x0 ) {
      printf("\n\n\n *** didn't find POI in ModelConfig.\n\n\n") ;
      return ;
   } else {
      printf(" Found POI in ModelConfig.\n") ;
   }
   myPOI -> setRange( 0., 100. ) ;

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

  SamplingDistribution * s = res_asymptotic_calculator->GetExpectedPValueDist(0);
  if ( s==0x0 ) {
     printf("\n\n\n *** Problem getting SamplingDistribution from res_asymptotic_calculator.\n\n") ;
     return ;
  } else {
     printf("\n\n Have SamplingDistribution from res_asymptotic_calculator.\n") ;
  }
  const std::vector<double> & values = s->GetSamplingDistribution();

  double maxSigma = res_asymptotic_calculator->fgAsymptoticMaxSigma;
  double dsig = 2* maxSigma/ (values.size() -1) ;

  int  im = (int) TMath::Floor ( ( -1 + maxSigma )/dsig + 0.5);
  int  i0 = (int) TMath::Floor ( ( maxSigma )/dsig + 0.5);
  int  ip = (int) TMath::Floor ( ( +1 + maxSigma )/dsig + 0.5);

  double CLs_expM = values[im];
  double CLs_exp0 = values[i0];
  double CLs_expP = values[ip];

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



} // ws_run_cls_asymptotic.



