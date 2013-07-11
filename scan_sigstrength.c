
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "TAxis.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;

   void scan_sigstrength( const char* wsfile = "outputfiles/ws.root", double scanLow = 0., double scanHigh = 2.0, int nScanPoints = 40, double ymax = 9.  ) {

      gStyle->SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;
      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      ws->Print() ;

      RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;
      cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }

      RooAbsPdf* likelihood = ws->pdf("likelihood") ;
      if ( likelihood == 0x0 ) { printf("\n\n *** can't find likelihood in workspace.\n\n" ) ; return ; }
      printf("\n\n Likelihood:\n") ;
      likelihood -> Print() ;



      rv_sig_strength -> setConstant( kFALSE ) ;

      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0), Hesse(true), Strategy(1) ) ;
      double minNllSusyFloat = fitResult->minNll() ;
      double susy_ss_atMinNll = rv_sig_strength -> getVal() ;

      RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
      RooMsgService::instance().getStream(1).removeTopic(Fitting) ;



      double ssVals[200] ;
      double testStatVals[200] ;

      for ( int spi=0; spi<nScanPoints; spi++ ) {

         ssVals[spi] = scanLow + spi*(scanHigh-scanLow)/(nScanPoints-1.) ;

         rv_sig_strength -> setVal( ssVals[spi] ) ;
         rv_sig_strength -> setConstant( kTRUE ) ;

         //// RooFitResult* fitResult_sp = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         RooFitResult* fitResult_sp = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;

         double minNll_sp = fitResult_sp->minNll() ;

         testStatVals[spi] = 2.*( minNll_sp - minNllSusyFloat ) ;
         delete fitResult_sp ;

         printf("  %3d : %s = %6.1f,  test stat = %7.3f\n", spi, rv_sig_strength->GetName(), ssVals[spi], testStatVals[spi] ) ; fflush(stdout) ;

      } // spi.




      TGraph *ScanPlot = new TGraph( nScanPoints, ssVals, testStatVals ) ;
      ScanPlot -> SetName("sig_strength_PL_scan") ;

      ScanPlot->GetXaxis()->SetTitle("signal strength");
      ScanPlot->GetYaxis()->SetTitle("test statistic");

      ScanPlot->SetMarkerStyle(20);
      ScanPlot->SetMarkerColor(2);
      ScanPlot->SetLineColor(2);
      ScanPlot->SetLineWidth(3);

      ScanPlot->SetTitle("Profile likelihood scan of signal strength") ;

      TLine* line = new TLine() ;
      line->SetLineColor(4) ;

      //TCanvas *c0 = new TCanvas("c0","simple PL scan",700,600);

      //gPad->SetGridy(1) ;

      if ( ymax>0 ) {
         TH1F* hdummy = new TH1F("hdummy","Profile likelihood scan of signal strength",2,ssVals[0], ssVals[nScanPoints-1]) ;
         hdummy->GetXaxis()->SetTitle("signal strength");
         hdummy->GetYaxis()->SetTitle("test statistic");
         hdummy->SetMinimum(0.) ;
         hdummy->SetMaximum(ymax) ;
         hdummy->Draw() ;
         ScanPlot->Draw("C") ;
      } else {
         ScanPlot->Draw("AC");
      }



    //--- Use interpolation to find +/- 1 sigma and UL points.

      TString parName( rv_sig_strength->GetName() ) ;

      double tsAtPoiMin(99.) ;
      double poiMin(-1) ;
      double poiMinusOneSigma(-1.) ;
      double poiPlusOneSigma(-1.) ;
      double poiUL(-1.) ;
      int nsearch(1000) ;
      printf("\n\n") ;
      for ( int i=0; i<nsearch; i++ ) {
         double poi = scanLow + i*(scanHigh-scanLow)/(nsearch-1.) ;
         double tsVal = ScanPlot->Eval( poi, 0, "S" ) ;
         if ( poi < susy_ss_atMinNll ) {
            if ( tsVal <= 1.0 && poiMinusOneSigma < 0. ) {
               poiMinusOneSigma = poi ;
               printf(" -1 sigma   : %s = %7.2f\n", parName.Data(), poi) ;
            }
         } else {
            if ( tsVal >= 1.0 && poiPlusOneSigma < 0. ) {
               poiPlusOneSigma = poi ;
               printf(" +1 sigma   : %s = %7.2f\n", parName.Data(), poi) ;
            }
            if ( tsVal >= 2.71 && poiUL < 0. ) {
               poiUL = poi ;
               printf(" Upper lim. : %s = %7.2f\n", parName.Data(), poi) ;
            }
         }
         if ( tsVal < tsAtPoiMin ) {
            tsAtPoiMin = tsVal ;
            poiMin = poi ;
         }
         /// printf("  %4d : %s = %6.1f,  test stat = %7.3f\n", i, parName.Data(), poi, tsVal ) ;
      }
      printf(" Minimum    : %s = %7.2f\n", parName.Data(), poiMin ) ;
      double poiSignif(-1.) ;
      if ( scanLow <= 0 ) {
         double ts = ScanPlot->Eval( 0., 0, "S" ) ;
         if ( ts >=0. ) { poiSignif = sqrt( ts ) ; }
         printf(" Signif     : %5.2f\n", poiSignif ) ;
      }
      printf("\n\n") ;

      line->DrawLine(scanLow, 1.0, poiPlusOneSigma, 1.0) ;
      line->DrawLine(scanLow, 2.71, poiUL, 2.71) ;
      line->DrawLine(poiMinusOneSigma,0.,poiMinusOneSigma,1.0) ;
      line->DrawLine(poiPlusOneSigma,0.,poiPlusOneSigma,1.0) ;
      line->DrawLine(poiUL,0.,poiUL,2.71) ;

      TText* text = new TText() ;
      text->SetTextSize(0.05) ;
      char tstring[1000] ;

      sprintf( tstring, "SUSY sig. strength = %5.2f +%5.2f, -%5.2f", poiMin, poiPlusOneSigma-poiMin, poiMin-poiMinusOneSigma ) ;
      text->DrawTextNDC( 0.20, 0.80, tstring ) ;
      sprintf( tstring, "SUSY sig. strength 95%% CL 1-sided UL = %5.2f ", poiUL ) ;
      text->DrawTextNDC( 0.20, 0.70, tstring ) ;
      if ( scanLow == 0. && poiSignif>-1. ) {
         sprintf( tstring, "Significance = %5.2f", poiSignif ) ;
         text->DrawTextNDC( 0.20, 0.60, tstring ) ;
      }




   } // scan_sigstrength.





