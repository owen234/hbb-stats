
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
#include "RooMinuit.h"
#include "TAxis.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;

   void constrained_scan( const char* wsfile = "outputfiles/ws.root",
                          const char* new_poi_name="mu_bg_4b_msig_met1",
                          double constraintWidth=1.5,
                          int npoiPoints = 20,
                          double poiMinVal = 0.,
                          double poiMaxVal = 10.0,
                          double ymax = 9.,
                          int verbLevel=1  ) {

      TString outputdir("outputfiles") ;

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



      /////rv_sig_strength -> setConstant( kFALSE ) ;
      rv_sig_strength -> setVal(0.) ;
      rv_sig_strength -> setConstant( kTRUE ) ;

      likelihood->fitTo( *rds, Save(false), PrintLevel(0), Hesse(true), Strategy(1) ) ;
      //RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0), Hesse(true), Strategy(1) ) ;
      //double minNllSusyFloat = fitResult->minNll() ;
      //double susy_ss_atMinNll = rv_sig_strength -> getVal() ;

      RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
      RooMsgService::instance().getStream(1).removeTopic(Fitting) ;



     //-- Construct the new POI parameter.
      RooAbsReal* new_poi_rar(0x0) ;

      new_poi_rar = ws->var( new_poi_name ) ;
      if ( new_poi_rar == 0x0 ) {
         printf("\n\n New POI %s is not a variable.  Trying function.\n\n", new_poi_name ) ;
         new_poi_rar = ws->function( new_poi_name ) ;
         if ( new_poi_rar == 0x0 ) {
            printf("\n\n New POI %s is not a function.  I quit.\n\n", new_poi_name ) ;
            return ;
         } else {
            printf("\n Found it.\n\n") ;
         }
      } else {
         printf("\n\n     New POI %s is a variable with current value %.1f.\n\n", new_poi_name, new_poi_rar->getVal() ) ;
      }

       double startPoiVal = new_poi_rar->getVal() ;


       RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;

       RooRealVar* rrv_poiValue = new RooRealVar( "poiValue", "poiValue", 0., -10000., 10000. ) ;

       RooRealVar* rrv_constraintWidth = new RooRealVar("constraintWidth","constraintWidth", 0.1, 0.1, 1000. ) ;
       rrv_constraintWidth -> setVal( constraintWidth ) ;
       rrv_constraintWidth -> setConstant(kTRUE) ;


       RooMinuit* rminuit( 0x0 ) ;


       RooMinuit* rminuit_uc = new RooMinuit( *nll  ) ;

       rminuit_uc->setPrintLevel(verbLevel-1) ;
       rminuit_uc->setNoWarn() ;

       rminuit_uc->migrad() ;
       rminuit_uc->hesse() ;

       RooFitResult* rfr_uc = rminuit_uc->fit("mr") ;

       double floatParInitVal[10000] ;
       char   floatParName[10000][100] ;
       int nFloatParInitVal(0) ;
       RooArgList ral_floats = rfr_uc->floatParsFinal() ;
       TIterator* floatParIter = ral_floats.createIterator() ;
       {
          RooRealVar* par ;
          while ( (par = (RooRealVar*) floatParIter->Next()) ) {
             sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
             floatParInitVal[nFloatParInitVal] = par->getVal() ;
             nFloatParInitVal++ ;
          }
       }


       printf("\n\n Unbiased best value for new POI %s is : %7.1f\n\n", new_poi_rar->GetName(), new_poi_rar->getVal() ) ;
       double best_poi_val = new_poi_rar->getVal() ;

       char minuit_formula[10000] ;
       sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)",
         nll->GetName(),
         rrv_constraintWidth->GetName(),
         new_poi_rar->GetName(), rrv_poiValue->GetName(),
         new_poi_rar->GetName(), rrv_poiValue->GetName()
          ) ;

       printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
       RooFormulaVar* new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
           RooArgList( *nll,
                       *rrv_constraintWidth,
                       *new_poi_rar, *rrv_poiValue,
                       *new_poi_rar, *rrv_poiValue
                       ) ) ;

       printf("\n\n Current value is %.2f\n\n",
            new_minuit_var->getVal() ) ;

       rminuit = new RooMinuit( *new_minuit_var ) ;


       RooAbsReal* plot_var = nll ;

       printf("\n\n Current value is %.2f\n\n",
            plot_var->getVal() ) ;


       rminuit->setPrintLevel(verbLevel-1) ;
       if ( verbLevel <=0 ) { rminuit->setNoWarn() ; }


       if ( poiMinVal < 0. && poiMaxVal < 0. ) {

          printf("\n\n Automatic determination of scan range.\n\n") ;

          if ( startPoiVal <= 0. ) {
             printf("\n\n *** POI starting value zero or negative %g.  Quit.\n\n\n", startPoiVal ) ;
             return ;
          }

          poiMinVal = startPoiVal - 3.5 * sqrt(startPoiVal) ;
          poiMaxVal = startPoiVal + 6.0 * sqrt(startPoiVal) ;

          if ( poiMinVal < 0. ) { poiMinVal = 0. ; }

          printf("    Start val = %g.   Scan range:   %g  to  %g\n\n", startPoiVal, poiMinVal, poiMaxVal ) ;


       }



    //----------------------------------------------------------------------------------------------


       double poiVals_scanDown[1000] ;
       double nllVals_scanDown[1000] ;

       //-- Do scan down from best value.

       printf("\n\n +++++ Starting scan down from best value.\n\n") ;

       double minNllVal(1.e9) ;

       for ( int poivi=0; poivi < npoiPoints/2 ; poivi++ ) {

          ////double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*(npoiPoints-1)) ;
          double poiValue = best_poi_val - poivi*(best_poi_val-poiMinVal)/(1.*(npoiPoints/2-1)) ;

          rrv_poiValue -> setVal( poiValue ) ;
          rrv_poiValue -> setConstant( kTRUE ) ;


       //+++++++++++++++++++++++++++++++++++

          rminuit->migrad() ;
          rminuit->hesse() ;
          RooFitResult* rfr = rminuit->save() ;

       //+++++++++++++++++++++++++++++++++++


          if ( verbLevel > 0 ) { rfr->Print("v") ; }


          float fit_minuit_var_val = rfr->minNll() ;

          printf(" %02d : poi constraint = %.2f : allvars : MinuitVar, createNLL, PV, POI :    %.5f   %.5f   %.5f   %.5f\n",
                poivi, rrv_poiValue->getVal(), fit_minuit_var_val, nll->getVal(), plot_var->getVal(), new_poi_rar->getVal() ) ;
          cout << flush ;



          poiVals_scanDown[poivi] = new_poi_rar->getVal() ;
          nllVals_scanDown[poivi] = plot_var->getVal() ;

          if ( nllVals_scanDown[poivi] < minNllVal ) { minNllVal = nllVals_scanDown[poivi] ; }

          delete rfr ;


       } // poivi


       printf("\n\n +++++ Resetting floats to best fit values.\n\n") ;

       for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
          RooRealVar* par = ws->var( floatParName[pi] ) ;
          par->setVal( floatParInitVal[pi] ) ;
       } // pi.

       printf("\n\n +++++ Starting scan up from best value.\n\n") ;

      //-- Now do scan up.

       double poiVals_scanUp[1000] ;
       double nllVals_scanUp[1000] ;

       for ( int poivi=0; poivi < npoiPoints/2 ; poivi++ ) {

          double poiValue = best_poi_val + poivi*(poiMaxVal-best_poi_val)/(1.*(npoiPoints/2-1)) ;

          rrv_poiValue -> setVal( poiValue ) ;
          rrv_poiValue -> setConstant( kTRUE ) ;


       //+++++++++++++++++++++++++++++++++++

          rminuit->migrad() ;
          rminuit->hesse() ;
          RooFitResult* rfr = rminuit->save() ;

       //+++++++++++++++++++++++++++++++++++


          if ( verbLevel > 0 ) { rfr->Print("v") ; }


          float fit_minuit_var_val = rfr->minNll() ;

          printf(" %02d : poi constraint = %.2f : allvars : MinuitVar, createNLL, PV, POI :    %.5f   %.5f   %.5f   %.5f\n",
                poivi, rrv_poiValue->getVal(), fit_minuit_var_val, nll->getVal(), plot_var->getVal(), new_poi_rar->getVal() ) ;
          cout << flush ;

          poiVals_scanUp[poivi] = new_poi_rar->getVal() ;
          nllVals_scanUp[poivi] = plot_var->getVal() ;

          if ( nllVals_scanUp[poivi] < minNllVal ) { minNllVal = nllVals_scanUp[poivi] ; }

          delete rfr ;


       } // poivi





       double poiVals[1000] ;
       double nllVals[1000] ;

       int pointCount(0) ;
       for ( int pi=0; pi<npoiPoints/2; pi++ ) {
          poiVals[pi] = poiVals_scanDown[(npoiPoints/2-1)-pi] ;
          nllVals[pi] = nllVals_scanDown[(npoiPoints/2-1)-pi] ;
          pointCount++ ;
       }
       for ( int pi=1; pi<npoiPoints/2; pi++ ) {
          poiVals[pointCount] = poiVals_scanUp[pi] ;
          nllVals[pointCount] = nllVals_scanUp[pi] ;
          pointCount++ ;
       }
       npoiPoints = pointCount ;

       printf("\n\n --- TGraph arrays:\n") ;
       for ( int i=0; i<npoiPoints; i++ ) {
          printf("  %2d : poi = %6.1f, nll = %g\n", i, poiVals[i], nllVals[i] ) ;
       }
       printf("\n\n") ;

       double nllDiffVals[1000] ;

       double poiAtMinlnL(-1.) ;
       double poiAtMinusDelta2(-1.) ;
       double poiAtPlusDelta2(-1.) ;
       for ( int poivi=0; poivi < npoiPoints ; poivi++ ) {
          nllDiffVals[poivi] = 2.*(nllVals[poivi] - minNllVal) ;
          double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*npoiPoints) ;
          if ( nllDiffVals[poivi] < 0.01 ) { poiAtMinlnL = poiValue ; }
          if ( poiAtMinusDelta2 < 0. && nllDiffVals[poivi] < 2.5 ) { poiAtMinusDelta2 = poiValue ; }
          if ( poiAtMinlnL > 0. && poiAtPlusDelta2 < 0. && nllDiffVals[poivi] > 2.0 ) { poiAtPlusDelta2 = poiValue ; }
       } // poivi

       printf("\n\n Estimates for poi at delta ln L = -2, 0, +2:  %g ,   %g ,   %g\n\n", poiAtMinusDelta2, poiAtMinlnL, poiAtPlusDelta2 ) ;




      //--- Main canvas

       TCanvas* cscan = (TCanvas*) gDirectory->FindObject("cscan") ;
       if ( cscan == 0x0 ) {
          printf("\n Creating canvas.\n\n") ;
          cscan = new TCanvas("cscan","Delta nll") ;
       }


       char gname[1000] ;

       TGraph* graph = new TGraph( npoiPoints, poiVals, nllDiffVals ) ;
       sprintf( gname, "scan_%s", new_poi_name ) ;
       graph->SetName( gname ) ;

       double poiBest(-1.) ;
       double poiMinus1stdv(-1.) ;
       double poiPlus1stdv(-1.) ;
       double poiMinus2stdv(-1.) ;
       double poiPlus2stdv(-1.) ;
       double twoDeltalnLMin(1e9) ;

       int nscan(1000) ;
       for ( int xi=0; xi<nscan; xi++ ) {

          double x = poiVals[0] + xi*(poiVals[npoiPoints-1]-poiVals[0])/(nscan-1) ;

          double twoDeltalnL = graph -> Eval( x, 0, "S" ) ;

          if ( poiMinus1stdv < 0. && twoDeltalnL < 1.0 ) { poiMinus1stdv = x ; printf(" set m1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( poiMinus2stdv < 0. && twoDeltalnL < 4.0 ) { poiMinus2stdv = x ; printf(" set m2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( twoDeltalnL < twoDeltalnLMin ) { poiBest = x ; twoDeltalnLMin = twoDeltalnL ; }
          if ( twoDeltalnLMin < 0.3 && poiPlus1stdv < 0. && twoDeltalnL > 1.0 ) { poiPlus1stdv = x ; printf(" set p1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( twoDeltalnLMin < 0.3 && poiPlus2stdv < 0. && twoDeltalnL > 4.0 ) { poiPlus2stdv = x ; printf(" set p2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}

          if ( xi%100 == 0 ) { printf( " %4d : poi=%6.2f,  2DeltalnL = %6.2f\n", xi, x, twoDeltalnL ) ; }

       }
       printf("\n\n POI estimate :  %g  +%g  -%g    [%g,%g],   two sigma errors: +%g  -%g   [%g,%g]\n\n",
               poiBest,
               (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv,
               (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv), poiMinus2stdv, poiPlus2stdv
               ) ;

       printf(" %s val,pm1sig,pm2sig: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
          new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv) ) ;

       char htitle[1000] ;
       sprintf(htitle, "%s profile likelihood scan: -2ln(L/Lm)", new_poi_name ) ;
       TH1F* hscan = new TH1F("hscan", htitle, 10, poiMinVal, poiMaxVal ) ;
       hscan->SetMinimum(0.) ;
       hscan->SetMaximum(ymax) ;


       hscan->DrawCopy() ;
       graph->SetLineColor(4) ;
       graph->SetLineWidth(3) ;
       graph->Draw("CP") ;
       gPad->SetGridx(1) ;
       gPad->SetGridy(1) ;
       cscan->Update() ;

       TLine* line = new TLine() ;
       line->SetLineColor(2) ;
       line->DrawLine(poiMinVal, 1., poiPlus1stdv, 1.) ;
       line->DrawLine(poiMinus1stdv,0., poiMinus1stdv, 1.) ;
       line->DrawLine(poiPlus1stdv ,0., poiPlus1stdv , 1.) ;

       TText* text = new TText() ;
       text->SetTextSize(0.04) ;
       char tstring[1000] ;

       sprintf( tstring, "%s = %.1f +%.1f -%.1f", new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
       text -> DrawTextNDC( 0.15, 0.85, tstring ) ;

       sprintf( tstring, "68%% interval [%.1f,  %.1f]", poiMinus1stdv, poiPlus1stdv ) ;
       text -> DrawTextNDC( 0.15, 0.78, tstring ) ;


       char hname[1000] ;
       sprintf( hname, "hscanout_%s", new_poi_name ) ;
       TH1F* hsout = new TH1F( hname,"scan results",4,0.,4.) ;
       double obsVal(-1.) ;
       hsout->SetBinContent(1, obsVal ) ;
       hsout->SetBinContent(2, poiPlus1stdv ) ;
       hsout->SetBinContent(3, poiBest ) ;
       hsout->SetBinContent(4, poiMinus1stdv ) ;
       TAxis* xaxis = hsout->GetXaxis() ;
       xaxis->SetBinLabel(1,"Observed val.") ;
       xaxis->SetBinLabel(2,"Model+1sd") ;
       xaxis->SetBinLabel(3,"Model") ;
       xaxis->SetBinLabel(4,"Model-1sd") ;

       char outrootfile[10000] ;
       sprintf( outrootfile, "%s/scan-ff-%s.root", outputdir.Data(), new_poi_name ) ;

       char outpdffile[10000] ;
       sprintf( outpdffile, "%s/scan-ff-%s.pdf", outputdir.Data(), new_poi_name ) ;

       cscan->Update() ; cscan->Draw() ;

       printf("\n Saving %s\n", outpdffile ) ;
       cscan->SaveAs( outpdffile ) ;



     //--- save in root file

       printf("\n Saving %s\n", outrootfile ) ;
       TFile fout(outrootfile,"recreate") ;
       graph->Write() ;
       hsout->Write() ;
       fout.Close() ;

       delete ws ;
       wstf->Close() ;




   } // constrained_scan.





