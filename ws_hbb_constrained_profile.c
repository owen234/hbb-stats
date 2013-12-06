
#include "TRandom2.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"
#include "RooMinuit.h"

#include <string.h>

  using namespace RooFit;
  using namespace RooStats;

  //==============================================================================================

   void ws_hbb_constrained_profile( const char* wsfile = "outputfiles/ws-data-unblind-sigmass-250.root",
                                   const char* new_poi_name = "mu_bg_3b_msig_met4",
                                   const char* ignore_observable_name = "N_3b_msig_met4",
                                   int npoiPoints = 20,
                                   double poiMinVal = 0.1,
                                   double poiMaxVal = 5.0,
                                   double constraintWidth = 150.,
                                   double ymax = 5.,
                                   int verbLevel=0,
                                   bool floatSusy = false,
                                   bool do_stat_only = false ) {


     gStyle->SetOptStat(0) ;

     //--- make output directory.

     char command[10000] ;
     sprintf( command, "basename %s", wsfile ) ;
     TString wsfilenopath = gSystem->GetFromPipe( command ) ;
     wsfilenopath.ReplaceAll(".root","") ;
     char outputdirstr[1000] ;
     sprintf( outputdirstr, "outputfiles/scans-%s", wsfilenopath.Data() ) ;
     TString outputdir( outputdirstr ) ;

     printf("\n\n Creating output directory: %s\n\n", outputdir.Data() ) ;
     sprintf(command, "mkdir -p %s", outputdir.Data() ) ;
     gSystem->Exec( command ) ;


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;



       if ( verbLevel > 0 ) { printf("\n\n Verbose level : %d\n\n", verbLevel) ; }


       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

       if ( verbLevel > 0 ) { ws->Print() ; }






       RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;

       if ( verbLevel > 0 ) {
          printf("\n\n\n  ===== RooDataSet ====================\n\n") ;
          rds->Print() ;
          rds->printMultiline(cout, 1, kTRUE, "") ;
       }





       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;
       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
       if ( rv_sig_strength == 0x0 ) {
          printf("\n\n\n *** can't find sig_strength in workspace.  Quitting.\n\n\n") ;
          return ;
       }


       if ( floatSusy ) {
          rv_sig_strength->setVal(0.) ;
          rv_sig_strength->setConstant( kFALSE ) ;
       } else {
          rv_sig_strength->setVal(0.) ;
          rv_sig_strength->setConstant( kTRUE ) ;
       }










       //-- do a prefit.

       printf("\n\n\n ====== Pre fit with unmodified nll var.\n\n") ;

       RooFitResult* dataFitResultSusyFixed = likelihood->fitTo(*rds, Save(true),Hesse(false),Minos(false),Strategy(1),PrintLevel(verbLevel));
       int dataSusyFixedFitCovQual = dataFitResultSusyFixed->covQual() ;
       if ( dataSusyFixedFitCovQual < 2 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFixedFitCovQual ) ; return ; }
       double dataFitSusyFixedNll = dataFitResultSusyFixed->minNll() ;

       if ( verbLevel > 0 ) {
          dataFitResultSusyFixed->Print("v") ;
       }

       printf("\n\n Nll value, from fit result : %.3f\n\n", dataFitSusyFixedNll ) ;

       delete dataFitResultSusyFixed ;





       //-- If doing stat-only, fix all nuisance parameters (which is probably all the ones that start with prim_).
       if ( do_stat_only ) {
          const RooArgSet* np_as = modelConfig -> GetNuisanceParameters() ;
          printf("\n\n ===== List of nuisance parameters.  Fixing them.\n\n" ) ;
          TIterator* np_iter = np_as -> createIterator() ;
          while ( RooRealVar* par = (RooRealVar*) np_iter -> Next() ) {
             RooRealVar* rrv = ws -> var( par -> GetName() ) ;
             printf( " %25s : val = %.5f\n", rrv -> GetName(), rrv -> getVal() ) ;
             rrv -> setConstant( kTRUE ) ;
          } // par
          printf("\n\n") ;
       }






       //-- find the new POI parameter.
       RooAbsReal* new_poi_rar(0x0) ;

       new_poi_rar = ws->var( new_poi_name ) ;
       if ( new_poi_rar == 0x0 ) {
          printf("\n\n New POI %s is not a variable.  Trying function.\n\n", new_poi_name ) ;
          new_poi_rar = ws->function( new_poi_name ) ;
          if ( new_poi_rar == 0x0 ) {
             printf("\n\n New POI %s is not a function.  I quit.\n\n", new_poi_name ) ;
             return ;
          }
       } else {
          printf("\n\n     New POI %s is a variable with current value %.1f.\n\n", new_poi_name, new_poi_rar->getVal() ) ;
       }

       double startPoiVal = new_poi_rar->getVal() ;

       if ( npoiPoints <=0 ) {
          printf("\n\n Quitting now.\n\n" ) ;
          return ;
       }


      //--- output of createNLL is what minuit uses, so use that.
       RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;

       RooRealVar* rrv_poiValue = new RooRealVar( "poiValue", "poiValue", 0., -10000., 10000. ) ;

       RooRealVar* rrv_constraintWidth = new RooRealVar("constraintWidth","constraintWidth", 0.1, 0.1, 1000. ) ;
       rrv_constraintWidth -> setVal( constraintWidth ) ;
       rrv_constraintWidth -> setConstant(kTRUE) ;




       if ( verbLevel > 0 ) {
          printf("\n\n ======= debug likelihood print\n\n") ;
          likelihood->Print("v") ;
          printf("\n\n ======= debug nll print\n\n") ;
          nll->Print("v") ;
       }






    //----------------------------------------------------------------------------------------------

       RooMinuit* rminuit( 0x0 ) ;

       RooFormulaVar* plot_var( 0x0 ) ;


       RooRealVar* rrv_obs(0x0) ;

       if ( strcmp( ignore_observable_name, "none" ) != 0 ) {

          TString pdfNameTS( ignore_observable_name ) ;
          pdfNameTS.ReplaceAll("N_","pdf_") ;
          RooAbsPdf* pdf_ignore = ws -> pdf( pdfNameTS ) ;

          if ( pdf_ignore == 0x0 ) {
             printf("\n\n *** Told to ignore %s but can't find pdf %s.\n\n", ignore_observable_name, pdfNameTS.Data() ) ;
             return ;
          }

          rrv_obs = ws -> var( ignore_observable_name ) ;
          if ( rrv_obs == 0x0 ) {
             printf("\n\n *** Can't find RooRealVar for observable %s\n\n", ignore_observable_name ) ;
             return ;
          }
          rrv_obs->setConstant(kTRUE) ;


          printf("\n\n Ignoring observable %s by removing pdf %s\n\n", ignore_observable_name, pdfNameTS.Data() ) ;

          char minuit_formula[10000] ;
          sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)+log(%s)",
                 nll->GetName(),
                 rrv_constraintWidth->GetName(),
                 new_poi_rar->GetName(), rrv_poiValue->GetName(),
                 new_poi_rar->GetName(), rrv_poiValue->GetName(),
                 pdf_ignore->GetName()
                 ) ;

          printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
          RooFormulaVar* new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
              RooArgList( *nll,
                          *rrv_constraintWidth,
                          *new_poi_rar, *rrv_poiValue,
                          *new_poi_rar, *rrv_poiValue,
                          *pdf_ignore
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               new_minuit_var->getVal() ) ;

          rminuit = new RooMinuit( *new_minuit_var ) ;




          char plotvar_formula[10000] ;
          sprintf( plotvar_formula, "%s+log(%s)",
                 nll->GetName(),
                 pdf_ignore->GetName()
                 ) ;

          printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
          plot_var = new RooFormulaVar("plot_var", plotvar_formula,
              RooArgList( *nll,
                          *pdf_ignore
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               plot_var->getVal() ) ;


       } else {


          RooFormulaVar* new_minuit_var(0x0) ;

          char minuit_formula[10000] ;
          sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)",
              nll->GetName(),
              rrv_constraintWidth->GetName(),
              new_poi_rar->GetName(), rrv_poiValue->GetName(),
              new_poi_rar->GetName(), rrv_poiValue->GetName()
              ) ;
          printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
          new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
              RooArgList( *nll,
                       *rrv_constraintWidth,
                       *new_poi_rar, *rrv_poiValue,
                       *new_poi_rar, *rrv_poiValue
                       ) ) ;


          printf("\n\n Current value is %.2f\n\n",
               new_minuit_var->getVal() ) ;

          rminuit = new RooMinuit( *new_minuit_var ) ;




          char plotvar_formula[10000] ;
          sprintf( plotvar_formula, "%s",
                 nll->GetName()
                 ) ;

          printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
          plot_var = new RooFormulaVar("plot_var", plotvar_formula,
              RooArgList( *nll
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               plot_var->getVal() ) ;


       }


       rminuit->setPrintLevel(verbLevel-1) ;
       if ( verbLevel <=0 ) { rminuit->setNoWarn() ; }

    //----------------------------------------------------------------------------------------------

       //-- If POI range is -1 to -1, automatically determine the range using the set value.

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


       double poiVals[1000] ;
       double nllVals[1000] ;
       double minNllVal(1.e9) ;


       for ( int poivi=0; poivi < npoiPoints ; poivi++ ) {

          double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*(npoiPoints-1)) ;

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





          poiVals[poivi] = new_poi_rar->getVal() ;
          nllVals[poivi] = plot_var->getVal() ;

          if ( nllVals[poivi] < minNllVal ) { minNllVal = nllVals[poivi] ; }

          delete rfr ;


       } // poivi

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
    // printf("\n\n POI estimate :  %g  +%g  -%g    [%g,%g]\n\n",
    //         poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv ) ;
       printf("\n\n POI estimate :  %6.2f  + %6.2f  - %6.2f    [ %6.2f , %6.2f ]\n\n",
               poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv ) ;

       printf(" %s val,pm1sig,pm2sig: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
          new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv) ) ;

       char htitle[1000] ;
       sprintf(htitle, "%s profile likelihood scan: -2ln(L/Lm)", new_poi_name ) ;
       TH1F* hscan = new TH1F("hscan", htitle, 10, poiMinVal, poiMaxVal ) ;
       hscan->SetMinimum(0.) ;
       hscan->SetMaximum(ymax) ;


       hscan->Draw() ;
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

       if ( poiBest > 4 ) {
          sprintf( tstring, "%s = %.1f +%.1f -%.1f", new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
       } else {
          sprintf( tstring, "%s = %.2f +%.2f -%.2f", new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
       }
       text -> DrawTextNDC( 0.15, 0.85, tstring ) ;

       if ( poiBest > 4 ) {
          sprintf( tstring, "68%% interval [%.1f,  %.1f]", poiMinus1stdv, poiPlus1stdv ) ;
       } else {
          sprintf( tstring, "68%% interval [%.2f,  %.2f]", poiMinus1stdv, poiPlus1stdv ) ;
       }
       text -> DrawTextNDC( 0.15, 0.78, tstring ) ;


       char hname[1000] ;
       sprintf( hname, "hscanout_%s", new_poi_name ) ;
       TH1F* hsout = new TH1F( hname,"scan results",4,0.,4.) ;
       double obsVal(-1.) ;
       if ( rrv_obs != 0 ) { obsVal = rrv_obs->getVal() ; }
       hsout->SetBinContent(1, obsVal ) ;
       hsout->SetBinContent(2, poiPlus1stdv ) ;
       hsout->SetBinContent(3, poiBest ) ;
       hsout->SetBinContent(4, poiMinus1stdv ) ;
       TAxis* xaxis = hsout->GetXaxis() ;
       xaxis->SetBinLabel(1,"Observed val.") ;
       xaxis->SetBinLabel(2,"Model+1sd") ;
       xaxis->SetBinLabel(3,"Model") ;
       xaxis->SetBinLabel(4,"Model-1sd") ;

       char susyfloatstring[100] ;
       if ( floatSusy ) {
          sprintf( susyfloatstring, "susyfloat" ) ;
       } else {
          sprintf( susyfloatstring, "susyfixed" ) ;
       }
       char stat_syst_string[100] ;
       if ( do_stat_only ) {
          sprintf( stat_syst_string, "stat-only" ) ;
       } else {
          sprintf( stat_syst_string, "stat-syst" ) ;
       }

       char outrootfile[10000] ;
       sprintf( outrootfile, "%s/scan-%s-%s-%s.root", outputdir.Data(), susyfloatstring, new_poi_name, stat_syst_string ) ;

       char outpdffile[10000] ;
       sprintf( outpdffile, "%s/scan-%s-%s-%s.pdf", outputdir.Data(), susyfloatstring, new_poi_name, stat_syst_string ) ;

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

   }

