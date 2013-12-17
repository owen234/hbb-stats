
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
#include "TGraph2D.h"
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

   void ws_halfblind_2Dscan_profile( const char* wsfile = "outputfiles/ws-data-unblind.root",
                                   const char* new_poi_name_x = "mu_bg_4b_msig_met3",
                                   int npoiPoints_x = 21, // use odd values only
                                   double poiMinVal_x = 0.2,
                                   double poiMaxVal_x = 1.2,
                                   const char* new_poi_name_y = "mu_bg_3b_msig_met3",
                                   int npoiPoints_y = 21, // use odd values only
                                   double poiMinVal_y = 1.3,
                                   double poiMaxVal_y = 3.5,
                                   double constraintWidth = 50.,
                                   double ymax = 5.,
                                   bool do_stat_only = false,
                                   int verbLevel=0 ) {

     if ( npoiPoints_x%2 != 1 ) { printf("\n\n Try again with an odd value for npoiPoints_x.  You used %d\n\n", npoiPoints_x ) ; return ; }
     if ( npoiPoints_y%2 != 1 ) { printf("\n\n Try again with an odd value for npoiPoints_y.  You used %d\n\n", npoiPoints_y ) ; return ; }


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

       RooRealVar* rrv_sig_strength = ws->var("sig_strength") ;
       if ( rrv_sig_strength == 0x0 ) {
          printf("\n\n\n *** can't find sig_strength in workspace.  Quitting.\n\n\n") ;
          return ;
       }





       //-- do BG only.
       rrv_sig_strength->setVal(0.) ;
       rrv_sig_strength->setConstant( kTRUE ) ;









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




       //-- Construct the new POI parameter.
       RooAbsReal* new_poi_rar_x(0x0) ;

       new_poi_rar_x = ws->var( new_poi_name_x ) ;
       if ( new_poi_rar_x == 0x0 ) {
          printf("\n\n New POI_x %s is not a variable.  Trying function.\n\n", new_poi_name_x ) ;
          new_poi_rar_x = ws->function( new_poi_name_x ) ;
          if ( new_poi_rar_x == 0x0 ) {
             printf("\n\n New POI_x %s is not a function.  I quit.\n\n", new_poi_name_x ) ;
             return ;
          }
       } else {
          printf("\n\n     New POI_x %s is a variable with current value %.1f.\n\n", new_poi_name_x, new_poi_rar_x->getVal() ) ;
       }


       RooAbsReal* new_poi_rar_y(0x0) ;

       new_poi_rar_y = ws->var( new_poi_name_y ) ;
       if ( new_poi_rar_y == 0x0 ) {
          printf("\n\n New POI_y %s is not a variable.  Trying function.\n\n", new_poi_name_y ) ;
          new_poi_rar_y = ws->function( new_poi_name_y ) ;
          if ( new_poi_rar_y == 0x0 ) {
             printf("\n\n New POI_y %s is not a function.  I quit.\n\n", new_poi_name_y ) ;
             return ;
          }
       } else {
          printf("\n\n     New POI_y %s is a variable with current value %.1f.\n\n", new_poi_name_y, new_poi_rar_y->getVal() ) ;
       }







       if ( npoiPoints_x <=0 ) {
          printf("\n\n Quitting now.\n\n" ) ;
          return ;
       }


   //  double startPoiVal_x = new_poi_rar_x->getVal() ;
   //  double startPoiVal_y = new_poi_rar_y->getVal() ;



      //--- The RooNLLVar is NOT equivalent to what minuit uses.
  //   RooNLLVar* nll = new RooNLLVar("nll","nll", *likelihood, *rds ) ;
  //   printf("\n\n Nll value, from construction : %.3f\n\n", nll->getVal() ) ;

      //--- output of createNLL IS what minuit uses, so use that.
       RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;

       RooRealVar* rrv_poiValue_x = new RooRealVar( "poiValue_x", "poiValue_x", 0., -10000., 10000. ) ;
       RooRealVar* rrv_poiValue_y = new RooRealVar( "poiValue_y", "poiValue_y", 0., -10000., 10000. ) ;
   /// rrv_poiValue->setVal( poiMinVal ) ;
   /// rrv_poiValue->setConstant(kTRUE) ;

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

       char ignore_observable_name[50][100] ;

       sprintf( ignore_observable_name[0], "4b_msig_met1" ) ;
       sprintf( ignore_observable_name[1], "4b_msig_met2" ) ;
       sprintf( ignore_observable_name[2], "4b_msig_met3" ) ;
       sprintf( ignore_observable_name[3], "4b_msig_met4" ) ;

       sprintf( ignore_observable_name[4], "3b_msig_met1" ) ;
       sprintf( ignore_observable_name[5], "3b_msig_met2" ) ;
       sprintf( ignore_observable_name[6], "3b_msig_met3" ) ;
       sprintf( ignore_observable_name[7], "3b_msig_met4" ) ;

       int n_ignore_observable(8) ;

       RooAbsReal* rar_ignore_pdf[100] ;
       RooAbsReal* rar_ignore_obs[100] ;

       char ignoreTermFormula[10000] ;
       RooArgSet  ignorePdfList ;

       for ( int ii=0; ii < n_ignore_observable; ii++ ) {

          char name[100] ;

          sprintf( name, "pdf_%s", ignore_observable_name[ii] ) ;
          rar_ignore_pdf[ii] = ws -> pdf( name ) ;
          if ( rar_ignore_pdf[ii] == 0x0 ) {
             printf("\n\n\n *** Told to ignore %s but can't find %s\n\n", ignore_observable_name[ii], name ) ;
             return ;
          }

          sprintf( name, "N_%s", ignore_observable_name[ii] ) ;
          rar_ignore_obs[ii] = ws -> var( name ) ;
          ( (RooRealVar*) rar_ignore_obs[ii] ) -> setConstant(kTRUE) ; // probably not necessary, but can't hurt.

          if ( ii==0 ) {
             sprintf( ignoreTermFormula, "log(@%d)", ii ) ;
          } else {
             char buffer[10000] ;
             sprintf( buffer, "%s+log(@%d)", ignoreTermFormula, ii ) ;
             sprintf( ignoreTermFormula, "%s", buffer ) ;
          }

          ignorePdfList.add( *rar_ignore_pdf[ii] ) ;

       } // ii

       printf("\n\n Creating ignore formula var with : %s\n", ignoreTermFormula ) ;

       RooFormulaVar* rfv_ignorePdfTerm = new RooFormulaVar("ignorePdfTerm", ignoreTermFormula, ignorePdfList ) ;
       rfv_ignorePdfTerm -> Print() ;


       char minuit_formula_unbiased_unconstrained[100000] ;
       sprintf( minuit_formula_unbiased_unconstrained, "%s+%s", nll->GetName(), rfv_ignorePdfTerm->GetName() ) ;
       RooFormulaVar* rfv_minuitvar_unbiased_unconstrained = new RooFormulaVar( "minuitvar_unbiased_unconstrained",
         minuit_formula_unbiased_unconstrained,
         RooArgList( *nll, *rfv_ignorePdfTerm ) ) ;

       RooMinuit* rminuit_ub_uc = new RooMinuit( *rfv_minuitvar_unbiased_unconstrained  ) ;

       rminuit_ub_uc->setPrintLevel(verbLevel-1) ;
       rminuit_ub_uc->setNoWarn() ;

       // rminuit_ub_uc->migrad() ;
       // rminuit_ub_uc->hesse() ;

       RooFitResult* rfr_ub_uc = rminuit_ub_uc->fit("mr") ;

       double floatParInitVal[10000] ;
       char   floatParName[10000][100] ;
       int nFloatParInitVal(0) ;
       RooArgList ral_floats = rfr_ub_uc->floatParsFinal() ;
       TIterator* floatParIter = ral_floats.createIterator() ;
       {
          RooRealVar* par ;
          while ( (par = (RooRealVar*) floatParIter->Next()) ) {
             sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
             floatParInitVal[nFloatParInitVal] = par->getVal() ;
             nFloatParInitVal++ ;
          }
       }


       printf("\n\n Unbiased best value for new POI_x %s is : %7.1f\n\n", new_poi_rar_x->GetName(), new_poi_rar_x->getVal() ) ;
       printf("\n\n Unbiased best value for new POI_y %s is : %7.1f\n\n", new_poi_rar_y->GetName(), new_poi_rar_y->getVal() ) ;

       double best_unbiased_poi_val_x = new_poi_rar_x->getVal() ;
       double best_unbiased_poi_val_y = new_poi_rar_y->getVal() ;

     //-------


       char constraint_x_formula[1000] ;
       sprintf( constraint_x_formula, "%s*(%s-%s)*(%s-%s)",
          rrv_constraintWidth->GetName(),
          new_poi_rar_x->GetName(),
          rrv_poiValue_x->GetName(),
          new_poi_rar_x->GetName(),
          rrv_poiValue_x->GetName() ) ;
       printf("\n\n poi_x constraint formula: %s\n\n", constraint_x_formula ) ;
       RooFormulaVar* rfv_constraint_x = new RooFormulaVar( "rfv_constraint_x", constraint_x_formula,
            RooArgList( *rrv_constraintWidth,
                        *new_poi_rar_x, *rrv_poiValue_x,
                        *new_poi_rar_x, *rrv_poiValue_x ) ) ;

       char constraint_y_formula[1000] ;
       sprintf( constraint_y_formula, "%s*(%s-%s)*(%s-%s)",
          rrv_constraintWidth->GetName(),
          new_poi_rar_y->GetName(),
          rrv_poiValue_y->GetName(),
          new_poi_rar_y->GetName(),
          rrv_poiValue_y->GetName() ) ;
       printf("\n\n poi_y constraint formula: %s\n\n", constraint_y_formula ) ;
       RooFormulaVar* rfv_constraint_y = new RooFormulaVar( "rfv_constraint_y", constraint_y_formula,
            RooArgList( *rrv_constraintWidth,
                        *new_poi_rar_y, *rrv_poiValue_y,
                        *new_poi_rar_y, *rrv_poiValue_y ) ) ;


       char minuit_formula[10000] ;
       sprintf( minuit_formula, "%s+%s+%s+%s",
         nll->GetName(),
         rfv_constraint_x->GetName(),
         rfv_constraint_y->GetName(),
         rfv_ignorePdfTerm->GetName() ) ;

       printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;

       RooFormulaVar* new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
           RooArgList( *nll,
                       *rfv_constraint_x,
                       *rfv_constraint_y,
                       *rfv_ignorePdfTerm
                       ) ) ;

       printf("\n\n Current value is %.2f\n\n",
            new_minuit_var->getVal() ) ;

       rminuit = new RooMinuit( *new_minuit_var ) ;

       char plotvar_formula[10000] ;
       sprintf( plotvar_formula, "%s+%s",
              nll->GetName(),
              rfv_ignorePdfTerm->GetName()
              ) ;

       printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
       plot_var = new RooFormulaVar("plot_var", plotvar_formula,
           RooArgList( *nll,
                       *rfv_ignorePdfTerm
                       ) ) ;

       printf("\n\n Current value is %.2f\n\n",
            plot_var->getVal() ) ;




       rminuit->setPrintLevel(verbLevel-1) ;
       if ( verbLevel <=0 ) { rminuit->setNoWarn() ; }

    //----------------------------------------------------------------------------------------------

   //  //-- If POI range is -1 to -1, automatically determine the range using the set value.

   //  if ( poiMinVal < 0. && poiMaxVal < 0. ) {

   //     printf("\n\n Automatic determination of scan range.\n\n") ;

   //     if ( startPoiVal <= 0. ) {
   //        printf("\n\n *** POI starting value zero or negative %g.  Quit.\n\n\n", startPoiVal ) ;
   //        return ;
   //     }

   //     poiMinVal = startPoiVal - 3.5 * sqrt(startPoiVal) ;
   //     poiMaxVal = startPoiVal + 6.0 * sqrt(startPoiVal) ;

   //     if ( poiMinVal < 0. ) { poiMinVal = 0. ; }

   //     printf("    Start val = %g.   Scan range:   %g  to  %g\n\n", startPoiVal, poiMinVal, poiMaxVal ) ;


   //  }



    //----------------------------------------------------------------------------------------------


       double poiVals_x_scan_xDown_yDown[100][100] ;
       double poiVals_y_scan_xDown_yDown[100][100] ;
       double nllVals_scan_xDown_yDown[100][100] ;

       double poiVals_x_scan_xDown_yUp[100][100] ;
       double poiVals_y_scan_xDown_yUp[100][100] ;
       double nllVals_scan_xDown_yUp[100][100] ;

       double poiVals_x_scan_xUp_yDown[100][100] ;
       double poiVals_y_scan_xUp_yDown[100][100] ;
       double nllVals_scan_xUp_yDown[100][100] ;

       double poiVals_x_scan_xUp_yUp[100][100] ;
       double poiVals_y_scan_xUp_yUp[100][100] ;
       double nllVals_scan_xUp_yUp[100][100] ;

       //-- Do scan down from best value.


       double minNllVal(1.e9) ;




       //--- scan down in y

       for ( int poivi_y=0; poivi_y < (npoiPoints_y+1)/2 ; poivi_y++ ) {

          double poiValue_y = best_unbiased_poi_val_y - poivi_y*(best_unbiased_poi_val_y-poiMinVal_y)/(1.*((npoiPoints_y+1)/2-1)) ;

          rrv_poiValue_y -> setVal( poiValue_y ) ;
          rrv_poiValue_y -> setConstant( kTRUE ) ;


          //--- scan down in x

          printf("\n\n +++++ Starting scan down x for down y point %d.\n\n", poivi_y ) ;

          for ( int poivi_x=0; poivi_x < (npoiPoints_x+1)/2 ; poivi_x++ ) {

             double poiValue_x = best_unbiased_poi_val_x - poivi_x*(best_unbiased_poi_val_x-poiMinVal_x)/(1.*((npoiPoints_x+1)/2-1)) ;

             rrv_poiValue_x -> setVal( poiValue_x ) ;
             rrv_poiValue_x -> setConstant( kTRUE ) ;


          //+++++++++++++++++++++++++++++++++++

             rminuit->migrad() ;
             rminuit->hesse() ;
             RooFitResult* rfr = rminuit->save() ;

          //+++++++++++++++++++++++++++++++++++


             if ( verbLevel > 0 ) { rfr->Print("v") ; }


             float fit_minuit_var_val = rfr->minNll() ;

             printf(" %02d,%02d : poi_x, poi_y constraint = %.2f, %.2f : allvars : MinuitVar, createNLL, PV, POI_x, POI_y :    %.5f   %.5f   %.5f   %.5f , %.5f\n",
                   poivi_x, poivi_y, rrv_poiValue_x->getVal(), rrv_poiValue_y->getVal(),
                   fit_minuit_var_val, nll->getVal(), plot_var->getVal(),
                   new_poi_rar_x->getVal(), new_poi_rar_y->getVal()
                   ) ;
             cout << flush ;

             poiVals_x_scan_xDown_yDown[poivi_x][poivi_y] = new_poi_rar_x->getVal() ;
             poiVals_y_scan_xDown_yDown[poivi_x][poivi_y] = new_poi_rar_y->getVal() ;
             nllVals_scan_xDown_yDown[poivi_x][poivi_y] = plot_var->getVal() ;

             if ( nllVals_scan_xDown_yDown[poivi_x][poivi_y] < minNllVal ) { minNllVal = nllVals_scan_xDown_yDown[poivi_x][poivi_y] ; }

             delete rfr ;


          } // poivi_x, scan down




          //-- Refit for best unbiased value.

          printf("\n\n +++++ Resetting floats to best unbiased fit values.\n\n") ;

          for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
             RooRealVar* par = ws->var( floatParName[pi] ) ;
             par->setVal( floatParInitVal[pi] ) ;
          } // pi.




          printf("\n\n +++++ Starting scan up x for down y point %d.\n\n", poivi_y ) ;

          //--- scan up in x

          for ( int poivi_x=0; poivi_x < (npoiPoints_x+1)/2 ; poivi_x++ ) {

             double poiValue_x = best_unbiased_poi_val_x + poivi_x*(poiMaxVal_x - best_unbiased_poi_val_x)/(1.*((npoiPoints_x+1)/2-1)) ;

             rrv_poiValue_x -> setVal( poiValue_x ) ;
             rrv_poiValue_x -> setConstant( kTRUE ) ;


          //+++++++++++++++++++++++++++++++++++

             rminuit->migrad() ;
             rminuit->hesse() ;
             RooFitResult* rfr = rminuit->save() ;

          //+++++++++++++++++++++++++++++++++++


             if ( verbLevel > 0 ) { rfr->Print("v") ; }


             float fit_minuit_var_val = rfr->minNll() ;

             printf(" %02d,%02d : poi_x, poi_y constraint = %.2f, %.2f : allvars : MinuitVar, createNLL, PV, POI_x, POI_y :    %.5f   %.5f   %.5f   %.5f , %.5f\n",
                   poivi_x, poivi_y, rrv_poiValue_x->getVal(), rrv_poiValue_y->getVal(),
                   fit_minuit_var_val, nll->getVal(), plot_var->getVal(),
                   new_poi_rar_x->getVal(), new_poi_rar_y->getVal()
                   ) ;
             cout << flush ;

             poiVals_x_scan_xUp_yDown[poivi_x][poivi_y] = new_poi_rar_x->getVal() ;
             poiVals_y_scan_xUp_yDown[poivi_x][poivi_y] = new_poi_rar_y->getVal() ;
             nllVals_scan_xUp_yDown[poivi_x][poivi_y] = plot_var->getVal() ;

             if ( nllVals_scan_xUp_yDown[poivi_x][poivi_y] < minNllVal ) { minNllVal = nllVals_scan_xUp_yDown[poivi_x][poivi_y] ; }

             delete rfr ;


          } // poivi_x, scan up


       } // poivi_y, scan down



     //-----------------

          //-- Refit for best unbiased value.

          printf("\n\n +++++ Resetting floats to best unbiased fit values.\n\n") ;

          for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
             RooRealVar* par = ws->var( floatParName[pi] ) ;
             par->setVal( floatParInitVal[pi] ) ;
          } // pi.




       //--- scan up in y

       for ( int poivi_y=0; poivi_y < (npoiPoints_y+1)/2 ; poivi_y++ ) {

          double poiValue_y = best_unbiased_poi_val_y + poivi_y*(poiMaxVal_y - best_unbiased_poi_val_y)/(1.*((npoiPoints_y+1)/2-1)) ;

          rrv_poiValue_y -> setVal( poiValue_y ) ;
          rrv_poiValue_y -> setConstant( kTRUE ) ;


          printf("\n\n +++++ Starting scan down x for up y point %d.\n\n", poivi_y ) ;

          //--- scan down in x

          for ( int poivi_x=0; poivi_x < (npoiPoints_x+1)/2 ; poivi_x++ ) {

             double poiValue_x = best_unbiased_poi_val_x - poivi_x*(best_unbiased_poi_val_x-poiMinVal_x)/(1.*((npoiPoints_x+1)/2-1)) ;

             rrv_poiValue_x -> setVal( poiValue_x ) ;
             rrv_poiValue_x -> setConstant( kTRUE ) ;


          //+++++++++++++++++++++++++++++++++++

             rminuit->migrad() ;
             rminuit->hesse() ;
             RooFitResult* rfr = rminuit->save() ;

          //+++++++++++++++++++++++++++++++++++


             if ( verbLevel > 0 ) { rfr->Print("v") ; }


             float fit_minuit_var_val = rfr->minNll() ;

             printf(" %02d,%02d : poi_x, poi_y constraint = %.2f, %.2f : allvars : MinuitVar, createNLL, PV, POI_x, POI_y :    %.5f   %.5f   %.5f   %.5f , %.5f\n",
                   poivi_x, poivi_y, rrv_poiValue_x->getVal(), rrv_poiValue_y->getVal(),
                   fit_minuit_var_val, nll->getVal(), plot_var->getVal(),
                   new_poi_rar_x->getVal(), new_poi_rar_y->getVal()
                   ) ;
             cout << flush ;

             poiVals_x_scan_xDown_yUp[poivi_x][poivi_y] = new_poi_rar_x->getVal() ;
             poiVals_y_scan_xDown_yUp[poivi_x][poivi_y] = new_poi_rar_y->getVal() ;
             nllVals_scan_xDown_yUp[poivi_x][poivi_y] = plot_var->getVal() ;

             if ( nllVals_scan_xDown_yUp[poivi_x][poivi_y] < minNllVal ) { minNllVal = nllVals_scan_xDown_yUp[poivi_x][poivi_y] ; }

             delete rfr ;


          } // poivi_x, scan down




          //-- Refit for best unbiased value.

          printf("\n\n +++++ Resetting floats to best unbiased fit values.\n\n") ;

          for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
             RooRealVar* par = ws->var( floatParName[pi] ) ;
             par->setVal( floatParInitVal[pi] ) ;
          } // pi.




          printf("\n\n +++++ Starting scan up x for up y point %d.\n\n", poivi_y ) ;

          //--- scan up in x

          for ( int poivi_x=0; poivi_x < (npoiPoints_x+1)/2 ; poivi_x++ ) {

             double poiValue_x = best_unbiased_poi_val_x + poivi_x*(poiMaxVal_x - best_unbiased_poi_val_x)/(1.*((npoiPoints_x+1)/2-1)) ;

             rrv_poiValue_x -> setVal( poiValue_x ) ;
             rrv_poiValue_x -> setConstant( kTRUE ) ;


          //+++++++++++++++++++++++++++++++++++

             rminuit->migrad() ;
             rminuit->hesse() ;
             RooFitResult* rfr = rminuit->save() ;

          //+++++++++++++++++++++++++++++++++++


             if ( verbLevel > 0 ) { rfr->Print("v") ; }


             float fit_minuit_var_val = rfr->minNll() ;

             printf(" %02d,%02d : poi_x, poi_y constraint = %.2f, %.2f : allvars : MinuitVar, createNLL, PV, POI_x, POI_y :    %.5f   %.5f   %.5f   %.5f , %.5f\n",
                   poivi_x, poivi_y, rrv_poiValue_x->getVal(), rrv_poiValue_y->getVal(),
                   fit_minuit_var_val, nll->getVal(), plot_var->getVal(),
                   new_poi_rar_x->getVal(), new_poi_rar_y->getVal()
                   ) ;
             cout << flush ;

             poiVals_x_scan_xUp_yUp[poivi_x][poivi_y] = new_poi_rar_x->getVal() ;
             poiVals_y_scan_xUp_yUp[poivi_x][poivi_y] = new_poi_rar_y->getVal() ;
             nllVals_scan_xUp_yUp[poivi_x][poivi_y] = plot_var->getVal() ;

             if ( nllVals_scan_xUp_yUp[poivi_x][poivi_y] < minNllVal ) { minNllVal = nllVals_scan_xUp_yUp[poivi_x][poivi_y] ; }

             delete rfr ;


          } // poivi_x, scan up


       } // poivi_y, scan up









       double poiVals_x[10000] ;
       double poiVals_y[10000] ;
       double nllVals[10000] ;

       int npoiPoints(0) ;

       for ( int pyi = (npoiPoints_y+1)/2-1; pyi>=0 ; pyi-- ) {
          for ( int pxi= (npoiPoints_x+1)/2-1; pxi>=0 ; pxi-- ) {
             poiVals_x[npoiPoints] = poiVals_x_scan_xDown_yDown[pxi][pyi] ;
             poiVals_y[npoiPoints] = poiVals_y_scan_xDown_yDown[pxi][pyi] ;
             nllVals[npoiPoints] = nllVals_scan_xDown_yDown[pxi][pyi] ;
             npoiPoints++ ;
          }
          for ( int pxi=1 ; pxi < (npoiPoints_x+1)/2; pxi++ ) {
             poiVals_x[npoiPoints] = poiVals_x_scan_xUp_yDown[pxi][pyi] ;
             poiVals_y[npoiPoints] = poiVals_y_scan_xUp_yDown[pxi][pyi] ;
             nllVals[npoiPoints] = nllVals_scan_xUp_yDown[pxi][pyi] ;
             npoiPoints++ ;
          }
       }
       for ( int pyi=1;  pyi < (npoiPoints_y+1)/2; pyi++ ) {
          for ( int pxi= (npoiPoints_x+1)/2-1; pxi>=0 ; pxi-- ) {
             poiVals_x[npoiPoints] = poiVals_x_scan_xDown_yUp[pxi][pyi] ;
             poiVals_y[npoiPoints] = poiVals_y_scan_xDown_yUp[pxi][pyi] ;
             nllVals[npoiPoints] = nllVals_scan_xDown_yUp[pxi][pyi] ;
             npoiPoints++ ;
          }
          for ( int pxi=1 ; pxi < (npoiPoints_x+1)/2; pxi++ ) {
             poiVals_x[npoiPoints] = poiVals_x_scan_xUp_yUp[pxi][pyi] ;
             poiVals_y[npoiPoints] = poiVals_y_scan_xUp_yUp[pxi][pyi] ;
             nllVals[npoiPoints] = nllVals_scan_xUp_yUp[pxi][pyi] ;
             npoiPoints++ ;
          }
       }





       double nllDiffVals[10000] ;

  ///  double poiAtMinlnL(-1.) ;
  ///  double poiAtMinusDelta2(-1.) ;
  ///  double poiAtPlusDelta2(-1.) ;
       for ( int poivi=0; poivi < npoiPoints ; poivi++ ) {
          nllDiffVals[poivi] = 2.*(nllVals[poivi] - minNllVal) ;
       // double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*npoiPoints) ;
       // if ( nllDiffVals[poivi] < 0.01 ) { poiAtMinlnL = poiValue ; }
       // if ( poiAtMinusDelta2 < 0. && nllDiffVals[poivi] < 2.5 ) { poiAtMinusDelta2 = poiValue ; }
       // if ( poiAtMinlnL > 0. && poiAtPlusDelta2 < 0. && nllDiffVals[poivi] > 2.0 ) { poiAtPlusDelta2 = poiValue ; }
       } // poivi

  ///  printf("\n\n Estimates for poi at delta ln L = -2, 0, +2:  %g ,   %g ,   %g\n\n", poiAtMinusDelta2, poiAtMinlnL, poiAtPlusDelta2 ) ;

       printf("\n\n --- TGraph arrays:\n") ;
       for ( int i=0; i<npoiPoints; i++ ) {
          printf("  %2d : poi_x, poi_y = %6.1f, %6.1f   dnll = %g\n", i, poiVals_x[i], poiVals_y[i], nllDiffVals[i] ) ;
       }
       printf("\n\n") ;



      //--- Main canvas

       TCanvas* cscan = (TCanvas*) gDirectory->FindObject("cscan") ;
       if ( cscan == 0x0 ) {
          printf("\n Creating canvas.\n\n") ;
          cscan = new TCanvas("cscan","Delta nll") ;
       }


       char gname[1000] ;

       TGraph2D* graph = new TGraph2D( npoiPoints, poiVals_x, poiVals_y, nllDiffVals ) ;
       sprintf( gname, "scan_%s_vs_%s", new_poi_name_y, new_poi_name_x ) ;
       graph->SetName( gname ) ;

 //    double poiBest(-1.) ;
 //    double poiMinus1stdv(-1.) ;
 //    double poiPlus1stdv(-1.) ;
 //    double poiMinus2stdv(-1.) ;
 //    double poiPlus2stdv(-1.) ;
 //    double twoDeltalnLMin(1e9) ;

 //    int nscan(1000) ;
 //    for ( int xi=0; xi<nscan; xi++ ) {

 //       double x = poiVals[0] + xi*(poiVals[npoiPoints-1]-poiVals[0])/(nscan-1) ;

 //       double twoDeltalnL = graph -> Eval( x, 0, "S" ) ;

 //       if ( poiMinus1stdv < 0. && twoDeltalnL < 1.0 ) { poiMinus1stdv = x ; printf(" set m1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
 //       if ( poiMinus2stdv < 0. && twoDeltalnL < 4.0 ) { poiMinus2stdv = x ; printf(" set m2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
 //       if ( twoDeltalnL < twoDeltalnLMin ) { poiBest = x ; twoDeltalnLMin = twoDeltalnL ; }
 //       if ( twoDeltalnLMin < 0.3 && poiPlus1stdv < 0. && twoDeltalnL > 1.0 ) { poiPlus1stdv = x ; printf(" set p1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
 //       if ( twoDeltalnLMin < 0.3 && poiPlus2stdv < 0. && twoDeltalnL > 4.0 ) { poiPlus2stdv = x ; printf(" set p2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}

 //       if ( xi%100 == 0 ) { printf( " %4d : poi=%6.2f,  2DeltalnL = %6.2f\n", xi, x, twoDeltalnL ) ; }

 //    }
 //    printf("\n\n POI estimate :  %g  +%g  -%g    [%g,%g],   two sigma errors: +%g  -%g   [%g,%g]\n\n",
 //            poiBest,
 //            (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv,
 //            (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv), poiMinus2stdv, poiPlus2stdv
 //            ) ;

 //    printf(" %s val,pm1sig,pm2sig: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
 //       new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv) ) ;

 //    char htitle[1000] ;
 //    sprintf(htitle, "%s profile likelihood scan: -2ln(L/Lm)", new_poi_name ) ;
 //    TH1F* hscan = new TH1F("hscan", htitle, 10, poiMinVal, poiMaxVal ) ;
 //    hscan->SetMinimum(0.) ;
 //    hscan->SetMaximum(ymax) ;


 //    hscan->DrawCopy() ;
 //    graph->SetLineColor(4) ;
 //    graph->SetLineWidth(3) ;
 //    graph->Draw("CP") ;
 //    gPad->SetGridx(1) ;
 //    gPad->SetGridy(1) ;
 //    cscan->Update() ;

       graph -> Draw("surf1") ;

 //    TLine* line = new TLine() ;
 //    line->SetLineColor(2) ;
 //    line->DrawLine(poiMinVal, 1., poiPlus1stdv, 1.) ;
 //    line->DrawLine(poiMinus1stdv,0., poiMinus1stdv, 1.) ;
 //    line->DrawLine(poiPlus1stdv ,0., poiPlus1stdv , 1.) ;

 //    TText* text = new TText() ;
 //    text->SetTextSize(0.04) ;
 //    char tstring[1000] ;

 //    sprintf( tstring, "%s = %.1f +%.1f -%.1f", new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
 //    text -> DrawTextNDC( 0.15, 0.85, tstring ) ;

 //    sprintf( tstring, "68%% interval [%.1f,  %.1f]", poiMinus1stdv, poiPlus1stdv ) ;
 //    text -> DrawTextNDC( 0.15, 0.78, tstring ) ;


 //    char hname[1000] ;
 //    sprintf( hname, "hscanout_%s", new_poi_name ) ;
 //    TH1F* hsout = new TH1F( hname,"scan results",4,0.,4.) ;
 //    double obsVal(-1.) ;
 //    hsout->SetBinContent(1, obsVal ) ;
 //    hsout->SetBinContent(2, poiPlus1stdv ) ;
 //    hsout->SetBinContent(3, poiBest ) ;
 //    hsout->SetBinContent(4, poiMinus1stdv ) ;
 //    TAxis* xaxis = hsout->GetXaxis() ;
 //    xaxis->SetBinLabel(1,"Observed val.") ;
 //    xaxis->SetBinLabel(2,"Model+1sd") ;
 //    xaxis->SetBinLabel(3,"Model") ;
 //    xaxis->SetBinLabel(4,"Model-1sd") ;

       char outrootfile[10000] ;
       char outpdffile[10000] ;
       if ( !do_stat_only ) {
          sprintf( outrootfile, "%s/scan-hb-%s-vs-%s.root", outputdir.Data(), new_poi_name_y, new_poi_name_x ) ;
          sprintf( outpdffile, "%s/scan-hb-%s-vs-%s.pdf", outputdir.Data(), new_poi_name_y, new_poi_name_x ) ;
       } else {
          sprintf( outrootfile, "%s/scan-hb-%s-vs-%s-stat-only.root", outputdir.Data(), new_poi_name_y, new_poi_name_x ) ;
          sprintf( outpdffile, "%s/scan-hb-%s-vs-%s-stat-only.pdf", outputdir.Data(), new_poi_name_y, new_poi_name_x ) ;
       }


       cscan->Update() ; cscan->Draw() ;

       printf("\n Saving %s\n", outpdffile ) ;
       cscan->SaveAs( outpdffile ) ;



     //--- save in root file

       printf("\n Saving %s\n", outrootfile ) ;
       TFile fout(outrootfile,"recreate") ;
       graph->Write() ;
   //  hsout->Write() ;
       fout.Close() ;

       delete ws ;
       wstf->Close() ;

   }

