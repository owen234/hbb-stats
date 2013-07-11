//
//
//  This is a simple setup for one METsig bin, where signal is ignored everywhere except in
//  the 4b mass signal box and the 3b sample is not used (only 4b and 2b).
//
//
//
#include "TMath.h"
#include "TSystem.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooStats/ModelConfig.h"


#include "getFileValue.c"

   using namespace RooFit ;
   using namespace RooStats ;

   const int bins_of_nb(2) ;
   const int max_bins_of_met(50) ;

   int bin_nb_count[2] = { 2, 4 } ;
   int met_bin_index[1] ;

   RooArgSet* globalObservables ;
   RooArgSet* allNuisances ;
   RooArgSet* allNuisancePdfs ;

  //--- prototypes here.

   RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) ;



  //===========================================================================================

   void build_hbb_workspace_onebin( const char* infile = "outputfiles/input-file.txt", const char* outfile = "outputfiles/ws.root", int the_met_bin=1 ) {


    //-------------------------------------------------------------------------

     //-- Create workspace and other RooStats things.

      printf("\n\n Creating workspace.\n\n") ;

      RooWorkspace workspace("ws") ;
      workspace.autoImportClassCode(true) ;

      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      RooArgSet* observedParametersList = new RooArgSet("observables") ;




    //-------------------------------------------------------------------------

      printf("\n\n Reading input file: %s\n\n", infile ) ;

      float fileVal ;
      char pname[1000] ;
      char formula[1000] ;


     //-- only look at one met bin.
      sprintf( pname, "bins_of_met" ) ;
      int bins_of_met = 1 ;
      met_bin_index[0] = the_met_bin ;

      //-- save bins_of_met in the workspace for convenience.
      RooRealVar bom( "bins_of_met", "bins_of_met", bins_of_met, 0., 1000. ) ;
      bom.setConstant(kTRUE) ;
      workspace.import(bom) ;


      //-- save bins_of_nb in the workspace for convenience.
      RooRealVar bonb( "bins_of_nb", "bins_of_nb", bins_of_nb, 0., 1000. ) ;
      bonb.setConstant(kTRUE) ;
      workspace.import(bonb) ;


      RooRealVar* rv_N_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooRealVar* rv_N_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

      RooRealVar* rv_smc_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooRealVar* rv_smc_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

      RooAbsReal* rv_Rsigsb_corr[bins_of_nb][max_bins_of_met]  ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "N_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_N_msig[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_N_msig[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_N_msig[nbi][mbi] -> setConstant( kTRUE ) ;
            observedParametersList -> add( *rv_N_msig[nbi][mbi] ) ;

            sprintf( pname, "N_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_N_msb[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_N_msb[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_N_msb[nbi][mbi] -> setConstant( kTRUE ) ;
            observedParametersList -> add( *rv_N_msb[nbi][mbi] ) ;

            sprintf( pname, "smc_%db_msig_met%d", bin_nb_count[nbi], mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_smc_msig[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_smc_msig[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_smc_msig[nbi][mbi] -> setConstant( kTRUE ) ;

            sprintf( pname, "smc_%db_msb_met%d", bin_nb_count[nbi], mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_smc_msb[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_smc_msb[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_smc_msb[nbi][mbi] -> setConstant( kTRUE ) ;

            float corrVal, corrSyst ;
            sprintf( pname, "Rsigsb_syst_%db_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            if ( !getFileValue( infile, pname, corrSyst ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            sprintf( pname, "Rsigsb_corr_%db_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            if ( !getFileValue( infile, pname, corrVal  ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }

            rv_Rsigsb_corr[nbi][mbi] = makeLognormalConstraint( pname, corrVal, corrSyst ) ;


         } // mbi.

      } // nbi.

     //-- Finished reading input from file.

    //-------------------------------------------------------------------------

      printf("\n\n Creating and importing dataset into workspace.\n\n") ;

      RooDataSet* dsObserved = new RooDataSet("hbb_observed_rds", "hbb observed data values", *observedParametersList ) ;
      dsObserved -> add( *observedParametersList ) ;
      workspace.import( *dsObserved ) ;

    //-------------------------------------------------------------------------

     //-- Define all floats.

      printf("\n\n Defining all unconstrained floats (Ratios, signal strength).\n\n") ;

      double R_msigmsb_initialval(0.15) ;

      RooRealVar* rv_R_msigmsb[50] ;

      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

         sprintf( pname, "R_msigmsb_met%d", met_bin_index[mbi] ) ;
         printf( "  %s\n", pname ) ;
         rv_R_msigmsb[mbi] = new RooRealVar( pname, pname, R_msigmsb_initialval, 0., 3. ) ;
         rv_R_msigmsb[mbi] -> setConstant( kFALSE ) ;
         rv_R_msigmsb[mbi] -> Print() ;

      } // mbi.

      printf("\n") ;

      sprintf( pname, "sig_strength" ) ;
      RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 1.0, 0., 10. ) ;
      rv_sig_strength -> setConstant(kFALSE) ;
      rv_sig_strength -> Print() ;
      printf("  %s\n\n", pname ) ;

    //-------------------------------------------------------------------------

     //-- Define all mu parameters.

      printf("\n\n Defining mu parameters.\n\n") ;

      RooAbsReal* rv_mu_bg_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
      RooAbsReal* rv_mu_bg_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      RooAbsReal* rv_mu_sig_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooAbsReal* rv_mu_sig_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "mu_bg_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_bg_msb[nbi][mbi] = new RooRealVar( pname, pname, rv_N_msb[nbi][mbi] -> getVal(), 0., 1.e6 ) ;
            rv_mu_bg_msb[nbi][mbi] -> Print() ;

            sprintf( formula, "@0 * @1 * @2" ) ;
            sprintf( pname, "mu_bg_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_bg_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_Rsigsb_corr[nbi][mbi], *rv_R_msigmsb[mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
            rv_mu_bg_msig[nbi][mbi] -> Print() ;

            sprintf( formula, "@0 * @1" ) ;
            sprintf( pname, "mu_sig_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_sig_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_sig_strength, *rv_smc_msig[nbi][mbi] ) ) ;
            rv_mu_sig_msig[nbi][mbi] -> Print() ;

            sprintf( formula, "@0 * @1" ) ;
            sprintf( pname, "mu_sig_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_sig_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_sig_strength, *rv_smc_msb[nbi][mbi] ) ) ;
            rv_mu_sig_msb[nbi][mbi] -> Print() ;


         } // mbi.

      } // nbi.

     //-- Finished defining mu parameters.

    //-------------------------------------------------------------------------

     //-- Defining small n's

     printf("\n\n Defining small n's.\n\n") ; fflush(stdout) ;

     RooAbsReal* rv_n_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
     RooAbsReal* rv_n_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( formula, "@0 + @1" ) ;

            sprintf( pname, "n_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ; fflush(stdout) ;
            rv_n_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msig[nbi][mbi], *rv_mu_bg_msig[nbi][mbi] ) ) ;
            rv_n_msig[nbi][mbi] -> Print() ;
            workspace.import( *rv_n_msig[nbi][mbi] ) ;

            sprintf( pname, "n_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ; fflush(stdout) ;
            rv_n_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msb[nbi][mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
            rv_n_msb[nbi][mbi] -> Print() ;
            workspace.import( *rv_n_msb[nbi][mbi] ) ;

     ////// sprintf( formula, "@0" ) ;

     ////// sprintf( pname, "n_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
     ////// printf( "  %s\n", pname ) ; fflush(stdout) ;
     ////// rv_n_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_bg_msig[nbi][mbi] ) ) ;
     ////// rv_n_msig[nbi][mbi] -> Print() ;
     ////// workspace.import( *rv_n_msig[nbi][mbi] ) ;

     ////// sprintf( pname, "n_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
     ////// printf( "  %s\n", pname ) ; fflush(stdout) ;
     ////// rv_n_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_bg_msb[nbi][mbi] ) ) ;
     ////// rv_n_msb[nbi][mbi] -> Print() ;
     ////// workspace.import( *rv_n_msb[nbi][mbi] ) ;

         } // mbi.

      } // nbi.

    //-------------------------------------------------------------------------

     //-- Define the Poisson pdfs for the observables.

      printf("\n\n Defining Poisson pdfs for the observables.\n\n") ;

      RooAbsReal* rv_pdf_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
      RooAbsReal* rv_pdf_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      RooArgSet pdflist ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "pdf_%db_msig_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_pdf_msig[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msig[nbi][mbi], *rv_n_msig[nbi][mbi] ) ;
            rv_pdf_msig[nbi][mbi] -> Print() ;

            pdflist.add( *rv_pdf_msig[nbi][mbi] ) ;

            sprintf( pname, "pdf_%db_msb_met%d", bin_nb_count[nbi], met_bin_index[mbi] ) ;
            printf( "  %s\n", pname ) ;
            rv_pdf_msb[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msb[nbi][mbi], *rv_n_msb[nbi][mbi] ) ;
            rv_pdf_msb[nbi][mbi] -> Print() ;

            pdflist.add( *rv_pdf_msb[nbi][mbi] ) ;

         } // mbi.

      } // nbi.

    //-------------------------------------------------------------------------

     //-- Build the likelihood.

      printf("\n\n Building the likelihood.\n\n") ;

      pdflist.add( *allNuisancePdfs ) ;

      pdflist.Print() ;
      printf("\n") ;

      RooProdPdf* likelihood = new RooProdPdf( "likelihood", "hbb likelihood", pdflist ) ;
      likelihood->Print() ;


    //-------------------------------------------------------------------------


      printf("\n\n Running a test fit.\n\n") ;


      dsObserved -> Print() ;
      dsObserved -> printMultiline(cout, 1, kTRUE, "") ;


      printf("\n\n =============================================\n\n") ;
      likelihood -> fitTo( *dsObserved, PrintLevel(3), Hesse(0), Minos(0) ) ;
      printf("\n\n =============================================\n\n") ;


   // printf("\n\n Estimated background in signal bin %s is %.1f \n\n\n ", rv_mu_bg_msig[1][0] -> GetName(), rv_mu_bg_msig[1][0] -> getVal()  ) ;



   /////   workspace.import( *likelihood ) ;



     //-- Set up RooStats models.

      printf("\n\n Setting up S+B model.\n\n") ;

      RooArgSet poi( *rv_sig_strength, "poi" ) ;
      RooUniform signal_prior( "signal_prior", "signal_prior", *rv_sig_strength ) ;

      ModelConfig sbModel ("SbModel");
      sbModel.SetWorkspace( workspace ) ;
      sbModel.SetPdf( *likelihood ) ;
      sbModel.SetParametersOfInterest( poi );
      sbModel.SetPriorPdf(signal_prior);
      sbModel.SetObservables( *observedParametersList );
      sbModel.SetNuisanceParameters( *allNuisances );
      sbModel.SetGlobalObservables( *globalObservables );

      workspace.Print() ;

      printf("\n\n Doing fit for S+B model.\n" ) ; fflush(stdout) ;

      RooAbsReal* pNll = sbModel.GetPdf()->createNLL(*dsObserved);
      RooAbsReal* pProfile = pNll->createProfile(RooArgSet());
      pProfile->getVal();
      RooArgSet* pPoiAndNuisance = new RooArgSet();
      pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
      if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
      printf("\n\n Will save these parameter points that correspond to the fit to data.\n\n") ; fflush(stdout) ;
      pPoiAndNuisance->Print("v");
      sbModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (sbModel);

      delete pProfile ;
      delete pNll ;
      delete pPoiAndNuisance ;

      printf("\n\n Setting up BG-only model.\n\n") ;

      ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
      bModel.SetName("BModel");
      bModel.SetWorkspace(workspace);

      printf("\n\n Doing fit for BG-only model.\n" ) ; fflush(stdout) ;
      pNll = bModel.GetPdf()->createNLL(*dsObserved);
      pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
      ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
      pProfile->getVal();
      pPoiAndNuisance = new RooArgSet();
      pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
      if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
      printf("\n\n Should use these parameter points to generate pseudo data for bkg only.\n\n") ; fflush(stdout) ;
      pPoiAndNuisance->Print("v");
      bModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (bModel);

      delete pProfile ;
      delete pNll ;
      delete pPoiAndNuisance ;

      workspace.Print() ;

      printf("\n\n Saving workspace in : %s\n\n", outfile ) ;

      gSystem->Exec(" mkdir -p outputfiles " ) ;

      workspace.writeToFile( outfile ) ;




   } // build_hbb_workspace1.

  //==============================================================================================

    RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) {

       if ( NP_err <= 0. ) {
          printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char pname[1000] ;
       sprintf( pname, "prim_%s", NP_name ) ;

       printf(" makeLognormalConstraint : creating primary log-normal variable %s\n", pname ) ;
       RooRealVar* np_prim_rrv = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_rrv -> setVal( 0. ) ;
       np_prim_rrv -> setConstant( kFALSE ) ;

       sprintf( pname, "prim_mean_%s", NP_name ) ;
       RooRealVar* np_prim_mean = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_mean->setConstant(kTRUE) ;

       sprintf( pname, "prim_sigma_%s", NP_name ) ;
       RooConstVar* np_prim_sigma = new RooConstVar( pname, pname, 1. ) ;


       char pdfname[1000] ;
       sprintf( pdfname, "pdf_prim_%s", NP_name ) ;
       RooGaussian* np_prim_pdf = new RooGaussian( pdfname, pdfname, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

       allNuisances -> add( *np_prim_rrv ) ;
       allNuisancePdfs -> add( *np_prim_pdf ) ;
       globalObservables -> add( *np_prim_mean ) ;


       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       //-- compute the log-normal-distributed parameter from the primary parameter.

       //--- This is the new way.  RMS of lognormal is much closer to sigma when sigma is
       //    large, doing it this way.  When sigma/mean is small, they are about the same.
       //    That is, exp(sigma/mean) is close to (sigma/mean + 1).  This one is better when
       //    sigma/mean is not small.  The high-side tail is not as strong.
       //
        RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( ( @1/@0 + 1. ), @2)",
                  RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
       //------------------------------------------------------------------------------------------


       printf("  makeLognormalConstraint : created log-normal nuisance parameter %s : val = %g\n", NP_name, np_rfv -> getVal() ) ;

       return np_rfv ;


    } // makeLognormalConstraint.

  //==============================================================================================







