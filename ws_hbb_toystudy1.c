
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TRandom3.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"

#include <string.h>

  using std::cout ;

  using namespace RooFit;
  using namespace RooStats;

  float inputval_n_msig[3][4] ; // first index is nb, second is metsig
  float inputval_n_msb[3][4]  ; // first index is nb, second is metsig

  float susyval_n_msig[3][4] ; // first index is nb, second is metsig
  float susyval_n_msb[3][4]  ; // first index is nb, second is metsig


  TRandom3* tran(0x0) ;

  RooDataSet* genToyData( RooWorkspace* ws, float gen_sig_strength ) ;


  //------------

   void ws_hbb_toystudy1( const char* wsfile = "outputfiles/ws-metsig-4metbin-w3b-wpu-csyst5-nosignal-sigmass-350-withMSB1.root",
                          int nToys_per_ssval=10,
                          int n_ssvals = 1,
                          float gen_sig_strength[] = 0,
                          const char* output_root_file = "outputfiles/toy-study.root"
                          ) {


       tran = new TRandom3(12345) ;

       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
       ws->Print() ;


       printf("\n\n") ;
       for ( int nbi=0; nbi<3; nbi++ ) {
          for ( int mbi=0; mbi<4; mbi++ ) {

             char pname[100] ;
             RooConstVar* rcv(0x0) ;

             sprintf( pname, "N_%db_msig_met%d_realVal", nbi+2, mbi+1 ) ;
             rcv = (RooConstVar*) ws->obj( pname ) ;
             if ( rcv == 0x0 ) { printf("\n\n *** Can't find %s in workspace of file %s\n\n", pname, wsfile ) ; return ; }
             inputval_n_msig[nbi][mbi] = rcv -> getVal() ;

             sprintf( pname, "N_%db_msb_met%d_realVal", nbi+2, mbi+1 ) ;
             rcv = (RooConstVar*) ws->obj( pname ) ;
             if ( rcv == 0x0 ) { printf("\n\n *** Can't find %s in workspace of file %s\n\n", pname, wsfile ) ; return ; }
             inputval_n_msb[nbi][mbi] = rcv -> getVal() ;

             printf(  "Input observable vals, %db, METsig%d :  SIG = %6.2f,  SB = %6.2f\n", nbi+2, mbi+1, inputval_n_msig[nbi][mbi], inputval_n_msb[nbi][mbi] ) ;

          } // mbi.
       } // nbi.
       printf("\n\n") ;





       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

       printf("\n\n\n  ===== SbModel ====================\n\n") ;
       modelConfig->Print() ;







       RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;
       printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

       rds->Print() ;
       rds->printMultiline(cout, 1, kTRUE, "") ;





       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       /// const RooArgSet* observables = modelConfig->GetObservables() ;

       /// const RooArgSet* nuisanceParameters = modelConfig->GetNuisanceParameters() ;

       RooRealVar* rrv_sig_strength = ws->var("sig_strength") ;
       if ( rrv_sig_strength == 0x0 ) {
          printf("\n\n\n *** can't find sig_strength in workspace.  Quitting.\n\n\n") ;
          return ;
       } else {
          printf(" current value is : %8.3f\n", rrv_sig_strength->getVal() ) ; fflush(stdout) ;
          //-------
          // printf(" fixing to zero.\n") ;
          // rrv_sig_strength->setVal(0.) ;
          // rrv_sig_strength->setConstant(kTRUE) ;
          //-------
             printf(" Letting it float.\n") ;
             rrv_sig_strength->setVal(0.) ;
             rrv_sig_strength->setConstant(kFALSE) ;
          //-------
       }

       printf("\n\n\n  ===== Doing a fit ====================\n\n") ;

       RooFitResult* preFitResult = likelihood->fitTo( *rds, Save(true) ) ;
       const RooArgList preFitFloatVals = preFitResult->floatParsFinal() ;
       {
         TIterator* parIter = preFitFloatVals.createIterator() ;
         while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
            printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
         }
       }



      //-- collect the nominal signal contributions to the observables at the theory Xsec.

       printf("\n\n") ;
       rrv_sig_strength->setVal(1.0) ;
       for ( int nbi=0; nbi<3; nbi++ ) {
          for ( int mbi=0; mbi<4; mbi++ ) {

             char pname[100] ;
             RooRealVar* rrv(0x0) ;

             sprintf( pname, "mu_sig_%db_msig_met%d", nbi+2, mbi+1 ) ;
             rrv = (RooRealVar*) ws -> obj( pname ) ;
             if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in workspace of file %s.\n\n", pname, wsfile ) ; return ; }
             susyval_n_msig[nbi][mbi] = rrv -> getVal() ;

             sprintf( pname, "mu_sig_%db_msb_met%d", nbi+2, mbi+1 ) ;
             rrv = (RooRealVar*) ws -> obj( pname ) ;
             if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in workspace of file %s.\n\n", pname, wsfile ) ; return ; }
             susyval_n_msb[nbi][mbi] = rrv -> getVal() ;

             printf(  "SUSY vals at theory Xsec, %db, METsig%d :  SIG = %6.2f,  SB = %6.2f\n", nbi+2, mbi+1, susyval_n_msig[nbi][mbi], susyval_n_msb[nbi][mbi] ) ;

          } // mbi.
       } // nbi.
       printf("\n\n") ;

       rrv_sig_strength->setVal(0.0) ;
       rrv_sig_strength->setConstant(kFALSE) ;










       printf("\n\n\n  ===== Begin Toy study ====================\n\n") ;

 ////  RooMCStudy mcs( *likelihood,
 ////                  *observables,
 ////                  Constrain( *nuisanceParameters ),
 ////                  Silence(),
 ////                  FitOptions(PrintLevel(-1),PrintEvalErrors(-1))
 ////                  ) ;
 ////
 //// mcs.generate( nToys, 1, true, "" ) ;

      int    fit_cov_qual ;
      int    true_ssi ;
      double true_sig_strength_val ;
      double fit_sig_strength_val ;
      double fit_sig_strength_err ;
      double fit_sig_strength_err_for_pull ;



      TTree* tt(0x0) ;

      tt = new TTree("toytt", "Toy TTree" ) ;

      tt->Branch( "fit_cov_qual",  &fit_cov_qual,  "fit_cov_qual/I"  ) ;
      tt->Branch( "true_ssi",  &true_ssi,  "true_ssi/I"  ) ;
      tt->Branch( "true_sig_strength_val",  &true_sig_strength_val,  "true_sig_strength_val/D"  ) ;
      tt->Branch( "fit_sig_strength_val",  &fit_sig_strength_val,  "fit_sig_strength_val/D"  ) ;
      tt->Branch( "fit_sig_strength_err",  &fit_sig_strength_err,  "fit_sig_strength_err/D"  ) ;
      tt->Branch( "fit_sig_strength_err_for_pull",  &fit_sig_strength_err_for_pull,  "fit_sig_strength_err_for_pull/D"  ) ;


      for ( int ssvi=0; ssvi<n_ssvals; ssvi++ ) {

         for ( int ti=0; ti<nToys_per_ssval; ti++ ) {

            true_sig_strength_val = gen_sig_strength[ssvi] ;
            true_ssi = ssvi ;

            //-- initialize all floating parameters to the values from the pre fit.
            TIterator* parIter = preFitFloatVals.createIterator() ;
            while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
               RooRealVar* rrv = ws->var( par->GetName() ) ;
               rrv->setVal( par->getVal() ) ;
               //printf(" %20s : %8.2f\n", rrv->GetName(), rrv->getVal() ) ;
            }
            //printf("\n\n") ;

            //-------
            //  RooAbsData* toyrds = (RooAbsData*) mcs.genData(ti) ;
            //-------
            RooAbsData* toyrds = genToyData( ws, gen_sig_strength[ssvi] ) ;
            //-------

        //  const RooArgSet* dsras = toyrds->get() ;
        //  TIterator* obsIter = dsras->createIterator() ;
        //  while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
        //     if ( strcmp( obs->GetName(), "Nsig" ) == 0 ) { 
        //        gen_nsig = obs->getVal() ;
        //        cout << "   found Nsig : " << gen_nsig << endl << flush ;
        //        break ;
        //     }
        //  }

            RooArgSet ras ;
            ras.add( *rrv_sig_strength ) ;

            RooFitResult* rfr = likelihood->fitTo( *toyrds, Save(true), Hesse(true), Minos(ras), Strategy(1), PrintLevel(-1), PrintEvalErrors(-1) ) ;


            fit_cov_qual = rfr->covQual() ;
            fit_sig_strength_val = rrv_sig_strength -> getVal() ;
            fit_sig_strength_err = rrv_sig_strength -> getError() ;
            if ( fit_sig_strength_val > gen_sig_strength[ssvi] ) {
               fit_sig_strength_err_for_pull = rrv_sig_strength -> getErrorHi() ;
            } else {
               fit_sig_strength_err_for_pull = rrv_sig_strength -> getErrorLo() ;
            }


            tt->Fill() ;

            printf("\n\n\n  ===== RooDataSet for toy %d ====================\n\n", ti) ;

            toyrds->Print() ;
            toyrds->printMultiline(cout, 1, kTRUE, "") ;
            printf("\n\n\n Sig str %.1f, Toy %4d Fit covariance quality : %d\n", gen_sig_strength[ssvi], ti, fit_cov_qual ) ;
            printf("  %d :  Sig str %.1f, Toy %4d Fit signal strength: %.2f\n", ssvi, gen_sig_strength[ssvi], ti, fit_sig_strength_val ) ;

            printf("\n\n") ;

            delete rfr ;
            delete toyrds ;


        } // ti.

     } // ssvi.

     gSystem->Exec("mkdir -p outputfiles") ;
     TFile f( output_root_file, "recreate") ;
     tt->Write() ;
     f.Close() ;



   } // ws_hbb_toystudy1

   //==================================================================================

   RooDataSet* genToyData( RooWorkspace* ws, float gen_sig_strength ) {

      RooArgSet obsList ;

      for ( int nbi=0; nbi<3; nbi++ ) {
         for ( int mbi=0; mbi<4; mbi++ ) {

            char pname[1000] ;
            RooRealVar* rrv(0x0) ;
            int nobs(0) ;

            sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
            rrv = (RooRealVar*) ws -> obj( pname ) ;
            if ( rrv == 0x0 ) { printf("\n\n *** genToyData : can't find %s in workspace.\n\n", pname ) ; return 0x0 ; }
            nobs = tran -> Poisson( inputval_n_msig[nbi][mbi] + gen_sig_strength * susyval_n_msig[nbi][mbi] ) ;
            rrv -> setVal( nobs ) ;
            obsList.add( *rrv ) ;

            sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
            rrv = (RooRealVar*) ws -> obj( pname ) ;
            if ( rrv == 0x0 ) { printf("\n\n *** genToyData : can't find %s in workspace.\n\n", pname ) ; return 0x0 ; }
            nobs = tran -> Poisson( inputval_n_msb[nbi][mbi] + gen_sig_strength * susyval_n_msb[nbi][mbi] ) ;
            rrv -> setVal( nobs ) ;
            obsList.add( *rrv ) ;

         } // mbi.
      } // nbi.

      RooDataSet* rds = new RooDataSet("toyfit_hbb_observed_rds", "hbb toy observed data values", obsList ) ;
      rds->add( obsList ) ;

      return rds ;

   } // genToyData.

   //==================================================================================

