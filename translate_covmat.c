
#include "TMath.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"


#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"

#include "histio.c"

#include "scan_sigstrength.c"
#include "cormat_plot.c"


  using namespace RooFit;
  using namespace RooStats;

  int first_met_bin_array_index(0) ;

  double compute_dBGdP_nbSIG( const char* bg_par_name, const char* fit_par_name, RooWorkspace* ws ) ;

  double compute_dBGdP_nbSB( const char* bg_par_name, const char* fit_par_name ) ;

  double compute_bg_cov( const char* bg_par_name_1, const char* bg_par_name_2, RooFitResult* rfr, RooWorkspace* ws ) ;

  //-------

   void translate_covmat( const char* wsfile = "outputfiles/ws-data-unblind-sigmass-500.root" ) {

      TText* tt_title = new TText() ;
      tt_title -> SetTextAlign(33) ;

      gStyle -> SetOptStat(0) ;

      gDirectory->Delete("h*") ;

      TFile* wstf = new TFile( wsfile ) ;

      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      ws->Print() ;

      int bins_of_met = TMath::Nint( ws->var("bins_of_met")->getVal()  ) ;
      printf("\n\n Bins of MET : %d\n\n", bins_of_met ) ;

      first_met_bin_array_index = TMath::Nint( ws->var("first_met_bin_array_index")->getVal()  ) ;
      if ( first_met_bin_array_index == 0 ) {
         printf(" Using first MET bin.\n" ) ;
      } else {
         printf(" *** NOT Using first MET bin.\n" ) ;
      }

      int bins_of_nb = TMath::Nint( ws->var("bins_of_nb")->getVal()  ) ;
      printf("\n\n Bins of nb : %d\n\n", bins_of_nb ) ;

      int nb_lookup[10] ;
      if ( bins_of_nb == 2 ) {
         nb_lookup[0] = 2 ;
         nb_lookup[1] = 4 ;
      } else if ( bins_of_nb == 3 ) {
         nb_lookup[0] = 2 ;
         nb_lookup[1] = 3 ;
         nb_lookup[2] = 4 ;
      }


      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }

      ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

      RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;

      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooAbsPdf* likelihood = modelConfig->GetPdf() ;

      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(3) ) ;
      fitResult->Print() ;


      char bg_par_name_list[100][100] ;
      int nbgpar(0) ;


      for ( int mbi=3; mbi>=0; mbi-- ) {
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_4b_msig_met%d", mbi+1 ) ;
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_4b_msb_met%d", mbi+1 ) ;
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_3b_msig_met%d", mbi+1 ) ;
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_3b_msb_met%d", mbi+1 ) ;
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_2b_msig_met%d", mbi+1 ) ;
         sprintf( bg_par_name_list[nbgpar++], "mu_bg_2b_msb_met%d", mbi+1 ) ;
      } // mbi.

      double bg_par_val[100] ;
      for ( int bgpi=0; bgpi<nbgpar; bgpi++ ) {
         double val(0.) ;
         RooRealVar* rrv = (RooRealVar*) ws -> obj( bg_par_name_list[bgpi] ) ;
         if ( rrv != 0x0 ) val = rrv -> getVal() ;
         bg_par_val[bgpi] = val ;
      } // bgpi.


      double bg_cov_mat[100][100] ;

      printf( "\n\n Computing BG par covariance matrix.\n\n" ) ;

      for ( int rbgpi=0; rbgpi<nbgpar; rbgpi++ ) {
         printf("  %s : ", bg_par_name_list[rbgpi] ) ;
         for ( int cbgpi=0; cbgpi<nbgpar; cbgpi++ ) {
            bg_cov_mat[rbgpi][cbgpi] = compute_bg_cov( bg_par_name_list[rbgpi], bg_par_name_list[cbgpi], fitResult, ws ) ;
            printf(".") ; fflush(stdout) ;
         } // rbgpi.
         printf("\n") ;
      } // cbgpi.

      printf("\n\n === Correlation matrix for BG pars.\n\n" ) ;

      TH2F* h_cormat = new TH2F( "h_cormat", "correlation matrix", nbgpar, 0.5, nbgpar+0.5, nbgpar, 0.5, nbgpar+0.5 ) ;
      h_cormat -> SetMinimum(-1.) ;
      h_cormat -> SetMaximum(+1.) ;
      h_cormat -> SetLabelSize( 0.03, "x" ) ;
      h_cormat -> SetLabelSize( 0.03, "y" ) ;

      for ( int rbgpi=0; rbgpi<nbgpar; rbgpi++ ) {
         printf( " %35s : ", bg_par_name_list[rbgpi] ) ;
         h_cormat -> GetXaxis() -> SetBinLabel( rbgpi+1, bg_par_name_list[rbgpi] ) ;
         for ( int cbgpi=0; cbgpi<nbgpar; cbgpi++ ) {
            double rsig = sqrt( bg_cov_mat[rbgpi][rbgpi] ) ;
            double csig = sqrt( bg_cov_mat[cbgpi][cbgpi] ) ;
            double rho = bg_cov_mat[rbgpi][cbgpi] / ( rsig * csig ) ;
            printf( " %7.4f ", rho ) ;
            h_cormat -> SetBinContent( rbgpi+1, nbgpar-cbgpi, rho ) ;
            h_cormat -> GetYaxis() -> SetBinLabel( nbgpar-cbgpi, bg_par_name_list[cbgpi] ) ;
         } // cbgpi.
         printf("\n") ;
      } // rbgpi.

      printf("\n\n") ;




      printf("\n\n ======== LandS correlated uncertainty parameters.\n\n" ) ;

      for ( int rbgpi=0; rbgpi<nbgpar; rbgpi++ ) {
         printf( " %35s : ", bg_par_name_list[rbgpi] ) ;
         for ( int cbgpi=0; cbgpi<nbgpar; cbgpi++ ) {
            double val(0.) ;
            if ( rbgpi == cbgpi ) {
               double sum(0.) ;
               for ( int sbgpi=0; sbgpi<nbgpar; sbgpi++ ) {
                  if ( sbgpi == rbgpi ) continue ;
                  if ( bg_par_val[rbgpi] != 0 && bg_par_val[sbgpi] != 0 ) {
                     sum += fabs( bg_cov_mat[rbgpi][sbgpi] ) / ( bg_par_val[rbgpi] * bg_par_val[sbgpi] ) ;
                  }
               } // bgpi.
               double diff(0.) ;
               if ( bg_par_val[rbgpi] != 0 ) {
                  diff = bg_cov_mat[rbgpi][cbgpi] / ( bg_par_val[rbgpi] * bg_par_val[cbgpi]) - sum ;
               }
               if ( diff > 0 ) val = sqrt( diff ) ;
            } else {
               if ( bg_par_val[rbgpi] != 0 && bg_par_val[cbgpi] != 0 ) {
                  if ( bg_cov_mat[rbgpi][cbgpi] > 0 ) {
                     val = sqrt( fabs(bg_cov_mat[rbgpi][cbgpi]) / ( bg_par_val[rbgpi] * bg_par_val[cbgpi]) ) ;
                  } else {
                     val = -1. * sqrt( fabs(bg_cov_mat[rbgpi][cbgpi]) / ( bg_par_val[rbgpi] * bg_par_val[cbgpi]) ) ;
                  }
               }
            }
            val += 1. ;
            printf( " %6.3f ", val ) ;
         } // cbgpi.
         printf("\n") ;
      } // rbgpi.

      printf("\n\n") ;


      printf(" ==== Values and derived errors from diagonal elements.\n\n" ) ;
      for ( int bgpi=0; bgpi<nbgpar; bgpi++ ) {
         printf( "  %35s :  %8.3f +/- %6.3f\n", bg_par_name_list[bgpi], bg_par_val[bgpi], sqrt( bg_cov_mat[bgpi][bgpi] ) ) ;
      } // bgpi.

      printf("\n\n") ;



      gStyle -> SetPadRightMargin( 0.15 ) ;
      gStyle -> SetPadTopMargin( 0.15 ) ;
      gStyle -> SetPadLeftMargin( 0.25 ) ;
      gStyle -> SetPadBottomMargin( 0.25 ) ;

      TCanvas* can_cormat = new TCanvas( "can_cormat", "correlation matrix", 800, 800 ) ;

      h_cormat -> GetXaxis() -> LabelsOption("v") ;
      h_cormat -> Draw("colz") ;
      can_cormat -> SaveAs("outputfiles/bg-cormat.pdf") ;

      TFile* tfcm = new TFile("outputfiles/bg-cormat.root","recreate") ;
      h_cormat -> Write() ;
      tfcm -> Close() ;

   } // translate_covmat

 //===========================================================================================

  double compute_dBGdP_nbSIG( const char* bg_par_name, const char* fit_par_name, RooWorkspace* ws ) {

  // bool verb(true) ;
     bool verb(false) ;

     int bg_nb(-1) ;
     int bg_mb(-1) ;

     sscanf( bg_par_name, "mu_bg_%db_msig_met%d", &bg_nb, &bg_mb ) ;

     if (verb) printf( " compute_dBGdP_nbSIG : computing derivative of %db, METsig %d with respect to %s \n", bg_nb, bg_mb, fit_par_name ) ;
     if ( bg_nb == -1 || bg_mb == -1 ) {
        printf("\n\n *** compute_dBGdP_nbSIG : unrecognized background string %s\n\n", bg_par_name ) ;
        return 0. ;
     }

     RooRealVar* bg_par_rrv = (RooRealVar*) ws -> obj( bg_par_name ) ;
     if ( bg_par_rrv == 0x0 ) {
        printf( "\n\n *** compute_dBGdP_nbSIG : can't find formula for %s in ws.\n\n", bg_par_name ) ;
        return 0. ;
     }
     double bg_par_val = bg_par_rrv->getVal() ;
     if ( verb ) { printf( " value of %s is %.4f\n", bg_par_name, bg_par_val ) ; }

     char pname[100] ;
     bool fp_found(false) ;
     char ws_name[100] ;
     bool is_prim_par(false) ;

    //-- R SIG/SB
     sprintf( pname, "R_msigmsb_met%d", bg_mb ) ;
     if ( strcmp( fit_par_name, pname ) == 0 ) {
        if (verb) printf( " Found R_sigsb.\n" ) ;
        sprintf( ws_name, "%s", pname ) ;
        fp_found = true ;
     }

    //-- SB yield
     sprintf( pname, "mu_bg_%db_msb_met%d", bg_nb, bg_mb ) ;
     if ( strcmp( fit_par_name, pname ) == 0 ) {
        if (verb) printf( " Found SB yield.\n" ) ;
        sprintf( ws_name, "%s", pname ) ;
        fp_found = true ;
     }

    //-- Closure syst
     sprintf( pname, "prim_Rsigsb_corr_%db_met%d", bg_nb, bg_mb ) ;
     if ( strcmp( fit_par_name, pname ) == 0 ) {
        if (verb) printf( " Found ABCD closure syst.\n" ) ;
        sprintf( ws_name, "Rsigsb_corr_%db_met%d", bg_nb, bg_mb ) ;
        fp_found = true ;
        is_prim_par = true ;
     }

    //-- BG sample comp syst.
     sprintf( pname, "prim_background_sample_comp" ) ;
     if ( strcmp( fit_par_name, pname ) == 0 ) {
        if (verb) printf( " Found BG sample comp syst.\n" ) ;
        sprintf( ws_name, "background_sample_comp" ) ;
        fp_found = true ;
        is_prim_par = true ;
     }

     if ( !fp_found ) {
        if (verb) printf(" compute_dBGdP_nbSIG : %s does not depend on %s.  Returning zero.\n", bg_par_name, fit_par_name ) ;
        return 0. ;
     }



     RooRealVar* rrv = (RooRealVar*) ws->obj( ws_name ) ;
     if ( rrv == 0x0 ) {
        printf(" \n\n *** compute_dBGdP_nbSIG : can't find %s in ws.\n\n", pname ) ;
        return 0. ;
     }

     if ( is_prim_par ) {

        RooRealVar* prim_rrv = (RooRealVar*) ws->obj( fit_par_name ) ;
        if ( prim_rrv == 0x0 ) { printf(" \n\n *** %s missing from ws.\n\n", fit_par_name ) ; return 0. ; }
        double original_prim_val = prim_rrv -> getVal() ;
        double original_lognorm_val = rrv -> getVal() ;
        prim_rrv -> setVal( original_prim_val + 1. ) ;
        double var_lognorm_val = rrv -> getVal() ;
        prim_rrv -> setVal( original_prim_val ) ;
        double dlognorm_dfitpar = var_lognorm_val - original_lognorm_val ;
        if ( verb ) printf(" %s lognorm, nom = %.4f, var = %.4f\n", fit_par_name, original_lognorm_val, var_lognorm_val ) ;

        double fit_par_val = rrv->getVal() ;
        if ( verb ) printf("  Value of %s is %.4f\n", ws_name, fit_par_val ) ;
        if ( fit_par_val > 0. && bg_par_val > 0. ) {
           return dlognorm_dfitpar * (bg_par_val / fit_par_val) ;
        } else {
           return 0. ;
        }

     } else {

        double fit_par_val = rrv->getVal() ;
        if ( verb ) printf("  Value of %s is %.4f\n", ws_name, fit_par_val ) ;
        if ( fit_par_val > 0. && bg_par_val > 0. ) {
           return bg_par_val / fit_par_val ;
        } else {
           return 0. ;
        }

     }


  } // compute_dBGdP_nbSIG

 //===========================================================================================

  double compute_dBGdP_nbSB( const char* bg_par_name, const char* fit_par_name ) {

     if ( strcmp( bg_par_name, fit_par_name ) == 0 ) {
        return 1. ;
     } else {
        return 0. ;
     }

  } // compute_dBGdP_nbSB

 //===========================================================================================

  double compute_bg_cov( const char* bg_par_name_1, const char* bg_par_name_2, RooFitResult* rfr, RooWorkspace* ws ) {

      bool verb(false) ;

      TString ts_bg_par_name_1( bg_par_name_1 ) ;
      TString ts_bg_par_name_2( bg_par_name_2 ) ;
      bool bgp1_is_sig(false) ;
      bool bgp2_is_sig(false) ;
      if ( ts_bg_par_name_1.Contains("_msig_") ) { bgp1_is_sig = true ; }
      if ( ts_bg_par_name_2.Contains("_msig_") ) { bgp2_is_sig = true ; }


      const RooArgList fitpars = rfr -> floatParsFinal() ;
      int nfp = fitpars.getSize() ;

      const TMatrixDSym covMat = rfr -> covarianceMatrix() ;

      double bg_cov(0.) ;

      for ( int fpri=0; fpri<nfp; fpri++ ) {
         for ( int fpci=0; fpci<nfp; fpci++ ) {

            RooRealVar* rpar = (RooRealVar*) fitpars.at(fpri) ;
            double rp_dBGdP(0.) ;
            if ( bgp1_is_sig ) {
               rp_dBGdP = compute_dBGdP_nbSIG( bg_par_name_1, rpar->GetName(), ws ) ;
            } else {
               rp_dBGdP = compute_dBGdP_nbSB( bg_par_name_1, rpar->GetName() ) ;
            }

            RooRealVar* cpar = (RooRealVar*) fitpars.at(fpci) ;
            double cp_dBGdP(0.) ;
            if ( bgp2_is_sig ) {
               cp_dBGdP = compute_dBGdP_nbSIG( bg_par_name_2, cpar->GetName(), ws ) ;
            } else {
               cp_dBGdP = compute_dBGdP_nbSB( bg_par_name_2, cpar->GetName() ) ;
            }

            if ( cp_dBGdP > 0 && rp_dBGdP > 0 ) {
               if (verb) {
                  printf(" %d,%d : value of derivative d %s / d %s is %.4f and d %s / d %s is %.4f\n",
                     fpri, fpci, bg_par_name_1, rpar->GetName(), rp_dBGdP, bg_par_name_2, cpar->GetName(), cp_dBGdP ) ;
                  printf(" %d,%d : fit par cov is %.6f.  Contribution to sum is %.6f\n",
                     fpri, fpci, covMat( fpri, fpci ), rp_dBGdP * cp_dBGdP * covMat( fpri, fpci ) ) ;
                  printf("\n") ;
               }
            }

            bg_cov += rp_dBGdP * cp_dBGdP * covMat( fpri, fpci ) ;

         } // fpci.
      } // fpri.


      return bg_cov ;


  } // compute_bg_cov

 //===========================================================================================








