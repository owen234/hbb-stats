
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"

  using namespace RooFit;
  using namespace RooStats;


   const int bins_of_nb(3) ;
   const int bins_of_met(4) ;
   const int first_met_bin_array_index(0) ;

   int nb_lookup[3] = { 2, 3, 4 } ;

   void gen_lands_from_ws( const char* ws_root_file, bool bg_at_ml = false ) {


      char pname[1000] ;


      TFile infile( ws_root_file, "read" ) ;
      if ( ! infile.IsOpen() ) { printf("\n\n *** Can't open input workspace root file : %s\n\n", ws_root_file ) ; return ; }

      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( infile.Get("ws") );

      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }

      ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

      RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;

      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooAbsPdf* likelihood = modelConfig->GetPdf() ;





     //-- unpack observables.

      int obs_N_msig[10][50] ; // first index is n btags, second is met bin.
      int obs_N_msb[10][50]  ; // first index is n btags, second is met bin.

      printf("\n\n") ;

      const RooArgSet* dsras = rds->get() ;
      TIterator* obsIter = dsras->createIterator() ;
      while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
         for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
            for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
               sprintf( pname, "N_%db_msig_met%d", nb_lookup[nbi], mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msig[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
               sprintf( pname, "N_%db_msb_met%d", nb_lookup[nbi], mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msb[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
            } // mbi.
         } // nbi.
      } // obs iterator.

      printf("\n\n") ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
            printf("    Observed events :  SIG = %4d ,  BG = %4d\n",  obs_N_msig[nbi][mbi], obs_N_msb[nbi][mbi] ) ;
         } // mbi.
      } // nbi.





     //-- signal counts at signal strength of 1.

      printf("\n\n") ;

      rv_sig_strength -> setVal( 1. ) ;
      float susy_msig[10][50] ;  // first index is n btags, second is met bin.
      float susy_msb[10][50]  ;  // first index is n btags, second is met bin.
      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "mu_sig_%db_msig_met%d", nbi+2, mbi+1 ) ;
            RooRealVar* rrv_susy_msig_obs = (RooRealVar*) ws->obj( pname ) ;
            if ( rrv_susy_msig_obs == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
            susy_msig[nbi][mbi] = rrv_susy_msig_obs -> getVal() ;

            sprintf( pname, "mu_sig_%db_msb_met%d", nbi+2, mbi+1 ) ;
            RooRealVar* rrv_susy_msb_obs = (RooRealVar*) ws->obj( pname ) ;
            if ( rrv_susy_msb_obs == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
            susy_msb[nbi][mbi] = rrv_susy_msb_obs -> getVal() ;

            printf( "   SUSY prediction for %db, METsig %d :   SIG = %5.2f,  SB = %5.2f\n", nbi+2, mbi+1, susy_msig[nbi][mbi], susy_msb[nbi][mbi] ) ;

         } // mbi.
      } // nbi.






     //-- Get BG values from fit with either signal strength set to zero or at ML value.

      rv_sig_strength -> setVal( 0. ) ;
      if ( bg_at_ml ) {
         printf("\n\n  ===== Fit with signal strength fixed to zero.\n\n") ;
         rv_sig_strength -> setConstant( kFALSE ) ;
      } else {
         printf("\n\n  ===== Fit with signal strength floating.\n\n") ;
         rv_sig_strength -> setConstant( kTRUE ) ;
      }

      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(3) ) ;
      fitResult->Print() ;

      printf("\n\n") ;



      float smbg_msig[10][50] ;  // first index is n btags, second is met bin.
      float smbg_msb[10][50]  ;  // first index is n btags, second is met bin.

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    //  This is the old way.  BG values are biased by observed counts in SIG bins.
    //
      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "mu_bg_%db_msig_met%d", nbi+2, mbi+1 ) ;
            RooRealVar* rrv_smbg_msig_obs = (RooRealVar*) ws->obj( pname ) ;
            if ( rrv_smbg_msig_obs == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
            smbg_msig[nbi][mbi] = rrv_smbg_msig_obs -> getVal() ;

            sprintf( pname, "mu_bg_%db_msb_met%d", nbi+2, mbi+1 ) ;
            RooRealVar* rrv_smbg_msb_obs = (RooRealVar*) ws->obj( pname ) ;
            if ( rrv_smbg_msb_obs == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
            smbg_msb[nbi][mbi] = rrv_smbg_msb_obs -> getVal() ;

            printf( "   BG value for %db, METsig %d :   SIG = %7.2f,  SB = %7.2f\n", nbi+2, mbi+1, smbg_msig[nbi][mbi], smbg_msb[nbi][mbi] ) ;

         } // mbi.
      } // nbi.
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    //  Instead, set BG to simple ABCD calc value for 3b and 3b SIG bins and the
    //  observed BG value for the rest.
    //
    //  Jan 27: this method gives worse agreement than method above.  I hate LandS.
    //
    //for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
    //   for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
    //      smbg_msb[nbi][mbi] = obs_N_msb[nbi][mbi] ;
    //      if ( nbi==0 ) { // 2b
    //         smbg_msig[nbi][mbi] = obs_N_msig[nbi][mbi] ;
    //      } else { // 3b or 4b
    //         float sig_over_sb_ratio_2b = 0. ;
    //         if ( obs_N_msb[0][mbi] > 0. ) { sig_over_sb_ratio_2b = (1.0 * obs_N_msig[0][mbi]) / (1.0 * obs_N_msb[0][mbi]) ; }
    //         smbg_msig[nbi][mbi] = sig_over_sb_ratio_2b * obs_N_msb[nbi][mbi] ;
    //      }
    //      printf( "   BG value for %db, METsig %d :   SIG = %7.2f,  SB = %7.2f\n", nbi+2, mbi+1, smbg_msig[nbi][mbi], smbg_msb[nbi][mbi] ) ;
    //   } // mbi.
    //} // nbi.

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      printf("\n\n") ;





     //-- ABCD closure

      float abcd_closure[10][50] ;  // first index is n btags, second is met bin.

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=first_met_bin_array_index; mbi<bins_of_met; mbi++ ) {
            if ( nbi==0 ) {
               abcd_closure[nbi][mbi] = 0. ;
            } else {
               sprintf( pname, "sigma_Rsigsb_corr_%db_met%d", nbi+2, mbi+1 ) ;
               RooRealVar* rrv = (RooRealVar*) ws->obj( pname ) ;
               if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
               abcd_closure[nbi][mbi] = rrv -> getVal() ;
               printf( "    ABCD closure syst, %db,  METsig %d :   %.3f\n", nbi+2, mbi+1, abcd_closure[nbi][mbi] ) ;
            }
         } // mbi.
      } // nbi.

      printf("\n\n") ;







     //-- BG sample comp
      float bg_sample_comp(0.) ;
      {
          sprintf( pname, "sigma_background_sample_comp" ) ;
          RooRealVar* rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          bg_sample_comp = rrv -> getVal() ;
          printf( "    background sample comp :   %.3f\n", bg_sample_comp ) ;
      }

      printf("\n\n") ;





     //-- Lumi
      float lumi_syst(0.) ;
      {
          sprintf( pname, "sigma_luminosity_uncertainty" ) ;
          RooRealVar* rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          lumi_syst = rrv -> getVal() ;
          printf( "    luminosity uncertainty :   %.3f\n", lumi_syst ) ;
      }

      printf("\n\n") ;





     //-- Trigger efficiency

      float trig_syst_metsig1(0.) ;
      float trig_syst_metsig2(0.) ;
      float trig_syst_metsig34(0.) ;
      {
          RooRealVar* rrv(0x0) ;
          float trig_eff_val(1.), trig_eff_err(0.) ;

          sprintf( pname, "sigma_trig_eff_corr_metsig1" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_err = rrv -> getVal() ;

          sprintf( pname, "mean_trig_eff_corr_metsig1" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_val = rrv -> getVal() ;

          if ( trig_eff_val > 0. ) { trig_syst_metsig1 = trig_eff_err / trig_eff_val ; }
          printf( "    trigger efficiency, METsig1  :   %.3f +/- %.3f   (%.3f)\n", trig_eff_val, trig_eff_err, trig_syst_metsig1 ) ;


          sprintf( pname, "sigma_trig_eff_corr_metsig2" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_err = rrv -> getVal() ;

          sprintf( pname, "mean_trig_eff_corr_metsig2" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_val = rrv -> getVal() ;

          if ( trig_eff_val > 0. ) { trig_syst_metsig2 = trig_eff_err / trig_eff_val ; }
          printf( "    trigger efficiency, METsig2  :   %.3f +/- %.3f   (%.3f)\n", trig_eff_val, trig_eff_err, trig_syst_metsig2 ) ;


          sprintf( pname, "sigma_trig_eff_corr_metsig34" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_err = rrv -> getVal() ;

          sprintf( pname, "mean_trig_eff_corr_metsig34" ) ;
          rrv = (RooRealVar*) ws->obj( pname ) ;
          if ( rrv == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
          trig_eff_val = rrv -> getVal() ;

          if ( trig_eff_val > 0. ) { trig_syst_metsig34 = trig_eff_err / trig_eff_val ; }
          printf( "    trigger efficiency, METsig34 :   %.3f +/- %.3f   (%.3f)\n", trig_eff_val, trig_eff_err, trig_syst_metsig34 ) ;

      }

      printf("\n\n") ;


     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     //-- Generate LandS output file.

      TString infile_ts( ws_root_file ) ;

      TString outfile_name = infile_ts.ReplaceAll(".root","-lands.txt") ;
      printf(" output file : %s\n", outfile_name.Data() ) ;

      FILE* outfile_lands = fopen( outfile_name.Data(), "w" ) ;

      printf("\n\n #----------------------------------------------------------------------------\n\n") ;
      printf("\n\n Creating LandS datacard: %s\n\n\n", outfile_name.Data() ) ;

      float par_width_NSD(6.) ;

      int nchan(0) ;
      int npars(0) ;

      fprintf( outfile_lands, "---\n" ) ;
      fprintf( outfile_lands, "## datacard for EW SUSY H(bb)H(bb)+MET ABCD using 4b, 3b, and 2b in %d bins of METsig\n\n\n", bins_of_met ) ;
      fprintf( outfile_lands, "## generated by gen_lands_from_ws.c from %s\n\n", ws_root_file ) ;
      nchan = 6 * bins_of_met ;
      npars = 59 ;
      fprintf( outfile_lands, "imax %d  number of channels (%d bins of met for 4b, 3b, and 2b, higgs mass SIG and SB)\n", nchan, bins_of_met ) ;
      fprintf( outfile_lands, "jmax 1  number of backgrounds\n" ) ;
      fprintf( outfile_lands, "kmax %d  number of nuisance parameters\n", npars ) ;
      fprintf( outfile_lands, "---\n\n\n" ) ;


      fprintf( outfile_lands, "-------------------" ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         fprintf( outfile_lands, "------------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
      }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "bin                " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "               N4bsig_met%d                 N4bsb_met%d                  N3bsig_met%d                 N3bsb_met%d                  N2bsig_met%d                 N2bsb_met%d   ", mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1 ) ; }
      fprintf( outfile_lands, "\n" ) ;
      fprintf( outfile_lands, "Observation        " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "                    %6d                     %6d                       %6d                     %6d                       %6d                     %6d   ",
          obs_N_msig[2][mbi] ,
          obs_N_msb[2][mbi]  ,
          obs_N_msig[1][mbi] ,
          obs_N_msb[1][mbi]  ,
          obs_N_msig[0][mbi] ,
          obs_N_msb[0][mbi]    ) ; }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "-------------------" ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         fprintf( outfile_lands, "------------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
      }
      fprintf( outfile_lands, "\n" ) ;



      fprintf( outfile_lands, "bin                " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "   N4bsig_met%d N4bsig_met%d     N4bsb_met%d  N4bsb_met%d      N3bsig_met%d N3bsig_met%d     N3bsb_met%d  N3bsb_met%d      N2bsig_met%d N2bsig_met%d     N2bsb_met%d  N2bsb_met%d   ", mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1, mbi+1 ) ; }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "process            " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "        signal        smbg         signal        smbg           signal        smbg         signal        smbg           signal        smbg         signal        smbg   " ) ; }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "process            " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "             0           1              0           1                0           1              0           1                0           1              0           1   " ) ; }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "rate               " ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile_lands, "        %6.2f     %7.1f         %6.2f     %7.1f           %6.2f     %7.1f         %6.2f     %7.1f           %6.2f     %7.1f         %6.2f     %7.1f   ",
           susy_msig[2][mbi],  smbg_msig[2][mbi],
           susy_msb[2][mbi] ,  smbg_msb[2][mbi] ,
           susy_msig[1][mbi],  smbg_msig[1][mbi],
           susy_msb[1][mbi] ,  smbg_msb[1][mbi] ,
           susy_msig[0][mbi],  smbg_msig[0][mbi],
           susy_msb[0][mbi] ,  smbg_msb[0][mbi]
           ) ; }
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "-------------------" ) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         fprintf( outfile_lands, "------------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
      }
      fprintf( outfile_lands, "\n" ) ;







     //--- flat correlation parameters to propagate the statistical errors on the ABCD method.

      for ( int rmbi=1; rmbi<=bins_of_met; rmbi++ ) {

         float lands_par_error ;
         float nobs ;

         fprintf( outfile_lands, "c_m%d_sig_allnb   lnU    ", rmbi ) ;
         //nobs = n_4b_msb[rmbi] ;
         nobs = obs_N_msb[2][rmbi-1] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -        %4.2f              -           -                -        %4.2f              -           -                -        %4.2f              -           -        ", lands_par_error, lands_par_error, lands_par_error ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;



         fprintf( outfile_lands, "c_m%d_R_4b        lnU    ", rmbi ) ;
         //nobs = n_4b_msb[rmbi] ;
         nobs = obs_N_msb[2][rmbi-1] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -        %4.2f              -        %4.2f                -           -              -           -                -           -              -           -        ", lands_par_error, lands_par_error ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;




         fprintf( outfile_lands, "c_m%d_R_3b        lnU    ", rmbi ) ;
         //nobs = n_3b_msb[rmbi] ;
         nobs = obs_N_msb[1][rmbi-1] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -           -              -           -                -        %4.2f              -        %4.2f                -           -              -           -        ", lands_par_error, lands_par_error ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;



         fprintf( outfile_lands, "c_m%d_R_2b        lnU    ", rmbi ) ;
         //nobs = n_2b_msb[rmbi] ;
         nobs = obs_N_msb[0][rmbi-1] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -        %4.2f              -        %4.2f        ", lands_par_error, lands_par_error ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;


      } // rmbi.







     //-- ABCD closure syst.

      for ( int rmbi=1; rmbi<=bins_of_met; rmbi++ ) {
         fprintf( outfile_lands, "closure_4b_met%d  lnN    ", rmbi ) ;
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -       %5.3f              -           -                -           -              -           -                -           -              -           -        ", (1. + abcd_closure[2][rmbi-1]) ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;
      } // rmbi.

      for ( int rmbi=1; rmbi<=bins_of_met; rmbi++ ) {
         fprintf( outfile_lands, "closure_3b_met%d  lnN    ", rmbi ) ;
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               fprintf( outfile_lands, "        -           -              -           -                -       %5.3f              -           -                -           -              -           -        ", (1. + abcd_closure[1][rmbi-1]) ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -        " ) ;
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;
      } // rmbi.






     //-- BG sample comp.

      fprintf( outfile_lands, "bg_sample_comp   lnN    " ) ;
      for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
          fprintf( outfile_lands, "        -       %5.3f              -           -                -       %5.3f              -           -                -           -              -           -        ", (1.+bg_sample_comp), (1.+bg_sample_comp) ) ;
      } // cmbi.
      fprintf( outfile_lands, "\n" ) ;








     //-- Lumi syst.

      fprintf( outfile_lands, "lumi_syst        lnN    " ) ;
      for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
         for ( int nbi=0; nbi<3; nbi++ ) {
            fprintf( outfile_lands, "    %5.3f           -          %5.3f           -        ", (1.+lumi_syst), (1.+lumi_syst) ) ;
         } // nbi.
      } // cmbi.
      fprintf( outfile_lands, "\n" ) ;






     //-- Trigger syst

      fprintf( outfile_lands, "trig_eff_met1    lnN    " ) ;
      for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
         for ( int nbi=0; nbi<3; nbi++ ) {
            if ( cmbi == 1 ) {
               fprintf( outfile_lands, "    %5.3f           -          %5.3f           -        ", (1.+trig_syst_metsig1), (1.+trig_syst_metsig1) ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -        " ) ;
            }
         } // nbi.
      } // cmbi.
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "trig_eff_met2    lnN    " ) ;
      for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
         for ( int nbi=0; nbi<3; nbi++ ) {
            if ( cmbi == 2 ) {
               fprintf( outfile_lands, "    %5.3f           -          %5.3f           -        ", (1.+trig_syst_metsig2), (1.+trig_syst_metsig2) ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -        " ) ;
            }
         } // nbi.
      } // cmbi.
      fprintf( outfile_lands, "\n" ) ;

      fprintf( outfile_lands, "trig_eff_met34   lnN    " ) ;
      for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
         for ( int nbi=0; nbi<3; nbi++ ) {
            if ( cmbi > 2 ) {
               fprintf( outfile_lands, "    %5.3f           -          %5.3f           -        ", (1.+trig_syst_metsig34), (1.+trig_syst_metsig34) ) ;
            } else {
               fprintf( outfile_lands, "        -           -              -           -        " ) ;
            }
         } // nbi.
      } // cmbi.
      fprintf( outfile_lands, "\n" ) ;




     //-- All shape systematics, one line each.

      int nsyst(6) ;
      char syst_name[6][10] = { "btagSF", "ISR", "JER", "JES", "PDF", "PU" } ;


      for ( int si=0; si<nsyst; si++ ) {

         printf("\n\n ==== %s\n\n", syst_name[si] ) ;

         fprintf( outfile_lands, " %14s  lnN  ", syst_name[si] ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            printf(" METsig bin %d\n", mbi+1 ) ;

            char hname[1000] ;

            sprintf( hname, "h_syst_%s_msig_met%d_nom", syst_name[si], mbi+1 ) ;
            TH1F* h_msig_nom = (TH1F*) infile.Get( hname ) ;
            if ( h_msig_nom == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

            sprintf( hname, "h_syst_%s_msig_met%d_p1s", syst_name[si], mbi+1 ) ;
            TH1F* h_msig_p1s = (TH1F*) infile.Get( hname ) ;
            if ( h_msig_p1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

            sprintf( hname, "h_syst_%s_msig_met%d_m1s", syst_name[si], mbi+1 ) ;
            TH1F* h_msig_m1s = (TH1F*) infile.Get( hname ) ;
            if ( h_msig_m1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }


            sprintf( hname, "h_syst_%s_msb_met%d_nom", syst_name[si], mbi+1 ) ;
            TH1F* h_msb_nom = (TH1F*) infile.Get( hname ) ;
            if ( h_msb_nom == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

            sprintf( hname, "h_syst_%s_msb_met%d_p1s", syst_name[si], mbi+1 ) ;
            TH1F* h_msb_p1s = (TH1F*) infile.Get( hname ) ;
            if ( h_msb_p1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

            sprintf( hname, "h_syst_%s_msb_met%d_m1s", syst_name[si], mbi+1 ) ;
            TH1F* h_msb_m1s = (TH1F*) infile.Get( hname ) ;
            if ( h_msb_m1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

            for ( int nbi=4; nbi>1; nbi-- ) {

               printf("     nb = %d\n", nbi ) ;

               float msig_nom = h_msig_nom -> GetBinContent( nbi ) ;
               float msig_p1s = h_msig_p1s -> GetBinContent( nbi ) ;
               float msig_m1s = h_msig_m1s -> GetBinContent( nbi ) ;

               float msb_nom  = h_msb_nom -> GetBinContent( nbi ) ;
               float msb_p1s  = h_msb_p1s -> GetBinContent( nbi ) ;
               float msb_m1s  = h_msb_m1s -> GetBinContent( nbi ) ;

               float msig_ave_frac_var = 1. ;
               if ( msig_nom > 0. ) { msig_ave_frac_var = 0.5* ( msig_p1s - msig_m1s ) / msig_nom ; }
               float msig_lands_par = 1. + msig_ave_frac_var ;

               float msb_ave_frac_var = 1. ;
               if ( msb_nom > 0. ) { msb_ave_frac_var = 0.5* ( msb_p1s - msb_m1s ) / msb_nom ; }
               float msb_lands_par = 1. + msb_ave_frac_var ;

               printf( "          msig :  m1s, nom, p1s = %.2f, %.2f, %.2f,   ave frac var = %5.2f\n",
                         msig_m1s, msig_nom, msig_p1s, msig_ave_frac_var ) ;
               printf( "          msb  :  m1s, nom, p1s = %.2f, %.2f, %.2f,   ave frac var = %5.2f\n",
                         msb_m1s, msb_nom, msb_p1s, msb_ave_frac_var ) ;

               fprintf( outfile_lands, "      %5.3f           -          %5.3f           -      ", msig_lands_par, msb_lands_par ) ;

            } // nbi.


         } // mbi.

         fprintf( outfile_lands, "\n" ) ;

      } // si.







     //-- Signal MC stat err, one line per observable (each is statistically independent).

      printf("\n\n ==== MC stat errors \n\n" ) ;

      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

         char hname[1000] ;

         sprintf( hname, "h_syst_mcstat_msig_met%d_nom", mbi+1 ) ;
         TH1F* h_msig_nom = (TH1F*) infile.Get( hname ) ;
         if ( h_msig_nom == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

         sprintf( hname, "h_syst_mcstat_msig_met%d_p1s", mbi+1 ) ;
         TH1F* h_msig_p1s = (TH1F*) infile.Get( hname ) ;
         if ( h_msig_p1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

         sprintf( hname, "h_syst_mcstat_msig_met%d_m1s", mbi+1 ) ;
         TH1F* h_msig_m1s = (TH1F*) infile.Get( hname ) ;
         if ( h_msig_m1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }


         sprintf( hname, "h_syst_mcstat_msb_met%d_nom", mbi+1 ) ;
         TH1F* h_msb_nom = (TH1F*) infile.Get( hname ) ;
         if ( h_msb_nom == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

         sprintf( hname, "h_syst_mcstat_msb_met%d_p1s", mbi+1 ) ;
         TH1F* h_msb_p1s = (TH1F*) infile.Get( hname ) ;
         if ( h_msb_p1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

         sprintf( hname, "h_syst_mcstat_msb_met%d_m1s", mbi+1 ) ;
         TH1F* h_msb_m1s = (TH1F*) infile.Get( hname ) ;
         if ( h_msb_m1s == 0x0 ) { printf("\n\n *** Can't find %s in %s.\n\n", hname, ws_root_file ) ; return ; }

         for ( int nbi=4; nbi>1; nbi-- ) {

            float msig_nom = h_msig_nom -> GetBinContent( nbi ) ;
            float msig_p1s = h_msig_p1s -> GetBinContent( nbi ) ;

            float msb_nom  = h_msb_nom -> GetBinContent( nbi ) ;
            float msb_p1s  = h_msb_p1s -> GetBinContent( nbi ) ;

            float msig_stat_err_lands_par = 1. ;
            if ( msig_nom > 0 ) { msig_stat_err_lands_par = msig_p1s / msig_nom ; }

            float msb_stat_err_lands_par = 1. ;
            if ( msb_nom > 0 ) { msb_stat_err_lands_par = msb_p1s / msb_nom ; }

            printf( "           MC stat err lands par, met%d, %db,    SIG = %5.3f     SB = %5.3f\n", mbi+1, nbi, msig_stat_err_lands_par, msb_stat_err_lands_par ) ;



            fprintf( outfile_lands, "MC_st_m%d_%db_SIG  lnN  ", mbi+1, nbi ) ;
            for ( int i=0; i<(mbi*3+(4-nbi)); i++ )     { fprintf( outfile_lands, "          -           -              -           -      " ) ; }
            fprintf( outfile_lands, "      %5.3f           -              -           -         ", msig_stat_err_lands_par ) ;
            for ( int i=0; i<((3-mbi)*3+(nbi-2)); i++ ) { fprintf( outfile_lands, "       -           -              -           -         " ) ; }
            fprintf( outfile_lands, "\n" ) ;


            fprintf( outfile_lands, "MC_st_m%d_%db_SB   lnN  ", mbi+1, nbi ) ;
            for ( int i=0; i<(mbi*3+(4-nbi)); i++ )     { fprintf( outfile_lands, "          -           -              -           -      " ) ; }
            fprintf( outfile_lands, "          -           -          %5.3f           -         ", msb_stat_err_lands_par ) ;
            for ( int i=0; i<((3-mbi)*3+(nbi-2)); i++ ) { fprintf( outfile_lands, "       -           -              -           -         " ) ; }
            fprintf( outfile_lands, "\n" ) ;


         } // nbi.

      } // mbi.




      fclose( outfile_lands ) ;


   } // gen_lands_from_ws

