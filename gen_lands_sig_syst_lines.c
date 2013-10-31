
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"

   void gen_lands_sig_syst_lines( const char* ws_root_file ) {


      TFile infile( ws_root_file, "read" ) ;
      if ( ! infile.IsOpen() ) { printf("\n\n *** Can't open input workspace root file : %s\n\n", ws_root_file ) ; return ; }




      TString infile_ts( ws_root_file ) ;

      TString outfile_name = infile_ts.ReplaceAll(".root","-lands-sig-syst.txt") ;
      printf(" output file : %s\n", outfile_name.Data() ) ;

      FILE* outfile = fopen( outfile_name, "w" ) ;


     //-- All shape systematics, one line each.

      int bins_of_met(4) ;

      int nsyst(6) ;
      char syst_name[6][10] = { "btagSF", "ISR", "JER", "JES", "PDF", "PU" } ;


      for ( int si=0; si<nsyst; si++ ) {

         printf("\n\n ==== %s\n\n", syst_name[si] ) ;

         fprintf( outfile, " %10s  lnN   ", syst_name[si] ) ;

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

               fprintf( outfile, "  %5.3f   -    %5.3f   -   ", msig_lands_par, msb_lands_par ) ;

            } // nbi.


         } // mbi.

         fprintf( outfile, "\n" ) ;

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



            fprintf( outfile, "   MC_stat_met%d_%db_SIG   lnN  ", mbi+1, nbi ) ;
            for ( int i=0; i<(mbi*3+(4-nbi)); i++ )     { fprintf( outfile, "    -        -        -        -    " ) ; }
            fprintf( outfile, "  %5.3f      -        -        -    ", msig_stat_err_lands_par ) ;
            for ( int i=0; i<((3-mbi)*3+(nbi-2)); i++ ) { fprintf( outfile, "    -        -        -        -    " ) ; }
            fprintf( outfile, "\n" ) ;


            fprintf( outfile, "   MC_stat_met%d_%db_SB    lnN  ", mbi+1, nbi ) ;
            for ( int i=0; i<(mbi*3+(4-nbi)); i++ )     { fprintf( outfile, "    -        -        -        -    " ) ; }
            fprintf( outfile, "    -        -       %5.3f     -    ", msb_stat_err_lands_par ) ;
            for ( int i=0; i<((3-mbi)*3+(nbi-2)); i++ ) { fprintf( outfile, "    -        -        -        -    " ) ; }
            fprintf( outfile, "\n" ) ;


         } // nbi.

      } // mbi.




      fclose( outfile ) ;


   } // gen_lands_sig_syst_lines

