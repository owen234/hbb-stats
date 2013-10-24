
#include "TMath.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCanvas.h"


#include "getFileValue.c"
#include "updateFileValue.c"

#include <fstream>


   const int bins_of_nb(3) ;
   const int bins_of_nb_sigmc(4) ;
   const int max_bins_of_met(50) ;
   int       bins_of_met(4) ;

   float smc_msig_val[bins_of_nb_sigmc][max_bins_of_met] ;
   float smc_msb_val[bins_of_nb_sigmc][max_bins_of_met] ;

   float smc_msig_err[bins_of_nb_sigmc][max_bins_of_met] ;
   float smc_msb_err[bins_of_nb_sigmc][max_bins_of_met] ;

   char btag_catname[4][10] = { "NT", "2b", "3b", "4b" } ;

  //--- prototypes here.

   bool readSignalCounts( const char* susy_counts_file, float sig_mass ) ;

  //===========================================================================================

   void draw_susy_event_counts( const char* susy_counts_file = "outputfiles/susy-signal-counts-4metbin-w3b-wpu.txt", float sig_mass = 250. ) {

      gDirectory -> Delete( "h*" ) ;

    //-------------------------------------------------------------------------

      printf("\n\n Reading input file: %s\n\n", susy_counts_file ) ;

      float fileVal ;
      char pname[1000] ;
      char command[10000] ;


      if ( !readSignalCounts( susy_counts_file, sig_mass ) ) {
         printf("\n\n *** Can't find signal mass of %.0f in %s\n\n", sig_mass, susy_counts_file ) ;
         return ;
      }

      TH1F* h_msig_nom[bins_of_nb_sigmc] ;
      TH1F* h_msb_nom[bins_of_nb_sigmc] ;

      TH1F* h_msig_p1s[bins_of_nb_sigmc] ;
      TH1F* h_msb_p1s[bins_of_nb_sigmc] ;

      TH1F* h_msig_m1s[bins_of_nb_sigmc] ;
      TH1F* h_msb_m1s[bins_of_nb_sigmc] ;

      TCanvas* can1 = (TCanvas*) gDirectory -> FindObject("can1") ;
      if ( can1 == 0x0 ) {
         can1 = new TCanvas( "can1", "SUSY signal counts", 600, 800 ) ;
      }
      can1 -> Clear() ;
      can1 -> Divide( 2, 4 ) ;

      int canind(1) ;

      for ( int nbi=0; nbi<bins_of_nb_sigmc; nbi++ ) {

         char hname[1000] ;
         char htitle[1000] ;

         sprintf( hname, "h_msig_%s_nom", btag_catname[nbi] ) ;
         sprintf( htitle, "SIG, %s", btag_catname[nbi] ) ;
         h_msig_nom[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_msig_%s_p1s", btag_catname[nbi] ) ;
         sprintf( htitle, "SIG, %s, +1 sigma", btag_catname[nbi] ) ;
         h_msig_p1s[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_msig_%s_m1s", btag_catname[nbi] ) ;
         sprintf( htitle, "SIG, %s, -1 sigma", btag_catname[nbi] ) ;
         h_msig_m1s[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;


         sprintf( hname, "h_msb_%s_nom", btag_catname[nbi] ) ;
         sprintf( htitle, "SB, %s", btag_catname[nbi] ) ;
         h_msb_nom[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_msb_%s_p1s, +1 sigma", btag_catname[nbi] ) ;
         sprintf( htitle, "SB, %s", btag_catname[nbi] ) ;
         h_msb_p1s[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_msb_%s_m1s, -1 sigma", btag_catname[nbi] ) ;
         sprintf( htitle, "SB, %s", btag_catname[nbi] ) ;
         h_msb_m1s[nbi] = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            h_msig_nom[nbi] -> SetBinContent( mbi+1, smc_msig_val[nbi][mbi] ) ;
            h_msig_p1s[nbi] -> SetBinContent( mbi+1, smc_msig_val[nbi][mbi] + smc_msig_err[nbi][mbi] ) ;
            h_msig_m1s[nbi] -> SetBinContent( mbi+1, smc_msig_val[nbi][mbi] - smc_msig_err[nbi][mbi] ) ;

            h_msb_nom[nbi] -> SetBinContent( mbi+1, smc_msb_val[nbi][mbi] ) ;
            h_msb_p1s[nbi] -> SetBinContent( mbi+1, smc_msb_val[nbi][mbi] + smc_msb_err[nbi][mbi] ) ;
            h_msb_m1s[nbi] -> SetBinContent( mbi+1, smc_msb_val[nbi][mbi] - smc_msb_err[nbi][mbi] ) ;

            char binlabel[100] ;
            sprintf( binlabel, "S bin %d", mbi+1 ) ;

            h_msig_nom[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;
            h_msig_p1s[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;
            h_msig_m1s[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;

            h_msb_nom[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;
            h_msb_p1s[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;
            h_msb_m1s[nbi] -> GetXaxis() -> SetBinLabel( mbi+1, binlabel ) ;

         } // mbi.

         h_msig_nom[nbi] -> SetLineWidth(2) ;
         h_msig_p1s[nbi] -> SetLineWidth(2) ;
         h_msig_m1s[nbi] -> SetLineWidth(2) ;

         h_msb_nom[nbi] -> SetLineWidth(2) ;
         h_msb_p1s[nbi] -> SetLineWidth(2) ;
         h_msb_m1s[nbi] -> SetLineWidth(2) ;

         h_msig_nom[nbi] -> SetLineColor(1) ;
         h_msig_p1s[nbi] -> SetLineColor(2) ;
         h_msig_m1s[nbi] -> SetLineColor(4) ;

         h_msb_nom[nbi] -> SetLineColor(1) ;
         h_msb_p1s[nbi] -> SetLineColor(2) ;
         h_msb_m1s[nbi] -> SetLineColor(4) ;

         can1 -> cd( canind++ ) ;
         h_msb_p1s[nbi] -> Draw() ;
         h_msb_m1s[nbi] -> Draw( "same" ) ;
         h_msb_nom[nbi] -> Draw( "same" ) ;

         can1 -> cd( canind++ ) ;
         h_msig_p1s[nbi] -> Draw() ;
         h_msig_m1s[nbi] -> Draw( "same" ) ;
         h_msig_nom[nbi] -> Draw( "same" ) ;

      } // nbi.


   } // draw_susy_event_counts







  //==============================================================================================


   bool readSignalCounts( const char* susy_counts_file, float sig_mass ) {

      ifstream infp ;
      infp.open( susy_counts_file ) ;
      if ( !infp.good() ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", susy_counts_file ) ;
         return false ;
      }

      int ArraySize = 4 + bins_of_met * bins_of_nb_sigmc * 2 * 2 ; // 2 (sig vs sb) * 2 (vals and errs).

      char command[10000] ;
      sprintf(command, "tail -1 %s | awk '{print NF}' | grep -q %d", susy_counts_file, ArraySize ) ;
      int returnStat = gSystem->Exec(command ) ;
      if ( returnStat !=0 ) {
         printf("\n\n\n *** setSusyScanPoint : expecting %d fields per line in input file %s.  Found ", ArraySize, susy_counts_file ) ; cout << flush ;
         sprintf( command, "tail -1 %s | awk '{print NF}'", susy_counts_file ) ;
         gSystem->Exec(command ) ; cout << flush ;
         printf("\n\n") ;
         return false ;
      }

      float ArrayContent[ArraySize] ;

      bool found(false) ;

      while ( infp.good() ) {

         TString line ;
         line.ReadLine( infp ) ;
         TObjArray* tokens = line.Tokenize(" ") ;
         printf(" number of fields : %d\n", tokens->GetEntries() ) ;
         if ( tokens->GetEntries() != ArraySize ) continue ;

         for ( int i = 0 ; infp && i < ArraySize; ++ i ) {
            TObjString* str = (TObjString*) (tokens->At(i)) ;
            sscanf( (str->GetString()).Data(), "%f", &(ArrayContent[i]) )  ;
         }
         printf( " 1st column: %.0f\n", ArrayContent[0] ) ;

         if ( TMath::Nint( fabs( ArrayContent[0] - sig_mass ) ) == 0 ) {

            printf("  found mass point %.0f in file %s\n", sig_mass, susy_counts_file ) ;

            found = true ;

            break ;

         }

      } // still reading input file?

      if ( !found ) { printf("\n\n *** Did not find mass point %.0f in %s.\n\n", sig_mass, susy_counts_file ) ; return false ; }

      for ( int nbi=0; nbi<bins_of_nb_sigmc; nbi++ ) {
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
            smc_msig_val[nbi][mbi] = ArrayContent[ 4 + (                mbi)*(bins_of_nb_sigmc) + nbi ] ;
            smc_msb_val[nbi][mbi]  = ArrayContent[ 4 + (1*bins_of_met + mbi)*(bins_of_nb_sigmc) + nbi ] ;
            smc_msig_err[nbi][mbi] = ArrayContent[ 4 + (2*bins_of_met + mbi)*(bins_of_nb_sigmc) + nbi ] ;
            smc_msb_err[nbi][mbi]  = ArrayContent[ 4 + (3*bins_of_met + mbi)*(bins_of_nb_sigmc) + nbi ] ;
         } // mbi.
      } // nbi.

      printf("\n\n\n") ;
      printf("=============================================================================================================================================================================\n") ;
      printf("  METsig   |        4bSB              4bSIG         |       3bSB              3bSIG         |        2bSB              2bSIG        |        NTSB              NTSIG        |\n") ;
      printf("=============================================================================================================================================================================\n") ;
      fflush(stdout) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         printf(" met bin %d : ", mbi+1 ) ;
         for ( int nbi=(bins_of_nb_sigmc-1); nbi>=0; nbi-- ) {
            printf( "  %6.1f +/- %4.1f,  %6.1f +/- %4.1f    |",
                 smc_msb_val[nbi][mbi] , smc_msb_err[nbi][mbi],
                 smc_msig_val[nbi][mbi], smc_msig_err[nbi][mbi] ) ;
         } // nbi.
         printf("\n") ;
      } // mbi.
      printf("=====================================================================================================================================\n") ;



      return true ;


   } // readSignalCounts.

 //===============================================================================================================







