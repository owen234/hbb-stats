

//
//  This takes a background only input file as input and
//  produces a new input file where the observables are set to BG + SUSY values.
//  Takes the susy counts file from the input file to ensure consistency.
//  Adds a tag line to the new input file that warns user that signal was embedded and what signal mass was used.
//

#include "TMath.h"
#include "TSystem.h"

#include "RooPosDefCorrGauss.h"

#include "getFileValue.c"
#include "updateFileValue.c"

#include <fstream>


   const int bins_of_nb(3) ;
   const int max_bins_of_met(50) ;
   int       bins_of_met ;

   float smc_msig_val[bins_of_nb][max_bins_of_met] ;
   float smc_msb_val[bins_of_nb][max_bins_of_met] ;

   float smc_msig_err[bins_of_nb][max_bins_of_met] ;
   float smc_msb_err[bins_of_nb][max_bins_of_met] ;


  //--- prototypes here.

   bool readSignalCounts( const char* susy_counts_file, float sig_mass ) ;

  //===========================================================================================

   void add_susy_to_obs( const char* infile = "outputfiles/input-file.txt",
                         float sig_mass = 250.,
                         float susy_signal_strength = 1.,
                         const char* outfile = ""
                        ) {


    //-------------------------------------------------------------------------

      printf("\n\n Reading input file: %s\n\n", infile ) ;

      float fileVal ;
      char pname[1000] ;
      char command[10000] ;


      sprintf( pname, "bins_of_met" ) ;
      if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
      bins_of_met = TMath::Nint( fileVal ) ;

     //-- get signal input file and look for requested signal mass.
      char susy_counts_filename[10000] ;
      if ( !getFileStringValue( infile, "signal_counts_file", susy_counts_filename ) ) {
         printf("\n\n *** Can't find input susy counts file: signal_counts_file line of %s.\n\n", infile ) ;
         return ;
      }
      if ( !readSignalCounts( susy_counts_filename, sig_mass ) ) {
         printf("\n\n *** Can't find signal mass of %.0f in %s\n\n", sig_mass, susy_counts_filename ) ;
         return ;
      }


      char newoutfilename[10000] ;
      if ( strlen(outfile) == 0 ) {
         sprintf( command, "basename %s .txt", infile ) ;
         const char* filebase = gSystem -> GetFromPipe( command ) ;
         sprintf( newoutfilename, "outputfiles/%s-susy-%04d-embedded-sigstrength-%.1f.txt", filebase, TMath::Nint(sig_mass), susy_signal_strength ) ;
         printf( "\n\n New output filename: %s\n\n", newoutfilename ) ;
      } else {
         sprintf( newoutfilename, "%s", outfile ) ;
         printf( "\n\n Output filename: %s\n\n", newoutfilename ) ;
      }

      sprintf( command, "cp %s %s\n", infile, newoutfilename ) ;
      int returnStat = gSystem -> Exec( command ) ;
      if ( returnStat != 0 ) { printf("\n\n *** Problem creating output file %s\n\n", newoutfilename ) ; return ; }



      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            float inval, newval ;

            sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
            getFileValue( infile, pname, inval ) ;
            newval = inval + susy_signal_strength * smc_msig_val[nbi][mbi] ;
            printf("  %s : adding %6.1f susy counts to %6.1f BG counts.\n", pname, susy_signal_strength * smc_msig_val[nbi][mbi], inval ) ;
            updateFileValue( newoutfilename, pname, newval ) ;

            sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
            getFileValue( infile, pname, inval ) ;
            newval = inval + susy_signal_strength * smc_msb_val[nbi][mbi] ;
            printf("  %s : adding %6.1f susy counts to %6.1f BG counts.\n", pname, susy_signal_strength * smc_msb_val[nbi][mbi], inval ) ;
            updateFileValue( newoutfilename, pname, newval ) ;

         } // mbi.
      } // nbi.


      sprintf( command, "echo \"EMBEDDED_SUSY_SIGNAL %.0f , signal-strength %.1f , %s\" >> %s", 
                sig_mass, susy_signal_strength, susy_counts_filename, newoutfilename ) ;
      gSystem -> Exec( command ) ;

      printf("\n\n Created embed input file: %s\n\n\n", newoutfilename ) ;


   } // add_susy_to_obs.







  //==============================================================================================


   bool readSignalCounts( const char* susy_counts_file, float sig_mass ) {

      ifstream infp ;
      infp.open( susy_counts_file ) ;
      if ( !infp.good() ) {
         printf("\n\n *** Problem opening input file: %s.\n\n", susy_counts_file ) ;
         return false ;
      }

      int ArraySize = 2 + bins_of_met * bins_of_nb * 2 * 2 ; // 2 (sig vs sb) * 2 (vals and errs).

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

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
            smc_msig_val[nbi][mbi] = ArrayContent[ 2 + (               nbi)*(bins_of_met) + mbi ] ;
            smc_msb_val[nbi][mbi]  = ArrayContent[ 2 + (1*bins_of_nb + nbi)*(bins_of_met) + mbi ] ;
            smc_msig_err[nbi][mbi] = ArrayContent[ 2 + (2*bins_of_nb + nbi)*(bins_of_met) + mbi ] ;
            smc_msb_err[nbi][mbi]  = ArrayContent[ 2 + (3*bins_of_nb + nbi)*(bins_of_met) + mbi ] ;
         } // mbi.
      } // nbi.

      printf("\n\n\n") ;
      printf("=====================================================================================================================================\n") ;
      printf("  METsig   |        4bSB              4bSIG         |       3bSB              3bSIG         |        2bSB              2bSIG        |\n") ;
      printf("=====================================================================================================================================\n") ;
      fflush(stdout) ;
      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         printf(" met bin %d : ", mbi+1 ) ;
         for ( int nbi=(bins_of_nb-1); nbi>=0; nbi-- ) {
            printf( "  %6.1f +/- %4.1f,  %6.1f +/- %4.1f    |",
                 smc_msb_val[nbi][mbi] , smc_msb_val[nbi][mbi]  * smc_msb_err[nbi][mbi],
                 smc_msig_val[nbi][mbi], smc_msig_val[nbi][mbi] * smc_msig_err[nbi][mbi] ) ;
         } // nbi.
         printf("\n") ;
      } // mbi.
      printf("=====================================================================================================================================\n") ;



      return true ;


   } // readSignalCounts.

 //===============================================================================================================







