
#include "TSystem.h"
#include <fstream>

   bool updateFileValue( const char* inFile,
                         const char* parameterName,
                         double newValue,
                         const char* ufvname = "ufv-output.txt" ) {


      ifstream infp ;
      infp.open( inFile ) ;
      if ( !infp.good() ) {
         printf("\n\n *** updateFileValue: Problem opening input file: %s.\n\n", inFile ) ;
         return false ;
      }

      FILE* outfp ;
      if ( (outfp=fopen( ufvname,"w"))==NULL ) {
         printf("\n\n *** updateFileValue: Problem opening output file.\n\n" ) ;
         return false ;
      }

      while ( infp.good() ) {

         TString line ;
         line.ReadLine( infp ) ;
         TObjArray* tokens = line.Tokenize(" ") ;
         //// printf(" number of fields : %d\n", tokens->GetEntries() ) ;
         TObjString* first_field = (TObjString*) (tokens->At(0)) ;
         if ( tokens->GetEntries() > 1 && (first_field->GetString()).CompareTo( parameterName ) == 0 ) {
            float value ;
            TObjString* second_field = (TObjString*) (tokens->At(1)) ;
            sscanf( (second_field->GetString()).Data(), "%f", &value ) ;
            printf(" Found %s.  Changing value from %g to %g\n", parameterName, value, newValue ) ;
            fprintf( outfp, "%s  %f\n", (first_field->GetString()).Data(), newValue ) ;
         } else {
            fprintf( outfp, "%s\n", line.Data() ) ;
         }


      } // still reading?


      fflush( outfp ) ;
      fclose( outfp ) ;

      char command[10000] ;
      sprintf( command, "mv %s  %s\n", ufvname, inFile ) ;
      gSystem->Exec( command ) ;

      return true ;

   }


