
#include "TSystem.h"
#include <iostream>

//#include "getFileValue.h"

//=====================================================================================================

   bool getFileValue( const char* inFile,
                      const char* parameterName,
                      float& returnValue ) {


      returnValue = 1.0 ;

    // Include a blank space at the end to avoid multiple
    //  matches in grep when target is a substring found in other lines.

      char command[10000] ;
      sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
      TString commandOutput = gSystem->GetFromPipe( command ) ;

      /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

      char label[1000] ;
      float value ;
      sscanf( commandOutput.Data(), "%s %g", label, &value ) ;
      if ( strcmp( label, parameterName ) == 0 ) {
         printf(" Found %s.  Value is %g\n", parameterName, value ) ;
         returnValue = value ;
         return true ;
      }

      printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

      return false ;

   }


//=====================================================================================================

   bool getFileValueWithError( const char* inFile,
                      const char* parameterName,
                      float& returnValue, float& returnError ) {


      returnValue = 1.0 ;
      returnError = 0.0 ;

      char command[10000] ;
      sprintf( command, "grep \"%s \" %s\n", parameterName, inFile ) ;
      TString commandOutput = gSystem->GetFromPipe( command ) ;

      /// printf( " Output of command is : %s\n", commandOutput.Data() ) ;

      char label[1000] ;
      float value, err ;
      sscanf( commandOutput.Data(), "%s %g %g", label, &value, &err ) ;
      if ( strcmp( label, parameterName ) == 0 ) {
         printf(" Found %s.  Value is %g +/- %g\n", parameterName, value, err ) ;
         returnValue = value ;
         returnError = err ;
         return true ;
      }

      printf("\n\n *** Could not find parameter %s in file %s.\n\n", parameterName, inFile ) ;

      return false ;

   }

//=====================================================================================================

