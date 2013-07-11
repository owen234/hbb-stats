


void rundemosignif( const char* wsfile = "outputfiles/ws-nosig-250-1metbin.root", int calculatorType=0, int nToys=5000 ) {

   gROOT->LoadMacro("StandardHypoTestDemo.C+") ;

   StandardHypoTestDemo( wsfile,
                            "ws",
                            "SbModel",
                            "BModel",
                            "hbb_observed_rds",
                            calculatorType,
                            3,
                            nToys ) ;


}


