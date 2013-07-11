


void rundemo( const char* wsfile = "outputfiles/ws-nosig-250-1metbin.root", int calculatorType=3, float poiMax=1., float poiMin=0., int npoints=10, int ntoys=1000 ) {

   gROOT->LoadMacro("StandardHypoTestInvDemo.C+") ;

   StandardHypoTestInvDemo( wsfile,
                            "ws",
                            "SbModel",
                            "BModel",
                            "hbb_observed_rds",
                            calculatorType,
                            3,
                            true,
                            npoints,
                            poiMin,
                            poiMax,
                            ntoys ) ;


}


