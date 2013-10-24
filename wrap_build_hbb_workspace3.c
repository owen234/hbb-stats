

  void wrap_build_hbb_workspace3(
                              const char* infile = "outputfiles/input-file.txt",
                              const char* outfile = "outputfiles/ws.root",
                              float sig_mass = 250.,
                              bool use3b = true,
                              bool combine_top_metbins = false,
                              int arg_syst_type = 2, // 1 = Gaussian, 2 = log-normal
                              bool drop_first_met_bin = false
                             ) {

     gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
     gROOT->LoadMacro("RooAsymAbsProd.cxx+") ;
     gROOT->LoadMacro("build_hbb_workspace3.c+") ;

     build_hbb_workspace3( infile, outfile, sig_mass, use3b, combine_top_metbins, arg_syst_type, drop_first_met_bin ) ;


  }

