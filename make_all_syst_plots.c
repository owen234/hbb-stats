

   void make_all_syst_plots() {

      gROOT->LoadMacro("draw_syst_hists.c+") ;

      int nsigmass(14) ;
      int sig_mass[14] = { 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500 } ;

      int nsyst(7) ;
      char syst[7][10] = { "btagSF", "ISR", "JER", "JES", "PDF", "PU", "mcstat" } ;

      for ( int smi=0; smi<nsigmass; smi++ ) {

         char wsfile[10000] ;

         sprintf( wsfile, "outputfiles/ws-metsig-4metbin-w3b-wpu-csyst5-nosignal-sigmass-%d-withMSB1.root", sig_mass[smi] ) ;

         for ( int si=0; si<nsyst; si++ ) {
            draw_syst_hists( wsfile, syst[si], sig_mass[smi] ) ;
         } // si.

      } // smi.

   }


