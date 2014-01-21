
#include "TChain.h"
#include "TH3D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"


   void gen_btag_sf_syst_files( ) {

      int lsp_mass = 1 ;

   // int bins_of_met = 4 ;
   // double met_bin_edges[5] = { 30., 50., 100., 150., 10000. } ;

      int bins_of_susymass = 1000 ;
      double susymass_bin_edges[1001] ;
      for ( int i=0; i<=1000; i++ ) { susymass_bin_edges[i] = i + 0.5 ; } // centered on integers between 1 and 1000.

      int bins_of_nbtag = 5 ;
      double nbtag_bin_edges[6] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5 } ;

      gDirectory -> Delete( "h*" ) ;

      TChain* sigmc_chain = new TChain( "reducedTree" ) ;

      //////sigmc_chain -> Add( "../reduced-trees-slim-oct22-v71-1-pj-s/*TChiHH_2b2b*.root" ) ;
      sigmc_chain -> Add( "../reduced-trees-jan08-v71-7-pj/reducedTree.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to500_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1992_v71-slim.root" ) ;


      printf("\n\n Setting up cuts.\n\n") ;

      char basiccuts[10000] ;
      sprintf( basiccuts, "cutPV==1&&passCleaning==1&&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40" ) ;

      char triggercuts[10000] ;
      sprintf( triggercuts, "(passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1||passMC_DiCentralPFJet30_PFMHT80==1)" ) ;
      char triggercuts_data[10000] ;
      sprintf( triggercuts_data, "(pass_DiCentralPFJet30_PFMET80_BTagCSV07==1||pass_PFMET150==1||pass_DiCentralPFJet30_PFMHT80==1)" ) ;

      char njetcuts[10000] ;
      sprintf( njetcuts, "njets20>=4&&njets20<=5" ) ;

      char skimcuts[10000] ;
      //--- It's safe both in data and MC to require both triggercuts and triggercuts_data.
      //    In data, passMC_* = pass_*.  In MC, pass_* is always true.
      sprintf( skimcuts, "((%s)&&(%s)&&(%s)&&(%s))", basiccuts, triggercuts, triggercuts_data, njetcuts ) ;


      char leptonveto[10000] ;
      sprintf( leptonveto, "%s", "nMuons==0&&nElectrons==0&&nIsoPFcands10_010==0&&nTausLoose==0" ) ;

      char maxdrcut[10000] ;
      sprintf( maxdrcut, "%s", "deltaRmax_hh<2.2" ) ;

      char mindphicut[10000] ;
      sprintf( mindphicut, "%s", "((METsig>50&&minDeltaPhi20_eta5_noIdAll_nobeta>0.3)||(METsig<50&&minDeltaPhi20_eta5_noIdAll_nobeta>0.5))" ) ;

      char jet2ptcut[10000] ;
      sprintf( jet2ptcut, "%s", "jetpt2>50" ) ;

      char masssigcuts[10000] ;
      sprintf( masssigcuts, "%s", "deltam<20&&avem>100&&avem<140" ) ;

      char masssbcuts[10000] ;
      sprintf( masssbcuts, "%s", "!(deltam<30&&avem>90&&avem<150)" ) ;

      char btag4bcuts[10000] ;
      sprintf( btag4bcuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4>0.244" ) ;

      char btag3bcuts[10000] ;
      sprintf( btag3bcuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4<0.244" ) ;

      char btag2bcuts[10000] ;
      sprintf( btag2bcuts, "%s", "CSVbest2>0.898&&CSVbest3<0.679" ) ;

      char allcommoncuts[10000] ;
      sprintf( allcommoncuts, "(%s)&&(%s)&&(%s)&&(%s)&&(%s)", skimcuts, leptonveto, maxdrcut, jet2ptcut, mindphicut ) ;

      char nomindphicuts[10000] ;
      sprintf( nomindphicuts, "(%s)&&(%s)&&(%s)&&(%s)", skimcuts, leptonveto, maxdrcut, jet2ptcut ) ;





     //--- fill histograms.

      TCanvas* can = new TCanvas("can","fill window") ;

      char cuts[100000] ;

      TH3D* h_msig_p1s = new TH3D( "h_msig_p1s", "higgs mass SIG, SF+1sig",  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;
      TH3D* h_msb_p1s  = new TH3D( "h_msb_p1s" , "higgs mass SB, SF+1sig" ,  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;

      sprintf( cuts, "(%s&&%s&&METsig>30)&&(m12==%d)", allcommoncuts, masssigcuts, lsp_mass ) ;
      printf(" h_msig_p1s : %s\n\n", cuts ) ;
      sigmc_chain -> Draw( "m0:nbtag_SFp1sig:nbtag_nomSF>>h_msig_p1s", cuts ) ;
      h_msig_p1s -> Draw("box") ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "(%s&&%s&&METsig>30)&&(m12==%d)", allcommoncuts, masssbcuts, lsp_mass ) ;
      printf(" h_msb_p1s : %s\n\n", cuts ) ;
      sigmc_chain -> Draw( "m0:nbtag_SFp1sig:nbtag_nomSF>>h_msb_p1s", cuts ) ;
      h_msb_p1s -> Draw("box") ;
      can -> Update() ; can -> Draw() ;



      TH3D* h_msig_m1s = new TH3D( "h_msig_m1s", "higgs mass SIG, SF-1sig",  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;
      TH3D* h_msb_m1s  = new TH3D( "h_msb_m1s" , "higgs mass SB, SF-1sig" ,  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;

      sprintf( cuts, "(%s&&%s&&METsig>30)&&(m12==%d)", allcommoncuts, masssigcuts, lsp_mass ) ;
      printf(" h_msig_m1s : %s\n\n", cuts ) ;
      sigmc_chain -> Draw( "m0:nbtag_SFm1sig:nbtag_nomSF>>h_msig_m1s", cuts ) ;
      h_msig_m1s -> Draw("box") ;
      can -> Update() ; can -> Draw() ;

      sprintf( cuts, "(%s&&%s&&METsig>30)&&(m12==%d)", allcommoncuts, masssbcuts, lsp_mass ) ;
      printf(" h_msb_m1s : %s\n\n", cuts ) ;
      sigmc_chain -> Draw( "m0:nbtag_SFm1sig:nbtag_nomSF>>h_msb_m1s", cuts ) ;
      h_msb_m1s -> Draw("box") ;
      can -> Update() ; can -> Draw() ;





     //--- compute fractions from histograms.


      TH3D* h_frac_msig_p1s = new TH3D( "h_frac_msig_p1s", "Transfer fraction, higgs mass SIG, SF+1sig",  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;
      TH3D* h_frac_msb_p1s  = new TH3D( "h_frac_msb_p1s" , "Transfer fraction, higgs mass SB, SF+1sig" ,  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;

      TH3D* h_frac_msig_m1s = new TH3D( "h_frac_msig_m1s", "Transfer fraction, higgs mass SIG, SF-1sig",  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;
      TH3D* h_frac_msb_m1s  = new TH3D( "h_frac_msb_m1s" , "Transfer fraction, higgs mass SB, SF-1sig" ,  bins_of_nbtag, nbtag_bin_edges,     bins_of_nbtag, nbtag_bin_edges,   bins_of_susymass, susymass_bin_edges ) ;


      for ( int hmi=1; hmi<=bins_of_susymass; hmi++ ) {

         for ( int fci=0; fci<bins_of_nbtag; fci++ ) {

            double sum_over_to_cat ;

           //--- SIG, +1 sigma
            sum_over_to_cat = 0. ;
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) { sum_over_to_cat += h_msig_p1s -> GetBinContent( fci+1, tci+1, hmi ) ; }
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
               double val(0.) ;
               if ( sum_over_to_cat > 0. ) { val = ( h_msig_p1s -> GetBinContent( fci+1, tci+1, hmi ) ) / sum_over_to_cat ; }
               h_frac_msig_p1s -> SetBinContent( fci+1, tci+1, hmi, val ) ;
            } // tci.

           //--- SB, +1 sigma
            sum_over_to_cat = 0. ;
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) { sum_over_to_cat += h_msb_p1s -> GetBinContent( fci+1, tci+1, hmi ) ; }
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
               double val(0.) ;
               if ( sum_over_to_cat > 0. ) { val = ( h_msb_p1s -> GetBinContent( fci+1, tci+1, hmi ) ) / sum_over_to_cat ; }
               h_frac_msb_p1s -> SetBinContent( fci+1, tci+1, hmi, val ) ;
            } // tci.


           //--- SIG, -1 sigma
            sum_over_to_cat = 0. ;
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) { sum_over_to_cat += h_msig_m1s -> GetBinContent( fci+1, tci+1, hmi ) ; }
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
               double val(0.) ;
               if ( sum_over_to_cat > 0. ) { val = ( h_msig_m1s -> GetBinContent( fci+1, tci+1, hmi ) ) / sum_over_to_cat ; }
               h_frac_msig_m1s -> SetBinContent( fci+1, tci+1, hmi, val ) ;
            } // tci.

           //--- SB, -1 sigma
            sum_over_to_cat = 0. ;
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) { sum_over_to_cat += h_msb_m1s -> GetBinContent( fci+1, tci+1, hmi ) ; }
            for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
               double val(0.) ;
               if ( sum_over_to_cat > 0. ) { val = ( h_msb_m1s -> GetBinContent( fci+1, tci+1, hmi ) ) / sum_over_to_cat ; }
               h_frac_msb_m1s -> SetBinContent( fci+1, tci+1, hmi, val ) ;
            } // tci.


         } // fci.

      } // hmi.



      gSystem -> Exec( "mkdir -p outputfiles" ) ;

      char outfile[10000] ;
      FILE* ofp ;

     //--- generate text output file, SIG.
      sprintf( outfile, "outputfiles/syst-btag-sf-SIG.txt" ) ;
      printf( "\n\n Creating output file : %s\n\n", outfile ) ;
      ofp = fopen( outfile, "w" ) ;

      for ( int hmi=1; hmi<=bins_of_susymass; hmi++ ) {

         if ( h_frac_msig_p1s -> Integral( 1, bins_of_nbtag,   1, bins_of_nbtag,    hmi, hmi ) > 0 ) {

            fprintf( ofp, " %5d  %5d    ", hmi, lsp_mass ) ;

           //--- SIG, +1 sigma
            for ( int fci=0; fci<bins_of_nbtag; fci++ ) {
               if ( fci == 1 ) continue ;
               for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
                  if ( tci == 1 ) continue ;
                  fprintf( ofp, " %.3f ", h_frac_msig_p1s -> GetBinContent( fci+1, tci+1, hmi ) ) ;
               }
               fprintf( ofp, "    " ) ;
            } // fci.

            fprintf( ofp, "            " ) ;

           //--- SIG, -1 sigma
            for ( int fci=0; fci<bins_of_nbtag; fci++ ) {
               if ( fci == 1 ) continue ;
               for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
                  if ( tci == 1 ) continue ;
                  fprintf( ofp, " %.3f ", h_frac_msig_m1s -> GetBinContent( fci+1, tci+1, hmi ) ) ; 
               }
               if ( fci<(bins_of_nbtag-1) ) { fprintf( ofp, "    " ) ; }
            } // fci.

            fprintf( ofp, "\n" ) ;

         } // any events for this susy mass value?

      } // hmi.

      fclose( ofp ) ;



     //--- generate text output file, SB.
      sprintf( outfile, "outputfiles/syst-btag-sf-SB.txt" ) ;
      printf( "\n\n Creating output file : %s\n\n", outfile ) ;
      ofp = fopen( outfile, "w" ) ;

      for ( int hmi=1; hmi<=bins_of_susymass; hmi++ ) {

         if ( h_frac_msb_p1s -> Integral( 1, bins_of_nbtag,   1, bins_of_nbtag,    hmi, hmi ) > 0 ) {

            fprintf( ofp, " %5d  %5d    ", hmi, lsp_mass ) ;

           //--- SB, +1 sigma
            for ( int fci=0; fci<bins_of_nbtag; fci++ ) {
               if ( fci == 1 ) continue ;
               for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
                  if ( tci == 1 ) continue ;
                  fprintf( ofp, " %.3f ", h_frac_msb_p1s -> GetBinContent( fci+1, tci+1, hmi ) ) ;
               }
               fprintf( ofp, "    " ) ;
            } // fci.

            fprintf( ofp, "            " ) ;

           //--- SB, -1 sigma
            for ( int fci=0; fci<bins_of_nbtag; fci++ ) {
               if ( fci == 1 ) continue ;
               for ( int tci=0; tci<bins_of_nbtag; tci++ ) {
                  if ( tci == 1 ) continue ;
                  fprintf( ofp, " %.3f ", h_frac_msb_m1s -> GetBinContent( fci+1, tci+1, hmi ) ) ; 
               }
               if ( fci<(bins_of_nbtag-1) ) { fprintf( ofp, "    " ) ; }
            } // fci.

            fprintf( ofp, "\n" ) ;

         } // any events for this susy mass value?

      } // hmi.

      fclose( ofp ) ;



      sprintf( outfile, "outputfiles/btag-sf-histograms.root" ) ;

      printf(" Creating %s\n", outfile ) ;

      TFile outroot( outfile, "recreate" ) ;

      h_msig_p1s -> Write() ;
      h_msb_p1s -> Write() ;
      h_msig_m1s -> Write() ;
      h_msb_m1s -> Write() ;

      h_frac_msig_p1s -> Write() ;
      h_frac_msb_p1s -> Write() ;
      h_frac_msig_m1s -> Write() ;
      h_frac_msb_m1s -> Write() ;

      outroot.Close() ;



   } // gen_btag_sf_syst_files

