
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
#include "histio.c"
#include "TMath.h"
#include "TString.h"



   TChain* sigchain ;

   float dataIntLumiIPB(20000.) ;

   float met_bin_edges_4bins[5] ;

   float met_bin_edges[100] ;




  //================================================================================================

   void gen_sig_input_file( const char* outfilename = "outputfiles/susy-signal-counts.txt",
                            const char* metvarname = "METsig",
                            bool usePUweight = false
                           ) {

      const int bins_of_met = 4 ;

      int lsp_mass(1) ;

      TString metvarname_nospecial(metvarname) ;
      metvarname_nospecial.ReplaceAll("/","_over_") ;
      metvarname_nospecial.ReplaceAll("(","_") ;
      metvarname_nospecial.ReplaceAll(")","") ;
      printf("\n\n %s\n\n", metvarname_nospecial.Data() ) ;


      ///bool fill_noweight_histograms(false) ;
      bool fill_noweight_histograms(true) ;

      if ( strcmp( metvarname, "METsig" ) == 0 ) {
         met_bin_edges_4bins[0] =  30. ;
         met_bin_edges_4bins[1] =  50. ;
         met_bin_edges_4bins[2] = 100. ;
         met_bin_edges_4bins[3] = 150. ;
         met_bin_edges_4bins[4] = 10000. ;
      } else if ( strcmp( metvarname, "METsig_2012" ) == 0 ) {
         met_bin_edges_4bins[0] =  25. ;
         met_bin_edges_4bins[1] =  40. ;
         met_bin_edges_4bins[2] =  85. ;
         met_bin_edges_4bins[3] = 120. ;
         met_bin_edges_4bins[4] = 10000. ;
      } else if ( strcmp( metvarname, "MET" ) == 0 ) {
         met_bin_edges_4bins[0] = 106. ;
         met_bin_edges_4bins[1] = 133. ;
         met_bin_edges_4bins[2] = 190. ;
         met_bin_edges_4bins[3] = 250. ;
         met_bin_edges_4bins[4] = 10000. ;
      } else if ( strcmp( metvarname, "MET/sqrt(HT30)" ) == 0 ) {
         met_bin_edges_4bins[0] = 5.9 ;
         met_bin_edges_4bins[1] = 7.7 ;
         met_bin_edges_4bins[2] = 10.8 ;
         met_bin_edges_4bins[3] = 13.0 ;
         met_bin_edges_4bins[4] = 10000. ;
      } else {
         printf("\n\n\n *** unrecognized met variable name : %s\n\n", metvarname ) ;
      }

      printf("\n\n  %s bins: ", metvarname ) ;
      for ( int bi=0; bi<=bins_of_met; bi++ ) { printf("  %.1f  ", met_bin_edges_4bins[bi]) ; }
      printf("\n\n") ;

      if ( bins_of_met == 4 ) {
         for ( int bi=0; bi<=bins_of_met; bi++ ) { met_bin_edges[bi] = met_bin_edges_4bins[bi] ; }
      } else {
         printf("\n\n *** Don't know how to do %d met bins.\n\n", bins_of_met ) ;
      }


      gDirectory -> Delete( "h*" ) ;




    //--- setup chains of reduced trees for the comps.

      printf("\n\n Setting up reduced tree chains.\n\n" ) ;


      //// char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-july11-2013-pt20" ;
      //// char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-skim-sept17-2013-v71-1s" ;
      char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-slim-oct08-2013-v71-5b" ;

      printf("\n\n\n   Reduced tree directory: %s\n\n\n", rtdir ) ;

      char pathandfile[10000] ;

      sigchain = new TChain("reducedTree") ;

     //--- old (private prod) signal MC.
 ///  sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.HiggsinoNLSP_chargino130_to_500_bino1-PU_S10-TChihh_v69-slimskim.root", rtdir ) ;
 ///  sigchain -> Add( pathandfile ) ;

     //--- new (official) signal MC.
      sprintf( pathandfile, "%s/reducedTree.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71_test-slimskim.root", rtdir ) ;
      sigchain -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-350to500_mLSP-1to370_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_test-slimskim.root", rtdir ) ;
      sigchain -> Add( pathandfile ) ;

      const int max_sig_points(100) ;
      float signal_weight[max_sig_points] ;
      float sigmass[max_sig_points] ;

      int nsigpoints(0) ;

    //-- See our twiki for an explanation of the weight factors.
    //    https://twiki.cern.ch/twiki/bin/view/CMS/SusyEWHHbbbb2013
    //
    //+++ This table is for the private production signal MC in the HiggsinoNLSP file
    //
  /// sigmass[nsigpoints] = 130. ;  signal_weight[nsigpoints] = 2.3045e-05 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 150. ;  signal_weight[nsigpoints] = 1.2294e-06 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 175. ;  signal_weight[nsigpoints] = 7.8243e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 200. ;  signal_weight[nsigpoints] = 4.6847e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 225. ;  signal_weight[nsigpoints] = 3.2248e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 250. ;  signal_weight[nsigpoints] = 7.1769e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 275. ;  signal_weight[nsigpoints] = 5.6463e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 300. ;  signal_weight[nsigpoints] = 2.9629e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 325. ;  signal_weight[nsigpoints] = 2.7941e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 350. ;  signal_weight[nsigpoints] = 2.4135e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 375. ;  signal_weight[nsigpoints] = 2.9794e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 400. ;  signal_weight[nsigpoints] = 1.8760e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 425. ;  signal_weight[nsigpoints] = 1.4117e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 450. ;  signal_weight[nsigpoints] = 1.0555e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 475. ;  signal_weight[nsigpoints] = 1.8952e-07 ;  nsigpoints++ ;
  /// sigmass[nsigpoints] = 500. ;  signal_weight[nsigpoints] = 1.4492e-07 ;  nsigpoints++ ;




    //+++ This table is for the official production files in SMS-TChiHH_2b2b files.
      sigmass[nsigpoints] = 130. ;  signal_weight[nsigpoints] = 2.1319e-05 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 150. ;  signal_weight[nsigpoints] = 7.1267e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 175. ;  signal_weight[nsigpoints] = 8.0680e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 200. ;  signal_weight[nsigpoints] = 4.9270e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 225. ;  signal_weight[nsigpoints] = 3.1318e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 250. ;  signal_weight[nsigpoints] = 5.0604e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 275. ;  signal_weight[nsigpoints] = 3.4774e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 300. ;  signal_weight[nsigpoints] = 2.3687e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 325. ;  signal_weight[nsigpoints] = 1.6850e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 350. ;  signal_weight[nsigpoints] = 2.1626e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 375. ;  signal_weight[nsigpoints] = 1.6090e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 400. ;  signal_weight[nsigpoints] = 1.1710e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 425. ;  signal_weight[nsigpoints] = 0.8738e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 450. ;  signal_weight[nsigpoints] = 0.6714e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 475. ;  signal_weight[nsigpoints] = 0.5033e-07 ;  nsigpoints++ ;
      sigmass[nsigpoints] = 500. ;  signal_weight[nsigpoints] = 0.3873e-07 ;  nsigpoints++ ;


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




     //--- Define cuts.

      printf("\n\n Setting up cuts.\n\n") ;

     //--- These are included in the skim definition of doSlimSkim.c.  Doesn't hurt to apply them here.
     //
      char basiccuts[10000] ;
      sprintf( basiccuts, "cutPV==1&&passCleaning==1&&buggyEvent==0&& MET/caloMET<2 && maxTOBTECjetDeltaMult<40" ) ;

      char triggercuts[10000] ;
      sprintf( triggercuts, "(passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1)" ) ;

      char njetcuts[10000] ;
      sprintf( njetcuts, "njets20>=4&&njets20<=5" ) ;

      char skimcuts[10000] ;
      sprintf( skimcuts, "((%s)&&(%s)&&(%s))", basiccuts, triggercuts, njetcuts ) ;


     //--- These are beyond the skim selection.

      char leptonveto[10000] ;
      /// sprintf( leptonveto, "%s", "nMuons==0&&nElectrons==0&&nIsoTracks15_005_03==0&&nTausLoose==0" ) ;
      sprintf( leptonveto, "%s", "nMuons==0&&nElectrons==0&&nIsoPFcands10_010==0&&nTausLoose==0" ) ;

      char drmaxcut[10000] ;
      sprintf( drmaxcut, "%s", "deltaRmax_hh<2.2" ) ;

      char mindphicut[10000] ;
      //// sprintf( mindphicut, "%s", "((METsig>50&&minDeltaPhi20>0.3)||(METsig<50&&minDeltaPhi20>0.5))" ) ;
      sprintf( mindphicut, "%s", "((METsig>50&&minDeltaPhi20_eta5_noIdAll_nobeta>0.3)||(METsig<50&&minDeltaPhi20_eta5_noIdAll_nobeta>0.5))" ) ;

      char jet2ptcut[10000] ;
      sprintf( jet2ptcut, "%s", "jetpt2>50" ) ;



      char masssigcuts[10000] ;
      sprintf( masssigcuts, "%s", "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140))" ) ;

      char masssbcuts[10000] ;
      sprintf( masssbcuts, "%s", "!(abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<30&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>90)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<150)))" ) ;

      char btag4bcuts[10000] ;
      sprintf( btag4bcuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4>0.244" ) ;

      char btag3bcuts[10000] ;
      sprintf( btag3bcuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4<0.244" ) ;

      char btag2bcuts[10000] ;
      sprintf( btag2bcuts, "%s", "CSVbest2>0.898&&CSVbest3<0.679" ) ;

      char allcommoncuts[10000] ;
      sprintf( allcommoncuts, "(%s)&&(%s)&&(%s)&&(%s)&&(%s)", skimcuts, leptonveto, drmaxcut, jet2ptcut, mindphicut ) ;





     //--- fill the histograms.

      TCanvas* can = new TCanvas("can","plots") ;

      printf("\n\n Filling the histograms.\n\n" ) ; fflush(stdout) ;

      char arg1[1000] ;
      char allcuts[10000] ;

      char puweight[100] ;
      if ( usePUweight ) {
         sprintf( puweight, "*PUweight" ) ;
      } else {
         sprintf( puweight, "" ) ;
      }

      TH1F* h_4b_msig[max_sig_points] ;
      TH1F* h_4b_msb[max_sig_points] ;
      TH1F* h_3b_msig[max_sig_points] ;
      TH1F* h_3b_msb[max_sig_points] ;
      TH1F* h_2b_msig[max_sig_points] ;
      TH1F* h_2b_msb[max_sig_points] ;

      TH1F* h_4b_msig_nw[max_sig_points] ;
      TH1F* h_4b_msb_nw[max_sig_points] ;
      TH1F* h_3b_msig_nw[max_sig_points] ;
      TH1F* h_3b_msb_nw[max_sig_points] ;
      TH1F* h_2b_msig_nw[max_sig_points] ;
      TH1F* h_2b_msb_nw[max_sig_points] ;

      for ( int spi=0; spi<nsigpoints; spi++ ) {

         char hname[1000] ;
         char htitle[1000] ;

         printf("\n\n ++++++++ Signal point, mgl = %.0f, sample weight = %.5f\n\n", sigmass[spi], signal_weight[spi]*dataIntLumiIPB ) ;

         sprintf( hname, "h_%s_4b_msig_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 4b, mass signal box, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_4b_msig[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msig[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssigcuts, btag4bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msig[spi] -> Print("all") ;



         sprintf( hname, "h_%s_3b_msig_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 3b, mass signal box, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_3b_msig[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msig[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssigcuts, btag3bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msig[spi] -> Print("all") ;



         sprintf( hname, "h_%s_2b_msig_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 2b, mass signal box, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_2b_msig[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msig[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssigcuts, btag2bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msig[spi] -> Print("all") ;






         sprintf( hname, "h_%s_4b_msb_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 4b, mass sideband, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_4b_msb[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msb[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssbcuts, btag4bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msb[spi] -> Print("all") ;



         sprintf( hname, "h_%s_3b_msb_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 3b, mass sideband, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_3b_msb[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msb[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssbcuts, btag3bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msb[spi] -> Print("all") ;



         sprintf( hname, "h_%s_2b_msb_smc_%.0f", metvarname_nospecial.Data(), sigmass[spi] ) ;
         sprintf( htitle, "%s, 2b, mass sideband, SUSY MC, mgl=%.0f", metvarname, sigmass[spi] ) ;
         h_2b_msb[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msb[spi] -> Sumw2() ;

         sprintf( arg1, "%s>>%s", metvarname, hname ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)*%g%s*%.0f", allcommoncuts, masssbcuts, btag2bcuts, sigmass[spi], lsp_mass, signal_weight[spi], puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         sigchain -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msb[spi] -> Print("all") ;







         if ( fill_noweight_histograms ) {

            sprintf( hname, "h_%s_4b_msig_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 4b, mass signal box, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_4b_msig_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_4b_msig_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssigcuts, btag4bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_4b_msig_nw[spi] -> Print("all") ;



            sprintf( hname, "h_%s_3b_msig_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 3b, mass signal box, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_3b_msig_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_3b_msig_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssigcuts, btag3bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_3b_msig_nw[spi] -> Print("all") ;



            sprintf( hname, "h_%s_2b_msig_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 2b, mass signal box, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_2b_msig_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_2b_msig_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssigcuts, btag2bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_2b_msig_nw[spi] -> Print("all") ;






            sprintf( hname, "h_%s_4b_msb_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 4b, mass sideband, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_4b_msb_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_4b_msb_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssbcuts, btag4bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_4b_msb_nw[spi] -> Print("all") ;



            sprintf( hname, "h_%s_3b_msb_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 3b, mass sideband, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_3b_msb_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_3b_msb_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssbcuts, btag3bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_3b_msb_nw[spi] -> Print("all") ;



            sprintf( hname, "h_%s_2b_msb_smc_%.0f_noweight", metvarname_nospecial.Data(), sigmass[spi] ) ;
            sprintf( htitle, "%s, 2b, mass sideband, SUSY MC, mgl=%.0f, no weight", metvarname, sigmass[spi] ) ;
            h_2b_msb_nw[spi] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_2b_msb_nw[spi] -> Sumw2() ;

            sprintf( arg1, "%s>>%s", metvarname, hname ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s)&&m0==%.0f&&m12==%d)", allcommoncuts, masssbcuts, btag2bcuts, sigmass[spi], lsp_mass ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            sigchain -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_2b_msb_nw[spi] -> Print("all") ;


         } // fill_noweight_histograms ?


      } // spi.





     //--- Print a nice table.


      float smc_4b_msig_val[max_sig_points][bins_of_met] ;
      float smc_3b_msig_val[max_sig_points][bins_of_met] ;
      float smc_2b_msig_val[max_sig_points][bins_of_met] ;

      float smc_4b_msb_val[max_sig_points][bins_of_met] ;
      float smc_3b_msb_val[max_sig_points][bins_of_met] ;
      float smc_2b_msb_val[max_sig_points][bins_of_met] ;

      float smc_4b_msig_err[max_sig_points][bins_of_met] ;
      float smc_3b_msig_err[max_sig_points][bins_of_met] ;
      float smc_2b_msig_err[max_sig_points][bins_of_met] ;

      float smc_4b_msb_err[max_sig_points][bins_of_met] ;
      float smc_3b_msb_err[max_sig_points][bins_of_met] ;
      float smc_2b_msb_err[max_sig_points][bins_of_met] ;



      printf("\n\n\n") ;
      printf("============================================================================================================================================\n") ;
      printf(" Mgl    METsig   |        4bSB              4bSIG        |        3bSB              3bSIG         |        2bSB              2bSIG         |\n") ;
      printf("============================================================================================================================================\n") ;
      fflush(stdout) ;

      for ( int spi=0; spi<nsigpoints; spi++ ) {

         printf(" Mgl = %.0f ----------------------------\n", sigmass[spi] ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            float metsiglow  = h_4b_msig[0]->GetXaxis()->GetBinLowEdge( mbi+1 ) ;
            float metsighigh = h_4b_msig[0]->GetXaxis()->GetBinLowEdge( mbi+2 ) ;
            if ( mbi==bins_of_met-1 ) metsighigh = 999. ;


            smc_4b_msig_val[spi][mbi] = h_4b_msig[spi]->GetBinContent(mbi+1) ;
            smc_3b_msig_val[spi][mbi] = h_3b_msig[spi]->GetBinContent(mbi+1) ;
            smc_2b_msig_val[spi][mbi] = h_2b_msig[spi]->GetBinContent(mbi+1) ;

            smc_4b_msb_val[spi][mbi] = h_4b_msb[spi]->GetBinContent(mbi+1) ;
            smc_3b_msb_val[spi][mbi] = h_3b_msb[spi]->GetBinContent(mbi+1) ;
            smc_2b_msb_val[spi][mbi] = h_2b_msb[spi]->GetBinContent(mbi+1) ;

           //-- define error as fractional error on yield.

            smc_4b_msig_err[spi][mbi] = 1. ;
            smc_3b_msig_err[spi][mbi] = 1. ;
            smc_2b_msig_err[spi][mbi] = 1. ;

            smc_4b_msb_err[spi][mbi]  = 1. ;
            smc_3b_msb_err[spi][mbi]  = 1. ;
            smc_2b_msb_err[spi][mbi]  = 1. ;

            if ( smc_4b_msig_val[spi][mbi] > 0. ) { smc_4b_msig_err[spi][mbi] = ( h_4b_msig[spi]->GetBinError(mbi+1) ) / ( smc_4b_msig_val[spi][mbi] )  ; }
            if ( smc_3b_msig_val[spi][mbi] > 0. ) { smc_3b_msig_err[spi][mbi] = ( h_3b_msig[spi]->GetBinError(mbi+1) ) / ( smc_3b_msig_val[spi][mbi] )  ; }
            if ( smc_2b_msig_val[spi][mbi] > 0. ) { smc_2b_msig_err[spi][mbi] = ( h_2b_msig[spi]->GetBinError(mbi+1) ) / ( smc_2b_msig_val[spi][mbi] )  ; }

            if ( smc_4b_msb_val[spi][mbi] > 0. ) { smc_4b_msb_err[spi][mbi] = ( h_4b_msb[spi]->GetBinError(mbi+1) ) / ( smc_4b_msb_val[spi][mbi] )  ; }
            if ( smc_3b_msb_val[spi][mbi] > 0. ) { smc_3b_msb_err[spi][mbi] = ( h_3b_msb[spi]->GetBinError(mbi+1) ) / ( smc_3b_msb_val[spi][mbi] )  ; }
            if ( smc_2b_msb_val[spi][mbi] > 0. ) { smc_2b_msb_err[spi][mbi] = ( h_2b_msb[spi]->GetBinError(mbi+1) ) / ( smc_2b_msb_val[spi][mbi] )  ; }


            printf( "%4.0f   [%3.0f,%3.0f] |  %6.1f +/- %4.1f,  %6.1f +/- %4.1f    |  %6.1f +/- %4.1f,  %6.1f +/- %4.1f     |  %6.1f +/- %4.1f,  %6.1f +/- %4.1f     |\n",
               sigmass[spi], metsiglow, metsighigh,
               smc_4b_msb_val[spi][mbi]  , smc_4b_msb_val[spi][mbi]  * smc_4b_msb_err[spi][mbi],
               smc_4b_msig_val[spi][mbi] , smc_4b_msig_val[spi][mbi] * smc_4b_msig_err[spi][mbi],
               smc_3b_msb_val[spi][mbi]  , smc_3b_msb_val[spi][mbi]  * smc_3b_msb_err[spi][mbi],
               smc_3b_msig_val[spi][mbi] , smc_3b_msig_val[spi][mbi] * smc_3b_msig_err[spi][mbi],
               smc_2b_msb_val[spi][mbi]  , smc_2b_msb_val[spi][mbi]  * smc_2b_msb_err[spi][mbi],
               smc_2b_msig_val[spi][mbi] , smc_2b_msig_val[spi][mbi] * smc_2b_msig_err[spi][mbi]
             ) ;

         } // mbi.

      } // spi.

      printf("\n\n\n") ;




     //-- print out no-weight tables, if the histograms were filled.
      if ( fill_noweight_histograms ) {

         printf("\n\n\n") ;
         printf("================================================================================\n") ;
         printf(" Mgl    METsig   |   4bSB   4bSIG    |   3bSB   3bSIG     |   2bSB   2bSIG     |\n") ;
         printf("================================================================================\n") ;
         fflush(stdout) ;

         for ( int spi=0; spi<nsigpoints; spi++ ) {

            printf(" Mgl = %.0f ----------------------------\n", sigmass[spi] ) ;

            for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {

               float metsiglow  = h_4b_msig[0]->GetXaxis()->GetBinLowEdge( hbi ) ;
               float metsighigh = h_4b_msig[0]->GetXaxis()->GetBinLowEdge( hbi+1 ) ;
               if ( hbi==bins_of_met ) metsighigh = 999. ;


               float nmsig_4b_val_smc = h_4b_msig_nw[spi]->GetBinContent(hbi) ;
               float nmsig_3b_val_smc = h_3b_msig_nw[spi]->GetBinContent(hbi) ;
               float nmsig_2b_val_smc = h_2b_msig_nw[spi]->GetBinContent(hbi) ;

               float nmsb_4b_val_smc = h_4b_msb_nw[spi]->GetBinContent(hbi) ;
               float nmsb_3b_val_smc = h_3b_msb_nw[spi]->GetBinContent(hbi) ;
               float nmsb_2b_val_smc = h_2b_msb_nw[spi]->GetBinContent(hbi) ;


               printf( "%4.0f   [%3.0f,%3.0f] |  %6.0f %6.0f    |  %6.0f %6.0f     |  %6.0f %6.0f     |\n",
                  sigmass[spi], metsiglow, metsighigh,
                  nmsb_4b_val_smc, nmsig_4b_val_smc,
                  nmsb_3b_val_smc, nmsig_3b_val_smc,
                  nmsb_2b_val_smc, nmsig_2b_val_smc
                ) ;

            } // hbi.

         } // spi.

         printf("\n\n\n") ;




      } // fill_noweight_histograms ?








     //--- Generate the likelihood input file.

      gSystem -> Exec("mkdir -p outputfiles") ;

      printf("\n\n output file with path: %s\n\n", outfilename ) ; fflush(stdout) ;




      char command[10000] ;
      sprintf( command, "ls %s >& /dev/null", outfilename ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat == 0 ) {
         char mvfile[10000] ;
         sprintf( mvfile, "%s-old", outfilename ) ;
         printf("\n\n *** Output file already exists.  Moving it to %s\n\n", mvfile ) ;
         sprintf( command, "mv %s %s", outfilename, mvfile ) ;
         gSystem->Exec( command ) ;
      }
      FILE* outfile = fopen( outfilename, "w" ) ;

      fprintf( outfile, "met_variable  %s\n", metvarname ) ;
      fprintf( outfile, "bins_of_met  %d  ", bins_of_met ) ;
      for ( int mbi=0; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "   %.2f   ", met_bin_edges[mbi] ) ; }
      fprintf( outfile, "\n" ) ;

      for ( int spi=0; spi<nsigpoints; spi++ ) {

         fprintf( outfile, "  %5.0f  %5.0f    ", sigmass[spi], 0. ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_2b_msig_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_3b_msig_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_4b_msig_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_2b_msb_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_3b_msb_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %7.2f ", smc_4b_msb_val[spi][mbi] ) ; } fprintf( outfile, "    " ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_2b_msig_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_3b_msig_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_4b_msig_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_2b_msb_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_3b_msb_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) { fprintf( outfile, " %5.3f ", smc_4b_msb_err[spi][mbi] ) ; } fprintf( outfile, "    " ) ;

         fprintf( outfile, "\n" ) ;

      } // spi.


      fclose( outfile ) ;


      printf("\n\n Created likelihood input file: %s\n\n", outfilename ) ;






     //--- save histograms.

      printf("\n\n Saving histograms to outputfiles/gen_sig_input.root\n\n") ;

      saveHist( "outputfiles/gen_sig_input.root", "h*" ) ;


   } // gen_sig_input_file







