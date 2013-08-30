
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
#include "histio.c"
#include "TMath.h"
#include "TString.h"


  //--- Note: the signal info is now done in gen_sig_input_file.c

 //----------------
 //int  nbgcomps(3) ;
 //char bgcompname[3][100] = { "tt", "znn", "qcd" } ;
 //TChain* bgcompchain[3] ;
 //----------------
   int  nbgcomps(2) ;
   char bgcompname[2][100] = { "tt", "znn" } ;
   TChain* bgcompchain[2] ;
 //----------------


   float dataIntLumiIPB(20000.) ;




   float met_bin_edges_1bins[2] = { 50., 10000. } ;

   float met_bin_edges_4bins[5] ;

   float met_bin_edges[100] ;




  //================================================================================================

   void gen_input_file2( const char* outfilename = "outputfiles/input-file.txt",
                        int bins_of_met = 4,
                        float min_met = 0.,
                        const char* metvarname = "METsig",
                        bool use3b = true,
                        bool usePUweight = true,
                        int  closure_syst_option = 1
                        ) {


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

      if ( bins_of_met == 1 ) {
         for ( int bi=0; bi<=bins_of_met; bi++ ) { met_bin_edges[bi] = met_bin_edges_1bins[bi] ; }
         met_bin_edges[0] = min_met ;
      } else if ( bins_of_met == 4 ) {
         for ( int bi=0; bi<=bins_of_met; bi++ ) { met_bin_edges[bi] = met_bin_edges_4bins[bi] ; }
      } else {
         printf("\n\n *** Don't know how to do %d met bins.\n\n", bins_of_met ) ;
      }


      gDirectory -> Delete( "h*" ) ;




    //--- setup chains of reduced trees for the comps.

      printf("\n\n Setting up reduced tree chains.\n\n" ) ;

      for ( int si=0; si<nbgcomps; si++ ) { bgcompchain[si] = new TChain("reducedTree") ; }

      /////char rtdir[10000] = "/data/cms/hadronic-susy-bjets/hbb/reduced-trees-may23-2013" ;
      ///char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-july08-2013" ;
      char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-july11-2013-pt20" ;
      //char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-july22-2013-nobeta" ;

      printf("\n\n\n   Reduced tree directory: %s\n\n\n", rtdir ) ;

      int compIndex(0) ;

      char pathandfile[10000] ;

     //--- ttbar, 1 and 2 lepton
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1800_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1799_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      compIndex++ ;

     //--- Znn
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1832_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1833_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1834_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1831_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1835_v69-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      compIndex++ ;

      if ( nbgcomps > 2 ) {
        //--- QCD
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1820_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1814_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1821_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1815_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1822_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1816_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1817_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1818_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1819_v69-slimskim.root", rtdir ) ;
         bgcompchain[compIndex] -> Add( pathandfile ) ;
         compIndex++ ;
      }




     //--- Define cuts.

      printf("\n\n Setting up cuts.\n\n") ;

     //--- These are included in the skim definition of doSlimSkim.c.  Doesn't hurt to apply them here.

      char basiccuts[10000] ;
      sprintf( basiccuts, "cutPV==1&&passCleaning==1&&buggyEvent==0" ) ;

      char triggercuts[10000] ;
      sprintf( triggercuts, "(passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_PFMET150==1)" ) ;

      char njetcuts[10000] ;
      sprintf( njetcuts, "njets20>=4&&njets20<=5" ) ;

      char skimcuts[10000] ;
      sprintf( skimcuts, "((%s)&&(%s)&&(%s))", basiccuts, triggercuts, njetcuts ) ;


     //--- These are beyond the skim selection.

      char masssigcuts[10000] ;
      sprintf( masssigcuts, "%s", "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140))" ) ;

      char masssbcuts[10000] ;
      sprintf( masssbcuts, "%s", "!(abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<30&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>90)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<150)))" ) ;

      char btag4cuts[10000] ;
      sprintf( btag4cuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4>0.244" ) ;

      char btag3cuts[10000] ;
      sprintf( btag3cuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4<0.244" ) ;

      char btag2cuts[10000] ;
      sprintf( btag2cuts, "%s", "CSVbest2>0.898&&CSVbest3<0.679" ) ;

      char leptonveto[10000] ;
      sprintf( leptonveto, "%s", "nMuons==0&&nElectrons==0&&nIsoTracks15_005_03==0&&nTausLoose==0" ) ;

      char drmaxcut[10000] ;
      sprintf( drmaxcut, "%s", "deltaRmax_hh<2.2" ) ;

      char mindphicut[10000] ;
      sprintf( mindphicut, "%s", "minDeltaPhi20>0.3" ) ;

      char jet2ptcut[10000] ;
      sprintf( jet2ptcut, "%s", "jetpt2>50" ) ;

      char allcommoncuts[10000] ;
      sprintf( allcommoncuts, "(%s)&&(%s)&&(%s)&&(%s)&&(%s)", skimcuts, leptonveto, drmaxcut, jet2ptcut, mindphicut ) ;




     //--- prepare histograms.

      printf("\n\n Booking histograms.\n\n") ;

      TH1F* h_4b_msig_bg[10] ;
      TH1F* h_3b_msig_bg[10] ;
      TH1F* h_2b_msig_bg[10] ;

      TH1F* h_4b_msb_bg[10] ;
      TH1F* h_3b_msb_bg[10] ;
      TH1F* h_2b_msb_bg[10] ;

      TH1F* h_4b_msig_bg_noweight[10] ;
      TH1F* h_3b_msig_bg_noweight[10] ;
      TH1F* h_2b_msig_bg_noweight[10] ;

      TH1F* h_4b_msb_bg_noweight[10] ;
      TH1F* h_3b_msb_bg_noweight[10] ;
      TH1F* h_2b_msb_bg_noweight[10] ;

      char hname[1000] ;
      char htitle[1000] ;

      for ( int si=0; si<nbgcomps; si++ ) {

         sprintf( hname, "h_%s_4b_msig_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 4b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_4b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_3b_msig_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 3b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_3b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_2b_msig_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 2b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_2b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_4b_msb_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 4b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_4b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msb_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_3b_msb_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 3b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_3b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msb_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_2b_msb_bg_%s", metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( htitle, "%s, 2b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_2b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msb_bg[si] -> Sumw2() ;

         if ( fill_noweight_histograms ) {
            sprintf( hname, "h_%s_4b_msig_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 4b, mass signal box, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_4b_msig_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_4b_msig_bg_noweight[si] -> Sumw2() ;

            sprintf( hname, "h_%s_3b_msig_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 3b, mass signal box, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_3b_msig_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_3b_msig_bg_noweight[si] -> Sumw2() ;

            sprintf( hname, "h_%s_2b_msig_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 2b, mass signal box, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_2b_msig_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_2b_msig_bg_noweight[si] -> Sumw2() ;

            sprintf( hname, "h_%s_4b_msb_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 4b, mass sideband, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_4b_msb_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_4b_msb_bg_noweight[si] -> Sumw2() ;

            sprintf( hname, "h_%s_3b_msb_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 3b, mass sideband, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_3b_msb_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_3b_msb_bg_noweight[si] -> Sumw2() ;

            sprintf( hname, "h_%s_2b_msb_bg_%s_noweight", metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( htitle, "%s, 2b, mass sideband, background, %s, no weighting", metvarname, bgcompname[si] ) ;
            h_2b_msb_bg_noweight[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
            h_2b_msb_bg_noweight[si] -> Sumw2() ;
         } // fill_noweight_histograms ?

      } // si.








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

      for ( int si=0; si<nbgcomps; si++ ) {

         printf("\n\n ++++++++ %s component\n\n", bgcompname[si] ) ;

         sprintf( arg1, "%s>>h_%s_4b_msig_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssigcuts, btag4cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msig_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_3b_msig_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssigcuts, btag3cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msig_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_2b_msig_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssigcuts, btag2cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msig_bg[si] -> Print("all") ;


         sprintf( arg1, "%s>>h_%s_4b_msb_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssbcuts, btag4cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msb_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_3b_msb_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssbcuts, btag3cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msb_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_2b_msb_bg_%s", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s))*weight3%s*%.0f", allcommoncuts, masssbcuts, btag2cuts, puweight, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msb_bg[si] -> Print("all") ;

         if ( fill_noweight_histograms ) {

            sprintf( arg1, "%s>>h_%s_4b_msig_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssigcuts, btag4cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_4b_msig_bg_noweight[si] -> Print("all") ;

            sprintf( arg1, "%s>>h_%s_3b_msig_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssigcuts, btag3cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_3b_msig_bg_noweight[si] -> Print("all") ;

            sprintf( arg1, "%s>>h_%s_2b_msig_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssigcuts, btag2cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_2b_msig_bg_noweight[si] -> Print("all") ;


            sprintf( arg1, "%s>>h_%s_4b_msb_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssbcuts, btag4cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_4b_msb_bg_noweight[si] -> Print("all") ;

            sprintf( arg1, "%s>>h_%s_3b_msb_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssbcuts, btag3cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_3b_msb_bg_noweight[si] -> Print("all") ;

            sprintf( arg1, "%s>>h_%s_2b_msb_bg_%s_noweight", metvarname, metvarname_nospecial.Data(), bgcompname[si] ) ;
            sprintf( allcuts, "((%s)&&(%s)&&(%s))", allcommoncuts, masssbcuts, btag2cuts ) ;
            printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
            bgcompchain[si] -> Draw( arg1, allcuts ) ;
            can->Update() ; can->Draw() ;
            h_2b_msb_bg_noweight[si] -> Print("all") ;

         } // fill_noweight_histograms ?


      } // si.






     //--- Print a nice table.


      float n_4b_msig[50] ;
      float n_3b_msig[50] ;
      float n_2b_msig[50] ;

      float n_4b_msb[50] ;
      float n_3b_msb[50] ;
      float n_2b_msb[50] ;

      float Rsigsb_4b_val[50] ;
      float Rsigsb_3b_val[50] ;
      float Rsigsb_2b_val[50] ;

      float Rsigsb_4b_err[50] ;
      float Rsigsb_3b_err[50] ;
      float Rsigsb_2b_err[50] ;

      float Rsigsb_4b_staterr2[50] ;
      float Rsigsb_3b_staterr2[50] ;
      float Rsigsb_2b_staterr2[50] ;

     //-- these metbinsum variables are for the top N-1 met bins.

      float bgsum_nmsig_4b_metbinsum_val(0.) ;
      float bgsum_nmsig_3b_metbinsum_val(0.) ;
      float bgsum_nmsig_2b_metbinsum_val(0.) ;

      float bgsum_nmsb_4b_metbinsum_val(0.) ;
      float bgsum_nmsb_3b_metbinsum_val(0.) ;
      float bgsum_nmsb_2b_metbinsum_val(0.) ;

      float bgsum_nmsig_4b_metbinsum_err2(0.) ;
      float bgsum_nmsig_3b_metbinsum_err2(0.) ;
      float bgsum_nmsig_2b_metbinsum_err2(0.) ;

      float bgsum_nmsb_4b_metbinsum_err2(0.) ;
      float bgsum_nmsb_3b_metbinsum_err2(0.) ;
      float bgsum_nmsb_2b_metbinsum_err2(0.) ;


      printf("\n\n\n") ;
      printf("============================================================================================================================================================================\n") ;
      printf(" METsig    comp  |   4bSB   4bSIG     4bSIG/4bSB    |   3bSB   3bSIG     3bSIG/3bSB     |   2bSB   2bSIG     2bSIG/2bSB     |    BG val.        3b pred.       2b pred.     |\n") ;
      printf("============================================================================================================================================================================\n") ;
      fflush(stdout) ;
      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {

         float metsiglow  = h_4b_msig_bg[0]->GetXaxis()->GetBinLowEdge( hbi ) ;
         float metsighigh = h_4b_msig_bg[0]->GetXaxis()->GetBinLowEdge( hbi+1 ) ;
         if ( hbi==bins_of_met ) metsighigh = 999. ;

         float bgsum_nmsig_4b_val(0.) ; float bgsum_nmsig_4b_err2(0.) ;
         float bgsum_nmsig_3b_val(0.) ; float bgsum_nmsig_3b_err2(0.) ;
         float bgsum_nmsig_2b_val(0.) ; float bgsum_nmsig_2b_err2(0.) ;

         float bgsum_nmsb_4b_val(0.) ; float bgsum_nmsb_4b_err2(0.) ;
         float bgsum_nmsb_3b_val(0.) ; float bgsum_nmsb_3b_err2(0.) ;
         float bgsum_nmsb_2b_val(0.) ; float bgsum_nmsb_2b_err2(0.) ;

         for ( int si=0; si<nbgcomps; si++ ) {

            float nmsig_4b_val = h_4b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_4b_err = h_4b_msig_bg[si]->GetBinError(hbi) ;
            float nmsig_3b_val = h_3b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_3b_err = h_3b_msig_bg[si]->GetBinError(hbi) ;
            float nmsig_2b_val = h_2b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_2b_err = h_2b_msig_bg[si]->GetBinError(hbi) ;

            float nmsb_4b_val = h_4b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_4b_err = h_4b_msb_bg[si]->GetBinError(hbi) ;
            float nmsb_3b_val = h_3b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_3b_err = h_3b_msb_bg[si]->GetBinError(hbi) ;
            float nmsb_2b_val = h_2b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_2b_err = h_2b_msb_bg[si]->GetBinError(hbi) ;


            printf( "[%3.0f,%3.0f] %6s |  %6.1f %6.1f                   |  %6.1f %6.1f                    |  %6.1f %6.1f                    |                                               |\n",
               metsiglow, metsighigh,
               bgcompname[si],
               nmsb_4b_val, nmsig_4b_val,
               nmsb_3b_val, nmsig_3b_val,
               nmsb_2b_val, nmsig_2b_val
             ) ;

            if ( si < nbgcomps ) {
               bgsum_nmsig_4b_val += nmsig_4b_val ;
               bgsum_nmsig_3b_val += nmsig_3b_val ;
               bgsum_nmsig_2b_val += nmsig_2b_val ;
               bgsum_nmsig_4b_err2 += ::pow( nmsig_4b_err, 2.) ;
               bgsum_nmsig_3b_err2 += ::pow( nmsig_3b_err, 2.) ;
               bgsum_nmsig_2b_err2 += ::pow( nmsig_2b_err, 2.) ;
               bgsum_nmsb_4b_val += nmsb_4b_val ;
               bgsum_nmsb_3b_val += nmsb_3b_val ;
               bgsum_nmsb_2b_val += nmsb_2b_val ;
               bgsum_nmsb_4b_err2 += ::pow( nmsb_4b_err, 2.) ;
               bgsum_nmsb_3b_err2 += ::pow( nmsb_3b_err, 2.) ;
               bgsum_nmsb_2b_err2 += ::pow( nmsb_2b_err, 2.) ;
               if ( hbi > 1 ) {
                  bgsum_nmsig_4b_metbinsum_val += nmsig_4b_val ;
                  bgsum_nmsig_3b_metbinsum_val += nmsig_3b_val ;
                  bgsum_nmsig_2b_metbinsum_val += nmsig_2b_val ;
                  bgsum_nmsig_4b_metbinsum_err2 += ::pow( nmsig_4b_err, 2.) ;
                  bgsum_nmsig_3b_metbinsum_err2 += ::pow( nmsig_3b_err, 2.) ;
                  bgsum_nmsig_2b_metbinsum_err2 += ::pow( nmsig_2b_err, 2.) ;
                  bgsum_nmsb_4b_metbinsum_val += nmsb_4b_val ;
                  bgsum_nmsb_3b_metbinsum_val += nmsb_3b_val ;
                  bgsum_nmsb_2b_metbinsum_val += nmsb_2b_val ;
                  bgsum_nmsb_4b_metbinsum_err2 += ::pow( nmsb_4b_err, 2.) ;
                  bgsum_nmsb_3b_metbinsum_err2 += ::pow( nmsb_3b_err, 2.) ;
                  bgsum_nmsb_2b_metbinsum_err2 += ::pow( nmsb_2b_err, 2.) ;
               }
            }

            if ( si == (nbgcomps-1) ) {

               Rsigsb_4b_val[hbi] =  bgsum_nmsig_4b_val / bgsum_nmsb_4b_val ;
               Rsigsb_4b_err[hbi] =  Rsigsb_4b_val[hbi] * sqrt( bgsum_nmsig_4b_err2 / ::pow(bgsum_nmsig_4b_val,2) + bgsum_nmsb_4b_err2 / ::pow(bgsum_nmsb_4b_val,2) ) ;
               Rsigsb_3b_val[hbi] =  bgsum_nmsig_3b_val / bgsum_nmsb_3b_val ;
               Rsigsb_3b_err[hbi] =  Rsigsb_3b_val[hbi] * sqrt( bgsum_nmsig_3b_err2 / ::pow(bgsum_nmsig_3b_val,2) + bgsum_nmsb_3b_err2 / ::pow(bgsum_nmsb_3b_val,2) ) ;
               Rsigsb_2b_val[hbi] =  bgsum_nmsig_2b_val / bgsum_nmsb_2b_val ;
               Rsigsb_2b_err[hbi] =  Rsigsb_2b_val[hbi] * sqrt( bgsum_nmsig_2b_err2 / ::pow(bgsum_nmsig_2b_val,2) + bgsum_nmsb_2b_err2 / ::pow(bgsum_nmsb_2b_val,2) ) ;

               Rsigsb_4b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_4b_val > 0. ) { Rsigsb_4b_staterr2[hbi] += Rsigsb_4b_val[hbi] * Rsigsb_4b_val[hbi] / bgsum_nmsig_4b_val ; }
               if ( bgsum_nmsb_4b_val  > 0. ) { Rsigsb_4b_staterr2[hbi] += Rsigsb_4b_val[hbi] * Rsigsb_4b_val[hbi] / bgsum_nmsb_4b_val ; }
               if ( Rsigsb_4b_staterr2[hbi] < 0.0001 ) { Rsigsb_4b_staterr2[hbi] = 0.15 ; }

               Rsigsb_3b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_3b_val > 0. ) { Rsigsb_3b_staterr2[hbi] += Rsigsb_3b_val[hbi] * Rsigsb_3b_val[hbi] / bgsum_nmsig_3b_val ; }
               if ( bgsum_nmsb_3b_val  > 0. ) { Rsigsb_3b_staterr2[hbi] += Rsigsb_3b_val[hbi] * Rsigsb_3b_val[hbi] / bgsum_nmsb_3b_val ; }
               if ( Rsigsb_3b_staterr2[hbi] < 0.0001 ) { Rsigsb_3b_staterr2[hbi] = 0.15 ; }

               Rsigsb_2b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_2b_val > 0. ) { Rsigsb_2b_staterr2[hbi] += Rsigsb_2b_val[hbi] * Rsigsb_2b_val[hbi] / bgsum_nmsig_2b_val ; }
               if ( bgsum_nmsb_2b_val  > 0. ) { Rsigsb_2b_staterr2[hbi] += Rsigsb_2b_val[hbi] * Rsigsb_2b_val[hbi] / bgsum_nmsb_2b_val ; }
               if ( Rsigsb_2b_staterr2[hbi] < 0.0001 ) { Rsigsb_2b_staterr2[hbi] = 0.15 ; }

               float bgpred_from3b_val = bgsum_nmsb_4b_val * Rsigsb_3b_val[hbi] ;
               float bgpred_from3b_err = bgpred_from3b_val * sqrt( bgsum_nmsb_4b_err2/::pow(bgsum_nmsb_4b_val,2.) + ::pow( Rsigsb_3b_err[hbi]/Rsigsb_3b_val[hbi],2.) ) ;

               float bgpred_from2b_val = bgsum_nmsb_4b_val * Rsigsb_2b_val[hbi] ;
               float bgpred_from2b_err = bgpred_from2b_val * sqrt( bgsum_nmsb_4b_err2/::pow(bgsum_nmsb_4b_val,2.) + ::pow( Rsigsb_2b_err[hbi]/Rsigsb_2b_val[hbi],2.) ) ;

               printf( "[%3.0f,%3.0f] %6s |  %6.1f %6.1f   %5.3f +/- %5.3f |  %6.1f %6.1f   %5.3f +/- %5.3f  |  %6.1f %6.1f   %5.3f +/- %5.3f  | %4.1f +/- %4.1f   %4.1f +/- %4.1f  %4.1f +/- %4.1f  |\n",
                  metsiglow, metsighigh,
                  "bg sum",
                  bgsum_nmsb_4b_val, bgsum_nmsig_4b_val,   Rsigsb_4b_val[hbi], Rsigsb_4b_err[hbi],
                  bgsum_nmsb_3b_val, bgsum_nmsig_3b_val,   Rsigsb_3b_val[hbi], Rsigsb_3b_err[hbi],
                  bgsum_nmsb_2b_val, bgsum_nmsig_2b_val,   Rsigsb_2b_val[hbi], Rsigsb_2b_err[hbi],
                  bgsum_nmsig_4b_val, sqrt(bgsum_nmsig_4b_err2),
                  bgpred_from3b_val, bgpred_from3b_err,
                  bgpred_from2b_val, bgpred_from2b_err
                ) ;
                printf("----------------------------------------------------------------------------------------------------------------------------|                                               |\n") ;
            }

         } // si.

         printf("============================================================================================================================================================================\n") ;

         n_4b_msig[hbi] = bgsum_nmsig_4b_val ;
         n_3b_msig[hbi] = bgsum_nmsig_3b_val ;
         n_2b_msig[hbi] = bgsum_nmsig_2b_val ;
         n_4b_msb[hbi]  = bgsum_nmsb_4b_val ;
         n_3b_msb[hbi]  = bgsum_nmsb_3b_val ;
         n_2b_msb[hbi]  = bgsum_nmsb_2b_val ;

      } // hbi.

      printf("\n\n\n") ;

     //-- Compute ratios for the combination of the top N-1 met bins.

      float Rsigsb_4b_metbinsum_val =  bgsum_nmsig_4b_metbinsum_val / bgsum_nmsb_4b_metbinsum_val ;
      float Rsigsb_4b_metbinsum_err =  Rsigsb_4b_metbinsum_val * sqrt( bgsum_nmsig_4b_metbinsum_err2 / ::pow(bgsum_nmsig_4b_metbinsum_val,2) + bgsum_nmsb_4b_metbinsum_err2 / ::pow(bgsum_nmsb_4b_metbinsum_val,2) ) ;
      float Rsigsb_3b_metbinsum_val =  bgsum_nmsig_3b_metbinsum_val / bgsum_nmsb_3b_metbinsum_val ;
      float Rsigsb_3b_metbinsum_err =  Rsigsb_3b_metbinsum_val * sqrt( bgsum_nmsig_3b_metbinsum_err2 / ::pow(bgsum_nmsig_3b_metbinsum_val,2) + bgsum_nmsb_3b_metbinsum_err2 / ::pow(bgsum_nmsb_3b_metbinsum_val,2) ) ;
      float Rsigsb_2b_metbinsum_val =  bgsum_nmsig_2b_metbinsum_val / bgsum_nmsb_2b_metbinsum_val ;
      float Rsigsb_2b_metbinsum_err =  Rsigsb_2b_metbinsum_val * sqrt( bgsum_nmsig_2b_metbinsum_err2 / ::pow(bgsum_nmsig_2b_metbinsum_val,2) + bgsum_nmsb_2b_metbinsum_err2 / ::pow(bgsum_nmsb_2b_metbinsum_val,2) ) ;

      float Rsigsb_4b_metbinsum_staterr2 = 0. ;
      if ( bgsum_nmsig_4b_metbinsum_val > 0. ) { Rsigsb_4b_metbinsum_staterr2 += Rsigsb_4b_metbinsum_val * Rsigsb_4b_metbinsum_val / bgsum_nmsig_4b_metbinsum_val ; }
      if ( bgsum_nmsb_4b_metbinsum_val  > 0. ) { Rsigsb_4b_metbinsum_staterr2 += Rsigsb_4b_metbinsum_val * Rsigsb_4b_metbinsum_val / bgsum_nmsb_4b_metbinsum_val ; }
      if ( Rsigsb_4b_metbinsum_staterr2 < 0.0001 ) { Rsigsb_4b_metbinsum_staterr2 = 0.15 ; }

      float Rsigsb_3b_metbinsum_staterr2 = 0. ;
      if ( bgsum_nmsig_3b_metbinsum_val > 0. ) { Rsigsb_3b_metbinsum_staterr2 += Rsigsb_3b_metbinsum_val * Rsigsb_3b_metbinsum_val / bgsum_nmsig_3b_metbinsum_val ; }
      if ( bgsum_nmsb_3b_metbinsum_val  > 0. ) { Rsigsb_3b_metbinsum_staterr2 += Rsigsb_3b_metbinsum_val * Rsigsb_3b_metbinsum_val / bgsum_nmsb_3b_metbinsum_val ; }
      if ( Rsigsb_3b_metbinsum_staterr2 < 0.0001 ) { Rsigsb_3b_metbinsum_staterr2 = 0.15 ; }

      float Rsigsb_2b_metbinsum_staterr2 = 0. ;
      if ( bgsum_nmsig_2b_metbinsum_val > 0. ) { Rsigsb_2b_metbinsum_staterr2 += Rsigsb_2b_metbinsum_val * Rsigsb_2b_metbinsum_val / bgsum_nmsig_2b_metbinsum_val ; }
      if ( bgsum_nmsb_2b_metbinsum_val  > 0. ) { Rsigsb_2b_metbinsum_staterr2 += Rsigsb_2b_metbinsum_val * Rsigsb_2b_metbinsum_val / bgsum_nmsb_2b_metbinsum_val ; }
      if ( Rsigsb_2b_metbinsum_staterr2 < 0.0001 ) { Rsigsb_2b_metbinsum_staterr2 = 0.15 ; }

      printf("\n\n\n") ;
      printf("  SIG/SB ratios from top %d %s bins.\n", bins_of_met-1, metvarname ) ;
      printf("   R SIG/SB, 4b : %.3f +/- %.3f MC stats (%.3f expected data stat err)\n", Rsigsb_4b_metbinsum_val, Rsigsb_4b_metbinsum_err, sqrt(Rsigsb_4b_metbinsum_staterr2) ) ;
      printf("   R SIG/SB, 3b : %.3f +/- %.3f MC stats (%.3f expected data stat err)\n", Rsigsb_3b_metbinsum_val, Rsigsb_3b_metbinsum_err, sqrt(Rsigsb_3b_metbinsum_staterr2) ) ;
      printf("   R SIG/SB, 2b : %.3f +/- %.3f MC stats (%.3f expected data stat err)\n", Rsigsb_2b_metbinsum_val, Rsigsb_2b_metbinsum_err, sqrt(Rsigsb_2b_metbinsum_staterr2) ) ;
      printf("\n\n\n") ;





     //-- print out no-weight tables, if the histograms were filled.
      if ( fill_noweight_histograms ) {

          for ( int si=0; si<nbgcomps; si++ ) {

             printf("\n\n\n") ;

             printf("\n\n\n") ;
             printf("=============================================================================================================================\n") ;
             printf(" METsig    comp  |   4bSB   4bSIG     4bSIG/4bSB    |   3bSB   3bSIG     3bSIG/3bSB     |   2bSB   2bSIG     2bSIG/2bSB     |\n") ;
             printf("=============================================================================================================================\n") ;
             fflush(stdout) ;
             for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {
                float metsiglow  = h_4b_msig_bg_noweight[0]->GetXaxis()->GetBinLowEdge( hbi ) ;
                float metsighigh = h_4b_msig_bg_noweight[0]->GetXaxis()->GetBinLowEdge( hbi+1 ) ;
                if ( hbi==bins_of_met ) metsighigh = 999. ;
                float n4bsb  = h_4b_msb_bg_noweight[si]  -> GetBinContent(hbi) ;
                float n4bsig = h_4b_msig_bg_noweight[si] -> GetBinContent(hbi) ;
                float n3bsb  = h_3b_msb_bg_noweight[si]  -> GetBinContent(hbi) ;
                float n3bsig = h_3b_msig_bg_noweight[si] -> GetBinContent(hbi) ;
                float n2bsb  = h_2b_msb_bg_noweight[si]  -> GetBinContent(hbi) ;
                float n2bsig = h_2b_msig_bg_noweight[si] -> GetBinContent(hbi) ;
                float r4bval = 0. ;  float r4berr = 0. ;
                float r3bval = 0. ;  float r3berr = 0. ;
                float r2bval = 0. ;  float r2berr = 0. ;
                if ( n4bsb > 0 ) { r4bval = n4bsig / n4bsb ;  if ( n4bsig > 0 ) { r4berr = r4bval * sqrt( 1./n4bsig + 1./n4bsb ) ; } }
                if ( n3bsb > 0 ) { r3bval = n3bsig / n3bsb ;  if ( n3bsig > 0 ) { r3berr = r3bval * sqrt( 1./n3bsig + 1./n3bsb ) ; } }
                if ( n2bsb > 0 ) { r2bval = n2bsig / n2bsb ;  if ( n2bsig > 0 ) { r2berr = r2bval * sqrt( 1./n2bsig + 1./n2bsb ) ; } }

                printf( "[%3.0f,%3.0f] %6s |  %6.0f %6.0f   %5.3f +/- %5.3f |  %6.0f %6.0f   %5.3f +/- %5.3f  |  %6.0f %6.0f   %5.3f +/- %5.3f  |\n",
                   metsiglow, metsighigh,
                   bgcompname[si],
                   n4bsb, n4bsig, r4bval, r4berr,
                   n3bsb, n3bsig, r3bval, r3berr,
                   n2bsb, n2bsig, r2bval, r2berr
                 ) ;

             } // hbi.
             printf("=============================================================================================================================\n") ;

             printf("\n\n\n") ;

          } // si.

      } // fill_noweight_histograms ?





      //-- Compute SIG/SB corrections and uncertainties.
      //
      //  The fit will do a weighted average of the 2b, 3b, and 4b samples, so corrections
      //  should be made with respect to that average.  In computing the weighted average
      //  below, use sqrt(N) for uncertainty, not MC stats error, since that's what the
      //  fit will do.
      //
      //  Taking the uncertainty on the correction as the (MC stat error / value) on the difference
      //  with respect to the weighted ave combined with half of the difference with respect
      //  to the weighted ave in quadrature.
      //
      //  If bins_of_met is set to 1, only use 2b and 4b (ignore 3b).

      printf("\n\n") ;
      float wave_Rsigsb_val[50] ;
      float wave_Rsigsb_err[50] ;
      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {
         if ( bins_of_met > 1 && use3b ) {
            wave_Rsigsb_val[hbi] =  ( Rsigsb_2b_val[hbi] / Rsigsb_2b_staterr2[hbi] + Rsigsb_3b_val[hbi] / Rsigsb_3b_staterr2[hbi] + Rsigsb_4b_val[hbi] / Rsigsb_4b_staterr2[hbi] ) /
                                    (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_3b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ;
            wave_Rsigsb_err[hbi] = sqrt( 1. / (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_3b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ) ;
         } else {
            wave_Rsigsb_val[hbi] =  ( Rsigsb_2b_val[hbi] / Rsigsb_2b_staterr2[hbi] + Rsigsb_4b_val[hbi] / Rsigsb_4b_staterr2[hbi] ) /
                                    (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ;
            wave_Rsigsb_err[hbi] = sqrt( 1. / (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ) ;
         }
         printf(" expected sqrt(N) weighted Rsig/sb ave for METsig bin %d :  %5.3f +/- %5.3f\n", hbi, wave_Rsigsb_val[hbi], wave_Rsigsb_err[hbi] ) ;
      } // hbi

      //
      //  Add combination of top N-1 met bins, in case we want to use that.
      //
      float wave_Rsigsb_metbinsum_val(0.) ;
      float wave_Rsigsb_metbinsum_err(0.) ;

      if ( bins_of_met > 1 && use3b ) {
         wave_Rsigsb_metbinsum_val =  (  Rsigsb_2b_metbinsum_val / Rsigsb_2b_metbinsum_staterr2 + Rsigsb_3b_metbinsum_val / Rsigsb_3b_metbinsum_staterr2 + Rsigsb_4b_metbinsum_val / Rsigsb_4b_metbinsum_staterr2 ) /
                                      (                      1.  / Rsigsb_2b_metbinsum_staterr2 +                     1.  / Rsigsb_3b_metbinsum_staterr2 +                     1.  / Rsigsb_4b_metbinsum_staterr2 ) ;
         wave_Rsigsb_metbinsum_err = sqrt( 1. / (                      1.  / Rsigsb_2b_metbinsum_staterr2 +                     1.  / Rsigsb_3b_metbinsum_staterr2 +                     1.  / Rsigsb_4b_metbinsum_staterr2 ) ) ;
      } else {
         wave_Rsigsb_metbinsum_val =  (  Rsigsb_2b_metbinsum_val / Rsigsb_2b_metbinsum_staterr2 + Rsigsb_4b_metbinsum_val / Rsigsb_4b_metbinsum_staterr2 ) /
                                      (                      1.  / Rsigsb_2b_metbinsum_staterr2 +                     1.  / Rsigsb_4b_metbinsum_staterr2 ) ;
         wave_Rsigsb_metbinsum_err = sqrt( 1. / (                      1.  / Rsigsb_2b_metbinsum_staterr2 +                     1.  / Rsigsb_4b_metbinsum_staterr2 ) ) ;
      }
      printf("\n") ;
      printf(" expected sqrt(N) weighted Rsig/sb ave for higest %d METsig bins:  %5.3f +/- %5.3f\n", bins_of_met-1, wave_Rsigsb_metbinsum_val, wave_Rsigsb_metbinsum_err ) ;
      printf("\n") ;



      float correction_Rsigsb_4b[50] ;
      float correction_Rsigsb_3b[50] ;
      float correction_Rsigsb_2b[50] ;

      float syst_Rsigsb_4b[50] ;
      float syst_Rsigsb_3b[50] ;
      float syst_Rsigsb_2b[50] ;

      printf("\n\n") ;

      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {

         correction_Rsigsb_4b[hbi] = Rsigsb_4b_val[hbi] / wave_Rsigsb_val[hbi] ;
         correction_Rsigsb_3b[hbi] = Rsigsb_3b_val[hbi] / wave_Rsigsb_val[hbi] ;
         correction_Rsigsb_2b[hbi] = Rsigsb_2b_val[hbi] / wave_Rsigsb_val[hbi] ;

         if ( closure_syst_option == 1 ) {

           //-- nominal: mcstat + 1/2 correction
            syst_Rsigsb_4b[hbi] = sqrt( ::pow( (Rsigsb_4b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_4b[hbi]-1.)/2.,2. ) ) ;
            syst_Rsigsb_3b[hbi] = sqrt( ::pow( (Rsigsb_3b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_3b[hbi]-1.)/2.,2. ) ) ;
            syst_Rsigsb_2b[hbi] = sqrt( ::pow( (Rsigsb_2b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_2b[hbi]-1.)/2.,2. ) ) ;

         } else if ( closure_syst_option == 2 ) {

           //-- no syst.
            syst_Rsigsb_4b[hbi] = 0. ;
            syst_Rsigsb_3b[hbi] = 0. ;
            syst_Rsigsb_2b[hbi] = 0. ;

         } else if ( closure_syst_option == 3 ) {

           //-- aggressive: mcstat only
            syst_Rsigsb_4b[hbi] =  (Rsigsb_4b_err[hbi]/Rsigsb_4b_val[hbi]) ;
            syst_Rsigsb_3b[hbi] =  (Rsigsb_3b_err[hbi]/Rsigsb_4b_val[hbi]) ;
            syst_Rsigsb_2b[hbi] =  (Rsigsb_2b_err[hbi]/Rsigsb_4b_val[hbi]) ;

         } else if ( closure_syst_option == 4 ) {

           //-- nominal: mcstat + full correction
            syst_Rsigsb_4b[hbi] = sqrt( ::pow( (Rsigsb_4b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_4b[hbi]-1.),2. ) ) ;
            syst_Rsigsb_3b[hbi] = sqrt( ::pow( (Rsigsb_3b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_3b[hbi]-1.),2. ) ) ;
            syst_Rsigsb_2b[hbi] = sqrt( ::pow( (Rsigsb_2b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  ::pow( (correction_Rsigsb_2b[hbi]-1.),2. ) ) ;

         } else {

            printf("\n\n\n *** Unknown value for closure_syst_option argument.  Supported values are:\n\n" ) ;
            printf("   1  MC stat error + 1/2 non-closure in quadrature.\n") ;
            printf("   2  No closure syst.\n") ;
            printf("   3  MC stat error only.\n") ;
            printf("   4  MC stat error + full non-closure in quadrature.\n") ;
            printf("\n\n\n") ;

            return ;

         }

         if ( bins_of_met > 1 && use3b ) {
            printf("  Rsig/sb corrections for METsig bin %d :     4b = %5.3f +/- %5.3f,       3b =  %5.3f +/- %5.3f,      2b =  %5.3f +/- %5.3f\n",
                  hbi,
                  correction_Rsigsb_4b[hbi], syst_Rsigsb_4b[hbi],
                  correction_Rsigsb_3b[hbi], syst_Rsigsb_3b[hbi],
                  correction_Rsigsb_2b[hbi], syst_Rsigsb_2b[hbi] ) ;
         } else {
            printf("  Rsig/sb corrections for METsig bin %d :     4b = %5.3f +/- %5.3f,       2b =  %5.3f +/- %5.3f\n",
                  hbi,
                  correction_Rsigsb_4b[hbi], syst_Rsigsb_4b[hbi],
                  correction_Rsigsb_2b[hbi], syst_Rsigsb_2b[hbi] ) ;
         }

      } // hbi.

      float correction_Rsigsb_4b_metbinsum  =  Rsigsb_4b_metbinsum_val / wave_Rsigsb_metbinsum_val ;
      float correction_Rsigsb_3b_metbinsum  =  Rsigsb_3b_metbinsum_val / wave_Rsigsb_metbinsum_val ;
      float correction_Rsigsb_2b_metbinsum  =  Rsigsb_2b_metbinsum_val / wave_Rsigsb_metbinsum_val ;

      float syst_Rsigsb_4b_metbinsum(0.) ;
      float syst_Rsigsb_3b_metbinsum(0.) ;
      float syst_Rsigsb_2b_metbinsum(0.) ;

      if ( closure_syst_option == 1 ) {

        //-- nominal: mcstat + 1/2 correction
         syst_Rsigsb_4b_metbinsum  =  sqrt( ::pow( (Rsigsb_4b_metbinsum_err/Rsigsb_4b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_4b_metbinsum-1.)/2., 2. ) ) ;
         syst_Rsigsb_3b_metbinsum  =  sqrt( ::pow( (Rsigsb_3b_metbinsum_err/Rsigsb_3b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_3b_metbinsum-1.)/2., 2. ) ) ;
         syst_Rsigsb_2b_metbinsum  =  sqrt( ::pow( (Rsigsb_2b_metbinsum_err/Rsigsb_2b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_2b_metbinsum-1.)/2., 2. ) ) ;

      } else if ( closure_syst_option == 2 ) {

        //-- no syst.
         syst_Rsigsb_4b_metbinsum  =  0. ;
         syst_Rsigsb_3b_metbinsum  =  0. ;
         syst_Rsigsb_2b_metbinsum  =  0. ;

      } else if ( closure_syst_option == 3 ) {

        //-- aggressive: mcstat only
         syst_Rsigsb_4b_metbinsum  =  Rsigsb_4b_metbinsum_err/Rsigsb_4b_metbinsum_val ;
         syst_Rsigsb_3b_metbinsum  =  Rsigsb_3b_metbinsum_err/Rsigsb_3b_metbinsum_val ;
         syst_Rsigsb_2b_metbinsum  =  Rsigsb_2b_metbinsum_err/Rsigsb_2b_metbinsum_val ;

      } else if ( closure_syst_option == 4 ) {

        //-- nominal: mcstat + full correction
         syst_Rsigsb_4b_metbinsum  =  sqrt( ::pow( (Rsigsb_4b_metbinsum_err/Rsigsb_4b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_4b_metbinsum-1.), 2. ) ) ;
         syst_Rsigsb_3b_metbinsum  =  sqrt( ::pow( (Rsigsb_3b_metbinsum_err/Rsigsb_3b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_3b_metbinsum-1.), 2. ) ) ;
         syst_Rsigsb_2b_metbinsum  =  sqrt( ::pow( (Rsigsb_2b_metbinsum_err/Rsigsb_2b_metbinsum_val), 2. ) + ::pow( (correction_Rsigsb_2b_metbinsum-1.), 2. ) ) ;

      } else {

         printf("\n\n\n *** Unknown value for closure_syst_option argument.  Supported values are:\n\n" ) ;
         printf("   1  MC stat error + 1/2 non-closure in quadrature.\n") ;
         printf("   2  No closure syst.\n") ;
         printf("   3  MC stat error only.\n") ;
         printf("   4  MC stat error + full non-closure in quadrature.\n") ;
         printf("\n\n\n") ;

         return ;

      }




      printf("\n") ;
      if ( bins_of_met > 1 && use3b ) {
         printf("  Rsig/sb corrections for top %d METsig bins :     4b = %5.3f +/- %5.3f,       3b =  %5.3f +/- %5.3f,      2b =  %5.3f +/- %5.3f\n",
               bins_of_met-1,
               correction_Rsigsb_4b_metbinsum, syst_Rsigsb_4b_metbinsum,
               correction_Rsigsb_3b_metbinsum, syst_Rsigsb_3b_metbinsum,
               correction_Rsigsb_2b_metbinsum, syst_Rsigsb_2b_metbinsum ) ;
      } else {
         printf("  Rsig/sb corrections for top %d METsig bins :     4b = %5.3f +/- %5.3f,       2b =  %5.3f +/- %5.3f\n",
               bins_of_met-1,
               correction_Rsigsb_4b_metbinsum, syst_Rsigsb_4b_metbinsum,
               correction_Rsigsb_2b_metbinsum, syst_Rsigsb_2b_metbinsum ) ;
      }
      printf("\n") ;



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

      //-- all observables, BG totals.
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_4b_msig_met%d   %8.2f\n", mbi, (n_4b_msig[mbi] )  ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_4b_msb_met%d    %8.2f\n", mbi, (n_4b_msb[mbi]  )  ) ; }
      if ( bins_of_met > 1 && use3b ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_3b_msig_met%d   %8.2f\n", mbi, (n_3b_msig[mbi] )  ) ; } }
      if ( bins_of_met > 1 && use3b ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_3b_msb_met%d    %8.2f\n", mbi, (n_3b_msb[mbi]  )  ) ; } }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_2b_msig_met%d   %8.2f\n", mbi, (n_2b_msig[mbi] )  ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_2b_msb_met%d    %8.2f\n", mbi, (n_2b_msb[mbi]  )  ) ; }


      //-- Rsig/sb correction factors and systematics.
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
         fprintf( outfile, "Rsigsb_corr_4b_met%d   %5.3f\n", mbi, correction_Rsigsb_4b[mbi] ) ;
         fprintf( outfile, "Rsigsb_syst_4b_met%d   %5.3f\n", mbi, syst_Rsigsb_4b[mbi] ) ;
      }
      if ( bins_of_met > 1 && use3b ) {
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile, "Rsigsb_corr_3b_met%d   %5.3f\n", mbi, correction_Rsigsb_3b[mbi] ) ;
            fprintf( outfile, "Rsigsb_syst_3b_met%d   %5.3f\n", mbi, syst_Rsigsb_3b[mbi] ) ;
         }
      }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
         fprintf( outfile, "Rsigsb_corr_2b_met%d   %5.3f\n", mbi, correction_Rsigsb_2b[mbi] ) ;
         fprintf( outfile, "Rsigsb_syst_2b_met%d   %5.3f\n", mbi, syst_Rsigsb_2b[mbi] ) ;
      }


      char metbinsumstring[100] ;
      if ( bins_of_met == 4 ) { sprintf( metbinsumstring, "234" ) ; }

      fprintf( outfile, "Rsigsb_corr_4b_metbins%s   %5.3f\n", metbinsumstring, correction_Rsigsb_4b_metbinsum ) ;
      fprintf( outfile, "Rsigsb_syst_4b_metbins%s   %5.3f\n", metbinsumstring, syst_Rsigsb_4b_metbinsum ) ;
      if ( bins_of_met > 1 && use3b ) {
         fprintf( outfile, "Rsigsb_corr_3b_metbins%s   %5.3f\n", metbinsumstring, correction_Rsigsb_3b_metbinsum ) ;
         fprintf( outfile, "Rsigsb_syst_3b_metbins%s   %5.3f\n", metbinsumstring, syst_Rsigsb_3b_metbinsum ) ;
      }
      fprintf( outfile, "Rsigsb_corr_2b_metbins%s   %5.3f\n", metbinsumstring, correction_Rsigsb_2b_metbinsum ) ;
      fprintf( outfile, "Rsigsb_syst_2b_metbins%s   %5.3f\n", metbinsumstring, syst_Rsigsb_2b_metbinsum ) ;


     //-- give location of signal counts file
      fprintf( outfile, "signal_counts_file  outputfiles/susy-signal-counts-4metbin-w3b-wpu.txt\n" ) ;

     //-- give list of shape systematic files here.
   // fprintf( outfile, "list_of_shape_systs   JES  PDF  btag\n" ) ;
   // fprintf( outfile, "shape_syst_JES    test-input-files1/test-syst-file1.txt\n") ;
   // fprintf( outfile, "shape_syst_PDF    test-input-files1/test-syst-file1.txt\n") ;
   // fprintf( outfile, "shape_syst_btag   test-input-files1/test-syst-file1.txt\n") ;
      fprintf( outfile, "list_of_shape_systs   ALL\n" ) ;
      fprintf( outfile, "shape_syst_ALL    test-input-files1/test-syst-file1.txt\n") ;




      fclose( outfile ) ;


      printf("\n\n Created likelihood input file: %s\n\n", outfilename ) ;


















    //--- Also prepare a LandS input datacard.  What the hell is a datacard?  Sounds like ancient technology.

      char outfilename_lands[10000] ;
      sprintf( outfilename_lands, "%s-lands", outfilename ) ;

      printf("\n\n #----------------------------------------------------------------------------\n\n") ;
      printf("\n\n Creating LandS datacard: %s\n\n\n", outfilename_lands ) ;

      FILE* outfile_lands = fopen( outfilename_lands, "w" ) ;

      float par_width_NSD(6.) ;

      int nchan(0) ;
      int npars(0) ;

      fprintf( outfile_lands, "---\n" ) ;
      if ( use3b ) {
         fprintf( outfile_lands, "## datacard for ABCD using 4b, 3b, and 2b in %d bins of %s\n\n\n", bins_of_met, metvarname ) ;
         nchan = 6 * bins_of_met ;
         npars = 5 * bins_of_met ;
         fprintf( outfile_lands, "imax %d  number of channels (%d bins of met for 4b, 3b, and 2b, higgs mass SIG and SB)\n", nchan, bins_of_met ) ;
         fprintf( outfile_lands, "jmax 1  number of backgrounds\n" ) ;
         fprintf( outfile_lands, "kmax %d  number of nuisance parameters (5 for each bin of met)\n", npars ) ;
      } else {
         fprintf( outfile_lands, "## datacard for ABCD using 4b and 2b in %d bins of %s\n\n\n", bins_of_met, metvarname ) ;
         nchan = 4 * bins_of_met ;
         npars = 4 * bins_of_met ;
         fprintf( outfile_lands, "imax %d  number of channels (%d bins of met for 4b and 2b, higgs mass SIG and SB)\n", nchan, bins_of_met ) ;
         fprintf( outfile_lands, "jmax 1  number of backgrounds\n" ) ;
         fprintf( outfile_lands, "kmax %d  number of nuisance parameters (4 for each bin of met)\n", npars ) ;
      }
      fprintf( outfile_lands, "---\n\n\n" ) ;

      if ( use3b ) {

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "bin                " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "               N4bsig_met%d                 N4bsb_met%d                  N3bsig_met%d                 N3bsb_met%d                  N2bsig_met%d                 N2bsb_met%d ", mbi, mbi, mbi, mbi, mbi, mbi ) ; }
         fprintf( outfile_lands, "\n" ) ;
         fprintf( outfile_lands, "Observation        " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "                    %6.1f                     %6.1f                       %6.1f                     %6.1f                       %6.1f                     %6.1f ",
             n_4b_msig[mbi] ,
             n_4b_msb[mbi]  ,
             n_3b_msig[mbi] ,
             n_3b_msb[mbi]  ,
             n_2b_msig[mbi] ,
             n_2b_msb[mbi]    ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;

      } else {

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "--------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "bin                " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "               N4bsig_met%d                 N4bsb_met%d                  N2bsig_met%d                 N2bsb_met%d ", mbi, mbi, mbi, mbi ) ; }
         fprintf( outfile_lands, "\n" ) ;
         fprintf( outfile_lands, "Observation        " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "                    %6.1f                     %6.1f                       %6.1f                     %6.1f ",
             n_4b_msig[mbi] ,
             n_4b_msb[mbi]  ,
             n_2b_msig[mbi] ,
             n_2b_msb[mbi]   ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "--------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;

      }

      if ( use3b ) {

         fprintf( outfile_lands, "bin                " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "   N4bsig_met%d N4bsig_met%d     N4bsb_met%d  N4bsb_met%d      N3bsig_met%d N3bsig_met%d     N3bsb_met%d  N3bsb_met%d      N2bsig_met%d N2bsig_met%d     N2bsb_met%d  N2bsb_met%d ", mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "process            " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "        signal        smbg         signal        smbg           signal        smbg         signal        smbg           signal        smbg         signal        smbg " ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "process            " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "             0           1              0           1                0           1              0           1                0           1              0           1 " ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "rate               " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "        %6.2f     %7.1f         %6.2f     %7.1f           %6.2f     %7.1f         %6.2f     %7.1f           %6.2f     %7.1f         %6.2f     %7.1f ",
              -1., n_4b_msig[mbi],
              -1., n_4b_msb[mbi],
              -1., n_3b_msig[mbi],
              -1., n_3b_msb[mbi],
              -1., n_2b_msig[mbi],
              -1., n_2b_msb[mbi]
              ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;


      } else {

         fprintf( outfile_lands, "bin                " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "   N4bsig_met%d N4bsig_met%d     N4bsb_met%d  N4bsb_met%d      N2bsig_met%d N2bsig_met%d     N2bsb_met%d  N2bsb_met%d ", mbi, mbi, mbi, mbi, mbi, mbi, mbi, mbi ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "process            " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "        signal        smbg         signal        smbg           signal        smbg         signal        smbg " ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "process            " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "             0           1              0           1                0           1              0           1 " ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "rate               " ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile_lands, "         %5.1f     %7.1f          %5.1f     %7.1f            %5.1f     %7.1f          %5.1f     %7.1f ",
              -1., n_4b_msig[mbi],
              -1., n_4b_msb[mbi],
              -1., n_2b_msig[mbi],
              -1., n_2b_msb[mbi]
              ) ; }
         fprintf( outfile_lands, "\n" ) ;

         fprintf( outfile_lands, "-------------------" ) ;
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "--------------------------------------------------------------------------------------------------------------" ) ;
         }
         fprintf( outfile_lands, "\n" ) ;
      }


      for ( int rmbi=1; rmbi<=bins_of_met; rmbi++ ) {

         float lands_par_error ;
         float nobs ;

         fprintf( outfile_lands, "c_m%d_sig_allnb   lnU    ", rmbi ) ;
         nobs = n_4b_msb[rmbi] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -        %4.2f              -           -                -        %4.2f              -           -                -        %4.2f              -           -      ", lands_par_error, lands_par_error, lands_par_error ) ;
               } else {
                  fprintf( outfile_lands, "        -        %4.2f              -           -                -        %4.2f              -           -      ", lands_par_error, lands_par_error ) ;
               }
            } else {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -      " ) ;
               } else {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -      " ) ;
               }
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;



         fprintf( outfile_lands, "c_m%d_R_4b        lnU    ", rmbi ) ;
         nobs = n_4b_msb[rmbi] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -        %4.2f              -        %4.2f                -           -              -           -                -           -              -           -      ", lands_par_error, lands_par_error ) ;
               } else {
                  fprintf( outfile_lands, "        -        %4.2f              -        %4.2f                -           -              -           -      ", lands_par_error, lands_par_error ) ;
               }
            } else {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -      " ) ;
               } else {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -      " ) ;
               }
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;




         if ( use3b ) {
            fprintf( outfile_lands, "c_m%d_R_3b        lnU    ", rmbi ) ;
            nobs = n_3b_msb[rmbi] ;
            lands_par_error = 1. + par_width_NSD ;
            if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
            for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
               if ( cmbi == rmbi ) {
                  fprintf( outfile_lands, "        -           -              -           -                -        %4.2f              -        %4.2f                -           -              -           -      ", lands_par_error, lands_par_error ) ;
               } else {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -      " ) ;
               }
            } // cmbi.
            fprintf( outfile_lands, "\n" ) ;
         }



         fprintf( outfile_lands, "c_m%d_R_2b        lnU    ", rmbi ) ;
         nobs = n_2b_msb[rmbi] ;
         lands_par_error = 1. + par_width_NSD ;
         if ( nobs > 0 ) { lands_par_error = 1. + par_width_NSD / sqrt( nobs ) ; }
         for ( int cmbi=1; cmbi<=bins_of_met; cmbi++ ) {
            if ( cmbi == rmbi ) {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -        %4.2f              -        %4.2f      ", lands_par_error, lands_par_error ) ;
               } else {
                  fprintf( outfile_lands, "        -           -              -           -                -        %4.2f              -        %4.2f      ", lands_par_error, lands_par_error ) ;
               }
            } else {
               if ( use3b ) {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -                -           -              -           -      " ) ;
               } else {
                  fprintf( outfile_lands, "        -           -              -           -                -           -              -           -      " ) ;
               }
            }
         } // cmbi.
         fprintf( outfile_lands, "\n" ) ;


      } // rmbi.

      char format1[10000] ;
      char format2[10000] ;
      if ( use3b ) {
         sprintf( format1, "        -           -              -        %%4.2f                -           -              -           -                -           -              -           -      " ) ;
         sprintf( format2, "        -           -              -           -                -           -              -           -                -           -              -           -      " ) ;
      } else {
         sprintf( format1, "        -           -              -        %%4.2f                -           -              -           -      " ) ;
         sprintf( format2, "        -           -              -           -                -           -              -           -      " ) ;
      }

      fprintf( outfile_lands, "closure_m1       lnN    " ) ;
         fprintf( outfile_lands, format1, (1. + syst_Rsigsb_4b[1]) ) ;
      for ( int mbi=2; mbi<=bins_of_met; mbi++ ) {
         fprintf( outfile_lands, "%s", format2 ) ;
      }
      fprintf( outfile_lands, "\n" ) ;

      if ( bins_of_met >= 2 ) {
         fprintf( outfile_lands, "closure_m2       lnN    " ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, format1, (1. + syst_Rsigsb_4b[2]) ) ;
         for ( int mbi=3; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "%s", format2 ) ;
         }
         fprintf( outfile_lands, "\n" ) ;
      }

      if ( bins_of_met >= 3 ) {
         fprintf( outfile_lands, "closure_m3       lnN    " ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, format1, (1. + syst_Rsigsb_4b[3]) ) ;
         for ( int mbi=4; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile_lands, "%s", format2 ) ;
         }
         fprintf( outfile_lands, "\n" ) ;
      }

      if ( bins_of_met >= 4 ) {
         fprintf( outfile_lands, "closure_m2       lnN    " ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, "%s", format2 ) ;
         fprintf( outfile_lands, format1, (1. + syst_Rsigsb_4b[4]) ) ;
         fprintf( outfile_lands, "\n" ) ;
      }



      fclose( outfile_lands ) ;

      printf("\n\n Created LandS datacard file: %s\n\n", outfilename_lands ) ;





     //--- save histograms.

      printf("\n\n Saving histograms to outputfiles/gen_input.root\n\n") ;

      saveHist( "outputfiles/gen_input.root", "h*" ) ;


   } // gen_input_file2







