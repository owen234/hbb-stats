
#include "TMath.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"


#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"

#include "histio.c"

#include "scan_sigstrength.c"


  using namespace RooFit;
  using namespace RooStats;


   void labelBins( TH1F* hist ) ;

   void fitqual_plots( const char* wsfile = "outputfiles/ws.root", const char* plottitle="" ) {

      TText* tt_title = new TText() ;
      tt_title -> SetTextAlign(33) ;

      gStyle -> SetOptStat(0) ;
      gStyle -> SetLabelSize( 0.06, "y" ) ;
      gStyle -> SetLabelSize( 0.08, "x" ) ;
      gStyle -> SetLabelOffset( 0.010, "y" ) ;
      gStyle -> SetLabelOffset( 0.010, "x" ) ;
      gStyle -> SetTitleSize( 0.07, "y" ) ;
      gStyle -> SetTitleSize( 0.05, "x" ) ;
      gStyle -> SetTitleOffset( 1.50, "x" ) ;
      gStyle -> SetTitleH( 0.07 ) ;
      gStyle -> SetPadLeftMargin( 0.15 ) ;
      gStyle -> SetPadBottomMargin( 0.15 ) ;
      gStyle -> SetTitleX( 0.10 ) ;

      gDirectory->Delete("h*") ;

      TFile* wstf = new TFile( wsfile ) ;

      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      ws->Print() ;

      int bins_of_met = TMath::Nint( ws->var("bins_of_met")->getVal()  ) ;
      printf("\n\n Bins of MET : %d\n\n", bins_of_met ) ;

      int bins_of_nb = TMath::Nint( ws->var("bins_of_nb")->getVal()  ) ;
      printf("\n\n Bins of nb : %d\n\n", bins_of_nb ) ;

      int nb_lookup[10] ;
      if ( bins_of_nb == 2 ) {
         nb_lookup[0] = 2 ;
         nb_lookup[1] = 4 ;
      } else if ( bins_of_nb == 3 ) {
         nb_lookup[0] = 2 ;
         nb_lookup[1] = 3 ;
         nb_lookup[2] = 4 ;
      }

      TCanvas* cfq1 = (TCanvas*) gDirectory->FindObject("cfq1") ;
      if ( cfq1 == 0x0 ) {
         if ( bins_of_nb == 3 ) {
            cfq1 = new TCanvas("cfq1","hbb fit", 700, 1000 ) ;
         } else if ( bins_of_nb == 2 ) {
            cfq1 = new TCanvas("cfq1","hbb fit", 700, 750 ) ;
         } else {
            return ;
         }
      }

      RooRealVar* rv_sig_strength = ws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }

      ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

      RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;

      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooAbsPdf* likelihood = modelConfig->GetPdf() ;

      ///RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(3) ) ;
      fitResult->Print() ;


      char hname[1000] ;
      char htitle[1000] ;
      char pname[1000] ;




     //-- unpack observables.

      int obs_N_msig[10][50] ; // first index is n btags, second is met bin.
      int obs_N_msb[10][50]  ; // first index is n btags, second is met bin.

      const RooArgSet* dsras = rds->get() ;
      TIterator* obsIter = dsras->createIterator() ;
      while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
         for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
            for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
               sprintf( pname, "N_%db_msig_met%d", nb_lookup[nbi], mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msig[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
               sprintf( pname, "N_%db_msb_met%d", nb_lookup[nbi], mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msb[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
            } // mbi.
         } // nbi.
      } // obs iterator.


      printf("\n\n") ;
      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         printf(" nb=%d :  ", nb_lookup[nbi] ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
            printf("  sig=%3d, sb=%3d  |", obs_N_msig[nbi][mbi], obs_N_msb[nbi][mbi] ) ;
         } // mbi.
         printf("\n") ;
      } // nbi.
      printf("\n\n") ;




      int pad(1) ;

      cfq1->Clear() ;
      cfq1->Divide( 2, bins_of_nb+1 ) ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {


         sprintf( hname, "h_bg_%db_msig_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_bg_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_bg_msig -> SetFillColor( kBlue-9 ) ;
         labelBins( hist_bg_msig ) ;

         sprintf( hname, "h_bg_%db_msb_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sb, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_bg_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_bg_msb -> SetFillColor( kBlue-9 ) ;
         labelBins( hist_bg_msb ) ;

         sprintf( hname, "h_sig_%db_msig_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_sig_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_sig_msig -> SetFillColor( kMagenta+2 ) ;
         labelBins( hist_sig_msig ) ;

         sprintf( hname, "h_sig_%db_msb_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sb, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_sig_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_sig_msb -> SetFillColor( kMagenta+2 ) ;
         labelBins( hist_sig_msb ) ;

         sprintf( hname, "h_all_%db_msig_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_all_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_all_%db_msb_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sb, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_all_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_data_%db_msig_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_data_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_data_msig -> SetLineWidth(2) ;
         hist_data_msig -> SetMarkerStyle(20) ;
         labelBins( hist_data_msig ) ;

         sprintf( hname, "h_data_%db_msb_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sb, %db, MET", nb_lookup[nbi] ) ;
         TH1F* hist_data_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_data_msb -> SetLineWidth(2) ;
         hist_data_msb -> SetMarkerStyle(20) ;
         labelBins( hist_data_msb ) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {



            sprintf( pname, "mu_bg_%db_msig_met%d", nb_lookup[nbi], mbi+1 ) ;
            RooAbsReal* mu_bg_msig = ws->function( pname ) ;
            if ( mu_bg_msig == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_bg_msig -> SetBinContent( mbi+1, mu_bg_msig->getVal() ) ;

            sprintf( pname, "mu_sig_%db_msig_met%d", nb_lookup[nbi], mbi+1 ) ;
            RooAbsReal* mu_sig_msig = ws->function( pname ) ;
            if ( mu_sig_msig == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_sig_msig -> SetBinContent( mbi+1, mu_sig_msig->getVal() ) ;

            hist_all_msig -> SetBinContent( mbi+1, mu_bg_msig->getVal() + mu_sig_msig->getVal() ) ;

            hist_data_msig -> SetBinContent( mbi+1, obs_N_msig[nbi][mbi] ) ;



            sprintf( pname, "mu_bg_%db_msb_met%d", nb_lookup[nbi], mbi+1 ) ;
            RooAbsReal* mu_bg_msb = ws->function( pname ) ;
            if ( mu_bg_msb == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_bg_msb -> SetBinContent( mbi+1, mu_bg_msb->getVal() ) ;

            sprintf( pname, "mu_sig_%db_msb_met%d", nb_lookup[nbi], mbi+1 ) ;
            RooAbsReal* mu_sig_msb = ws->function( pname ) ;
            if ( mu_sig_msb == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_sig_msb -> SetBinContent( mbi+1, mu_sig_msb->getVal() ) ;

            hist_all_msb -> SetBinContent( mbi+1, mu_bg_msb->getVal() + mu_sig_msb->getVal() ) ;

            hist_data_msb -> SetBinContent( mbi+1, obs_N_msb[nbi][mbi] ) ;



         } // mbi.

         cfq1->cd( pad ) ;

         sprintf( hname, "h_stack_%db_msig_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         THStack* hstack_msig = new THStack( hname, htitle ) ;
         hstack_msig -> Add( hist_bg_msig ) ;
         hstack_msig -> Add( hist_sig_msig ) ;

         hist_data_msig -> Draw("e") ;
         hstack_msig -> Draw("same") ;
         hist_data_msig -> Draw("same e") ;
         hist_data_msig -> Draw("same axis") ;

         tt_title -> DrawTextNDC( 0.85, 0.85, plottitle ) ;

         pad++ ;



         cfq1->cd( pad ) ;

         sprintf( hname, "h_stack_%db_msb_met", nb_lookup[nbi] ) ;
         sprintf( htitle, "mass sig, %db, MET", nb_lookup[nbi] ) ;
         THStack* hstack_msb = new THStack( hname, htitle ) ;
         hstack_msb -> Add( hist_bg_msb ) ;
         hstack_msb -> Add( hist_sig_msb ) ;

         hist_data_msb -> Draw("e") ;
         hstack_msb -> Draw("same") ;
         hist_data_msb -> Draw("same e") ;
         hist_data_msb -> Draw("same axis") ;

         tt_title -> DrawTextNDC( 0.85, 0.85, plottitle ) ;

         pad++ ;



      } // nbi.




      TH1F* hist_R_msigmsb = new TH1F( "h_R_msigmsb", "R msig/msb vs met bin", bins_of_met, 0.5, 0.5+bins_of_met ) ;
      hist_R_msigmsb -> SetLineWidth(2) ;
      hist_R_msigmsb -> SetMarkerStyle(20) ;
      hist_R_msigmsb -> SetYTitle("R msig/msb") ;
      labelBins( hist_R_msigmsb ) ;


      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         sprintf( pname, "R_msigmsb_met%d", mbi+1 ) ;
         RooRealVar* rrv_R = ws->var( pname ) ;
         if ( rrv_R == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
         hist_R_msigmsb -> SetBinContent( mbi+1, rrv_R -> getVal() ) ;
         hist_R_msigmsb -> SetBinError( mbi+1, rrv_R -> getError() ) ;
      } // mbi.

      cfq1->cd( pad ) ;

      gPad->SetGridy(1) ;

      hist_R_msigmsb -> SetMaximum(0.35) ;
      hist_R_msigmsb -> Draw("e") ;

      tt_title -> DrawTextNDC( 0.85, 0.85, plottitle ) ;

      pad++ ;



      cfq1->cd( pad ) ;

      scan_sigstrength( wsfile ) ;

      tt_title -> DrawTextNDC( 0.85, 0.25, plottitle ) ;



      TString pdffile( wsfile ) ;
      pdffile.ReplaceAll("ws-","fitqual-") ;
      pdffile.ReplaceAll("root","pdf") ;


      cfq1->SaveAs( pdffile ) ;



      TString histfile( wsfile ) ;
      histfile.ReplaceAll("ws-","fitqual-") ;

      saveHist( histfile, "h*" ) ;



   } // fitqual_plots

 //===========================================================================================


   void labelBins( TH1F* hist ) {
      TAxis* xaxis = hist->GetXaxis() ;
      for ( int mbi=0; mbi<hist->GetNbinsX(); mbi++ ) {
         char label[1000] ;
         sprintf( label, "MET%d", mbi+1 ) ;
         xaxis->SetBinLabel( mbi+1, label ) ;
      } // mbi.
   }


 //===========================================================================================













