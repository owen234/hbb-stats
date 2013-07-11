
#include "TChain.h"
#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TText.h"


   void plottoy1( const char* infile = "outputfiles/toytt-withsig-250-ss1.0.root", bool doUL = false, const char* plottitle="" ) {

      gStyle -> SetOptStat("emro") ;
      gStyle -> SetStatW(0.30) ;
      gStyle -> SetStatH(0.20) ;
      gStyle -> SetTitleH(0.06) ;
      gStyle -> SetLabelSize( 0.05, "y" ) ;
      gStyle -> SetLabelSize( 0.05, "x" ) ;
      gStyle -> SetLabelOffset( 0.010, "y" ) ;
      gStyle -> SetLabelOffset( 0.010, "x" ) ;

      TCanvas* ctoy = (TCanvas*) gDirectory->FindObject("ctoy") ;
      if ( ctoy == 0x0 ) {
         ctoy = new TCanvas("ctoy", "toy results", 1000, 750 ) ;
      }
      ctoy -> Clear() ;
      ctoy -> Divide(2,2) ;

      TChain* toytt = new TChain("toytt") ;
      toytt -> Add( infile ) ;

      TText* tt_title = new TText() ;


      gDirectory -> Delete( "h*" ) ;

      TH1F* h_str_fit    = new TH1F("h_str_fit"   , "signal strength"             , 100, 0., 3. ) ;
      TH1F* h_str_err    = new TH1F("h_str_err"   , "signal strength stat err."   , 100, 0., 1. ) ;
      TH1F* h_str_pull   = new TH1F("h_str_pull"  , "signal strength pull"        , 100, -10., 10. ) ;
      TH1F* h_str_signif = new TH1F("h_str_signif", "signal strength significance", 100, 0., 10. ) ;
      TH1F* h_str_ul     = new TH1F("h_str_ul"    , "signal strength upper limit" , 100, 0., 2. ) ;

      h_str_fit    -> SetFillColor(11) ;
      h_str_err    -> SetFillColor(11) ;
      h_str_pull   -> SetFillColor(11) ;
      h_str_signif -> SetFillColor(11) ;
      h_str_ul     -> SetFillColor(11) ;

      ctoy -> cd(1) ;
      toytt -> Draw( "fit_sig_strength>>h_str_fit", "" ) ;
      tt_title->DrawTextNDC(0.15,0.83, plottitle ) ;

      ctoy -> cd(2) ;
      toytt -> Draw( "fit_sig_strength_err>>h_str_err", "" ) ;
      tt_title->DrawTextNDC(0.15,0.83, plottitle ) ;

      if ( !doUL ) {

         ctoy -> cd(3) ;
         toytt -> Draw( "(fit_sig_strength-1.)/fit_sig_strength_err>>h_str_pull", "" ) ;
         tt_title->DrawTextNDC(0.15,0.83, plottitle ) ;

         ctoy -> cd(4) ;
         toytt -> Draw( "fit_sig_signif>>h_str_signif", "" ) ;
         tt_title->DrawTextNDC(0.15,0.83, plottitle ) ;

      } else {

         ctoy -> cd(3) ;
         toytt -> Draw( "fit_sig_ul>>h_str_ul", "" ) ;
         tt_title->DrawTextNDC(0.15,0.83, plottitle ) ;

      }

      TString pdffile( infile ) ;
      pdffile.ReplaceAll( "toytt-", "toyplots-" ) ;
      pdffile.ReplaceAll( "root", "pdf" ) ;

      ctoy -> SaveAs( pdffile ) ;

   }

