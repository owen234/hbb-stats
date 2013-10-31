
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "histio.c"

   void draw_syst_hists( const char* inwsfile, const char* syst_name, int sig_mass = 0, float max_msig_evts = 8., float max_msb_evts = 16. ) {

      gStyle -> SetOptStat(0) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetTitleW( 0.9 ) ;


      gDirectory -> Delete( "h*" ) ;

      loadHist( inwsfile ) ;

      TLine* line = new TLine() ;
      line -> SetLineStyle(2) ;

     //----

      gSystem -> Exec( "mkdir -p outputfiles/syst-plots" ) ;

      TString infile_ts( inwsfile ) ;
      TObjArray* tokens = infile_ts.Tokenize("/") ;

      TObjString* tos = (TObjString*) (tokens -> At( tokens->GetEntries() - 1 ) ) ;
      TString fname( tos->GetString() ) ;
      printf( "ws file name : %s\n", fname.Data() ) ;

      TString pdfbasename = fname.ReplaceAll(".root","") ;
      printf(" pdf base name : %s\n", pdfbasename.Data() ) ;


      char sigsb_str[2][10] = { "msig", "msb" } ;

      for ( int ssbi=0; ssbi<2; ssbi++ ) {

         char cname[100] ;

         sprintf( cname, "can_evts_%s", sigsb_str[ssbi] ) ;
         TCanvas* can_evts = (TCanvas*) gDirectory -> FindObject( cname ) ;
         if ( can_evts == 0x0 ) {
            if ( ssbi == 0 ) {
               can_evts = new TCanvas( cname, "SIG observables", 1100, 300 ) ;
            } else {
               can_evts = new TCanvas( cname, "SB observables", 1100, 300 ) ;
            }
         }

         sprintf( cname, "can_frac_%s", sigsb_str[ssbi] ) ;
         TCanvas* can_frac = (TCanvas*) gDirectory -> FindObject( cname ) ;
         if ( can_frac == 0x0 ) {
            if ( ssbi == 0 ) {
               can_frac = new TCanvas( cname, "SIG observables, syst (%)", 1100, 300 ) ;
            } else {
               can_frac = new TCanvas( cname, "SB observables, syst (%)", 1100, 300 ) ;
            }
         }


         can_evts -> Clear() ;
         can_evts -> Divide(4,1) ;

         can_frac -> Clear() ;
         can_frac -> Divide(4,1) ;

         for ( int ci=1; ci<=4; ci++ ) {


            char hname[1000] ;
            char hnamev[1000] ;
            char hnamenf[1000] ;

            sprintf( hname, "h_syst_%s_%s_met%d_nom", syst_name, sigsb_str[ssbi], ci ) ;
            TH1F* hist_nom = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( hist_nom == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }
            sprintf( hnamenf, "%s_nf", hname ) ;
            TH1F* hist_nom_nf = (TH1F*) hist_nom -> Clone( hnamenf ) ;

            sprintf( hname, "h_syst_%s_%s_met%d_m1s", syst_name, sigsb_str[ssbi], ci ) ;
            TH1F* hist_m1s = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( hist_m1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }
            sprintf( hnamev, "%s_var", hname ) ;
            TH1F* hist_m1s_var = (TH1F*) hist_m1s -> Clone( hnamev ) ;

            sprintf( hname, "h_syst_%s_%s_met%d_p1s", syst_name, sigsb_str[ssbi], ci ) ;
            TH1F* hist_p1s = (TH1F*) gDirectory -> FindObject( hname ) ;
            if ( hist_p1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }
            sprintf( hnamev, "%s_var", hname ) ;
            TH1F* hist_p1s_var = (TH1F*) hist_p1s -> Clone( hnamev ) ;

            TString htitle ;

            htitle = hist_nom -> GetTitle() ;
            if ( sig_mass > 0 ) {
               char sigmassstr[1000] ;
               sprintf( sigmassstr, ", higgsino mass = %d", sig_mass ) ;
               htitle.ReplaceAll( ", nominal", sigmassstr ) ;
            } else {
               htitle.ReplaceAll( ", nominal", "" ) ;
            }
            hist_nom -> SetTitle( htitle ) ;

            htitle = hist_m1s_var -> GetTitle() ;
            if ( sig_mass > 0 ) {
               char sigmassstr[1000] ;
               sprintf( sigmassstr, ", higgsino mass = %d", sig_mass ) ;
               htitle.ReplaceAll( ", -1 sigma", sigmassstr ) ;
            } else {
               htitle.ReplaceAll( ", -1 sigma", "" ) ;
            }
            hist_m1s_var -> SetTitle( htitle ) ;



            hist_nom -> SetLineWidth( 2 ) ;
            hist_nom_nf -> SetLineWidth( 2 ) ;
            hist_m1s -> SetLineWidth( 2 ) ;
            hist_p1s -> SetLineWidth( 2 ) ;

            hist_m1s -> SetLineColor( 4 ) ;
            hist_p1s -> SetLineColor( 2 ) ;

            hist_nom -> SetFillColor( 18 ) ;


            float max_evts(0.) ;
            if ( ssbi == 0 ) {
               max_evts = max_msig_evts ;
            } else {
               max_evts = max_msb_evts ;
            }

            float hmax(0.) ;
            if ( max_evts < 0 ) {
               if ( 1.2*(hist_nom->GetMaximum()) > hmax ) { hmax = 1.2*(hist_nom->GetMaximum()) ; }
               if ( 1.2*(hist_m1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_m1s->GetMaximum()) ; }
               if ( 1.2*(hist_p1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_p1s->GetMaximum()) ; }
            } else {
               hmax = max_evts ;
            }

            hist_nom -> SetMaximum( hmax ) ;

            for ( int hbi=1; hbi <= hist_nom -> GetNbinsX(); hbi++ ) {
               float nom_val, p1s_val, m1s_val ;
               nom_val = hist_nom -> GetBinContent( hbi ) ;
               p1s_val = hist_p1s -> GetBinContent( hbi ) ;
               m1s_val = hist_m1s -> GetBinContent( hbi ) ;
               hist_p1s_var -> SetBinContent( hbi, 0. ) ;
               hist_m1s_var -> SetBinContent( hbi, 0. ) ;
               if ( nom_val > 0 ) {
                  hist_p1s_var -> SetBinContent( hbi, (p1s_val - nom_val)/nom_val ) ;
                  hist_m1s_var -> SetBinContent( hbi, (m1s_val - nom_val)/nom_val ) ;
                  printf( " hbi=%d : p1s, nom, m1s : %.2f %.2f %.2f\n", hbi, p1s_val, nom_val, m1s_val ) ;
               }
            } // hbi

            hist_m1s_var -> SetMinimum( -0.3 ) ;
            hist_m1s_var -> SetMaximum(  0.3 ) ;

            hist_m1s_var -> SetLineColor( 4 ) ;
            hist_p1s_var -> SetLineColor( 2 ) ;
            hist_m1s_var -> SetLineWidth( 2 ) ;
            hist_p1s_var -> SetLineWidth( 2 ) ;

            hist_m1s_var -> SetFillColor( 4 ) ;
            hist_m1s_var -> SetFillStyle( 3354 ) ;

            hist_p1s_var -> SetFillColor( 2 ) ;
            hist_p1s_var -> SetFillStyle( 3345 ) ;

            hist_nom     -> SetTitleOffset( 1.5, "y" ) ;
            hist_m1s_var -> SetTitleOffset( 1.5, "y" ) ;
            hist_nom     -> SetTitleSize( 0.05, "y" ) ;
            hist_m1s_var -> SetTitleSize( 0.05, "y" ) ;
            hist_nom     -> SetLabelSize( 0.07, "x" ) ;
            hist_m1s_var -> SetLabelSize( 0.07, "x" ) ;
            hist_nom     -> SetLabelSize( 0.05, "y" ) ;
            hist_m1s_var -> SetLabelSize( 0.05, "y" ) ;
            hist_nom     -> SetLabelOffset( 0.01, "x" ) ;
            hist_m1s_var -> SetLabelOffset( 0.01, "x" ) ;
            hist_nom     -> SetLabelOffset( 0.01, "y" ) ;
            hist_m1s_var -> SetLabelOffset( 0.01, "y" ) ;

            hist_nom -> SetYTitle( "Events at theory Xsec" ) ;
            hist_m1s_var -> SetYTitle( "Systematic (var-nom)/nom" ) ;


            can_evts -> cd( ci ) ;
            hist_nom -> Draw() ;
            hist_m1s -> Draw("same" ) ;
            hist_p1s -> Draw("same" ) ;
            hist_nom_nf -> Draw("same" ) ;
            hist_nom -> Draw("same axis" ) ;

            can_frac -> cd( ci ) ;
            hist_m1s_var -> Draw( ) ;
            hist_p1s_var -> Draw("same" ) ;
            line -> DrawLine( hist_m1s_var->GetBinLowEdge(1), 0., hist_m1s_var->GetBinLowEdge( hist_nom -> GetNbinsX() + 1 ), 0. ) ;



         } // ci

         char pdfname[10000] ;

         sprintf( pdfname, "outputfiles/syst-plots/%s-syst-%s-events-%s.pdf", pdfbasename.Data(), syst_name, sigsb_str[ssbi] ) ;
         can_evts -> SaveAs( pdfname ) ;
         sprintf( pdfname, "outputfiles/syst-plots/%s-syst-%s-frac-%s.pdf", pdfbasename.Data(), syst_name, sigsb_str[ssbi] ) ;
         can_frac -> SaveAs( pdfname ) ;

      } // ssbi.





   } // draw_syst_hists



