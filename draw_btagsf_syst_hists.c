
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TStyle.h"

#include "histio.c"

   void draw_btagsf_syst_hists( const char* inwsfile ) {

      gStyle -> SetOptStat(0) ;

      gDirectory -> Delete( "h*" ) ;

      loadHist( inwsfile ) ;

     //----

      TCanvas* can_sb = (TCanvas*) gDirectory -> FindObject( "can_sb" ) ;
      if ( can_sb == 0x0 ) {
         can_sb = new TCanvas( "can_sb", "SB observables", 1100, 300 ) ;
      }

      can_sb -> Clear() ;
      can_sb -> Divide(4,1) ;

      for ( int ci=1; ci<=4; ci++ ) {

         can_sb -> cd( ci ) ;

         char hname[1000] ;

         sprintf( hname, "h_btagsfsyst_msb_met%d_nom", ci ) ;
         TH1F* hist_nom = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_nom == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         sprintf( hname, "h_btagsfsyst_msb_met%d_m1s", ci ) ;
         TH1F* hist_m1s = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_m1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         sprintf( hname, "h_btagsfsyst_msb_met%d_p1s", ci ) ;
         TH1F* hist_p1s = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_p1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         hist_nom -> SetLineWidth( 2 ) ;
         hist_m1s -> SetLineWidth( 2 ) ;
         hist_p1s -> SetLineWidth( 2 ) ;

         hist_m1s -> SetLineColor( 4 ) ;
         hist_p1s -> SetLineColor( 2 ) ;

         float hmax(0.) ;
         if ( 1.2*(hist_nom->GetMaximum()) > hmax ) { hmax = 1.2*(hist_nom->GetMaximum()) ; }
         if ( 1.2*(hist_m1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_m1s->GetMaximum()) ; }
         if ( 1.2*(hist_p1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_p1s->GetMaximum()) ; }

         hist_nom -> SetMaximum( hmax ) ;
         hist_nom -> Draw() ;
         hist_m1s -> Draw("same" ) ;
         hist_p1s -> Draw("same" ) ;
         hist_nom -> Draw("same" ) ;

      } // ci


     //----

      TCanvas* can_sig = (TCanvas*) gDirectory -> FindObject( "can_sig" ) ;
      if ( can_sig == 0x0 ) {
         can_sig = new TCanvas( "can_sig", "SIG observables", 1100, 300 ) ;
      }

      can_sig -> Clear() ;
      can_sig -> Divide(4,1) ;

      for ( int ci=1; ci<=4; ci++ ) {

         can_sig -> cd( ci ) ;

         char hname[1000] ;

         sprintf( hname, "h_btagsfsyst_msig_met%d_nom", ci ) ;
         TH1F* hist_nom = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_nom == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         sprintf( hname, "h_btagsfsyst_msig_met%d_m1s", ci ) ;
         TH1F* hist_m1s = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_m1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         sprintf( hname, "h_btagsfsyst_msig_met%d_p1s", ci ) ;
         TH1F* hist_p1s = (TH1F*) gDirectory -> FindObject( hname ) ;
         if ( hist_p1s == 0x0 ) { printf("\n\n *** Can't find %s\n\n", hname ) ; return ; }

         hist_nom -> SetLineWidth( 2 ) ;
         hist_m1s -> SetLineWidth( 2 ) ;
         hist_p1s -> SetLineWidth( 2 ) ;

         hist_m1s -> SetLineColor( 4 ) ;
         hist_p1s -> SetLineColor( 2 ) ;

         float hmax(0.) ;
         if ( 1.2*(hist_nom->GetMaximum()) > hmax ) { hmax = 1.2*(hist_nom->GetMaximum()) ; }
         if ( 1.2*(hist_m1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_m1s->GetMaximum()) ; }
         if ( 1.2*(hist_p1s->GetMaximum()) > hmax ) { hmax = 1.2*(hist_p1s->GetMaximum()) ; }

         hist_nom -> SetMaximum( hmax ) ;
         hist_nom -> Draw() ;
         hist_m1s -> Draw("same" ) ;
         hist_p1s -> Draw("same" ) ;
         hist_nom -> Draw("same" ) ;

      } // ci



   } // draw_btagsf_syst_hists



