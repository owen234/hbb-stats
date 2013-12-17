
#include "TROOT.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TPad.h"
#include "TObjArray.h"
#include "TList.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TString.h"


   void draw_2dscan( const char* rootfile = "outputfiles/scans-ws-data-unblind/scan-hb-mu_bg_3b_msig_met3-vs-mu_bg_4b_msig_met3.root",
                     const char* gr_name = "scan_mu_bg_3b_msig_met3_vs_mu_bg_4b_msig_met3" ) {

      TCanvas* c_scan2d = (TCanvas*) gDirectory -> FindObject("c_scan2d") ;
      if ( c_scan2d == 0x0 ) {
         c_scan2d = new TCanvas( "c_scan2d", "2d scan", 700, 700 ) ;
      }


      TFile* rf = new TFile( rootfile, "read" ) ;
      TGraph2D* gr2d = (TGraph2D*) rf->Get( gr_name ) ;

      if ( gr2d == 0x0 ) {
         printf("\n\n *** can't find %s in %s.\n\n", gr_name, rootfile ) ;
         return ;
      }

      int np = gr2d -> GetN() ;
      double* xval = gr2d -> GetX() ;
      double* yval = gr2d -> GetY() ;
      double* zval = gr2d -> GetZ() ;
      double minz(1e9) ;
      int minz_pi(-1) ;
      for ( int pi=0; pi<np; pi++ ) {
         printf( "%3d : %.4f, %.4f, %.4f\n", pi, xval[pi], yval[pi], zval[pi] ) ;
         if ( zval[pi] < minz ) {
            minz = zval[pi] ;
            minz_pi = pi ;
         }
      } // pi
      if ( minz_pi < 0 ) { printf("error\n\n") ; return ; }
      double bestx = xval[minz_pi] ;
      double besty = yval[minz_pi] ;
      printf("\n\n Minimum at x=%.3f, y=%.3f\n\n", bestx, besty ) ;

      TGraph* gr_bestpoint = new TGraph(1,&bestx,&besty) ;
      gr_bestpoint -> SetMarkerStyle(20) ;
      gr_bestpoint -> SetMarkerSize(1.5) ;







      TH2* h2 = gr2d -> GetHistogram() ;
      h2 -> SetContour(3) ;
      h2 -> SetContourLevel(0,0.5) ;
      h2 -> SetContourLevel(1,1.0) ;
      h2 -> SetContourLevel(2,2.0) ;

      h2 -> Draw( "cont list" ) ;

      gPad -> Update() ;
      TObjArray *contours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours") ;
      TList *list         = (TList*)contours->At(1);
      //-- glitch in the contours near the edges.  Isn't always first in list.
      TGraph *gr1 = (TGraph*)list->First(); ;
      ////TGraph *gr1 = (TGraph*)list->At(1); ;
      if ( gr1 == 0x0 ) { printf("\n\n *** missing contour.\n\n" ) ; return ; }
      gr1->SetLineWidth(2) ;
      gr1->SetLineColor(1) ;





      h2 -> Draw( "cont1" ) ;
      gr1->Draw("c") ;
      gr_bestpoint -> Draw("p") ;


      int gr1_np = gr1 -> GetN() ;
      double* gr1_x = gr1->GetX() ;
      double* gr1_y = gr1->GetY() ;
      double gr1_xmax(-1e9) ;
      double gr1_y_at_xmax(0.) ;
      double gr1_xmin(1e9) ;
      double gr1_y_at_xmin(0.) ;
      double gr1_ymax(-1e9) ;
      double gr1_x_at_ymax(0.) ;
      double gr1_ymin(1e9) ;
      double gr1_x_at_ymin(0.) ;
      for ( int pi=0; pi<gr1_np; pi++ ) {
         if ( gr1_x[pi] > gr1_xmax ) { gr1_xmax = gr1_x[pi] ; gr1_y_at_xmax = gr1_y[pi] ; }
         if ( gr1_y[pi] > gr1_ymax ) { gr1_ymax = gr1_y[pi] ; gr1_x_at_ymax = gr1_x[pi] ; }
         if ( gr1_x[pi] < gr1_xmin ) { gr1_xmin = gr1_x[pi] ; gr1_y_at_xmin = gr1_y[pi] ; }
         if ( gr1_y[pi] < gr1_ymin ) { gr1_ymin = gr1_y[pi] ; gr1_x_at_ymin = gr1_x[pi] ; }
      } // pi

      double x_m1sigma = bestx - gr1_xmin ;
      double x_p1sigma = gr1_xmax - bestx ;
      double y_m1sigma = besty - gr1_ymin ;
      double y_p1sigma = gr1_ymax - besty ;
      printf(" x errors : +%.3f, -%.3f\n", x_p1sigma, x_m1sigma ) ;
      printf(" y errors : +%.3f, -%.3f\n", y_p1sigma, y_m1sigma ) ;
      printf(" x at ymin, ybest, ymax : %.3f, %.3f, %.3f\n", gr1_x_at_ymin, bestx, gr1_x_at_ymax ) ;
      printf(" y at xmin, xbest, xmax : %.3f, %.3f, %.3f\n", gr1_y_at_xmin, besty, gr1_y_at_xmax ) ;
      double rho_from_ymax = (gr1_x_at_ymax - bestx) / x_p1sigma ;
      double rho_from_ymin = (bestx - gr1_x_at_ymin) / x_m1sigma ;
      double rho_from_xmax = (gr1_y_at_xmax - besty) / y_p1sigma ;
      double rho_from_xmin = (besty - gr1_y_at_xmin) / y_m1sigma ;
      double ave_rho = 0.25 * (rho_from_ymax + rho_from_ymin + rho_from_xmax + rho_from_xmin ) ;
      printf(" rho from point at ymax, ymin, xmax, xmin (ave):  %.3f, %.3f, %.3f, %.3f  (%.3f)\n\n", rho_from_ymax, rho_from_ymin, rho_from_xmax, rho_from_xmin, ave_rho ) ;

      TLine* line = new TLine() ;
      line -> SetLineColor(2) ;
      line -> DrawLine( bestx, gr1_ymin, bestx, gr1_ymax ) ;
      line -> DrawLine( gr1_xmin, besty, gr1_xmax, besty ) ;

      double endbarfrac = 0.05 ;
      double endbarsizex = endbarfrac * ( gr2d -> GetXmax() - gr2d -> GetXmin() ) ;
      double endbarsizey = endbarfrac * ( gr2d -> GetYmax() - gr2d -> GetYmin() ) ;

      line -> DrawLine( bestx - endbarsizex, gr1_ymin, bestx + endbarsizex, gr1_ymin ) ;
      line -> DrawLine( bestx - endbarsizex, gr1_ymax, bestx + endbarsizex, gr1_ymax ) ;
      line -> DrawLine( gr1_xmin, besty - endbarsizey, gr1_xmin, besty + endbarsizey ) ;
      line -> DrawLine( gr1_xmax, besty - endbarsizey, gr1_xmax, besty + endbarsizey ) ;

      line -> SetLineColor(4) ;
      line -> DrawLine( gr1_x_at_ymax, gr1_ymax, gr1_x_at_ymax, gr1_ymax - endbarsizey ) ;
      line -> DrawLine( gr1_x_at_ymin, gr1_ymin, gr1_x_at_ymin, gr1_ymin + endbarsizey ) ;
      line -> DrawLine( gr1_xmin, gr1_y_at_xmin, gr1_xmin + endbarsizex, gr1_y_at_xmin ) ;
      line -> DrawLine( gr1_xmax, gr1_y_at_xmax, gr1_xmax - endbarsizex, gr1_y_at_xmax ) ;


      TString savename( rootfile ) ;
      savename.ReplaceAll( ".root", "-nice.pdf" ) ;
      c_scan2d -> SaveAs( savename ) ;

   } // draw_2dscan



