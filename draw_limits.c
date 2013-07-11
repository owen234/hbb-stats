
#include "TCanvas.h"

#include "limit_vs_minmet.c"


   void draw_limits( const char* number_file = "all-results.txt" ) {

       gStyle -> SetPadGridX(1) ;
       gStyle -> SetPadGridY(1) ;
       gStyle -> SetTitleW(0.95) ;

       int cw(800), ch(1200) ;

       int ci ;


       //----

       TCanvas* can_sig250 = new TCanvas("limits_sig250", "Upper limits, sig250", cw, ch ) ;

       can_sig250 -> Divide(2,3) ;

       ci = 1 ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(5, number_file) ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(15, number_file) ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(7, number_file) ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(17, number_file) ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(1, number_file) ;

       can_sig250 -> cd(ci++) ;
       limit_vs_minmet(3, number_file) ;

       can_sig250 -> SaveAs("outputfiles/limits_sig250.pdf") ;


       //----

       TCanvas* can_sig400 = new TCanvas("limits_sig400", "Upper limits, sig400", cw, ch ) ;

       can_sig400 -> Divide(2,3) ;

       ci = 1 ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(6, number_file) ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(16, number_file) ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(8, number_file) ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(18, number_file) ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(2, number_file) ;

       can_sig400 -> cd(ci++) ;
       limit_vs_minmet(4, number_file) ;

       can_sig400 -> SaveAs("outputfiles/limits_sig400.pdf") ;


   } // draw_limits



