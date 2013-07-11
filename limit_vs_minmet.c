
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"

    const int nmetcut(4) ;

    double ul_m2s[nmetcut] ;
    double ul_m1s[nmetcut] ;
    double ul_med[nmetcut] ;
    double ul_p1s[nmetcut] ;
    double ul_p2s[nmetcut] ;

    double min_metsig[nmetcut] ;
    double dex[nmetcut] ;

    double ymax ;

    TH2F* hdummy(0x0) ;

    bool read_numbers( int plotIndex = 0, const char* number_file = "all-results.txt" ) ;

    void limit_vs_minmet( int plotIndex = 0, const char* number_file = "all-results.txt" ) {

        gStyle -> SetOptStat(0) ;

        gStyle -> SetPadGridX(1) ;
        gStyle -> SetPadGridY(1) ;
        gStyle -> SetTitleOffset( 1.1, "x" ) ;
        gStyle -> SetTitleOffset( 1.3, "y" ) ;


        if ( ! read_numbers( plotIndex, number_file ) ) return ;

       int npoints ;
       if ( plotIndex < 10 ) { npoints = 4 ; } else { npoints = 2 ; }

       TGraphAsymmErrors* gr1s = new TGraphAsymmErrors( npoints, min_metsig, ul_med, dex, dex, ul_m1s, ul_p1s ) ;
       TGraphAsymmErrors* gr2s = new TGraphAsymmErrors( npoints, min_metsig, ul_med, dex, dex, ul_m2s, ul_p2s ) ;

       gr2s->SetFillColor(5) ;

       gr1s->SetFillColor(3) ;
       gr1s->SetLineColor(1) ;
       gr1s->SetLineWidth(2) ;
       gr1s->SetMarkerStyle(20) ;





       hdummy -> Draw() ;

       if ( plotIndex < 10 ) {
          gr2s -> Draw( "3" ) ;
          gr1s -> Draw( "3" ) ;
          gr1s -> Draw( "PL" ) ;
       } else {
          gr2s -> Draw( "3" ) ;
          gr1s -> Draw( "3" ) ;
          TLine* medline = new TLine() ;
          medline -> SetLineWidth(2) ;
          medline -> DrawLine( min_metsig[0], ul_med[0], min_metsig[1], ul_med[0] ) ;
       }

       hdummy -> Draw("axis same") ;

       TLine* line = new TLine() ;
       line -> SetLineWidth(2) ;
       line -> SetLineColor(2) ;
       line -> DrawLine(0., 1., 175., 1. ) ;


    }

  //=============================================


    bool read_numbers( int plotIndex, const char* number_file  ) {

       char target_label[1000] ;
       int  target_label_n(0) ;
       int  target_lines(0) ;

       if ( plotIndex == 0 ) {

          printf("\n\n\n") ;
          printf("  1 = abcd-nosig-250\n" ) ;
          printf("  2 = abcd-nosig-400\n" ) ;
          printf("  3 = onebin-nosig-250\n" ) ;
          printf("  4 = onebin-nosig-250\n" ) ;
          printf("  5 = asymptotic-nosig-250\n" ) ;
          printf("  6 = asymptotic-nosig-400\n" ) ;
          printf("  7 = toys-nosig-250\n" ) ;
          printf("  8 = toys-nosig-400\n" ) ;
          printf("\n") ;
          printf(" 15 = asymptotic-nosig-250, 4 met bins\n") ;
          printf(" 16 = asymptotic-nosig-400, 4 met bins\n") ;
          printf(" 17 = toys-nosig-250, 4 met bins\n") ;
          printf(" 18 = toys-nosig-400, 4 met bins\n") ;
          printf("\n\n\n") ;
          return false ;

       }


       if ( plotIndex == 1 ) {
          sprintf( target_label, "post-abcd-nosig-250" ) ; target_label_n = 19 ; target_lines = 4 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy1", "Upper Limit vs min METsig, LandS ABCD, sig250", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 2 ) {
          sprintf( target_label, "post-abcd-nosig-400" ) ; target_label_n = 19 ; target_lines = 4 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy2", "Upper Limit vs min METsig, LandS ABCD, sig400", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 3 ) {
          sprintf( target_label, "post-onebin-nosig-250" ) ; target_label_n = 21 ; target_lines = 4 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy3", "Upper Limit vs min METsig, LandS one observable, sig250", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 4 ) {
          sprintf( target_label, "post-onebin-nosig-400" ) ; target_label_n = 21 ; target_lines = 4 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy4", "Upper Limit vs min METsig, LandS one observable, sig400", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 5 ) {
          sprintf( target_label, "limit-asymptotic-nosig-250-1metbin" ) ; target_label_n = 34 ; target_lines = 4 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy5", "Upper Limit vs min METsig, RooStats asymptotic, sig250", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 6 ) {
          sprintf( target_label, "limit-asymptotic-nosig-400-1metbin" ) ; target_label_n = 34 ; target_lines = 4 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy6", "Upper Limit vs min METsig, RooStats asymptotic, sig400", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 7 ) {
          sprintf( target_label, "limit-toys-nosig-250-1metbin" ) ; target_label_n = 28 ; target_lines = 4 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy7", "Upper Limit vs min METsig, RooStats toys, sig250", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 8 ) {
          sprintf( target_label, "limit-toys-nosig-400-1metbin" ) ; target_label_n = 28 ; target_lines = 4 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy8", "Upper Limit vs min METsig, RooStats toys, sig400", 2, 0., 175., 2., 0., ymax ) ;
       }

       if ( plotIndex == 15 ) {
          sprintf( target_label, "limit-asymptotic-nosig-250-4metbin" ) ; target_label_n = 34 ; target_lines = 1 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy15", "Upper Limit vs min METsig, RooStats asymptotic, sig250, 4 met bins", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 16 ) {
          sprintf( target_label, "limit-asymptotic-nosig-400-4metbin" ) ; target_label_n = 34 ; target_lines = 1 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy16", "Upper Limit vs min METsig, RooStats asymptotic, sig400, 4 met bins", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 17 ) {
          sprintf( target_label, "limit-toys-nosig-250-4metbin" ) ; target_label_n = 28 ; target_lines = 1 ;
          ymax = 1.2 ;
          hdummy = new TH2F("hdummy17", "Upper Limit vs min METsig, RooStats toys, sig250, 4 met bins", 2, 0., 175., 2., 0., ymax ) ;
       }
       if ( plotIndex == 18 ) {
          sprintf( target_label, "limit-toys-nosig-400-4metbin" ) ; target_label_n = 28 ; target_lines = 1 ;
          ymax = 3.0 ;
          hdummy = new TH2F("hdummy18", "Upper Limit vs min METsig, RooStats toys, sig400, 4 met bins", 2, 0., 175., 2., 0., ymax ) ;
       }



       if ( target_lines == 4 ) {
          min_metsig[0] =  30. ;
          min_metsig[1] =  50. ;
          min_metsig[2] = 100. ;
          min_metsig[3] = 150. ;
          hdummy -> SetXTitle("Minimum METsig") ;
       } else {
          min_metsig[0] =  30. ;
          min_metsig[1] = 150. ;
       }
       hdummy -> SetYTitle("Signal strength UL, 95%") ;


       FILE* infp ;
       if ( (infp=fopen( number_file, "r" )) == NULL ) {
          printf("\n\n *** problems opening all-results.txt.\n\n") ;
          return false ;
       }

       int linecount(0) ;
       while ( !feof(infp) ) {
          char label[1000] ; char s1[100] ; char s2[100] ;
          float m2s, m1s, med, p1s, p2s ;
          fscanf( infp, "%s %s %s %f %f %f %f %f", label, s1, s2,   &m2s, &m1s, &med, &p1s, &p2s ) ;
          if ( feof(infp) ) break ;
          if ( strncmp( label, target_label, target_label_n ) == 0 ) {
             ul_m2s[linecount] = med - m2s ;
             ul_m1s[linecount] = med - m1s ;
             ul_med[linecount] = med ;
             ul_p1s[linecount] = p1s - med ;
             ul_p2s[linecount] = p2s - med ;
             linecount++ ;
          }
       }
       printf("\n found %d lines matching %s.\n\n", linecount, target_label ) ;

       if ( plotIndex > 10 ) {
          ul_m2s[1] = ul_m2s[0] ;
          ul_m1s[1] = ul_m1s[0] ;
          ul_med[1] = ul_med[0] ;
          ul_p1s[1] = ul_p1s[0] ;
          ul_p2s[1] = ul_p2s[0] ;
       }

       if ( linecount != target_lines ) return false ;

       return true ;

    }


  //=============================================



