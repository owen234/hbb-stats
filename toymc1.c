
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
#include "TString.h"
#include "TRandom.h"


#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"

  using namespace RooFit;
  using namespace RooStats;
  using std::cout ;

#include "build_hbb_workspace1.c"

   double toy_mean_N_msig[bins_of_nb][max_bins_of_met] ;  // bins_of_nb and max_bins_of_met are defined in build_hbb_workspace1.c (included above).
   double toy_mean_N_msb[bins_of_nb][max_bins_of_met]  ;  // bins_of_nb and max_bins_of_met are defined in build_hbb_workspace1.c (included above).

   RooAbsPdf* likelihood ;
   RooRealVar* rv_sig_strength ;
   RooDataSet* rds ;
   RooWorkspace* toyws ;
   int bins_of_met ;
   TRandom* tran ;
   TFile* ttfile ;
   TTree* toytt ;

   bool getMeanValsFromInput( const char* infile ) ;
   float findUL( RooDataSet* toyds ) ;

   RooDataSet* genToyData() ;




  //===================================================================================

   void toymc1( const char* infile = "outputfiles/input-file.txt", int nToy = 100 ) {

      TString toywsfile( infile ) ;
      toywsfile.ReplaceAll("input-","toyws-") ;
      toywsfile.ReplaceAll("txt","root") ;

      build_hbb_workspace1( infile, toywsfile ) ;

      TFile* wstf = new TFile( toywsfile ) ;

      toyws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      toyws->Print() ;

      bins_of_met = TMath::Nint( toyws->var("bins_of_met")->getVal()  ) ;
      printf("\n\n Bins of MET : %d\n\n", bins_of_met ) ;

      rds = (RooDataSet*) toyws->obj( "hbb_observed_rds" ) ;
      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      likelihood = toyws->pdf("likelihood") ;

      rv_sig_strength = toyws->var("sig_strength") ;
      if ( rv_sig_strength == 0x0 ) { printf("\n\n *** can't find sig_strength in workspace.\n\n" ) ; return ; }

      getMeanValsFromInput( infile ) ;

      tran = new TRandom(12345) ;

      TH1F* hss = new TH1F("hss","signal strength", 100, 0., 5. ) ;

      TString outfile( infile ) ;
      outfile.ReplaceAll( "input-", "toytt-" ) ;
      outfile.ReplaceAll( "txt", "root" ) ;


      ttfile = new TFile( outfile, "recreate" ) ;
      toytt = new TTree( "toytt", "Toy MC study" ) ;
      float fit_sig_strength ;
      float fit_sig_strength_err ;
      float fit_sig_signif ;
      float fit_sig_ul ;
      toytt -> Branch( "fit_sig_strength", &fit_sig_strength, "fit_sig_strength/F" ) ;
      toytt -> Branch( "fit_sig_strength_err", &fit_sig_strength_err, "fit_sig_strength_err/F" ) ;
      toytt -> Branch( "fit_sig_signif", &fit_sig_signif, "fit_sig_signif/F" ) ;
      toytt -> Branch( "fit_sig_ul", &fit_sig_ul, "fit_sig_ul/F" ) ;


      RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
      RooMsgService::instance().getStream(1).removeTopic(Fitting) ;

      for ( int ti=0; ti<nToy; ti++ ) {

         RooDataSet* toyds = genToyData() ;
         if ( toyds == 0x0 ) { printf("\n\n *** generation of toy ds failed!\n\n") ; return ; }
         //toyds -> printMultiline( cout, 1, kTRUE, "" ) ;

        //-- fit with signal floating.

         rv_sig_strength -> setConstant( kFALSE ) ;

         RooFitResult* fitResult = likelihood->fitTo( *toyds, Save(true), PrintLevel(-1), Hesse(true), Strategy(1) ) ;

         double minNllFloat = fitResult -> minNll() ;

         fit_sig_strength = rv_sig_strength->getVal() ;
         fit_sig_strength_err = rv_sig_strength->getError() ;


        //-- significance

         rv_sig_strength -> setVal(0.00000001) ;
         rv_sig_strength -> setConstant( kTRUE ) ;
         RooFitResult* fitResult0 = likelihood->fitTo( *toyds, Save(true), PrintLevel(-1), Hesse(true), Strategy(1) ) ;

         float minNll0 = fitResult0 -> minNll() ;
         float testStat = 2.*( minNll0 - minNllFloat ) ;
         fit_sig_signif = 0. ;
         if ( testStat > 0 ) { fit_sig_signif = sqrt( testStat ) ; }
         delete fitResult0 ;


        //-- upper limit.

         fit_sig_ul = findUL( toyds ) ;




         printf("  toy %4d , signal strength : %5.2f +/- %5.2f,  signif = %5.2f, 1-sided UL = %5.2f\n",
             ti, fit_sig_strength, fit_sig_strength_err, fit_sig_signif, fit_sig_ul ) ;

         hss->Fill( fit_sig_strength ) ;

         toytt -> Fill() ;

         delete toyds ;
         delete fitResult ;

      } // ti.

      gStyle->SetOptStat("emr") ;
      hss->Draw() ;

      ttfile -> cd() ;
      toytt -> Write() ;
      ttfile -> Close() ;

   } // toymc1

  //===================================================================================

   bool getMeanValsFromInput( const char* infile ) {

      char pname[1000] ;
      float fileVal ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return false ; }
            toy_mean_N_msig[nbi][mbi] = fileVal ;

            sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return false ; }
            toy_mean_N_msb[nbi][mbi] = fileVal ;

         } // mbi.
      } // nbi.

      return true ;

   } // getMeanValsFromInput.

  //===================================================================================


   RooDataSet* genToyData() {

      RooArgSet obsList ;

      char pname[1000] ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            int nobs ;
            RooRealVar* rrv ;

            sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
            rrv = toyws -> var( pname ) ;
            if ( rrv == 0x0 ) { printf("\n\n\n *** genToyData: can't find %s in workspace.\n\n", pname ) ; fflush(stdout) ; return 0x0 ; }
            nobs = tran->Poisson( toy_mean_N_msig[nbi][mbi] ) ;
            rrv -> setVal( nobs ) ;
            obsList.add( *rrv ) ;

            sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
            rrv = toyws -> var( pname ) ;
            if ( rrv == 0x0 ) { printf("\n\n\n *** genToyData: can't find %s in workspace.\n\n", pname ) ; fflush(stdout) ; return 0x0 ; }
            nobs = tran->Poisson( toy_mean_N_msb[nbi][mbi] ) ;
            rrv -> setVal( nobs ) ;
            obsList.add( *rrv ) ;

         } // mbi.
      } // nbi.

      RooDataSet* trds = new RooDataSet("toyfit_hbb_observed_rds", "toy dataset", obsList ) ;
      trds -> add( obsList ) ;

      return trds ;


   } // genToyData.

  //===================================================================================

   float findUL( RooDataSet* toyds ) {

      rv_sig_strength -> setConstant( kFALSE ) ;
      RooFitResult* rfr_float = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;
      double minNll_float = rfr_float -> minNll() ;
      double sig_strength_float = rv_sig_strength -> getVal() ;
      double sig_strength_err  = rv_sig_strength -> getError() ;
      delete rfr_float ;

      float sig_strength_a = sig_strength_float + 1.5 * sig_strength_err ;
      float sig_strength_b = sig_strength_float + 2.0 * sig_strength_err ;

      rv_sig_strength -> setConstant( kTRUE ) ;

      rv_sig_strength -> setVal( sig_strength_a ) ;
      RooFitResult* rfr_a = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;
      double minNll_a = rfr_a -> minNll() ;
      double ts_a = 2.*( minNll_a - minNll_float ) ;
      delete rfr_a ;

      rv_sig_strength -> setVal( sig_strength_b ) ;
      RooFitResult* rfr_b = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;
      double minNll_b = rfr_b -> minNll() ;
      double ts_b = 2.*( minNll_b - minNll_float ) ;
      delete rfr_b ;


      float sig_strength_g1 = sig_strength_a + (2.70 - ts_a) * ( sig_strength_b - sig_strength_a ) / ( ts_b - ts_a ) ;

      rv_sig_strength -> setVal( sig_strength_g1 ) ;
      RooFitResult* rfr_g1 = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;
      double minNll_g1 = rfr_g1 -> minNll() ;
      double ts_g1 = 2.*( minNll_g1 - minNll_float ) ;
      delete rfr_g1 ;

      float sig_strength_g2 = sig_strength_a + (2.70 - ts_a) * ( sig_strength_g1 - sig_strength_a ) / ( ts_g1 - ts_a ) ;

      //rv_sig_strength -> setVal( sig_strength_g2 ) ;
      //RooFitResult* rfr_g2 = likelihood -> fitTo( *toyds, Save(true), Hesse(false), Strategy(1), PrintLevel(-1) ) ;
      //double minNll_g2 = rfr_g2 -> minNll() ;
      //double ts_g2 = 2.*( minNll_g2 - minNll_float ) ;
      //delete rfr_g2 ;

      return sig_strength_g2 ;

   } // findUL.

  //===================================================================================





