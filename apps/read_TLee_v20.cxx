#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "WCPLEEANA/TLee.h"

#include "WCPLEEANA/Configure_Lee.h"

#include "TApplication.h"

//#include <chrono>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*  
  usage:
  
  make clean
  make  
  ./read_TLee_v20 -f 1 -p 1
  
  ---> README:
  ---> Makefile: comment the line "ROOTSYS=/home/xji/data0/software/root_build", if you have your own "ROOTSYS"
  ---> minuit2 is in the ROOT
*/

int main(int argc, char** argv)
{
  TString roostr = "";
  
  double scaleF_POT = 1;
  int ifile = 1;
  
  for(int i=1; i<argc; i++) {    
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
  }
  
  cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;

  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);
  
  if( !config_Lee::flag_display_graphics ) {
    gROOT->SetBatch( 1 );
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  TApplication theApp("theApp",&argc,argv);
  
  TLee *Lee_test = new TLee();
    
  ////////// just do it one time in the whole procedure

  Lee_test->channels_observation   = config_Lee::channels_observation;
  Lee_test->syst_cov_flux_Xs_begin = config_Lee::syst_cov_flux_Xs_begin;
  Lee_test->syst_cov_flux_Xs_end   = config_Lee::syst_cov_flux_Xs_end;
  Lee_test->syst_cov_mc_stat_begin = config_Lee::syst_cov_mc_stat_begin;
  Lee_test->syst_cov_mc_stat_end   = config_Lee::syst_cov_mc_stat_end;  
 
  Lee_test->flag_lookelsewhere     = config_Lee::flag_lookelsewhere;
  
  Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
                                      config_Lee::detector_directory, config_Lee::mc_directory);

  int size_array_LEE_ch = sizeof(config_Lee::array_LEE_ch)/sizeof(config_Lee::array_LEE_ch[0]);
  for(int idx=0; idx<size_array_LEE_ch; idx++) {
    if( config_Lee::array_LEE_ch[idx]!=0 ) Lee_test->map_Lee_ch[config_Lee::array_LEE_ch[idx]] = 1;
  }
  
  Lee_test->scaleF_POT = scaleF_POT;
  
  Lee_test->Set_Spectra_MatrixCov();
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();

  ////////// can do any times
  
  Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
  Lee_test->flag_syst_detector   = config_Lee::flag_syst_detector;
  Lee_test->flag_syst_additional = config_Lee::flag_syst_additional;
  Lee_test->flag_syst_mc_stat    = config_Lee::flag_syst_mc_stat;
  
  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  Lee_test->Set_Collapse();

  /////////////////////////////
  /////////////////////////////
  
  //file_collapsed_covariance_matrix->Close();
  
  //////////

  if( config_Lee::flag_plotting_systematics ) Lee_test->Plotting_systematics();
  
  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit
  
  //Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  //Lee_test->Set_Collapse();

  if( config_Lee::flag_GoF_output2file_default_0 ) {

    TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");
  
    TTree *tree_config = new TTree("tree", "configure information");
    int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
    int flag_syst_detector = config_Lee::flag_syst_detector;
    int flag_syst_additional = config_Lee::flag_syst_additional;
    int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
    double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
    double user_scaleF_POT = scaleF_POT;
    vector<double>vc_val_GOF;
    vector<int>vc_val_GOF_NDF;
    tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
    tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
    tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
    tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
    tree_config->Branch("user_Lee_strength_for_output_covariance_matrix", &user_Lee_strength_for_output_covariance_matrix,
			"user_Lee_strength_for_output_covariance_matrix/D" );
    tree_config->Branch("user_scaleF_POT", &user_scaleF_POT, "user_scaleF_POT/D" );
    tree_config->Branch("vc_val_GOF", &vc_val_GOF);
    tree_config->Branch("vc_val_GOF_NDF", &vc_val_GOF_NDF);
    file_collapsed_covariance_matrix->cd();

    Lee_test->matrix_absolute_cov_newworld.Write("matrix_absolute_cov_newworld");// (bins, bins)
    Lee_test->matrix_absolute_flux_cov_newworld.Write("matrix_absolute_flux_cov_newworld");
    Lee_test->matrix_absolute_Xs_cov_newworld.Write("matrix_absolute_Xs_cov_newworld");
    Lee_test->matrix_absolute_detector_cov_newworld.Write("matrix_absolute_detector_cov_newworld");
    Lee_test->matrix_absolute_mc_stat_cov_newworld.Write("matrix_absolute_mc_stat_cov_newworld");
    Lee_test->matrix_absolute_additional_cov_newworld.Write("matrix_absolute_additional_cov_newworld");
                 
    for(auto it=Lee_test->matrix_input_cov_detector_sub.begin(); it!=Lee_test->matrix_input_cov_detector_sub.end(); it++) {
      int idx = it->first;
      roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
      Lee_test->matrix_absolute_detector_sub_cov_newworld[idx].Write(roostr);
    }
     
    Lee_test->matrix_pred_newworld.Write("matrix_pred_newworld");// (1, bins)
    Lee_test->matrix_data_newworld.Write("matrix_data_newworld");// (1, bins)

    for(auto it_sub=Lee_test->matrix_sub_flux_geant4_Xs_newworld.begin(); it_sub!=Lee_test->matrix_sub_flux_geant4_Xs_newworld.end(); it_sub++) {
      int index = it_sub->first;
      roostr = TString::Format("matrix_sub_flux_geant4_Xs_newworld_%d", index);
      Lee_test->matrix_sub_flux_geant4_Xs_newworld[index].Write(roostr);
    }
  
    ///////
    
    file_collapsed_covariance_matrix->cd();

    for(auto it=Lee_test->map_data_spectrum_ch_bin.begin(); it!=Lee_test->map_data_spectrum_ch_bin.end(); it++) {
      int val_ch = it->first;
      int size_map = it->second.size();
      int size_before = 0;
      for(int idx=1; idx<val_ch; idx++) {
	int size_current = Lee_test->map_data_spectrum_ch_bin[idx].size();
	size_before += size_current;
      }
      
      vector<int>vc_target_chs;
      for(int ibin=1; ibin<size_map; ibin++) {
	vc_target_chs.push_back( size_before + ibin -1 );
      }
      
      vector<int>vc_support_chs;

      Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 100 + val_ch );

      vc_val_GOF.push_back( Lee_test->val_GOF_noConstrain );
      vc_val_GOF_NDF.push_back( Lee_test->val_GOF_NDF );
      //cout<<" ---> "<<Lee_test->val_GOF_NDF<<"\t"<<Lee_test->val_GOF_noConstrain<<endl;
    }
    
    tree_config->Fill();
    tree_config->Write();
    file_collapsed_covariance_matrix->Close();
  }

  ///////////////////////////////////////////////////////////////////////////////////  
  ///////////////////////////////////////////////////////////////////////////////////
  
  // bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
  // bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
  // bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
  // bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
  // bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
  // bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
  // bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
  // bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  if( 0 ) {
    if( 1 ) {
      vector<int>vc_target_chs;
      for(int ibin=1; ibin<=25; ibin++) vc_target_chs.push_back( 26*2 +ibin -1 );      
      vector<int>vc_support_chs;      
      Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 101 );
    }
    if( 1 ) {
      vector<int>vc_target_chs;
      for(int ibin=1; ibin<=25; ibin++) vc_target_chs.push_back( 26*3 +ibin -1 );      
      vector<int>vc_support_chs;      
      Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 102 );
    }
      
  }
  
  ///////////////////////// test
  
  if( 0 ) {// 8a
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 25 );// oldworld, newworld
    for(int idx=1; idx<=25; idx++) matrix_gof_trans(idx-1, idx-1) = 1;
    for(int idx=1; idx<=25; idx++) matrix_gof_trans(26+idx-1, idx-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 25, matrix_gof_trans.GetNcols()-25, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 801);
  }
  
  if( 0 ) {
    vector<int>vc_target_chs;
    for(int ibin=1; ibin<=25; ibin++) vc_target_chs.push_back( 26+ibin -1 );
      
    vector<int>vc_support_chs;
    //for(int ibin=1; ibin<=26*3 + 11*3; ibin++) vc_support_chs.push_back( 26+ibin -1 );
 
    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 100001 );
  }
  
  if( 0 ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 25 + 26*2+11*3 );// oldworld, newworld
    for(int idx=1; idx<=25; idx++) matrix_gof_trans(idx-1, idx-1) = 1;
    for(int idx=1; idx<=25; idx++) matrix_gof_trans(26+idx-1, idx-1) = 1;   
    for(int idx=1; idx<=26*2+11*3; idx++) matrix_gof_trans(26*2+idx-1, 25+idx-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 25, matrix_gof_trans.GetNcols()-25, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 100001);
  }
  
  if( 0 ) {
    vector<int>vc_target_chs;
    for(int ibin=1; ibin<=10; ibin++) vc_target_chs.push_back( 26*4 + 11*2  + ibin -1 );
      
    vector<int>vc_support_chs;
    for(int ibin=1; ibin<=26*2; ibin++) vc_support_chs.push_back( 26*2 +ibin -1 );
    
    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 101 );
  }
    
  if( 0 ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back(1);
    vc_target_chs.push_back(2);
 
    vector<int>vc_support_chs;  
    vc_support_chs.push_back(3);
    vc_support_chs.push_back(4);

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 3002 );
  }

  //////////////////////////////////////////////////////////////////////////////////////// LEE strength fitting

  if( config_Lee::flag_Lee_strength_data ) {

    cout<<endl;
    cout<<" -----------------------------------"<<endl;
    cout<<" LEE strength fitting"<<endl;
    cout<<" -----------------------------------"<<endl;
    
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(1, 0);// (initial value, fix or not)

    // cout<<endl<<TString::Format(" ---> Best-fit Lee strength: LEEx = %6.4f, chi2/ndf = %6.3f/%d",
    // 				Lee_test->minimization_Lee_strength_val,
    // 				Lee_test->minimization_chi2,
    // 				Lee_test->minimization_NDF
    // 				)<<endl<<endl;    

    cout<<endl<<TString::Format(" ---> Best-fit Lee strength: LEEx = %6.4f, chi2 = %6.3f, status = %d warning",
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_chi2,
				Lee_test->minimization_status
				)<<endl;
    
    // double user_bestFit = Lee_test->minimization_Lee_strength_val;
    
    /////////////////////////////////////////
    
    double gmin = Lee_test->minimization_chi2;
    TGraph *gh_scan = new TGraph();
    double slow = 0;
    double shgh = 3;
    int nscan = 100;
    double val_max_dchi2 = 0;
    double step = (shgh-slow)/nscan;
    for(int idx=1; idx<=nscan; idx++) {
      //if( idx%(max(1, nscan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./nscan, idx)<<endl;
      double val_s = slow + (idx-1)*step;
      Lee_test->Minimization_Lee_strength_FullCov(val_s, 1);// (initial value, fix or not)
      double val_chi2 = Lee_test->minimization_chi2;
      gh_scan->SetPoint( gh_scan->GetN(), val_s, val_chi2 - gmin);
      if( val_max_dchi2<val_chi2 - gmin ) val_max_dchi2 = val_chi2 - gmin;
    }
    
    double val_dchi2at0 = gh_scan->Eval(0);
    double val_dchi2at1 = gh_scan->Eval(1);
    if( fabs(val_dchi2at0)<1e-6 ) val_dchi2at0 = 0;
    double val_dchi2at1_minus_dchi2at0 = val_dchi2at1 - val_dchi2at0;
    
    //cout<<endl<<Form(" ---> dchi2 at LEEx 0/1: %7.4f %7.4f", val_dchi2at0, val_dchi2at1 )<<endl<<endl;

    cout<<endl;
    cout<<Form(" ---> Dchi2 from @LEEx1 minus @LEExMIN: %7.4f", val_dchi2at1)<<endl<<endl;
    cout<<Form(" ---> Dchi2 from @LEEx0 minus @LEExMIN: %7.4f", val_dchi2at0)<<endl<<endl;   
    cout<<Form(" ---> Dchi2 from @LEEx1 minus @LEEx0:   %7.4f", val_dchi2at1_minus_dchi2at0)<<endl<<endl;
   

    // gh_scan->SetName("gh_dchi2");
    // gh_scan->SaveAs("file_dchi2_scan.root");
    
    
    // TCanvas *canv_gh_scan = new TCanvas("canv_gh_scan", "canv_gh_scan", 900, 650);
    // canv_gh_scan->SetLeftMargin(0.15); canv_gh_scan->SetRightMargin(0.1);
    // canv_gh_scan->SetTopMargin(0.1); canv_gh_scan->SetBottomMargin(0.15);    
    // gh_scan->Draw("al");
    // gh_scan->GetXaxis()->SetTitle("LEE strength"); gh_scan->GetYaxis()->SetTitle("#Delta#chi^{2}");    
    // gh_scan->GetXaxis()->SetLabelSize(0.05); gh_scan->GetXaxis()->SetTitleSize(0.05);
    // gh_scan->GetYaxis()->SetLabelSize(0.05); gh_scan->GetYaxis()->SetTitleSize(0.05);
    // gh_scan->GetXaxis()->CenterTitle(); gh_scan->GetYaxis()->CenterTitle();
    // gh_scan->GetXaxis()->SetTitleOffset(1.2);
    // //gh_scan->GetYaxis()->SetRangeUser(0, val_max_dchi2*1.1);

    // gh_scan->SetName("gh_scan");
    // gh_scan->SaveAs("gh_scan_1eNp.root");
    
    // // TLine *lineA_dchi2at1 = new TLine(1, 0, 1, val_dchi2at1);    
    // // lineA_dchi2at1->Draw("same");
    // // lineA_dchi2at1->SetLineWidth(2);
    // // lineA_dchi2at1->SetLineColor(kBlue);
    // // lineA_dchi2at1->SetLineStyle(7);
    // // TLine *lineB_dchi2at1 = new TLine(0, val_dchi2at1, 1, val_dchi2at1);    
    // // lineB_dchi2at1->Draw("same");
    // // lineB_dchi2at1->SetLineWidth(2);
    // // lineB_dchi2at1->SetLineColor(kBlue);
    // // lineB_dchi2at1->SetLineStyle(7);
    // // auto *tt_text_data = new TLatex( 0.2, val_dchi2at1*1.1, Form("#Delta#chi^{2} = %4.3f", val_dchi2at1) );
    // // tt_text_data->SetTextAlign(11); tt_text_data->SetTextSize(0.05); tt_text_data->SetTextAngle(0);
    // // tt_text_data->SetTextFont(42);  tt_text_data->Draw(); tt_text_data->SetTextColor(kBlue);

    // //canv_gh_scan->SaveAs("canv_gh_scan.png");

  }
 
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Advanced Statistics Analysis /////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  //
  // (*) Hypothesis test always reject the null
  //
  // [*] the realization of the functionality is a combination of the tools:
  //
  //         Lee_test->Set_measured_data();
  //
  //         Lee_test->scaleF_Lee = #;
  //         Lee_test->Set_Collapse();
  //
  //         Lee_test->Set_toy_Asimov();
  //
  //         Lee_test->Set_Variations( # );
  //         Lee_test->Set_toy_Variation( # );
  //
  //         Lee_test->Minimization_Lee_strength_FullCov(#, #);
  //
  
  /////////////////////////////////////////////////////// example: do fitting on Asimov sample

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
  
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }

  ////////////////////////////////////////////////////// example: do fitting on variation sample

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();

    Lee_test->Set_Variations( 10 );// generate 10 variation samples
    Lee_test->Set_toy_Variation( 4 );// use the 4th sample as the input data for the fitting
    
    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////// example: simple versus simple likelihood ratio test

  if( 0 ) {    
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// (initial value, fix or not)
    double val_chi2_Lee = Lee_test->minimization_chi2;

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);// (initial value, fix or not)
    double val_chi2_sm = Lee_test->minimization_chi2;

    double val_dchi2 = val_chi2_Lee - val_chi2_sm;
    
    cout<<endl<<TString::Format(" ---> dchi2 = Lee - sm: %7.4f, LEE %7.4f, sm %7.4f", val_dchi2, val_chi2_Lee, val_chi2_sm)<<endl<<endl;
  }

  ////////////////////////////////////////////////////////// sensitivity calcualtion by FC
  
  if( 0 ) {
    
    double chi2_null_null8sm_true8sm  = 0;
    double chi2_gmin_null8sm_true8sm  = 0;
    double chi2_null_null8Lee_true8Lee = 0;
    double chi2_gmin_null8Lee_true8Lee = 0;
    
    TFile *file_out = new TFile(TString::Format("file_out_%03d.root", ifile), "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, "chi2_null_null8sm_true8sm/D" );
    tree->Branch("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, "chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, "chi2_null_null8Lee_true8Lee/D" );
    tree->Branch("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, "chi2_gmin_null8Lee_true8Lee/D" );

    int N_toy = 1000;

    auto time_start = chrono::high_resolution_clock::now();
    
    for(int itoy=1; itoy<=N_toy; itoy++) {
            
      if( itoy%max(N_toy/10,1)==0 ) {
	cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./N_toy, itoy)<<endl;
      }
      cout<<Form(" running %6d", itoy)<<endl;
      
      int status_fit = 0;
          
      /////////////////////////////////// null8sm, true8sm
      
      Lee_test->scaleF_Lee = 0;
      Lee_test->Set_Collapse();    
      Lee_test->Set_Variations(1);
      Lee_test->Set_toy_Variation(1);
    
      Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      chi2_null_null8sm_true8sm = Lee_test->minimization_chi2;

      Lee_test->Minimization_Lee_strength_FullCov(1, 0);
      chi2_gmin_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      /////////////////////////////////// null8Lee, true8Lee
      
      Lee_test->scaleF_Lee = 1;
      Lee_test->Set_Collapse();    
      Lee_test->Set_Variations(1);
      Lee_test->Set_toy_Variation(1);
    
      Lee_test->Minimization_Lee_strength_FullCov(1, 1);
      chi2_null_null8Lee_true8Lee = Lee_test->minimization_chi2;

      Lee_test->Minimization_Lee_strength_FullCov(1, 0);
      chi2_gmin_null8Lee_true8Lee = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      ///////////////////////////////////
      
      if( status_fit!=0 ) continue;
      tree->Fill();
    }
    
    auto time_stop = chrono::high_resolution_clock::now();
    auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
    cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
    
    file_out->cd();
    tree->Write();
    file_out->Close();
    
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// nested PDF

if( 0 ) {    
  ofstream real_data;
  // append ifile to the name of the file
	real_data.open("/mnt/d/eLEE/plots/run1_to_run3/spectrum/plots/nested_PDF/old_TLee_improvement/real_data.txt");
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// (initial value, fix or not)
    double val_chi2_Lee = Lee_test->minimization_chi2;

    Lee_test->Minimization_Lee_strength_FullCov(1, 0);// (initial value, fix or not)
    double val_chi2_min = Lee_test->minimization_chi2;

    double val_dchi2 = val_chi2_Lee - val_chi2_min;
    
    real_data << "real data chi2_min, chi2_LEE, chi2_LEE_minus_min: " << val_chi2_min << ", " << val_chi2_Lee << ", " <<  val_dchi2 << "\n";
  }


    if( 0 ) {
    double chi2_null_null8Lee_true8Lee = 0; 
    double chi2_gmin_null8Lee_true8Lee = 0; 
    int status_fit = 0;                     
    int itoy = 0; 
    int success_toy = 0;                  

    TString file_name = TString::Format("/mnt/d/eLEE/plots/run1_to_run3/spectrum/plots/nested_PDF/old_TLee_improvement/file_out_%03d.root", ifile);
    TFile *file_out = new TFile(file_name, "RECREATE");
    if (!file_out || file_out->IsZombie()) {
        cerr << "[ERROR] 无法创建输出文件：" << file_name << endl;
        return 1;
    }

    TTree *tree = new TTree("tree", "Toy fit results (Lee strength)");
    tree->Branch("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, "chi2_null_null8Lee_true8Lee/D");
    tree->Branch("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, "chi2_gmin_null8Lee_true8Lee/D");
    tree->Branch("status_fit", &status_fit, "status_fit/I");

    int N_toy = 1000;
    for (itoy = 1; itoy <= N_toy; itoy++) {
        if (itoy % max(N_toy / 10, 1) == 0) {
            cout << TString::Format("Toy: %.2f%%, %d", itoy * 100.0 / N_toy, itoy) << endl;
        }

        status_fit = 0;
        chi2_null_null8Lee_true8Lee = 0;  
        chi2_gmin_null8Lee_true8Lee = 0;

        Lee_test->scaleF_Lee = 1;
        Lee_test->Set_Collapse();    
        Lee_test->Set_Variations(N_toy);
        Lee_test->Set_toy_Variation(itoy); 

        Lee_test->Minimization_Lee_strength_FullCov(1, 1);
        int status1 = Lee_test->minimization_status;
        if (status1 == 0) {
            chi2_null_null8Lee_true8Lee = Lee_test->minimization_chi2;
        }
        status_fit += status1;

        Lee_test->Minimization_Lee_strength_FullCov(1, 0);
        int status2 = Lee_test->minimization_status;
        if (status2 == 0) {
            chi2_gmin_null8Lee_true8Lee = Lee_test->minimization_chi2;
        }
        status_fit += status2;
        tree->Fill();

    }

    file_out->cd();          
    tree->Write();  
    file_out->Flush();       
    file_out->Close();       
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////// Sensitivity by Asimov sample

  if( 0 ) {

    cout<<endl;
    cout<<" -----------------------------------"<<endl;
    cout<<" Sensitivity by Asimov sample"<<endl;
    cout<<" -----------------------------------"<<endl;
    cout<<endl;
  
    ///////////////////////// reject SM
    
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
    
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting
    Lee_test->Minimization_Lee_strength_FullCov(0, 1);// (initial value, fix or not)

    double sigma_SM = sqrt( Lee_test->minimization_chi2 );
    cout<<TString::Format(" ---> Excluding  SM: %5.2f sigma", sigma_SM)<<endl;
    
    ///////////////////////// reject 1*LEE
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();
    
    Lee_test->Set_toy_Asimov();// use the Asimov sample as the input data for the fitting
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// (initial value, fix or not)

    double sigma_Lee = sqrt( Lee_test->minimization_chi2 );
    cout<<TString::Format(" ---> Excluding LEE: %5.2f sigma", sigma_Lee)<<endl<<endl;;
    
  }

  ////////////////////////////////////////////////  Feldman-Cousins approach --> heavy computation cost

  if( 0 ) {
    
    /////////////// range: [low, hgh] with step
    
    double Lee_true_low = 0;
    double Lee_true_hgh = 3;
    double Lee_step     = 0.02;
    
    /////////////// dchi2 distribution 
    
    //int num_toy = 2;    
    //Lee_test->Exe_Feldman_Cousins(Lee_true_low, Lee_true_hgh, Lee_step, num_toy, ifile);

    /////////////// dchi2 of Asimov sample
   
    //Lee_test->Exe_Fledman_Cousins_Asimov(Lee_true_low, Lee_true_hgh, Lee_step);

    /////////////// dchi2 of measured data

    if( config_Lee::flag_Lee_scan_data ) {
      Lee_test->Set_measured_data();    
      TMatrixD matrix_data_input_fc = Lee_test->matrix_data_newworld;    
      Lee_test->Exe_Fiedman_Cousins_Data( matrix_data_input_fc, Lee_true_low, Lee_true_hgh, Lee_step );
    }
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" -----------------------------------"<<endl;
  cout<<" Check the initialization at the end"<<endl;
  cout<<" -----------------------------------"<<endl;
  
  cout<<endl;
  cout<<" ---> flag_syst_flux_Xs    "<<Lee_test->flag_syst_flux_Xs<<endl;
  cout<<" ---> flag_syst_detector   "<<Lee_test->flag_syst_detector<<endl;
  cout<<" ---> flag_syst_additional "<<Lee_test->flag_syst_additional<<endl;
  cout<<" ---> flag_syst_mc_stat    "<<Lee_test->flag_syst_mc_stat<<endl;  
  cout<<endl;

  cout<<" ---> LEE channel size (set by array_LEE_ch in config): "<<Lee_test->map_Lee_ch.size()<<endl;
  if( (int)(Lee_test->map_Lee_ch.size()) ) {
    for(auto it_map_Lee=Lee_test->map_Lee_ch.begin(); it_map_Lee!=Lee_test->map_Lee_ch.end(); it_map_Lee++) {
      cout<<" ---> LEE channel: "<< it_map_Lee->first<<endl;
    }
  }
  cout<<endl;

  cout<<" ---> MC stat file: "<<Lee_test->syst_cov_mc_stat_begin<<".log - "<<Lee_test->syst_cov_mc_stat_end<<".log"<<endl;
  cout<<endl;
  
  cout<<" ---> check: scaleF_POT "<<scaleF_POT<<", file_flag "<<ifile<<endl<<endl;

  cout<<endl<<endl;
  cout<<" ---> Complete all the program"<<endl;
  cout<<endl<<endl;
  
  if( config_Lee::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;
    
    theApp.Run();
  }
  
  return 0;
}
