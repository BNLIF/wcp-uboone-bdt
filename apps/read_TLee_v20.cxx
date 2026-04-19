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

  Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
  Lee_test->flag_syst_reweight    = config_Lee::flag_syst_reweight;
  Lee_test->flag_syst_reweight_cor    = config_Lee::flag_syst_reweight_cor;
    
  ////////// just do it one time in the whole procedure

  Lee_test->channels_observation   = config_Lee::channels_observation;
  Lee_test->syst_cov_flux_Xs_begin = config_Lee::syst_cov_flux_Xs_begin;
  Lee_test->syst_cov_flux_Xs_end   = config_Lee::syst_cov_flux_Xs_end;
  Lee_test->syst_cov_mc_stat_begin = config_Lee::syst_cov_mc_stat_begin;
  Lee_test->syst_cov_mc_stat_end   = config_Lee::syst_cov_mc_stat_end;  
 
  Lee_test->flag_lookelsewhere     = config_Lee::flag_lookelsewhere;

  Lee_test->scaleF_POT = scaleF_POT;
  Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
                                      config_Lee::detector_directory, config_Lee::mc_directory);
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

  Lee_test->flag_Lee_minimization_after_constraint = config_Lee::flag_Lee_minimization_after_constraint;
  
  //////////

  if( 0 ) {// shape only covariance

    int nbins = Lee_test->bins_newworld;
    TMatrixD matrix_pred = Lee_test->matrix_pred_newworld;
    TMatrixD matrix_syst = Lee_test->matrix_absolute_cov_newworld;    
    TMatrixD matrix_shape(nbins, nbins);
    TMatrixD matrix_mixed(nbins, nbins);
    TMatrixD matrix_norm(nbins, nbins);
    
    ///
    double N_T = 0;
    for(int idx=0; idx<nbins; idx++) N_T += matrix_pred(0, idx);

    ///
    double M_kl = 0;
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {
	M_kl += matrix_syst(i,j);
      }
    }

    ///
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {      
	double N_i = matrix_pred(0, i);
	double N_j = matrix_pred(0, j);

	double M_ij = matrix_syst(i,j);      
	double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
	double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);

	matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
	matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;	
	matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
      }
    }

    Lee_test->matrix_absolute_cov_newworld = matrix_shape;
    
  }// shape only covariance

  if( 0 ) {// shape only covariance

    int nbins = Lee_test->bins_newworld;
    TMatrixD matrix_pred = Lee_test->matrix_pred_newworld;
    TMatrixD matrix_syst = Lee_test->matrix_absolute_flux_cov_newworld;    
    TMatrixD matrix_shape(nbins, nbins);
    TMatrixD matrix_mixed(nbins, nbins);
    TMatrixD matrix_norm(nbins, nbins);
    
    ///
    double N_T = 0;
    for(int idx=0; idx<nbins; idx++) N_T += matrix_pred(0, idx);

    ///
    double M_kl = 0;
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {
	M_kl += matrix_syst(i,j);
      }
    }

    ///
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {      
	double N_i = matrix_pred(0, i);
	double N_j = matrix_pred(0, j);

	double M_ij = matrix_syst(i,j);      
	double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
	double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);

	matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
	matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;	
	matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
      }
    }

    Lee_test->matrix_absolute_flux_cov_newworld = matrix_shape;
    
  }// shape only covariance

  /////////////////////////////
  /////////////////////////////
  
  TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");
  
  TTree *tree_config = new TTree("tree", "configure information");
  int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
  int flag_syst_detector = config_Lee::flag_syst_detector;
  int flag_syst_additional = config_Lee::flag_syst_additional;
  int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
  int flag_syst_reweight = config_Lee::flag_syst_reweight;
  int flag_syst_reweight_cor = config_Lee::flag_syst_reweight_cor;
  double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  double user_scaleF_POT = scaleF_POT;
  vector<double>vc_val_GOF;
  vector<int>vc_val_GOF_NDF;
  tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
  tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
  tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
  tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
  tree_config->Branch("flag_syst_reweight", &flag_syst_reweight, "flag_syst_reweight/I" );
  tree_config->Branch("flag_syst_reweight_cor", &flag_syst_reweight_cor, "flag_syst_reweight_cor/I" );
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
  
  //file_collapsed_covariance_matrix->Close();
  
  //////////

  if( config_Lee::flag_plotting_systematics ) Lee_test->Plotting_systematics();
  
  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit
  
  //Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  //Lee_test->Set_Collapse();

  if( config_Lee::flag_GoF_output2file_default_0 ) {
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
  
  bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
  bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
  bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
  bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
  bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8
  
  ///////////////////////// gof
  
  if( flag_nueCC_FC_by_all ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 1 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 2 );
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 8 );
  }

  ///////////////////////// gof
  
  if( flag_nueCC_LowE_FC_by_all ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
  }

  
  if( 0 ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*6 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*6 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
  }

  
  if( 0 ) {// first 6 bins--> 1 bin, constrained by others
    int nbins_first = 6;
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 1 + (26-nbins_first) + 26*3 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=nbins_first; ibin++) matrix_gof_trans(ibin-1, 0) = 1;
    for( int ibin=1; ibin<=26*4+11*3-nbins_first; ibin++) matrix_gof_trans(nbins_first+ibin-1, ibin) = 1;
    
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 1, matrix_gof_trans.GetNcols()-1, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 201);
  }  
  
  if( 0 ) {// first 6 bins, constrained by others
    int nbins_first = 6;
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( nbins_first, matrix_gof_trans.GetNcols()-nbins_first, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 202);
  }

  ///////////////////////// gof
  
  if( flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, (26-8) + 26*3 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=(26-8) + 26*3 + 11*3; ibin++) matrix_gof_trans(8+ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( (26-8), matrix_gof_trans.GetNcols()-(26-8), matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 6);    
  }

  ///////////////////////// gof
  
  if( flag_nueCC_PC_by_numuCC_pi0) {    
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 5 );
  }
  
  ///////////////////////// gof
  
  if( flag_both_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 3 );
    vc_target_chs.push_back( 4 );
    
    vector<int>vc_support_chs;
    
    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 1 );
  }

  if( 0 ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 52 );// oldworld, newworld
    for( int ibin=1; ibin<=52; ibin++) matrix_gof_trans(52+ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 52, 0, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 123);

    TFile *file_numu = new TFile("file_numu.root", "recreate");
    matrix_gof_pred.Write("matrix_gof_pred");
    matrix_gof_data.Write("matrix_gof_meas");
    matrix_gof_syst.Write("matrix_gof_syst");
    file_numu->Close();
  }
  
  ///////////////////////// gof
  
  if( flag_CCpi0_FC_by_numuCC ) { 
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 5 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 2 );
  }
  
  ///////////////////////// gof
  
  if( flag_CCpi0_PC_by_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 6 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 3 );
  }
  
  ///////////////////////// gof
  
  if( flag_NCpi0_by_numuCC ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 7 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 4 );
  }

  ////////////////////////////////////////////////////////////////////

  if( 0 ) {
    vector<int>vc_target_chs;
    vc_target_chs.push_back( 1 );
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 101 );
  }

  if( 0 ) {
    vector<int>vc_target_chs;
    //vc_target_chs.push_back( 1 );
    vc_target_chs.push_back( 2 );
    
    vector<int>vc_support_chs;
    vc_support_chs.push_back( 3 );
    vc_support_chs.push_back( 4 );
    vc_support_chs.push_back( 5 );
    vc_support_chs.push_back( 6 );
    vc_support_chs.push_back( 7 );

    Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 102 );
  }

  
  if( 0 ) {
    vector<int>vc_target_chs;
    for(int idx=1; idx<=8; idx++) vc_target_chs.push_back( idx-1 );
    
    vector<int>vc_support_chs;
    for(int idx=9; idx<=Lee_test->bins_newworld; idx++) {
      if( idx>=26+1 && idx<=26+8 ) continue;
      vc_support_chs.push_back( idx-1 );
    }

    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 101 );
  }

  

  
  ////////////////////////////////////////////////////////////////////

  bool flag_publicnote = 0;
  
  if( flag_publicnote ) {
    
    ///////////////////////// gof
    
    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 3 );
      vc_target_chs.push_back( 4 ); 
      vector<int>vc_support_chs;      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 1 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 5 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 2 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 6 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 3 );
    }

    if( 1 ) {
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 7 );
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 4 );
    }

    if( 0 ) {    
      vector<int>vc_target_chs;
      vc_target_chs.push_back( 2 );
      
      vector<int>vc_support_chs;
      vc_support_chs.push_back( 3 );
      vc_support_chs.push_back( 4 );
      vc_support_chs.push_back( 5 );
      vc_support_chs.push_back( 6 );
      vc_support_chs.push_back( 7 );
      
      Lee_test->Exe_Goodness_of_fit( vc_target_chs, vc_support_chs, 5 );
    }
    
    if( 0 ) {
      TMatrixD matrix_gof_trans( Lee_test->bins_newworld, (26-8) + 26*3 + 11*3 );// oldworld, newworld
      for( int ibin=1; ibin<=(26-8) + 26*3 + 11*3; ibin++) matrix_gof_trans(8+ibin-1, ibin-1) = 1;
    
      TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
      matrix_gof_trans_T.Transpose( matrix_gof_trans );

      TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

      Lee_test->Exe_Goodness_of_fit( (26-8), matrix_gof_trans.GetNcols()-(26-8), matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 6);    
    }

    if( 0 ) {
      TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
      for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
      TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
      matrix_gof_trans_T.Transpose( matrix_gof_trans );

      TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

      Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
    }

  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////// LEE strength fitting

  if( config_Lee::flag_Lee_strength_data ) {

    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// (initial value, fix or not)

    // cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.3f, %5.3f +/- %5.3f",
    //                             Lee_test->minimization_chi2,
    //                             Lee_test->minimization_Lee_strength_val,
    //                             Lee_test->minimization_Lee_strength_err
    //                             )<<endl<<endl;

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: %6.4f,  chi2 %6.3f",                          
                                Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_chi2
                                )<<endl<<endl;    

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: %6.4f +/- %6.4f,  chi2 %6.3f",                          
                                Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err,
				Lee_test->minimization_chi2
                                )<<endl<<endl;    

    /////////////////////////////////////////
    
    double gmin = Lee_test->minimization_chi2;
    TGraph *gh_scan = new TGraph();
    double slow = 0;
    double shgh = 3;
    int nscan = 100;
    double val_max_dchi2 = 0;
    double step = (shgh-slow)/nscan;
    for(int idx=1; idx<=nscan; idx++) {
      if( idx%(max(1, nscan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./nscan, idx)<<endl;
      double val_s = slow + (idx-1)*step;
      Lee_test->Minimization_Lee_strength_FullCov(val_s, 1);// (initial value, fix or not)
      double val_chi2 = Lee_test->minimization_chi2;
      gh_scan->SetPoint( gh_scan->GetN(), val_s, val_chi2 - gmin);
      if( val_max_dchi2<val_chi2 - gmin ) val_max_dchi2 = val_chi2 - gmin;
    }
    
    double val_dchi2at0 = gh_scan->Eval(0);
    double val_dchi2at1 = gh_scan->Eval(1);
    if( fabs(val_dchi2at0)<1e-6 ) val_dchi2at0 = 0;
    
    cout<<endl<<Form(" ---> dchi2 at LEE 0/1: %7.4f %7.4f", val_dchi2at0, val_dchi2at1 )<<endl<<endl;
    
    TCanvas *canv_gh_scan = new TCanvas("canv_gh_scan", "canv_gh_scan", 900, 650);
    canv_gh_scan->SetLeftMargin(0.15); canv_gh_scan->SetRightMargin(0.1);
    canv_gh_scan->SetTopMargin(0.1); canv_gh_scan->SetBottomMargin(0.15);    
    gh_scan->Draw("al");
    gh_scan->GetXaxis()->SetTitle("LEE strength"); gh_scan->GetYaxis()->SetTitle("#Delta#chi^{2}");    
    gh_scan->GetXaxis()->SetLabelSize(0.05); gh_scan->GetXaxis()->SetTitleSize(0.05);
    gh_scan->GetYaxis()->SetLabelSize(0.05); gh_scan->GetYaxis()->SetTitleSize(0.05);
    gh_scan->GetXaxis()->CenterTitle(); gh_scan->GetYaxis()->CenterTitle();
    gh_scan->GetXaxis()->SetTitleOffset(1.2);
    gh_scan->GetYaxis()->SetRangeUser(0, val_max_dchi2*1.1);
   
    TLine *lineA_dchi2at1 = new TLine(1, 0, 1, val_dchi2at1);    
    lineA_dchi2at1->Draw("same");
    lineA_dchi2at1->SetLineWidth(2);
    lineA_dchi2at1->SetLineColor(kBlue);
    lineA_dchi2at1->SetLineStyle(7);
    TLine *lineB_dchi2at1 = new TLine(0, val_dchi2at1, 1, val_dchi2at1);    
    lineB_dchi2at1->Draw("same");
    lineB_dchi2at1->SetLineWidth(2);
    lineB_dchi2at1->SetLineColor(kBlue);
    lineB_dchi2at1->SetLineStyle(7);
    auto *tt_text_data = new TLatex( 0.2, val_dchi2at1*1.1, Form("#Delta#chi^{2} = %4.3f", val_dchi2at1) );
    tt_text_data->SetTextAlign(11); tt_text_data->SetTextSize(0.05); tt_text_data->SetTextAngle(0);
    tt_text_data->SetTextFont(42);  tt_text_data->Draw(); tt_text_data->SetTextColor(kBlue);

    canv_gh_scan->SaveAs("canv_gh_scan.png");
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// MicroBooNE suggested /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  if( config_Lee::flag_chi2_data_H0 ) {

    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting    
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();
    
    Lee_test->Minimization_Lee_strength_FullCov(0, 1);
    
    double chi2 = Lee_test->minimization_chi2;
    int ndf = Lee_test->bins_newworld;
    double p_value = TMath::Prob( chi2, ndf );
    
    cout<<Form(" ---> flag_chi2_data_H0, chi2/ndf %8.2f %3d %8.4f, p-value %f", chi2,ndf, chi2/ndf, p_value)<<endl;
  }

  
  if( config_Lee::flag_dchi2_H0toH1 ) {
    
    Lee_test->Set_measured_data();// use the measured data as the input data for the fitting
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);    
    double chi2_H0 = Lee_test->minimization_chi2;
    
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);    
    double chi2_H1 = Lee_test->minimization_chi2;

    double dchi2 = chi2_H0 - chi2_H1;
    
    int ndf = 1;
    double p_value = TMath::Prob( dchi2, ndf );    
    cout<<Form(" ---> flag_dchi2_H0toH1, chi2/ndf %8.2f %3d %8.4f, p-value %f", dchi2, ndf, dchi2/ndf, p_value)<<endl;

    p_value = TMath::Prob( -dchi2, ndf );
    cout<<Form(" ---> flag_dchi2_H1toH0, chi2/ndf %8.2f %3d %8.4f, p-value %f", -dchi2, ndf, -dchi2/ndf, p_value)<<endl;
    
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

    int N_toy = 500;
        
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
    
    file_out->cd();
    tree->Write();
    file_out->Close();
    
  }

  //////////////////////////////////////////////// Sensitivity by Asimov sample

  if( 0 ) {

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
    
    // int num_toy = 2;    
    // Lee_test->Exe_Feldman_Cousins(Lee_true_low, Lee_true_hgh, Lee_step, num_toy, ifile);

    /////////////// dchi2 of Asimov sample
    
    Lee_test->Exe_Feldman_Cousins_Asimov(Lee_true_low, Lee_true_hgh, Lee_step);

    /////////////// dchi2 of measured data
    /*
    Lee_test->Set_measured_data();    
    TMatrixD matrix_data_input_fc = Lee_test->matrix_data_newworld;    
    Lee_test->Exe_Feldman_Cousins_Data( matrix_data_input_fc, Lee_true_low, Lee_true_hgh, Lee_step );
    */
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ---> flag_syst_flux_Xs    "<<Lee_test->flag_syst_flux_Xs<<endl;
  cout<<" ---> flag_syst_detector   "<<Lee_test->flag_syst_detector<<endl;
  cout<<" ---> flag_syst_additional "<<Lee_test->flag_syst_additional<<endl;
  cout<<" ---> flag_syst_mc_stat    "<<Lee_test->flag_syst_mc_stat<<endl; 
  cout<<" ---> flag_syst_reweight       "<<Lee_test->flag_syst_reweight<<endl;
  cout<<" ---> flag_syst_reweight_cor   "<<Lee_test->flag_syst_reweight_cor<<endl; 
  
  cout<<endl<<endl;
  cout<<" ---> Finish all the program"<<endl;
  cout<<endl<<endl;
  
  if( config_Lee::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;
    
    theApp.Run();
  }
  
  return 0;
}
