#pragma once

namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  /// paper used
  // TString spectra_file = "./TLee_input_normal_7chs_Jun29_nearside/merge.root";
  // TString flux_Xs_directory = "./TLee_input_normal_7chs_Jun29_nearside/flux_Xs/";
  // TString detector_directory = "./TLee_input_normal_7chs_Jun29_nearside/det/";
  // TString mc_directory = "./TLee_input_normal_7chs_Jun29_nearside/mc_stat/";

  TString spectra_file = "./new_TLee_input_opendata5e19/merge.root";
  TString flux_Xs_directory = "./new_TLee_input_opendata5e19/flux_Xs/";
  TString detector_directory = "./new_TLee_input_opendata5e19/det/";
  TString mc_directory = "./new_TLee_input_opendata5e19/mc_stat/";

  ///
  int array_LEE_ch[4] = {8,9,0,0};// element value "0" will note be set to the LEE_ch
  //int array_LEE_ch[4] = {12,13,14,15};// element value "0" will note be set to the LEE_ch
  //int array_LEE_ch[6] = {1,2,3,4,5,6};// element value "0" will note be set to the LEE_ch
  //int array_LEE_ch[8] = {2,3,5,6,7,8,0,0};// element value "0" will note be set to the LEE_ch

  int channels_observation = 7;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse
                               // NOTE: This value is not used in the lastest version
  
  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 98;
   
  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
   
  ////////// display graphics flag

  bool flag_display_graphics = 0;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;

  double Lee_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 0;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF         = 0;

  bool flag_GoF_output2file_default_0 = 0;

  bool flag_lookelsewhere             = 0;
  
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;
  bool flag_Lee_scan_data     = 0;

  bool flag_GOF = 0;
}
