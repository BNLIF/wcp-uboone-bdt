#pragma once

namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  /* TString spectra_file = "./data_framework_Doc33131/merge_all.root"; */
  /* TString flux_Xs_directory = "./data_framework_Doc33131/flux_Xs/"; */
  /* TString detector_directory = "./data_framework_Doc33131/det_both/"; */
  /* TString mc_directory = "./data_framework_Doc33131/mc_stat/"; */

  //TString detector_directory = "/home/xji/data0/work/505_TLee/TLee_for_data_5e19/data_framework/det_11stat/";
  //TString detector_directory = "/home/xji/data0/work/505_TLee/TLee_for_data_5e19/data_framework/det_norandom/";  
  
  /* TString spectra_file = "./new_new_TLee_input_fakeset5_myself/merge.root"; */
  /* TString flux_Xs_directory = "./new_new_TLee_input_fakeset5_myself/flux_Xs/"; */
  /* TString detector_directory = "./new_new_TLee_input_fakeset5_myself/det/"; */
  /* TString mc_directory = "./new_new_TLee_input_fakeset5_myself/mc_stat/"; */

  /* TString spectra_file = "./new_new_TLee_input_fakeset5_1mu0p1muNp/merge.root"; */
  /* TString flux_Xs_directory = "./new_new_TLee_input_fakeset5_1mu0p1muNp/flux_Xs/"; */
  /* TString detector_directory = "./new_new_TLee_input_fakeset5_1mu0p1muNp/det/"; */
  /* TString mc_directory = "./new_new_TLee_input_fakeset5_1mu0p1muNp/mc_stat/"; */

  /* TString spectra_file = "./new_new_TLee_input_fakeset7/merge.root"; */
  /* TString flux_Xs_directory = "./new_new_TLee_input_fakeset7/flux_Xs/"; */
  /* TString detector_directory = "./new_new_TLee_input_fakeset7/det/"; */
  /* TString mc_directory = "./new_new_TLee_input_fakeset7/mc_stat/"; */

  /* TString spectra_file = "./new_new_TLee_input_normal_highstatDetVar/merge.root"; */
  /* TString flux_Xs_directory = "./new_new_TLee_input_normal_highstatDetVar/flux_Xs/"; */
  /* TString detector_directory = "./new_new_TLee_input_normal_highstatDetVar/det/"; */
  /* TString mc_directory = "./new_new_TLee_input_normal_highstatDetVar/mc_stat/"; */
          
  /* TString spectra_file = "./TLee_input_NumiReinteractionIssue_cutbased/merge.root"; */
  /* TString flux_Xs_directory = "./TLee_input_NumiReinteractionIssue_cutbased/flux_Xs/"; */
  /* TString detector_directory = "./TLee_input_NumiReinteractionIssue_cutbased/det/"; */
  /* TString mc_directory = "./TLee_input_NumiReinteractionIssue_cutbased/mc_stat/"; */
  
  /* TString spectra_file = "./numi_result_syst_0_0/merge.root"; */
  /* TString flux_Xs_directory = "./numi_result_syst_0_0/XsFlux/"; */
  /* TString detector_directory = "./TLee_input_Mar25_final_16chs_standard/det/"; */
  /* TString mc_directory = "./numi_result_syst_0_0/mc_stat/"; */  
  
  /* TString spectra_file = "./TLee_input_eLEEnote1/merge.root"; */
  /* TString flux_Xs_directory = "./TLee_input_eLEEnote1/flux_Xs/"; */
  /* TString detector_directory = "./TLee_input_eLEEnote1/det/"; */
  /* TString mc_directory = "./TLee_input_eLEEnote1/mc_stat/"; */
  
  /*
  TString spectra_file = "./numi_summary_noosc/merge.root";
  TString flux_Xs_directory = "./numi_summary_noosc/XsFlux/";
  TString detector_directory = "./numi_summary_noosc/det/";
  TString mc_directory = "./numi_summary_noosc/mcstat/";
  */
  /*
  TString spectra_file = "./numi_summary_osc/merge.root";
  TString flux_Xs_directory = "./numi_summary_osc/XsFlux/";
  TString detector_directory = "./numi_summary_osc/det/";
  TString mc_directory = "./numi_summary_osc/mcstat/";
  */


  
  TString spectra_file = "./new_TLee_input_opendata5e19/merge.root";
  TString flux_Xs_directory = "./new_TLee_input_opendata5e19/flux_Xs/";
  TString detector_directory = "./new_TLee_input_opendata5e19/det/";
  TString mc_directory = "./new_TLee_input_opendata5e19/mc_stat/";


  
  int channels_observation = 7;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse

  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 19;//cov_18.root is uncorrelated reweighting and cov_19.root is correlated
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 99;
   

  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
 
  
  ////////// display graphics flag

  bool flag_display_graphics = 1;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;
  bool flag_syst_reweight        = 0;
  bool flag_syst_reweight_cor    = 0;


  double Lee_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 0;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF         = 0;

  bool flag_GoF_output2file_default_0 = 1;

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

  ////////// MicroBooNE suggested

  bool flag_chi2_data_H0 = 0;
  bool flag_dchi2_H0toH1 = 0;
  
  ////////// Advanced tools
  
  ///// void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
  ///// do the fitting on the spectra and cov_total after constraint ?
  bool flag_Lee_minimization_after_constraint = 0;// hardcoded, only for the standard 7-ch fitting
  
}
