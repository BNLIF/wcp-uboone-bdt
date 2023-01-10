#ifndef LEEANAC_MASTER_COV_MATRIX
#define LEEANAC_MASTER_COV_MATRIX

#include "TFile.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include <map>
#include <set>

#include <iostream>
#include <fstream>

// #include "WCPLEEANA/eval.h"
// #include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/weights.h"
#include "WCPLEEANA/pfeval.h"

namespace LEEana{
  
  class JointXsecHelper{
    public: 

      JointXsecHelper(TString xs_ch_filename = "");
     
      int chWgtbyHistName(TString histname) const;
                                        
      std::set<TString> getXsNames() { return fXsNames; }
      std::map<TString, int> getXsInfo() { return fXsInfo; }

    private:
      std::map<TString, int> fXsInfo;
      std::set<TString> fXsNames;
  };

  class CovMatrix{
  public:
    // CovMatrix(TString cov_filename = "./configurations/cov_input.txt", TString cv_filename = "./configurations/cv_input.txt", TString file_filename = "./configurations/file_ch.txt", TString rw_filename = "./configurations/rw_cv_input.txt");
    CovMatrix(TString cov_filename = "../configurations/test_7ch_xs/cov_input.txt", TString cv_filename = "../configurations/test_7ch_xs/cv_input.txt", TString file_filename = "../configurations/test_7ch_xs/file_ch.txt", TString rw_filename = "../configurations/test_7ch_xs/rw_cv_input.txt");
    ~CovMatrix();

    // void add_xs_config(TString xs_ch_filename = "../configurations/xs_ch.txt", TString xs_real_bin_filename = "../configurations/xs_real_bin.txt");
    void add_xs_config(TString xs_ch_filename = "../configurations/test_7ch_xs/xs_ch.txt", TString xs_real_bin_filename = "../configurations/test_7ch_xs/xs_real_bin.txt");
    void add_osc_config(TString osc_ch_filename = "./configurations/osc_ch.txt", TString osc_pars_filename = "./configurations/osc_parameter.txt");
    
    void add_rw_config(int run=1);   
    void print_rw(std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > rw_info);
 
    void print_ch_info();
    void print_filetype_info();
    void print_systematics();
    void print_matrix();
    
    void print_cvfile_info();

    // histogram ...
    TString get_ch_name(int ch);
    TString get_ch_var(int ch);
    int get_ch(TString name);
    
    std::tuple<int, float, float> get_ch_hist(int ch);

    // ... filetype related ...
    int get_ch_filetype(int ch);
    std::vector<int> get_filetype_chs(int filetype);

    //
    std::set<int> get_xfs_filetypes(){return xfs_filetypes;};
    std::set<int> get_det_filetypes(){return det_filetypes;};
    std::set<int> get_add_filetypes(){return add_filetypes;};

    bool get_sys_xs_flux(int ch);
    bool get_sys_det(int ch);
    std::pair<bool, float> get_sys_add(int ch);
    int get_sys_mc_same(int ch);
    std::map<int, std::vector<int> > get_mcstat_same_chs(){return map_mcstat_same_chs;};

    int get_obsch(int ch);
    int get_covch(int ch);
    int get_obsch_fcov(int covch);
    std::pair<int, int> get_obsch_info(int obsch);
    std::pair<int, int> get_covch_info(int covch);

    int get_covch_startbin(int covch){return map_covch_startbin[covch];}; 
    
    TMatrixD* get_mat_collapse(){return mat_collapse;};
    TMatrixD* get_add_cov_matrix(){return mat_add_cov;};
    
    float get_ext_pot(TString filename);
    std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > get_histograms(TString filename, int flag = 0);
    std::map<TString, std::tuple<int, int, TString, float, int, double, int> > get_map_inputfile_info(){return map_inputfile_info;};
    
    std::map<TString, std::pair<TString, int> > get_map_pred_histo_histo_err2_lee(){return map_pred_histo_histo_err2_lee;};

    // Now the cross uncertainty term
    std::map<std::pair<TString, TString>, std::pair<TString, int> > get_map_pair_hist_hists_cros(){return map_pair_histo_histos_cros;};

    std::map<int, std::set<std::set<std::pair<TString, int> > > > get_map_pred_obsch_histos(){return map_pred_obsch_histos;};

    int get_obsch_name(TString name);
    int get_covch_name(TString name);

    void fill_data_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram);
    void fill_pred_histograms(int run, std::map<int, std::vector<TH1F*> >& map_obsch_histos, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > >& map_obsch_bayes, std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > >& map_obsch_infos, std::map<TString, std::pair<TH1F*, double> >& map_name_histogram, float lee_strength, std::map<int, double> map_data_period_pot, bool flag_breakdown, std::map<int, std::vector<TH1F*> >& map_obsch_subhistos);

    void fill_pred_R_signal(int run, TMatrixD* mat_R, TVectorD* vec_signal,  std::map<int, double>& map_data_period_pot, std::map<TString, std::tuple<TH1F*, TH2F*, double> >& map_name_xs_hists);
    
    // Xs related 
    void gen_xs_cov_matrix(int run, std::map<int, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_covch_hists, std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists, TVectorD* vec_mean,  TMatrixD* cov_xf_mat, TVectorD* vec_signal, TMatrixD* mat_R);
    std::pair<std::vector<int>, std::vector<int> > get_events_weights_xs(TString input_filename, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos);

    void fill_xs_histograms(int num, int tot_num, int acc_no, int no, int tot_no, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists);
    void fill_xs_histograms(std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists);

    
    // XsFlux related
    void gen_xf_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean,  TMatrixD* cov_xf_mat);

    std::pair<std::vector<int>, std::vector<int> > get_events_weights(TString input_filename, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos);

    void fill_xf_histograms(int num, int tot_num, int acc_no, int no, int tot_no, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
     void fill_xf_histograms(std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);

     // prediction statistics
     void gen_pred_stat_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TMatrixD* cov_mat);
     void get_pred_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos);
     void fill_pred_stat_histograms(std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
     void fill_pred_stat_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist, std::map<TString, std::set<std::pair<TString, double> > >& map_filename_histoname_ratio);
     
     // data statistics
     void gen_data_stat_cov_matrix(int run, std::map<int, TH1F*>& map_obsch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TMatrixD* cov_mat_bootstrapping);
     void get_data_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int,  std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos);
     void fill_data_stat_histograms(std::map<TString, std::vector< std::tuple<int, int, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
     void fill_data_stat_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
     
     // detector ...
    void gen_det_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TVectorD* vec_mean_diff, TMatrixD* cov_mat_bootstrapping, TMatrixD* cov_det_mat, int flag_gp);
    void get_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > >&map_all_events,  std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos);

    void fill_det_histograms(std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
    void fill_det_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist);
    
    std::pair<double,double> get_bayes_errors(double num);

    void add_disabled_ch_name(TString name);
    std::set<TString>& get_disabled_ch_names(){return disabled_ch_names;};
    void remove_disabled_ch_name(TString name);
    // xs related
    int get_xs_nsignals();
    int get_xs_nmeas();
    bool is_xs_chname(TString name);
    int get_cut_file(){return cut_file;};

    std::map<TString, int>& get_map_cut_xs_bin(){return map_cut_xs_bin;};

    void init_spec_weights(int num, int num1 = 1000, double strength = 0.1);
    std::vector<float> get_spec_weight(LEEana::EvalInfo& eval, LEEana::PFevalInfo& pfeval);

    bool get_osc_flag(){return flag_osc;};
    bool is_osc_channel(TString ch_name);
    double get_osc_weight(EvalInfo& eval, PFevalInfo& pfeval);


    bool get_flag_reweight(){return flag_reweight;};
    std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > get_rw_info(bool ov=false){
      if(!(ov)) return rw_info;
      std::get<0>(rw_info) = 1;
      for(size_t rw=0; rw<std::get<1>(rw_info).size(); rw++){ std::get<0>( std::get<1>(rw_info)[rw]  ) = 1; }
      return rw_info;
    };

  private:
    TGraph *gl, *gh;
    int g_llimit, g_hlimit;
    
    TMatrixD* mat_collapse;
    TMatrixD* mat_add_cov;
    
    // basic information about the channels
    // name, var_name, bin, llmit, hlimit, weight, obs_no, lee
    std::map<int, std::tuple<TString, TString, int, float, float, TString, int, int> > map_ch_hist;
    std::map<TString, int> map_name_ch;
    
    // information regarding ch and their filetype
    std::map<int, int> map_ch_filetype;
    std::map<int, std::vector<int> > map_filetype_chs;

    // information regarding systematics (xs_flux, det, additional, mc stat...)
    // xs_flux, det, additional relative uncertainties, mc_stat
    std::map<int, std::tuple<int, int, float, int> > map_ch_systematics;
    std::set<int> xfs_filetypes;
    std::set<int> det_filetypes;
    std::set<int> add_filetypes;

    // MC statistics ... same channels  (diagnal to start with ...)
    std::map<int, std::vector<int> > map_mcstat_same_chs;
    std::map<int, std::set<int> > map_filetype_mcstats;
    
    // prepare covariance matrix structure ...
    std::map<int, int> map_ch_obsch;
    std::map<int, int> map_ch_covch;

    std::map<int, int> map_covch_obsch; // map ...
    
    std::map<int, int> map_obsch_nbin; // record the bin number + 1
    std::map<int, int> map_covch_nbin; // record the bin number + 1
    
    // covariance matrix internally ...
    std::map<int, int> map_covch_startbin; 
    std::map<int, int> map_obsch_startbin;
    
    // covariance matrix in observation
    std::map<int, int> map_covchbin_obschbin;

    // CV related input ...
    std::map<int, TString> map_filetype_name;
    std::map<int, std::vector<TString> > map_filetype_inputfiles;
    std::map<TString, int> map_inputfile_filetype;
    
    // filetype, period, outfile_name, external pot if any, file_no, norm_pot
    std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info;
    std::map<int, int> map_fileno_period;
    std::map<TString, std::vector<TString> > map_inputfile_cuts;

    // histogram infos (name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms;
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms_err2;
    std::map<TString, std::vector< std::tuple<TString, int, float, float, TString, TString, TString, TString > > > map_inputfile_histograms_cros;

    std::map<TString, TString> map_histogram_inputfile;
    std::map<TString, std::pair<int, double> > map_histogram_covch_add;
    
    // structure of summing histograms together for prediction ...
    std::map<int, std::set<int> > map_pred_obsch_covch; // OK
    std::map<int, std::set<int> > map_pred_covch_ch; // OK
    
    std::map<int, std::set<std::pair<TString, TString> > > map_pred_ch_subch; // OK
    std::map<std::pair<TString, TString> , std::set<std::pair<TString, int> > > map_pred_subch_histos; //OK
    
    
    
    std::map<TString, std::pair<TString, int> > map_pred_histo_histo_err2_lee; //OK

    // Now the cross uncertainty term
    std::map<std::pair<TString, TString>, std::pair<TString, int> > map_pair_histo_histos_cros; // OK

    // obsch --> subchannel --> histos
    std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_obsch_histos; // total ...
    // covch --> subchannel --> histos 
    std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_covch_histos; // total ...

    std::set<TString> disabled_ch_names;
    
    // Xs related
    int cut_file;
    
    JointXsecHelper xsechelper;
    std::set<TString> xs_signal_ch_names;
    std::map<TString, int> xs_signal_ch_info;
    
    std::map<int, TString> map_xs_bin_cut;
    std::map<TString, int> map_cut_xs_bin;
    std::map<int, double> map_xs_bin_constant;
    std::map<int, std::pair<double, double> > map_xs_bin_errs;

    // Osc related
    bool flag_osc;
    std::set<TString> osc_signal_ch_names;
    double osc_par_delta_m2_eV2;
    double osc_par_sin22theta_ee;
    double osc_par_delta_m2_41_eV2;
    double osc_par_sin2_2theta_14; // = ee
    double osc_par_sin2_theta_24;
    double osc_par_sin2_theta_34;

    //reweighting related
    bool flag_reweight = false;
    int rw_type = 0;
    std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > >  rw_info;  
 
    // special weights ...
    bool flag_spec_weights;
    std::vector<std::vector<float> > spec_weights;
    
  };
}

#endif
