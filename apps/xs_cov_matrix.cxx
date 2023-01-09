#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TVectorD.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./xf_cov_matrix -r[#sys 1-14]" << std::endl;
  }
  int run = 17; // run 1 ... xs ...
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    }
  }
  if (run !=17 && run!=0 && run!=18 && run!=19) {
    std::cout  << "Force run = 17, check configurations ..." << std::endl;
    run = 17;
  }

  TString xf_input_config_file = "../configurations/test_7ch_xs/xf_input.txt";
  if (run>17) xf_input_config_file = "../configurations/test_7ch_xs/rw_sys_input.txt";

  CovMatrix cov("../configurations/test_7ch_xs/cov_input.txt", xf_input_config_file, "../configurations/test_7ch_xs/xf_file_ch.txt", "../configurations/test_7ch_xs/rw_cv_input.txt");
  // cov.add_disabled_ch_name("BG_nueCC_FC_overlay");
  // cov.add_disabled_ch_name("BG_nueCC_PC_overlay");
  // cov.add_disabled_ch_name("BG_nueCC_FC_dirt");
  // cov.add_disabled_ch_name("BG_nueCC_PC_dirt");

  cov.add_xs_config();
  
  if(run==0) {
    cov.add_rw_config(run);
    std::cout<<"Including reweighting uncertainty in xs systematics"<<std::endl;
    run=17;
  }
  if(run==18 || run==19) {
    cov.add_rw_config(run);
    std::cout<<"Calculating reweighting uncertainty like xs systematics"<<std::endl;
  }


  cov.print_rw(cov.get_rw_info());

  // Get the file based on runno ...
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();
  // Construct the histogram ...

  // outfilename ...
  TString outfile_name;

  TH1F *htemp; TH1F *htemp1; TH1F *htemp2; TH2F *htemp3;
  
  // measure, signal, signal_nominal, R, number ...
  std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> > map_histoname_hists;
  //std::map<TString, TH1F* > map_histoname_hist;
  std::map<int, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> > map_covch_hists;
  //std::map<int, TH1F* > map_covch_hist;
  
  for (auto it = map_inputfile_info.begin(); it!=map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);

    
    if (period == run){
      outfile_name = out_filename;
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename, 0);
      
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	TString histoname = std::get<0>(*it1);
	Int_t nbin = std::get<1>(*it1);
	float llimit = std::get<2>(*it1);
	float hlimit = std::get<3>(*it1);
	TString var_name = std::get<4>(*it1);
	TString ch_name = std::get<5>(*it1);
	TString add_cut = std::get<6>(*it1);
	TString weight = std::get<7>(*it1);
	
	//	std::cout << std::get<0>(*it1)  << " " << std::get<1>(*it1) << " " << std::get<4>(*it1) << " " << std::get<5>(*it1) << " " << std::get<6>(*it1) << " " << std::get<7>(*it1) << std::endl;
	htemp = 0; htemp1 = 0; htemp2 = 0; htemp3 = 0;
	int htemp_n = 0;
	// create histograms ...
	htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
	
	bool flag_spec = cov.is_xs_chname(ch_name) ;
	//std::cout << histoname << " " << ch_name << " " << cov.is_xs_chname(ch_name) << std::endl;
	int covch = cov.get_covch_name(ch_name);
	
	if (flag_spec){
	  int nbins1 = cov.get_xs_nsignals();
	  TString temp_histoname = histoname + "_signal";
	  htemp1 = new TH1F(temp_histoname, temp_histoname,nbins1, 0.5,nbins1+0.5);
	  temp_histoname = histoname + "_signalbar";
	  htemp2 = new TH1F(temp_histoname, temp_histoname,nbins1, 0.5,nbins1+0.5);
	  temp_histoname = histoname + "_R";
	  htemp3 = new TH2F(temp_histoname, temp_histoname,nbin, llimit, hlimit,nbins1, 0.5,nbins1+0.5);
	  htemp_n = 4;
	}else{
	  htemp_n = 1;
	}
	
	map_histoname_hists[histoname] = std::make_tuple(htemp, htemp1, htemp2, htemp3, htemp_n);
	
	
	if (map_covch_hists.find(covch) == map_covch_hists.end()){
	  TH1F *htemp_1 = (TH1F*)htemp->Clone(Form("pred_covch_%d",covch));
	  if (flag_spec){
	    TH1F *htemp1_1 = (TH1F*)htemp1->Clone(Form("pred_covch_%d_signal",covch));
	    TH1F *htemp2_1 = (TH1F*)htemp2->Clone(Form("pred_covch_%d_signalbar",covch));
	    TH2F *htemp3_1 = (TH2F*)htemp3->Clone(Form("pred_covch_%d_R",covch));
	    map_covch_hists[covch] = std::make_tuple(htemp_1, htemp1_1, htemp2_1, htemp3_1, htemp_n);
	  }else{
	    map_covch_hists[covch] = std::make_tuple(htemp_1, htemp1, htemp2, htemp3, htemp_n);
	  }
	}
      }
      //  std::cout << input_filename << " " << filetype << " " << out_filename << std::endl; 
    }
  }
  
  // hack ...
  outfile_name = "../hist_rootfiles/XsFlux/cov_xs.root";
  if(run==18) outfile_name = "../hist_rootfiles/XsFlux/cov_rw.root";
  if(run==19) outfile_name = "../hist_rootfiles/XsFlux/cov_rw_cor.root";
  std::cout << outfile_name << std::endl;
  TMatrixD* cov_add_mat = cov.get_add_cov_matrix();
  // create a covariance matrix for bootstrapping ...
  TMatrixD* cov_xs_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  TVectorD* vec_mean = new TVectorD(cov_add_mat->GetNrows());

  //std::cout << cov.get_xs_nmeas() << " " << cov.get_xs_nsignals() << std::endl;
  // matrix and R ...
  TVectorD* vec_signal = new TVectorD(cov.get_xs_nsignals());

  TMatrixD *mat_R = new TMatrixD(cov_add_mat->GetNrows(),cov.get_xs_nsignals());

  
  cov.gen_xs_cov_matrix(run, map_covch_hists, map_histoname_hists, vec_mean, cov_xs_mat, vec_signal, mat_R);

  TMatrixD* frac_cov_xs_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  
  for (size_t i=0; i!= frac_cov_xs_mat->GetNrows(); i++){
    double val_1 = (*vec_mean)(i);
    for (size_t j=0; j!=frac_cov_xs_mat->GetNrows();j++){
      double val_2 = (*vec_mean)(j);
      double val = (*cov_xs_mat)(i,j);
      if (val_1 ==0 && val_2 == 0){
	(*frac_cov_xs_mat)(i,j) = 0;
      }else if (val_1 ==0 || val_2 ==0){
	if (val !=0){
	  if (i==j){
	    (*frac_cov_xs_mat)(i,j) = 0.;
	  }else{
	    (*frac_cov_xs_mat)(i,j) = 0;
	  }
	}else{
	  (*frac_cov_xs_mat)(i,j) = 0;
	}
      }else{
	(*frac_cov_xs_mat)(i,j) = val/val_1/val_2;
      }
    }
  }
  
  TFile *file = new TFile(outfile_name,"RECREATE");
  vec_mean->Write(Form("vec_mean_%d",run));
  cov_xs_mat->Write(Form("cov_xf_mat_%d",run));
  frac_cov_xs_mat->Write(Form("frac_cov_xf_mat_%d",run));
  vec_signal->Write(Form("vec_signal_%d",run));
  mat_R->Write(Form("mat_R_%d",run));
  
  // save central ... results ...
  // for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
  //  ((TH1F*)it->second)->SetDirectory(file);
  // }
  
  for (auto it = map_covch_hists.begin(); it != map_covch_hists.end(); it++){
    auto results = it->second;
    TH1F *h1 = std::get<0>(results);
    TH1F *h2 = std::get<1>(results);
    TH1F *h3 = std::get<2>(results);
    TH2F *h4 = std::get<3>(results);
    int num = std::get<4>(results);
    if (num == 1){
      h1->SetDirectory(file);
    }else{
      h1->SetDirectory(file);
      h2->SetDirectory(file);
      h3->SetDirectory(file);
      h4->SetDirectory(file);
    }
  }
  
  file->Write();
  file->Close();
  
  return 0;
}
