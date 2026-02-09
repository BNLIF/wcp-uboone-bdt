#include <iostream>
#include <sstream>
#include <memory>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"
#include "WCPLEEANA/bdt.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TPDF.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./merge_hist_xs -r[#run, 0 for all] -e[1 for standard, 2 for Bayesian] -L[LEE strength]" << std::endl;
  }
  
  int run = 1; // run 1 ...
  int flag_err = 1;// 1 for standard, 2 for Bayesian ...
  float lee_strength = 0; // no LEE strength ...
  int flag_display = 0;
  int flag_breakdown = 0;
  int flag_bdt = 0;
  TString bdt_varname = "";

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 'e':
      flag_err = atoi(&argv[i][2]); // error for plotting
      break;
    case 'l':
      lee_strength = atof(&argv[i][2]);
      break;
    case 'd':
      flag_display = atoi(&argv[i][2]);
      break;
    case 's':
      flag_breakdown = atoi(&argv[i][2]);
      break;
    case 'b':
      flag_bdt = 1;
      bdt_varname = &argv[i][2];  // load bdt model
      break;
    }
  }

  CovMatrix cov;
  cov.add_xs_config();
  std::shared_ptr<TMVA::Reader> reader = 0;
  if(flag_bdt){
    reader = std::shared_ptr<TMVA::Reader>(fetch_bdtreader(bdt_varname));
    cov.set_bdt_reader(reader);
  }

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();

  TFile *temp_file;
  TH1F *htemp = 0; TH1F *htemp1 =0; TH2F *htemp2 = 0; int num = 0;
  TTree *T;
  Double_t pot;
  std::map<TString, std::pair<TH1F*, double> > map_name_histogram;

  std::map<TString, std::tuple<TH1F*, TH2F*, double> > map_name_xs_hists;

  // data POT ...
  std::map<int, double> map_data_period_pot;
  //  std::vector<TH1F*> temp_histograms;
  
  // open all the histograms ...
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    temp_file = new TFile(out_filename);
    T = (TTree*)temp_file->Get("T");
    T->SetBranchAddress("pot",&pot);
    T->GetEntry(0);

    if (filetype==5 || filetype == 15){
      map_data_period_pot[period] = pot;
    }
    
    // std::cout << pot << std::endl;
    
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
    std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
    std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
    std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
    
    for (size_t i=0;i!=all_histo_infos.size();i++){
      htemp = (TH1F*)temp_file->Get(std::get<0>(all_histo_infos.at(i)));
      TString temp_name = std::get<0>(all_histo_infos.at(i));
      TString temp_name1 = temp_name + "_signal";
      TString temp_name2 = temp_name + "_R";
      TH1F *htemp1 = (TH1F*)temp_file->Get(temp_name1);
      TH2F *htemp2 = (TH2F*)temp_file->Get(temp_name2);
      // std::cout << temp_name << " " << htemp1 << " " << htemp2 << std::endl;
      //      std::cout << std::get<0>(all_histo_infos.at(i)) << " " << htemp->GetSum() << std::endl;
      //      temp_histograms.push_back(htemp);
      map_name_histogram[std::get<0>(all_histo_infos.at(i))] = std::make_pair(htemp, pot);
      if (htemp1 !=0 && htemp2 !=0){
	map_name_xs_hists[std::get<0>(all_histo_infos.at(i))] = std::make_tuple(htemp1, htemp2, pot);
      }
    }
  }

  
 
  
  // create histograms for data, create histograms for predictions
  // obsch --> histograms (data, prediction, prediction_err2
  std::map<int, std::vector<TH1F*> > map_obsch_histos;
  // obsch --> break down histograms (truth label: add_cut) prediction 
  std::map<int, std::vector<TH1F*> > map_obsch_subhistos;
  // Bayesian error needed ...
  // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_bayes;
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_infos;
  
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (filetype == 5 || filetype == 15){
      // name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	int obsch = cov.get_obsch_name(std::get<5>(*it1));	
	htemp = map_name_histogram[std::get<0>(*it1)].first;

	std::vector<TH1F*> vec_histos;
	
	TH1F *hdata = (TH1F*)htemp->Clone(Form("data_%d",obsch));
	hdata->Reset();
	TH1F *hpred = (TH1F*)htemp->Clone(Form("pred_%d",obsch));
	hpred->Reset();
	TH1F *hpred_err2 = (TH1F*)htemp->Clone(Form("pred_err2_%d",obsch));
	hpred_err2->Reset();

	vec_histos.push_back(hdata);
	vec_histos.push_back(hpred);
	vec_histos.push_back(hpred_err2);

	// get histograms ...
	map_obsch_histos[obsch] = vec_histos;
	//map_obsch_bayes[obsch].resize(htemp->GetNbinsX()+1);
	// for (Int_t i=0;i!=htemp->GetNbinsX()+1;i++){
	//   std::vector< std::tuple<double, double, double> > temp;
	//   map_obsch_bayes[obsch].push_back(temp);
	// }
	
	
	//std::cout << std::get<5>(*it1) << " " << obsch << " " << htemp->GetSum() << std::endl;
      }
      
      //      break;
    }
  }
  
  // get data histograms ...
  cov.fill_data_histograms(run, map_obsch_histos, map_name_histogram);
  
  // get predictions and its uncertainties ...,
  cov.fill_pred_histograms(run, map_obsch_histos, map_obsch_bayes, map_obsch_infos, map_name_histogram, lee_strength, map_data_period_pot, flag_breakdown, map_obsch_subhistos);

  
  // matrix and R ...
  TVectorD* vec_signal = new TVectorD(cov.get_xs_nsignals());
  // additional covariance matrix ...
  TMatrixD* mat_add_cov = cov.get_add_cov_matrix();
  TMatrixD *mat_R = new TMatrixD(mat_add_cov->GetNrows(),cov.get_xs_nsignals());
  cov.fill_pred_R_signal(run, mat_R, vec_signal,  map_data_period_pot,  map_name_xs_hists);
  
  
  // get Bayesian errrors ...
  if (flag_err==2){
    std::cout << lee_strength << " " << run << std::endl;
    
    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1);
      TH1F *h2 = it->second.at(2);
      int obsch = it->first;

      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = map_obsch_bayes[obsch];

      //std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << h1->GetNbinsX() << std::endl;

      //if (obsch !=1) continue;
      
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
	Bayes bayes;
	//	if (i!=0) continue;
	//double temp = 0, temp1=0;
	double zero_weight = 0, zero_count = 0;
	double nonzero_meas = 0, nonzero_sigma2 = 0, nonzero_weight = 0;
	for (auto it1 = bayes_inputs.begin(); it1!=bayes_inputs.end(); it1++){
	  //bayes.add_meas_component(std::get<0>((*it1).at(i)), std::get<1>((*it1).at(i)), std::get<2>((*it1).at(i)));
	  // temp += std::get<0>((*it1).at(i));
	  // temp1 += std::get<1>((*it1).at(i));
	  //std::cout << i << " " << std::get<0>((*it1).at(i)) << " " << std::get<1>((*it1).at(i)) << " " << std::get<2>((*it1).at(i)) << " " << std::endl;
	  // approximation for zeros ...
	  if (std::get<0>((*it1).at(i)) == 0){
	    zero_weight += pow(std::get<2>((*it1).at(i)),2);
	    zero_count += std::get<2>((*it1).at(i));
	  }else{
	    nonzero_meas += std::get<0>((*it1).at(i));
	    nonzero_sigma2 += std::get<1>((*it1).at(i));
	    nonzero_weight += std::get<2>((*it1).at(i));
	  }
	}
	if (zero_count != 0)	bayes.add_meas_component(0,0,zero_weight/zero_count);
	if (nonzero_meas != 0)  bayes.add_meas_component(nonzero_meas, nonzero_sigma2, nonzero_weight);
	
	bayes.do_convolution();
	
	double cov = bayes.get_covariance();
	//double cov1 = bayes.get_covariance_mc();
	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1) << std::endl;
	// std::cout << temp << " " << temp1 << " " << h1->GetBinContent(i+1) << " " << h2->GetBinContent(i+1) << std::endl;
        if(isnan(cov) || isinf(cov)) {
            cov = bayes.get_covariance_mc();
            //cov = h1->SetBinError(i+1, h1->GetBinError(i));
        }
	h1->SetBinError(i+1,sqrt(cov));
      }
      // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
      //std::map<int, std::vector< std::vector< std::tuple<double, double, double> > > > map_obsch_bayes;
    }
  }else if (flag_err == 1){
    // prediction ...
    TFile *file3 = new TFile("merge_xs.root","RECREATE");
    
    TMatrixD* mat_collapse = cov.get_mat_collapse();
    mat_collapse->Write("mat_collapse");
    
    
    
    std::map<int, TH1F*> map_covch_histo;
    
    for (auto it = map_obsch_infos.begin(); it != map_obsch_infos.end(); it++){
      int obsch = it->first;
      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = it->second;
      TH1F *htemp = (TH1F*)map_obsch_histos[obsch].at(0); // data histogram ...
      
      for (size_t i=0;i!=bayes_inputs.size(); i++){
	int covch = std::get<3>(bayes_inputs.at(i).front());

	
	// double add_sys = std::get<4>(bayes_inputs.at(i).front());
	int start_bin = cov.get_covch_startbin(covch);
	
	if (map_covch_histo.find(covch) == map_covch_histo.end()){
	  TH1F *hnew = (TH1F*)htemp->Clone(Form("histo_%d",covch));
	  hnew->Reset();
	  map_covch_histo[covch] = hnew;
	}
	TH1F *htemp1 = map_covch_histo[covch];

	//std::cout << covch << " " << obsch << " " << i << " " << htemp1->GetSum() << std::endl;
	for (size_t j=0;j!=bayes_inputs.at(i).size();j++){
	  htemp1->SetBinContent(j+1, htemp1->GetBinContent(j+1) + std::get<0>(bayes_inputs.at(i).at(j)));
	  htemp1->SetBinError(j+1, sqrt(pow(htemp1->GetBinError(j+1),2) + std::get<1>(bayes_inputs.at(i).at(j))));
	  
	  (*mat_add_cov)(start_bin + j, start_bin + j) += std::get<4>(bayes_inputs.at(i).at(j));
	  
	}
	
      }
    }
    
    for (auto it = map_obsch_histos.begin(); it != map_obsch_histos.end(); it++){
      int obsch = it->first;
      TH1F *hdata = it->second.at(0);
      TH1F *hmc = it->second.at(1);
      
      hdata->SetName(Form("hdata_obsch_%d",obsch));
      hmc->SetName(Form("hmc_obsch_%d",obsch));
      hdata->SetDirectory(file3);
      hmc->SetDirectory(file3);
    }
    
    mat_add_cov->Write("cov_mat_add");
    vec_signal->Write("vec_signal");
    mat_R->Write("mat_R");
    
    file3->Write();
    file3->Close();
  }

}
