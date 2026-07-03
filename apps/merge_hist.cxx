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
    std::cout << "./merge_hist -r[#run, 0 for all] -e[1 for standard, 2 for Bayesian] -L[LEE strength]" << std::endl;
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
  std::shared_ptr<TMVA::Reader> reader = 0;
  if(flag_bdt){
    reader = std::shared_ptr<TMVA::Reader>(fetch_bdtreader(bdt_varname));
    cov.set_bdt_reader(reader);
  }

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();

  
  TFile *temp_file;
  TH1F *htemp;
  TTree *T;
  Double_t pot;
  std::map<TString, std::pair<TH1F*, double> > map_name_histogram;

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

    if (filetype==150 || filetype == 15){
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
      
      //std::cout << out_filename << " " << std::get<0>(all_histo_infos.at(i)) << " " << htemp << std::endl;
      //if (htemp == 0) continue;
      //      temp_histograms.push_back(htemp);
      map_name_histogram[std::get<0>(all_histo_infos.at(i))] = std::make_pair(htemp, pot);
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

    if (filetype ==150 || filetype == 15){
      // name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	int obsch = cov.get_obsch_name(std::get<5>(*it1));	
	htemp = map_name_histogram[std::get<0>(*it1)].first;
	//	if (htemp ==0 ) continue;
	//	std::cout << filetype << " " << obsch << " " << htemp << std::endl;
	
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

    /* for (auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){ */
    /*         for(size_t i=0; i<it->second.size(); i++){ */
    /*             std::cout<<"DEBUG2: "<<it->first<<": "<<map_obsch_subhistos[it->first].at(i)->GetName()<<std::endl; */
    /*         } */
    /* } */ 

  
  
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
	//std::cout << bayes_inputs.size() << std::endl;

	double zero_weight = 0, zero_count = 0;
	double nonzero_meas = 0, nonzero_sigma2 = 0, nonzero_weight = 0;
	for (auto it1 = bayes_inputs.begin(); it1!=bayes_inputs.end(); it1++){
	  // approximation for zeros ...
	  if (std::get<0>((*it1).at(i)) == 0){
	    zero_weight += pow(std::get<2>((*it1).at(i)),2);
	    zero_count += std::get<2>((*it1).at(i));
	  }else{
	    nonzero_meas += std::get<0>((*it1).at(i));
	    nonzero_sigma2 += std::get<1>((*it1).at(i));
	    nonzero_weight += std::get<2>((*it1).at(i));
	  }
	  // bayes.add_meas_component(std::get<0>((*it1).at(i)), std::get<1>((*it1).at(i)), std::get<2>((*it1).at(i)));
	  // temp += std::get<0>((*it1).at(i));
	  // temp1 += std::get<1>((*it1).at(i));
	  // std::cout << i << " " << std::get<0>((*it1).at(i)) << " " << std::get<1>((*it1).at(i)) << " " << std::get<2>((*it1).at(i)) << " " << std::endl;
	}
	if (zero_count != 0)	bayes.add_meas_component(0,0,zero_weight/zero_count);
	if (nonzero_meas != 0)  bayes.add_meas_component(nonzero_meas, nonzero_sigma2, nonzero_weight);
	
	bayes.do_convolution();
	
	double cov = bayes.get_covariance();
	//	double cov1 = bayes.get_covariance_mc(); //cov = cov1;
	if(isnan(cov) || isinf(cov)) {
            cov = bayes.get_covariance_mc();
            //cov = h1->SetBinError(i+1, h1->GetBinError(i));
        }
	//	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1) << " " << cov1 << std::endl;

	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1)  << std::endl;

	
	// std::cout << temp << " " << temp1 << " " << h1->GetBinContent(i+1) << " " << h2->GetBinContent(i+1) << std::endl;
       
	h1->SetBinError(i+1,sqrt(cov));
      }
      // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2 
      //std::map<int, std::vector< std::vector< std::tuple<double, double, double> > > > map_obsch_bayes;
    }
  }

  if (flag_display == 1){
    // plotting ...
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    //gStyle->SetOptStat(0);

    if(flag_breakdown == 0){
    gROOT->ProcessLine(".x DrawOption.cc");
    TCanvas c1("ToyMC","ToyMC",2000,800);
    c1.Divide(4,2);
    c1.Draw();
    
    if (flag_err==1){
      
      for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
	TH1F *h1 = it->second.at(1);
	TH1F *h2 = it->second.at(2);
	for (int i=0;i!=h1->GetNbinsX()+1;i++){
	  h1->SetBinError(i+1,sqrt(h2->GetBinContent(i+1)));
	}
      }
      
      
      c1.cd(1);
      TGraphErrors *g10 = new TGraphErrors();
      TGraphErrors *g11 = new TGraphErrors();
      for (int i=0;i!=map_obsch_histos[1].at(0)->GetNbinsX()+1;i++){
	double x = map_obsch_histos[1].at(0)->GetBinCenter(i+1);
	double y = map_obsch_histos[1].at(0)->GetBinContent(i+1);
	double x_err = 0;
	double y_err = map_obsch_histos[1].at(0)->GetBinError(i+1);
	g10->SetPoint(i,x,y);
	g10->SetPointError(i,x_err,y_err);
	
	y = map_obsch_histos[1].at(1)->GetBinContent(i+1);
	y_err = map_obsch_histos[1].at(1)->GetBinError(i+1);
	
	g11->SetPoint(i,x,y);
	g11->SetPointError(i,x_err,y_err);
	
	//map_obsch_histos[1].at(0)->Draw();
	//map_obsch_histos[1].at(1)->Draw("same");
	//map_obsch_histos[1].at(1)->SetLineColor(2);
      }
      
      g10->Draw("A*"); g10->SetTitle("nueCC FC");
      g10->SetMarkerStyle(20);
      g11->Draw("*same");
      g11->SetMarkerStyle(21);
      g11->SetMarkerColor(2);
      g11->SetLineColor(2);
      
      // for (Int_t i=0;i!=map_obsch_histos[1].at(1)->GetNbinsX()+1;i++){
      //   std::cout << i << " " << map_obsch_histos[1].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[1].at(1)->GetBinError(i+1) << std::endl;
      // }
      
      c1.cd(2);
      // map_obsch_histos[3].at(0)->Draw();
      // map_obsch_histos[3].at(1)->Draw("same");
      // map_obsch_histos[3].at(1)->SetLineColor(2);
      
      // for (Int_t i=0;i!=map_obsch_histos[3].at(1)->GetNbinsX()+1;i++){
      //   std::cout << i << " " << map_obsch_histos[3].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[3].at(1)->GetBinError(i+1) << std::endl;
    // }

    TGraphErrors *g30 = new TGraphErrors();
    TGraphErrors *g31 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[3].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[3].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[3].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[3].at(0)->GetBinError(i+1);
      g30->SetPoint(i,x,y);
      g30->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[3].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[3].at(1)->GetBinError(i+1);

      g31->SetPoint(i,x,y);
      g31->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[3].at(0)->Draw();
      //map_obsch_histos[3].at(1)->Draw("same");
      //map_obsch_histos[3].at(1)->SetLineColor(2);
    }

    g30->Draw("A*");  g30->SetTitle("numuCC FC");
    g30->SetMarkerStyle(20);
    g31->Draw("*same");
    g31->SetMarkerStyle(21);
    g31->SetMarkerColor(2);
    g31->SetLineColor(2);

    
    c1.cd(3);
    // map_obsch_histos[5].at(0)->Draw();
    // map_obsch_histos[5].at(1)->Draw("same");
    // map_obsch_histos[5].at(1)->SetLineColor(2);
    
    // for (Int_t i=0;i!=map_obsch_histos[5].at(1)->GetNbinsX()+1;i++){
    //   std::cout << i << " " << map_obsch_histos[5].at(1)->GetBinContent(i+1) << " " << map_obsch_histos[5].at(1)->GetBinError(i+1) << std::endl;
    // }

    TGraphErrors *g50 = new TGraphErrors();
    TGraphErrors *g51 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[5].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[5].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[5].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[5].at(0)->GetBinError(i+1);
      g50->SetPoint(i,x,y);
      g50->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[5].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[5].at(1)->GetBinError(i+1);

      g51->SetPoint(i,x,y);
      g51->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[5].at(0)->Draw();
      //map_obsch_histos[5].at(1)->Draw("same");
      //map_obsch_histos[5].at(1)->SetLineColor(2);
    }

    g50->Draw("A*");  g50->SetTitle("CCpio FC");
    g50->SetMarkerStyle(20);
    g51->Draw("*same");
    g51->SetMarkerStyle(21);
    g51->SetMarkerColor(2);
    g51->SetLineColor(2);

    
    c1.cd(5);
    // map_obsch_histos[2].at(0)->Draw();
    // map_obsch_histos[2].at(1)->Draw("same");
    // map_obsch_histos[2].at(1)->SetLineColor(2);

    TGraphErrors *g20 = new TGraphErrors();
    TGraphErrors *g21 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[2].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[2].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[2].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[2].at(0)->GetBinError(i+1);
      g20->SetPoint(i,x,y);
      g20->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[2].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[2].at(1)->GetBinError(i+1);

      g21->SetPoint(i,x,y);
      g21->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[2].at(0)->Draw();
      //map_obsch_histos[2].at(1)->Draw("same");
      //map_obsch_histos[2].at(1)->SetLineColor(2);
    }

    g20->Draw("A*");g20->SetTitle("nueCC PC");
    g20->SetMarkerStyle(20);
    g21->Draw("*same");
    g21->SetMarkerStyle(21);
    g21->SetMarkerColor(2);
    g21->SetLineColor(2);
    

    
    c1.cd(6);
    // map_obsch_histos[4].at(0)->Draw();
    // map_obsch_histos[4].at(1)->Draw("same");
    // map_obsch_histos[4].at(1)->SetLineColor(2);

    TGraphErrors *g40 = new TGraphErrors();
    TGraphErrors *g41 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[4].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[4].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[4].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[4].at(0)->GetBinError(i+1);
      g40->SetPoint(i,x,y);
      g40->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[4].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[4].at(1)->GetBinError(i+1);

      g41->SetPoint(i,x,y);
      g41->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[4].at(0)->Draw();
      //map_obsch_histos[4].at(1)->Draw("same");
      //map_obsch_histos[4].at(1)->SetLineColor(2);
    }

    g40->Draw("A*"); g40->SetTitle("numuCC PC");
    g40->SetMarkerStyle(20);
    g41->Draw("*same");
    g41->SetMarkerStyle(21);
    g41->SetMarkerColor(2);
    g41->SetLineColor(2);
    
    c1.cd(7);
    // map_obsch_histos[6].at(0)->Draw();
    // map_obsch_histos[6].at(1)->Draw("same");
    // map_obsch_histos[6].at(1)->SetLineColor(2);

    TGraphErrors *g60 = new TGraphErrors();
    TGraphErrors *g61 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[6].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[6].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[6].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[6].at(0)->GetBinError(i+1);
      g60->SetPoint(i,x,y);
      g60->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[6].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[6].at(1)->GetBinError(i+1);

      g61->SetPoint(i,x,y);
      g61->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[6].at(0)->Draw();
      //map_obsch_histos[6].at(1)->Draw("same");
      //map_obsch_histos[6].at(1)->SetLineColor(2);
    }

    g60->Draw("A*");  g60->SetTitle("CCpio PC");
    g60->SetMarkerStyle(20);
    g61->Draw("*same");
    g61->SetMarkerStyle(21);
    g61->SetMarkerColor(2);
    g61->SetLineColor(2);

    
    c1.cd(4);
    // map_obsch_histos[7].at(0)->Draw();
    // map_obsch_histos[7].at(1)->Draw("same");
    // map_obsch_histos[7].at(1)->SetLineColor(2);

    TGraphErrors *g70 = new TGraphErrors();
    TGraphErrors *g71 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[7].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[7].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[7].at(0)->GetBinContent(i+1);
      double x_err = 0;
      double y_err = map_obsch_histos[7].at(0)->GetBinError(i+1);
      g70->SetPoint(i,x,y);
      g70->SetPointError(i,x_err,y_err);

      y = map_obsch_histos[7].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[7].at(1)->GetBinError(i+1);

      g71->SetPoint(i,x,y);
      g71->SetPointError(i,x_err,y_err);
      
      //map_obsch_histos[7].at(0)->Draw();
      //map_obsch_histos[7].at(1)->Draw("same");
      //map_obsch_histos[7].at(1)->SetLineColor(2);
    }

    g70->Draw("A*");  g70->SetTitle("NC pio");
    g70->SetMarkerStyle(20);
    g71->Draw("*same");
    g71->SetMarkerStyle(21);
    g71->SetMarkerColor(2);
    g71->SetLineColor(2);
    
  }else if (flag_err == 2){
    c1.cd(1);
    TGraphAsymmErrors *g10 = new TGraphAsymmErrors();
    TGraphErrors *g11 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[1].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[1].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[1].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g10->SetPoint(i,x,y);
      g10->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[1].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[1].at(1)->GetBinError(i+1);
      g11->SetPoint(i,x,y);
      g11->SetPointError(i,x_err,y_err);
    }
    g10->Draw("A*"); g10->SetTitle("nueCC FC");
    g10->SetMarkerStyle(20);
    g11->Draw("*same");
    g11->SetMarkerStyle(21);
    g11->SetMarkerColor(2);
    g11->SetLineColor(2);

    c1.cd(5);
    TGraphAsymmErrors *g20 = new TGraphAsymmErrors();
    TGraphErrors *g21 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[2].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[2].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[2].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g20->SetPoint(i,x,y);
      g20->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[2].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[2].at(1)->GetBinError(i+1);
      g21->SetPoint(i,x,y);
      g21->SetPointError(i,x_err,y_err);
    }
    g20->Draw("A*");  g20->SetTitle("nueCC PC");
    g20->SetMarkerStyle(20);
    g21->Draw("*same");
    g21->SetMarkerStyle(21);
    g21->SetMarkerColor(2);
    g21->SetLineColor(2);

    c1.cd(2);
    TGraphAsymmErrors *g30 = new TGraphAsymmErrors();
    TGraphErrors *g31 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[3].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[3].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[3].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g30->SetPoint(i,x,y);
      g30->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[3].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[3].at(1)->GetBinError(i+1);
      g31->SetPoint(i,x,y);
      g31->SetPointError(i,x_err,y_err);
    }
    g30->Draw("A*"); g30->SetTitle("numuCC FC");
    g30->SetMarkerStyle(20);
    g31->Draw("*same");
    g31->SetMarkerStyle(21);
    g31->SetMarkerColor(2);
    g31->SetLineColor(2);

    c1.cd(6);
    TGraphAsymmErrors *g40 = new TGraphAsymmErrors();
    TGraphErrors *g41 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[4].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[4].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[4].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g40->SetPoint(i,x,y);
      g40->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[4].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[4].at(1)->GetBinError(i+1);
      g41->SetPoint(i,x,y);
      g41->SetPointError(i,x_err,y_err);
    }
    g40->Draw("A*");  g40->SetTitle("numuCC PC");
    g40->SetMarkerStyle(20);
    g41->Draw("*same");
    g41->SetMarkerStyle(21);
    g41->SetMarkerColor(2);
    g41->SetLineColor(2);

    c1.cd(3);
    TGraphAsymmErrors *g50 = new TGraphAsymmErrors();
    TGraphErrors *g51 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[5].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[5].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[5].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g50->SetPoint(i,x,y);
      g50->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[5].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[5].at(1)->GetBinError(i+1);
      g51->SetPoint(i,x,y);
      g51->SetPointError(i,x_err,y_err);
    }
    g50->Draw("A*");  g50->SetTitle("CC pio FC");
    g50->SetMarkerStyle(20);
    g51->Draw("*same");
    g51->SetMarkerStyle(21);
    g51->SetMarkerColor(2);
    g51->SetLineColor(2);

    c1.cd(7);
    TGraphAsymmErrors *g60 = new TGraphAsymmErrors();
    TGraphErrors *g61 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[6].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[6].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[6].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g60->SetPoint(i,x,y);
      g60->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[6].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[6].at(1)->GetBinError(i+1);
      g61->SetPoint(i,x,y);
      g61->SetPointError(i,x_err,y_err);
    }
    g60->Draw("A*"); g60->SetTitle("CCpio PC");
    g60->SetMarkerStyle(20);
    g61->Draw("*same");
    g61->SetMarkerStyle(21);
    g61->SetMarkerColor(2);
    g61->SetLineColor(2);

    c1.cd(4);
    TGraphAsymmErrors *g70 = new TGraphAsymmErrors();
    TGraphErrors *g71 = new TGraphErrors();
    for (int i=0;i!=map_obsch_histos[7].at(0)->GetNbinsX()+1;i++){
      double x = map_obsch_histos[7].at(0)->GetBinCenter(i+1);
      double y = map_obsch_histos[7].at(0)->GetBinContent(i+1);
      double x_err = 0;
      auto result = cov.get_bayes_errors(y);
      double y_err = 0;
      g70->SetPoint(i,x,y);
      g70->SetPointError(i,0,0,result.first, result.second);
      y = map_obsch_histos[7].at(1)->GetBinContent(i+1);
      y_err = map_obsch_histos[7].at(1)->GetBinError(i+1);
      g71->SetPoint(i,x,y);
      g71->SetPointError(i,x_err,y_err);
    }
    g70->Draw("A*");g70->SetTitle("NCpio " );
    g70->SetMarkerStyle(20);
    g71->Draw("*same");
    g71->SetMarkerStyle(21);
    g71->SetMarkerColor(2);
    g71->SetLineColor(2);
    
  }
  theApp.Run();
  }
  else // flag_breakdown == true
  {
    gROOT->ProcessLine(".x DrawOption.cc");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(132);
    gStyle->SetLegendTextSize(0.04);

    int nchannels = map_obsch_subhistos.size();
    TCanvas *canvas[nchannels];
    TGraphAsymmErrors *gr[nchannels];
    TGraphAsymmErrors *gratio_mc[nchannels];
    TGraphAsymmErrors *gratio_data[nchannels];
    THStack *hstack[nchannels]; 
    TLegend *legend[nchannels]; 
    TLegend *legend2[nchannels]; 
    for(auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){
        int obschannel = it->first;
        std::cout<<"Channel: "<<obschannel<<std::endl;
        canvas[obschannel-1] = new TCanvas(Form("canvas%d", obschannel), Form("channel%d", obschannel), 1200, 900);
        TPad *pad1 = new TPad("pad1", "", 0.01,0.3,0.99,0.99,0,0,0);
        TPad *pad2 = new TPad("pad2", "", 0.01,0.01,0.99,0.3,0,0,0);
        pad1->SetBottomMargin(0);
        pad1->SetLeftMargin(0.12);
        pad1->SetRightMargin(0.1);
        pad2->SetTopMargin(0.05);
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.1);
        pad2->SetBottomMargin(0.20);
        pad1->Draw();
        pad2->Draw();
        hstack[obschannel-1] = new THStack(Form("hs%d", obschannel),"");
        legend[obschannel-1] = new TLegend(0.3, 0.65, 0.85, 0.92);
        TH1F* hdata = (TH1F*)map_obsch_histos[obschannel].at(0)->Clone("hdata");
        TH1F* hbadmatch = (TH1F*)hdata->Clone("hbadmatch");
        TH1F* hnumuCCinFV = (TH1F*)hdata->Clone("hnumuCCinFV");
        TH1F* hnueCCinFV = (TH1F*)hdata->Clone("hnueCCinFV");
        TH1F* hNCinFV = (TH1F*)hdata->Clone("hNCinFV");
        TH1F* hCCpi0inFV = (TH1F*)hdata->Clone("hCCpi0inFV");
        TH1F* hNCpi0inFV = (TH1F*)hdata->Clone("hNCpi0inFV");
        TH1F* houtFV = (TH1F*)hdata->Clone("houtFV");
        TH1F* hext = (TH1F*)hdata->Clone("hext");
        TH1F* hdirt = (TH1F*)hdata->Clone("hdirt");
        TH1F* hLEE = (TH1F*)hdata->Clone("hLEE");
        hbadmatch->Reset();
        hnumuCCinFV->Reset();
        hnueCCinFV->Reset();
        hNCinFV->Reset();
        hCCpi0inFV->Reset();
        hNCpi0inFV->Reset();
        houtFV->Reset();
        hext->Reset();
        hdirt->Reset();
        hLEE->Reset();
        bool flag_leeexist = false;
        for(size_t i=0; i<it->second.size(); i++){
            TH1F* htemp = map_obsch_subhistos[obschannel].at(i);
            std::string histname = htemp->GetName();
            std::istringstream sss(histname);
            for(std::string line; std::getline(sss, line, '_');){
                if(line == "badmatch") {
                    std::cout<<"badmatch"<<" "<<histname<<std::endl;
                    hbadmatch->Add(htemp);
                    break;
                }
                if(line == "numuCCinFV") {
                    std::cout<<"numuCCinFV"<<" "<<histname<<std::endl;
                    hnumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "nueCCinFV") {
                    std::cout<<"nueCCinFV"<<" "<<histname<<std::endl;
                    hnueCCinFV->Add(htemp);
                    break;
                }
                if(line == "NCinFV") {
                    std::cout<<"NCinFV"<<" "<<histname<<std::endl;
                    hNCinFV->Add(htemp);
                    break;
                }
                if(line == "CCpi0inFV") {
                    std::cout<<"CCpi0inFV"<<" "<<histname<<std::endl;
                    hCCpi0inFV->Add(htemp);
                    break;
                }
                if(line == "NCpi0inFV") {
                    std::cout<<"NCpi0inFV"<<" "<<histname<<std::endl;
                    hNCpi0inFV->Add(htemp);
                    break;
                }
                if(line == "outFV") {
                    std::cout<<"outFV"<<" "<<histname<<std::endl;
                    houtFV->Add(htemp);
                    break;
                }
                if(line == "ext") {
                    std::cout<<"ext"<<" "<<histname<<std::endl;
                    hext->Add(htemp);
                    break;
                }
                if(line == "dirt") {
                    std::cout<<"dirt"<<" "<<histname<<std::endl;
                    hdirt->Add(htemp);
                    break;
                }
                if(line == "LEE") {
                    std::cout<<"LEE"<<" "<<histname<<std::endl;
                    flag_leeexist = true;
                    hLEE->Add(htemp);
                    break;
                }
            }
        }
        pad1->cd(); 
        float datapot = 0;
        if(run == 0){
            for(auto it=map_data_period_pot.begin(); it!=map_data_period_pot.end(); it++)
            {
                datapot += it->second;
            }
        }else datapot = map_data_period_pot[run];

        gr[obschannel-1] = new TGraphAsymmErrors();
        legend[obschannel-1]->SetNColumns(2);
        legend[obschannel-1]->AddEntry((TObject*)0, Form("Data POT: %.3e", datapot), "");
        legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("BNB data, %.1f", hdata->Integral()), "lp");
        hstack[obschannel-1]->Add(hbadmatch); 
        legend[obschannel-1]->AddEntry(hbadmatch, Form("Cosmic, %.1f", hbadmatch->Integral()), "F"); 
        hbadmatch->SetFillStyle(3004);
        hbadmatch->SetFillColorAlpha(kRed+2, 0.5);
        hbadmatch->SetLineColor(kRed+2);
        hbadmatch->SetLineWidth(1);

        hstack[obschannel-1]->Add(hext); 
        legend[obschannel-1]->AddEntry(hext, Form("EXT, %.1f", hext->Integral()), "F"); 
        hext->SetFillStyle(3004);
        hext->SetFillColorAlpha(kOrange+3, 0.5);
        hext->SetLineColor(kOrange+3);
        hext->SetLineWidth(1);

        hstack[obschannel-1]->Add(hdirt); 
        legend[obschannel-1]->AddEntry(hdirt, Form("Dirt, %.1f", hdirt->Integral()), "F"); 
        hdirt->SetFillStyle(3224);
        hdirt->SetFillColorAlpha(kGray, 0.5);
        hdirt->SetLineColor(kGray+2);
        hdirt->SetLineWidth(1);

        hstack[obschannel-1]->Add(houtFV);
        legend[obschannel-1]->AddEntry(houtFV, Form("out FV, %.1f", houtFV->Integral()), "F"); 
        houtFV->SetFillStyle(3224);
        houtFV->SetFillColorAlpha(kOrange+1, 0.5);
        houtFV->SetLineColor(kOrange+1);
        houtFV->SetLineWidth(1);

        hstack[obschannel-1]->Add(hNCpi0inFV); 
        legend[obschannel-1]->AddEntry(hNCpi0inFV, Form("NC #pi^{0} in FV,  %.1f", hNCpi0inFV->Integral()), "F"); 
        hNCpi0inFV->SetFillStyle(1001);
        hNCpi0inFV->SetFillColorAlpha(38, 0.5);
        hNCpi0inFV->SetLineColor(38);
        hNCpi0inFV->SetLineWidth(1);

        hstack[obschannel-1]->Add(hCCpi0inFV); 
        legend[obschannel-1]->AddEntry(hCCpi0inFV, Form("CC #pi^{0} in FV, %.1f", hCCpi0inFV->Integral()), "F"); 
        hCCpi0inFV->SetFillStyle(1001);
        hCCpi0inFV->SetFillColorAlpha(30, 0.5);
        hCCpi0inFV->SetLineColor(30);
        hCCpi0inFV->SetLineWidth(1);

        hstack[obschannel-1]->Add(hNCinFV); 
        legend[obschannel-1]->AddEntry(hNCinFV, Form("NC in FV, %.1f", hNCinFV->Integral()), "F"); 
        hNCinFV->SetFillStyle(1001);
        hNCinFV->SetFillColorAlpha(kOrange+1, 0.5);
        hNCinFV->SetLineColor(kOrange+1);
        hNCinFV->SetLineWidth(1);

        hstack[obschannel-1]->Add(hnumuCCinFV); 
        legend[obschannel-1]->AddEntry(hnumuCCinFV, Form("#nu_{#mu} CC in FV, %.1f", hnumuCCinFV->Integral()), "F"); 
        hnumuCCinFV->SetFillStyle(1001);
        hnumuCCinFV->SetFillColorAlpha(kAzure+6, 0.5);
        hnumuCCinFV->SetLineColor(kAzure+6);
        hnumuCCinFV->SetLineWidth(1);

        hstack[obschannel-1]->Add(hnueCCinFV); 
        legend[obschannel-1]->AddEntry(hnueCCinFV, Form("#nu_{e} CC in FV, %.1f", hnueCCinFV->Integral()), "F"); 
        hnueCCinFV->SetFillStyle(1001);
        hnueCCinFV->SetFillColorAlpha(kGreen+1, 0.5);
        hnueCCinFV->SetLineColor(kGreen+1);
        hnueCCinFV->SetLineWidth(1);
       
        if(flag_leeexist){
        hstack[obschannel-1]->Add(hLEE); 
        legend[obschannel-1]->AddEntry(hLEE, Form("LEE, %.1f", hLEE->Integral()), "F");
        hLEE->SetFillStyle(1001);
        hLEE->SetFillColorAlpha(kMagenta, 0.5);
        hLEE->SetLineColor(kMagenta);
        hLEE->SetLineWidth(1);
        }

        TH1F* hmc = (TH1F*)map_obsch_histos[obschannel].at(1)->Clone("hmc");
        hmc->Draw("hist");
        hmc->GetYaxis()->SetTitle("Event counts");
        float mcymax = hmc->GetBinContent(hmc->GetMaximumBin());
        float dataymax = hdata->GetBinContent(hdata->GetMaximumBin());
        if(dataymax>mcymax) mcymax = dataymax;
        hmc->SetMaximum(2.0*mcymax);
        hmc->GetYaxis()->SetRangeUser(-0.02*mcymax, 1.6*mcymax);
        hmc->SetLineColor(kBlack);
        hmc->SetLineWidth(5);
        hstack[obschannel-1]->Draw("hist same");
        
        gratio_mc[obschannel-1] = new TGraphAsymmErrors();
        gratio_data[obschannel-1] = new TGraphAsymmErrors();
        float maxratio = 1.5;
        for(int i=0; i<hdata->GetNbinsX(); i++)
        {
            double x = hdata->GetBinCenter(i+1);
            double x_err = hdata->GetBinWidth(1)*0.5;
            double y = hdata->GetBinContent(i+1);
            double y_err = hdata->GetBinError(i+1);
            auto bayesError = cov.get_bayes_errors(y);
            double bayesError_low = bayesError.first;
            double bayesError_up = bayesError.second;
            double ymc = hmc->GetBinContent(i+1);
            double ymc_err = hmc->GetBinError(i+1);
            gr[obschannel-1]->SetPoint(i,x,y);
            gratio_mc[obschannel-1]->SetPoint(i,x,1);
            if(ymc!=0){ 
                gratio_data[obschannel-1]->SetPoint(i,x,y/ymc); 
                if(maxratio<y/ymc) maxratio = y/ymc;
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, ymc_err/ymc, ymc_err/ymc);
            }
            else { 
                gratio_data[obschannel-1]->SetPoint(i, x, 10); // invalid value 
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==2){ //update data point errors
                gr[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low, bayesError_up);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low/ymc, bayesError_up/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==1){
                gr[obschannel-1]->SetPointError(i, x_err, x_err, y_err, y_err);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, y_err/ymc, y_err/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
        }
        gr[obschannel-1]->Draw("P same"); 
        gr[obschannel-1]->SetLineWidth(2);
        gr[obschannel-1]->SetMarkerStyle(20);
        gr[obschannel-1]->SetMarkerSize(1.5);
        gr[obschannel-1]->SetLineColor(kBlack);

        //legend[obschannel-1]->SetFillStyle(0);
        legend[obschannel-1]->SetHeader(Form("#SigmaDATA/#Sigma(MC+EXT)=%.2f", hdata->Integral()/hmc->Integral()), "C");
        legend[obschannel-1]->Draw();
        pad2->cd();
        gratio_mc[obschannel-1]->Draw("a2");
        gratio_mc[obschannel-1]->SetFillColor(kGray);
        gratio_mc[obschannel-1]->GetYaxis()->SetRangeUser(0,int(1.5*maxratio));
        gratio_mc[obschannel-1]->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        gratio_mc[obschannel-1]->GetYaxis()->SetTitle("Data/Pred");
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleOffset(0.5);
        if(obschannel>=5) //hard coded at this moment
        {
            gratio_mc[obschannel-1]->GetXaxis()->SetTitle("Reco #pi^{0} mass [MeV]");
        }
        else gratio_mc[obschannel-1]->GetXaxis()->SetTitle("Reco neutrino energy [MeV]");
        gratio_mc[obschannel-1]->GetXaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetXaxis()->SetLabelSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleOffset(0.35);
        gratio_mc[obschannel-1]->GetYaxis()->SetLabelSize(0.1);
        gratio_data[obschannel-1]->Draw("P same");
        gratio_data[obschannel-1]->SetLineWidth(2);
        gratio_data[obschannel-1]->SetMarkerStyle(20);
        gratio_data[obschannel-1]->SetMarkerSize(1.5);
        gratio_data[obschannel-1]->SetLineColor(kBlack);

        TH1F* hist = (TH1F*)hdata->Clone("hist");
        hist->Reset();
        hist->Draw("axis same");

        TLine* line = new TLine(hmc->GetXaxis()->GetXmin(),1,hmc->GetXaxis()->GetXmax(),1);
        line->Draw();
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        legend2[obschannel-1] = new TLegend(0.2, 0.7, 0.8, 0.95);
        legend2[obschannel-1]->SetNColumns(2);
        legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat. uncertainty (Bayesian)", "F");
        legend2[obschannel-1]->AddEntry(gratio_data[obschannel-1],"Data stat. uncertainty (Bayesian)", "lp");
        legend2[obschannel-1]->SetTextSize(0.06);
        legend2[obschannel-1]->SetFillStyle(0);
        legend2[obschannel-1]->Draw();

        if(obschannel==1) canvas[obschannel-1]->Print("selection.pdf(");
        else if(obschannel==nchannels) canvas[obschannel-1]->Print("selection.pdf)");
        else canvas[obschannel-1]->Print("selection.pdf");
    } 
    theApp.Run();
  } // flag_breakdown end
  }


  
  
  // std::map<TString, std::pair<TString, int> > map_pred_histo_hist_err2_lee = cov.get_map_pred_histo_histo_err2_lee();
  // std::map<std::pair<TString, TString>, std::pair<TString, int> > map_pair_histo_histos_cross = cov.get_map_pair_hist_hists_cros();
  // std::map<int, std::set<std::set<std::pair<TString, int> > > > map_pred_obsch_histos = cov.get_map_pred_obsch_histos();

  if (flag_err == 1){
    // prediction ...
    TFile *file3 = new TFile("merge.root","RECREATE");
    
    TMatrixD* mat_collapse = cov.get_mat_collapse();
    mat_collapse->Write("mat_collapse");
    
    // additional covariance matrix ...
    TMatrixD* mat_add_cov = cov.get_add_cov_matrix();
    
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
	
	  //  if (covch == 8 && obsch == 1)
	  //  std::cout << "bin: " << j << " " << std::get<0>(bayes_inputs.at(i).at(j)) << " " << std::get<1>(bayes_inputs.at(i).at(j)) << " " << std::get<3>(bayes_inputs.at(i).at(j)) << std::endl;
	  //	  if (std::get<4>(bayes_inputs.at(i).at(j)) > 0)
	  //  std::cout << obsch << " " << i << " " << start_bin + j << " " <<  std::get<4>(bayes_inputs.at(i).at(j)) << std::endl;
	}
	
	//	std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << i << " " << covch << " " << start_bin << std::endl;
	//  break;
      }
    }
    
    //for (auto it = map_name_histogram.begin(); it != map_name_histogram.end(); it++){
    //  it->second.first->SetDirectory(file3);
    // }
    
    
    
    //    for (auto it = map_covch_histo.begin(); it != map_covch_histo.end(); it++){
    //  it->second->SetDirectory(file3);
    // }
    
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
    
    file3->Write();
    file3->Close();
  }  

}
