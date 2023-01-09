#include "WCPLEEANA/master_cov_matrix.h"



#include "TPrincipal.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TRandom3.h"

//#include "WCPLEEANA/cuts.h"
//#include "WCPLEEANA/pot.h"
//#include "WCPLEEANA/weights.h"
//#include "WCPLEEANA/pfeval.h"
//#include "WCPLEEANA/eval.h"
//#include "WCPLEEANA/kine.h"
//#include "WCPLEEANA/tagger.h"


#include <algorithm>


using namespace LEEana;


float leeweight(float Enu)
 {
         if(Enu<200 || Enu>=800) return 0.0;
         else if(Enu>=200 && Enu<250) return 6.37441;
         else if(Enu>=250 && Enu<300) return 5.64554;
         else if(Enu>=300 && Enu<350) return 3.73055;
         else if(Enu>=350 && Enu<400) return 1.50914;
         else if(Enu>=400 && Enu<450) return 1.07428;
         else if(Enu>=450 && Enu<500) return 0.754093;
         else if(Enu>=500 && Enu<600) return 0.476307;
         else if(Enu>=600 && Enu<800) return 0.152327;
         else return 0.0;
 }

#include "mcm_1.h"
#include "mcm_2.h"
#include "mcm_data_stat.h"
#include "mcm_pred_stat.h"

LEEana::CovMatrix::CovMatrix(TString cov_filename, TString cv_filename, TString file_filename, TString rw_filename){
  flag_osc = false;

  rw_type = 0;
  flag_reweight = false;

  std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > >  rw_info_vec;
  std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  >  rw_info_i;

  std::ifstream infile_rw(rw_filename);
  if(infile_rw.is_open()){
    while(!infile_rw.eof()){
      bool flag_reweight_i = 0;
      TString cut_str;
      TString var_str;
      double min_var;
      double max_var;
      bool use_underflow;
      bool use_overflow;
      bool equal_bins=true;
      std::vector<double> reweight;
      std::vector<double> bins;
      //reads the gerneral reweighig parameters
      std::string line;
      std::getline(infile_rw, line);
      while(line.empty()) {std::getline(infile_rw, line);}
      if(line == "end") break;
      std::istringstream iss(line);
      iss >> flag_reweight_i >> cut_str >> var_str >> min_var >> max_var >> use_underflow >> use_overflow >> equal_bins;
      if(flag_reweight_i) flag_reweight=1;

      if(!(equal_bins)){
        double tb;
        std::getline(infile_rw, line);
        while(line.empty()) {std::getline(infile_rw, line);} //skips spaces
        std::istringstream iss_rw(line);
        while(iss_rw>>tb){bins.push_back(tb);}
        if(bins[0]!=min_var) {
          min_var=bins[0];
          std::cerr<<"WARNING: bins and min_var do not match for "<<cut_str<<". Forcing min_var=bins[0]. Check reweighting configuration"<<std::endl;
        }
        if(bins.back()!=max_var){
          max_var=bins.back();
          std::cerr<<"WARNING: bins and max_var do not match for "<<cut_str<<". Forcing max_var=bins.back(). Check reweighting configuration"<<std::endl;
        }
      }

      //now get the weights
      double trw;
      std::getline(infile_rw, line);
      while(line.empty()) {std::getline(infile_rw, line);} //skips spaces
      std::istringstream iss_rw(line);
      while(iss_rw>>trw){reweight.push_back(trw);}
      rw_info_i = std::make_tuple(flag_reweight_i, cut_str, var_str, min_var, max_var, use_underflow, use_overflow, equal_bins, reweight, bins);
      rw_info_vec.push_back(rw_info_i);
      if(!(equal_bins)){
        if(use_underflow && use_overflow){
          if(bins.size()+1!=reweight.size()) std::cerr<<"WARNING: bins and reweighting function for "<<cut_str<<" do not match. Check reweighting configuration"<<std::endl;
        }else if(use_underflow || use_overflow){
           if(bins.size()!=reweight.size()) std::cerr<<"WARNING: bins and reweighting function for "<<cut_str<<" do not match. Check reweighting configuration"<<std::endl;
        }else if(bins.size()-1!=reweight.size()) std::cerr<<"WARNING: bins and reweighting function for "<<cut_str<<" do not match. Check reweighting configuration"<<std::endl;
      } 
    }
  }
  rw_info = std::make_tuple(flag_reweight, rw_info_vec); 
  


  std::ifstream infile(cov_filename);
  TString name, var_name;
  Int_t bin_num;
  Float_t low_limit, high_limit;
  Int_t obs_no;
  Int_t flag_xs_flux;
  Int_t flag_det;
  Float_t flag_add;
  Int_t flag_same_mc_stat;
  Int_t cov_sec_no;
  Int_t file_no;

  int ch_no = 0;
  
  int covbin = 0;
  int obscovbin = 0;
  TString weight;
  int lee_strength;
  
  while(!infile.eof()){
    infile >> name >> var_name >> bin_num >> low_limit >> high_limit >> obs_no >> flag_xs_flux >> flag_det >> flag_add >> flag_same_mc_stat >> cov_sec_no >> file_no >> weight >> lee_strength;
    // std::cout << name << " " << var_name << " " << low_limit << " " << bin_num << " " << file_no << std::endl;
    if (bin_num == -1) break;
    
    map_ch_hist[ch_no] = std::make_tuple(name, var_name, bin_num, low_limit, high_limit, weight, obs_no, lee_strength);
    map_name_ch[name] = ch_no;
    
    map_ch_filetype[ch_no] = file_no;
    map_filetype_chs[file_no].push_back(ch_no);

    if (file_no !=5 && file_no != 15) { // data
      map_ch_systematics[ch_no] = std::make_tuple(flag_xs_flux, flag_det, flag_add, flag_same_mc_stat);
      if (flag_xs_flux == 1) xfs_filetypes.insert(file_no);
      if (flag_det == 1) det_filetypes.insert(file_no);
      if (flag_add != 0) add_filetypes.insert(file_no);
      
      if (flag_same_mc_stat !=0) {
	map_mcstat_same_chs[flag_same_mc_stat].push_back(ch_no);
	map_filetype_mcstats[file_no].insert(flag_same_mc_stat);
      }
    
      // prepare for the matrix
      map_ch_obsch[ch_no] = obs_no;
      map_ch_covch[ch_no] = cov_sec_no;
      
      map_obsch_nbin[obs_no]     = bin_num + 1; // add the overflow bin
      map_covch_nbin[cov_sec_no] = bin_num + 1; // add the overflow bin
      
      map_covch_obsch[cov_sec_no] = obs_no; // ch map

      // relation among obs_no, cov_no, ch_no ...
      map_pred_obsch_covch[obs_no].insert(cov_sec_no);
      map_pred_covch_ch[cov_sec_no].insert(ch_no);
    }
    
    ch_no ++;
  }

  int total_obs_bin = 0;
  int start_bin = 0;
  for (auto it = map_obsch_nbin.begin(); it!= map_obsch_nbin.end(); it++){
    map_obsch_startbin[it->first] = start_bin;
    start_bin += it->second;
    total_obs_bin += it->second;
    //std::cout << it->first << " " << it->second << std::endl;
  }

  int total_cov_bin = 0;
  start_bin = 0;
  for (auto it = map_covch_nbin.begin(); it != map_covch_nbin.end(); it++){
    map_covch_startbin[it->first] = start_bin;
    start_bin += it->second;
    total_cov_bin += it->second;
  }
 
  // std::cout << "cov matrix : " << total_cov_bin << ", " << total_obs_bin << std::endl;
  mat_collapse = new TMatrixD(total_cov_bin,total_obs_bin);
  mat_add_cov = new TMatrixD(total_cov_bin, total_cov_bin);
  
  // form the large covariance matrix, and the start bin ...
  for (auto it = map_covch_obsch.begin(); it!= map_covch_obsch.end(); it++){
    //std::cout << it->first << " " << it->second << std::endl;
    int nbin = map_covch_nbin[it->first];
    int start_bin_cov = map_covch_startbin[it->first];
    int start_bin_obs = map_obsch_startbin[it->second];
    for (int i=0; i!=nbin; i++){
      map_covchbin_obschbin[start_bin_cov + i] =  start_bin_obs + i;
      (*mat_collapse)(start_bin_cov + i, start_bin_obs + i) = 1;
      //std::cout << start_bin_cov + i << " " << start_bin_obs + i << std::endl;
    }
  }

  int filetype;
  //  TString name;
  int period;
  TString input_filename;
  TString out_filename;
  float ext_pot;
  double norm_pot;
  int norm_period;
  //std::cout << cv_filename << std::endl;
  std::ifstream infile1(cv_filename);
  while(!infile1.eof()){
    infile1 >> filetype >> name >> period >> input_filename >> out_filename >> ext_pot >> file_no >> norm_pot >> norm_period;
    //std::cout << filetype << " " << out_filename << " " << file_no << std::endl;
    
    if (filetype == -1) break;
    
    map_filetype_name[filetype] = name;
    map_filetype_inputfiles[filetype].push_back(input_filename);
    map_inputfile_filetype[input_filename] = filetype;
    map_inputfile_info[input_filename] = std::make_tuple(filetype, period, out_filename, ext_pot, file_no, norm_pot, norm_period);
    map_fileno_period[file_no] = period;
  }

  std::ifstream infile2(file_filename);
  TString cut_name;
  while (!infile2.eof()){
    infile2 >> input_filename >> cut_name;
    if (input_filename == "end") break;
    map_inputfile_cuts[input_filename].push_back(cut_name);

    //    std::cout << input_filename << " " << cut_name << std::endl;
    
    // cut_name, input_filename --> filetype --> chs ...
    int filetype = map_inputfile_filetype[input_filename];
    std::vector<int> chs = map_filetype_chs[filetype];
    for (auto it = chs.begin(); it!= chs.end(); it++){
      map_pred_ch_subch[*it].insert(std::make_pair(std::get<0>(map_ch_hist[*it]),cut_name));
    }
  }


  // sort out the histograms ...
  for (auto it = map_inputfile_info.begin(); it!= map_inputfile_info.end(); it++){
    TString filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    int file_no = std::get<4>(it->second);
    for (auto it1 = map_inputfile_cuts[filename].begin(); it1 != map_inputfile_cuts[filename].end(); it1++){
      TString add_cut = *it1;

      //      std::cout << filename << " " << add_cut << std::endl;
      
      for (auto it2 = map_filetype_chs[filetype].begin(); it2 != map_filetype_chs[filetype].end(); it2++){
	int ch = *it2;
	auto it3 = map_ch_hist.find(ch);

	TString name = std::get<0>(it3->second);
	TString var_name = std::get<1>(it3->second);
	int nbin = std::get<2>(it3->second);
	float llimit = std::get<3>(it3->second);
	float hlimit = std::get<4>(it3->second);
	TString weight = std::get<5>(it3->second);
	int lee_strength = std::get<7>(it3->second);

	//std::cout << name << " " << lee_strength << std::endl;
	
	TString weight2 = weight + "_" + weight;
	TString histo_name = name + Form("_%d_",file_no) + var_name + "_" + add_cut;
	TString histo_name1 = histo_name + "_err2";

	map_histogram_inputfile[histo_name] = filename;
	map_histogram_inputfile[histo_name1] = filename;
	
	map_histogram_covch_add[histo_name] = std::make_pair(map_ch_covch[map_name_ch[name]], std::get<2>(map_ch_systematics[map_name_ch[name]]));
	map_histogram_covch_add[histo_name1] = std::make_pair(map_ch_covch[map_name_ch[name]], std::get<2>(map_ch_systematics[map_name_ch[name]]));
	
	map_inputfile_histograms[filename].push_back(std::make_tuple(histo_name, nbin, llimit, hlimit, var_name, name, add_cut, weight));
	map_inputfile_histograms_err2[filename].push_back(std::make_tuple(histo_name1, nbin, llimit, hlimit, var_name, name, add_cut, weight2));

	map_pred_subch_histos[std::make_pair(name,add_cut)].insert(std::make_pair(histo_name, period));
	map_pred_histo_histo_err2_lee[histo_name] = std::make_pair(histo_name1,lee_strength);
	//std::cout << histo_name << " " << " " << histo_name1 << " " << nbin << " " << llimit << " " << hlimit << " " << var_name << " " << name << " " << add_cut << std::endl;
	
	//	std::cout << filename << " " << add_cut << std::endl;
      }
      std::set<int> mc_stat = map_filetype_mcstats[filetype];
      for (auto it5 = mc_stat.begin(); it5 != mc_stat.end(); it5++){
	std::vector<int> correlated_chs = map_mcstat_same_chs[*it5];
	//std::cout << correlated_chs.size() << " " << filetype << " " << filename << std::endl;
	for (size_t i=0;i!=correlated_chs.size();i++){
	  int ch1 = correlated_chs.at(i);
	  int obsch1 = map_ch_obsch[ch1];
	  auto it3 = map_ch_hist.find(ch1);
	  TString name1 = std::get<0>(it3->second);
	  TString var_name1 = std::get<1>(it3->second);
	  int nbin1 = std::get<2>(it3->second);
	  float llimit1 = std::get<3>(it3->second);
	  float hlimit1 = std::get<4>(it3->second);
	  TString weight1 = std::get<5>(it3->second);
	  TString histo_name1 = name1 + Form("_%d_",file_no) + var_name1 + "_" + add_cut;
	  
	  for (size_t j=i+1;j<correlated_chs.size();j++){
	    int ch2 = correlated_chs.at(j);
	    auto it4 = map_ch_hist.find(ch2);
	    TString name2 = std::get<0>(it4->second);
	    TString var_name2 = std::get<1>(it4->second);
	    int nbin2 = std::get<2>(it4->second);
	    float llimit2 = std::get<3>(it4->second);
	    float hlimit2 = std::get<4>(it4->second);
	    TString weight2 = std::get<5>(it4->second);
	    TString histo_name2 = name2 + Form("_%d_",file_no) + var_name2 + "_" + add_cut;
	    TString histo_name = name1 + "_" + name2 + "_" + add_cut + Form("_%d",file_no);
	    TString weight = weight1 +"_" + weight2;

	    map_histogram_inputfile[histo_name] = filename;
	    map_histogram_covch_add[histo_name] = std::make_pair(-1,0);
	    
	    map_inputfile_histograms_cros[filename].push_back(std::make_tuple(histo_name, nbin1, llimit1, hlimit1, var_name1, name1, add_cut, weight));
	    
	    map_pair_histo_histos_cros[std::make_pair(histo_name1, histo_name2)] = std::make_pair(histo_name, obsch1);
	    
	  }
	}
      }
      
    }
  }


  // now form the final prediction map ...
  for (auto it = map_pred_obsch_covch.begin(); it!= map_pred_obsch_covch.end(); it++){
    int obsch = it->first;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      int covch = *it1;
      std::set<int> chs = map_pred_covch_ch[covch];
      for (auto it2 = chs.begin(); it2 != chs.end(); it2++){
	int ch = *it2;
	std::set<std::pair<TString, TString> > subchs = map_pred_ch_subch[ch];
	for (auto it3 = subchs.begin(); it3 != subchs.end(); it3++){
	  std::pair<TString, TString> subch = *it3;
	  std::set<std::pair<TString, int> > histos = map_pred_subch_histos[subch];

	  map_pred_obsch_histos[obsch].insert(histos);
	  map_pred_covch_histos[covch].insert(histos);
	  //for (auto it4 = histos.begin(); it4 != histos.end(); it4++){
	  //  TString histo = *it4;
	  //map_pred_obsch_histos[obsch].insert(histo);
	  //}
	  
	}
      }
    }
  }

  // for (auto it = map_pred_obsch_histos.begin(); it!=map_pred_obsch_histos.end();it++){
  //   std::cout << it->first << std::endl;
  //   for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
  //     std::cout << "sub: " << (*it1).size() << std::endl;
  //     for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
  // 	std::cout << it->first << " " << (*it2).first << " " << (*it2).second << std::endl;
  //     }
  //   }
  // }

  g_llimit = 0;
  g_hlimit = 100;
  
  Double_t x1[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		    31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
		    41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
		    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
		    61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
		    71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
		    81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
		    91, 92, 93, 94, 95, 96, 97, 98, 99, 100};

  
  Double_t yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117,
		    7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734,
		    16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803,
		    25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396,
		    34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975,
		    43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257,
		    53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075,
		    62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317,
		    71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907,
		    81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788}; 
  Double_t yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179,
		    14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133,
		    25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117,
		    36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549,
		    47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985,
		    58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713,
		    68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902,
		    79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665,
		    89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079,
		    100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};

  for (Int_t i=0;i!=101;i++){
    yl[i] = fabs(yl[i] - x1[i]);
    yh[i] = fabs(yh[i] - x1[i]);
  }
  
  gl = new TGraph(101,x1,yl);
  gh = new TGraph(101,x1,yh);


  flag_spec_weights = false;
  //init_spec_weights(2800,1000,0.15); // binning 25 neutrino energy, 7 cos theta, 20, muon energy ...
}


LEEana::CovMatrix::~CovMatrix(){

  //std::cout << "haha " << std::endl;
  delete mat_collapse;
  delete mat_add_cov;
  delete gl;
  delete gh;

  //  std::cout << "hehe " << std::endl;
}


void LEEana::CovMatrix::init_spec_weights(int num, int num1, double strength){
  flag_spec_weights = true;
  std::cout << "initialize special weights" << std::endl;
  spec_weights.resize(num);
  for (int i=0;i!=num;i++){
    spec_weights.at(i).resize(num1,0);
    for (int j=0;j!=num1;j++){
      float temp = gRandom->Gaus(0,strength);
      if (temp < -1) temp = -1;
      spec_weights.at(i).at(j) = temp;
    }
  } 
}

std::vector<float> LEEana::CovMatrix::get_spec_weight(LEEana::EvalInfo& eval, LEEana::PFevalInfo& pfeval){
  //std::cout << eval.truth_isCC << " " << eval.truth_nuPdg << " " << spec_weights.at(0).at(0) << " " << spec_weights.at(0).at(1) << std::endl;
  if (eval.truth_isCC==1 && fabs(eval.truth_nuPdg)==14 ){
    double Enu = eval.truth_nuEnergy;
    double Pmu = sqrt(pow(pfeval.truth_muonMomentum[0],2) + pow(pfeval.truth_muonMomentum[1],2) + pow(pfeval.truth_muonMomentum[2],2));
    double costheta = 1.;
    if (Pmu>0)
      costheta = pfeval.truth_muonMomentum[2]/Pmu;
    
    int index_Enu = 0, index_Pmu = 0, index_costheta = 0;
    if (Enu < 200){   // 17 bins ...
    }else if (Enu < 300){ index_Enu = 1;
    }else if (Enu < 400){ index_Enu = 2;
    }else if (Enu < 500){ index_Enu = 3;
    }else if (Enu < 600){ index_Enu = 4;
    }else if (Enu < 700){ index_Enu = 5;
    }else if (Enu < 800){ index_Enu = 6;
    }else if (Enu < 900){ index_Enu = 7;
    }else if (Enu < 1000){ index_Enu = 8;
    }else if (Enu < 1200){ index_Enu = 9;
    }else if (Enu < 1400){ index_Enu = 10;
    }else if (Enu < 1600){ index_Enu = 11;
    }else if (Enu < 1800){ index_Enu = 12;
    }else if (Enu < 2000){ index_Enu = 13;
    }else if (Enu < 2200){ index_Enu = 14;
    }else if (Enu < 2500){ index_Enu = 15;
    }else{ index_Enu = 16;
    }

    if (costheta < -0.5){  // 8 of them ...
    }else if (costheta < 0.){  index_costheta = 1;
    }else if (costheta < 0.2){ index_costheta = 2;
    }else if (costheta < 0.4){ index_costheta = 3;
    }else if (costheta < 0.6){ index_costheta = 4;
    }else if (costheta < 0.8){ index_costheta = 5;
    }else if (costheta < 0.9){ index_costheta = 6;
    }else{index_costheta = 7;
    }

    if (Pmu < 100){ // 14 ...
    }else if (Pmu < 200){ index_Pmu = 1;
    }else if (Pmu < 300){ index_Pmu = 2;
    }else if (Pmu < 400){ index_Pmu = 3;
    }else if (Pmu < 500){ index_Pmu = 4;
    }else if (Pmu < 600){ index_Pmu = 5;
    }else if (Pmu < 700){ index_Pmu = 6;
    }else if (Pmu < 800){ index_Pmu = 7;
    }else if (Pmu < 1000){ index_Pmu = 8;
    }else if (Pmu < 1200){ index_Pmu = 9;
    }else if (Pmu < 1400){ index_Pmu = 10;
    }else if (Pmu < 1600){ index_Pmu = 11;
    }else if (Pmu < 2000){ index_Pmu = 12;
    }else{ index_Pmu = 13;
    }

    int index = index_Enu * 8 * 14 + index_costheta * 14 + index_Pmu;

    //    std::cout << index << " " << spec_weights.at(index).at(0) << std::endl;
    
    return spec_weights.at(index);
    
  }else{
    std::vector<float> abc;
    abc.resize(spec_weights.front().size(),0);
    return abc;
  }
}


void LEEana::CovMatrix::add_osc_config(TString osc_ch_filename, TString osc_pars_filename ){
  flag_osc = true;
  std::ifstream infile1(osc_ch_filename);
  TString temp;
  while(!infile1.eof()){
    infile1 >> temp;
    if (temp == "End") break;
    osc_signal_ch_names.insert(temp);
  }
  std::ifstream infile2(osc_pars_filename);
  infile2 >> osc_par_delta_m2_eV2 >> osc_par_sin22theta_ee;
  infile2 >> osc_par_delta_m2_41_eV2 >> osc_par_sin2_2theta_14 >> osc_par_sin2_theta_24 >> osc_par_sin2_theta_34;

  std::cout << "===> Oscillation Mode: " << std::endl;
  std::cout << "===> Dm2: " << osc_par_delta_m2_eV2 << " eV2 " << std::endl;
  std::cout << "===> sin22theta_ee: " << osc_par_sin22theta_ee << std::endl;
  std::cout << "=====> Full oscillation mode: " << std::endl;
  std::cout << "=====> Dm2_41: " << osc_par_delta_m2_41_eV2 << " eV2 " << std::endl;
  std::cout << "=====> sin2_2theta_14: " << osc_par_sin2_2theta_14 << std::endl;
  std::cout << "=====> sin2_theta_24: " << osc_par_sin2_theta_24 << std::endl;
  std::cout << "=====> sin2_theta_34: " << osc_par_sin2_theta_34 << std::endl;
  
}

void LEEana::CovMatrix::add_rw_config( int run ){
  if (run==18) rw_type=1;
  if (run==19) rw_type=2;
  if (run==0) rw_type=3;
}

void LEEana::CovMatrix::print_rw(std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > rw_info ){
  std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double> > rw_info_i;

  if(std::get<0>(rw_info)){//Are you applying any reweighting?
    for(size_t rw=0; rw<std::get<1>(rw_info).size(); rw++){
      rw_info_i = std::get<1>(rw_info)[rw];
      if(std::get<0>(rw_info_i)) std::cout<<"Applying reweighting to "<<std::get<1>(rw_info_i)<<" in "<<std::get<2>(rw_info_i)<<std::endl;
    }
  }
}

bool LEEana::CovMatrix::is_osc_channel(TString ch_name){
  if (osc_signal_ch_names.find(ch_name) == osc_signal_ch_names.end())
    return false;
  else
    return true;
}

double LEEana::CovMatrix::get_osc_weight(EvalInfo& eval, PFevalInfo& pfeval){
  double weight = 1.0;
  // only support nueCC disappearance now ...
  if (fabs(eval.truth_nuPdg)==12){  
    weight = 1 - osc_par_sin22theta_ee * pow(TMath::Sin(1.267 * osc_par_delta_m2_eV2 * (pfeval.mcflux_gen2vtx + pfeval.mcflux_dk2gen)/pfeval.truth_nu_momentum[3]/1000.),2);
  }
  
  if(osc_par_delta_m2_41_eV2>=0){
  // nueCC to nueCC 
  if (fabs(eval.truth_nuPdg)==12 && fabs(pfeval.mcflux_ntype)==12 && eval.truth_isCC){ 
    weight = 0.0; 
  }
  // numuCC to numuCC
  else if (fabs(eval.truth_nuPdg)==14 && fabs(pfeval.mcflux_ntype)==14 && eval.truth_isCC){
    weight = 0.0;
  }
  // numuCC to nueCC 
  else if (fabs(eval.truth_nuPdg)==12 && fabs(pfeval.mcflux_ntype)==14 && eval.truth_isCC){ 
    weight = 0.006; 
  }
  // nueCC to numuCC 
  else if (fabs(eval.truth_nuPdg)==14 && fabs(pfeval.mcflux_ntype)==12 && eval.truth_isCC){ 
    weight = 0.; 
  }
  // 1 - (nue to sterile) NC  
  else if (fabs(eval.truth_nuPdg)==12 && fabs(pfeval.mcflux_ntype)==12 && !eval.truth_isCC){ 
    weight = 0.0; 
  }
  // 1 - (numu to sterile) NC  
  else if (fabs(eval.truth_nuPdg)==14 && fabs(pfeval.mcflux_ntype)==14 && !eval.truth_isCC){ 
    weight = 0.0; 
  }
  else{
    std::cout << "Error: unknown file/osc type: final nu flavor: "<<eval.truth_nuPdg<<" initial nu flavor: "<<pfeval.mcflux_ntype<<" CC interaction: "<<eval.truth_isCC<<std::endl;
    exit(1);
  }
  }
  return weight;
}


void LEEana::CovMatrix::add_xs_config(TString xs_ch_filename , TString xs_real_bin_filename ){
  std::ifstream infile1(xs_ch_filename);
  TString temp;
  while(!infile1.eof()){
    infile1 >> temp;
    if (temp == "End") break;
    xs_signal_ch_names.insert(temp);
  }
  //  std::cout << xs_signal_ch_names.size() << std::endl;

  std::ifstream infile2(xs_real_bin_filename);
  int bin_no;
  TString cut;
  double constant;
  double pot_err;
  double target_err;
  int tmp_cut_file;
  
  while(!infile2.eof()){
    infile2 >> bin_no >> tmp_cut_file >> cut >> constant >> pot_err >> target_err;
    if (bin_no==-1) break;
    cut_file = tmp_cut_file;
    map_xs_bin_cut[bin_no] = cut;
    map_cut_xs_bin[cut] = bin_no;
    map_xs_bin_constant[bin_no] = constant;
    map_xs_bin_errs[bin_no] = std::make_pair(pot_err, target_err);
  }
  
}

int LEEana::CovMatrix::get_xs_nsignals(){
  return map_xs_bin_cut.size();
}

int LEEana::CovMatrix::get_xs_nmeas(){
  int ncount = 0;
  for (auto it = xs_signal_ch_names.begin(); it != xs_signal_ch_names.end(); it++){
    int ch = map_name_ch[*it];
    ncount += std::get<2>(map_ch_hist[ch]) + 1; 
  }
  return ncount;
}

bool LEEana::CovMatrix::is_xs_chname(TString name){
  if (xs_signal_ch_names.find(name)==xs_signal_ch_names.end()){
    return false;
  }else{
    return true;
  }
}

void LEEana::CovMatrix::gen_xs_cov_matrix(int run, std::map<int, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_covch_hists, std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists, TVectorD* vec_mean,  TMatrixD* cov_xs_mat, TVectorD* vec_signal, TMatrixD* mat_R){
  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ; 
  std::map<int, TString> map_no_histoname; 
  
  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    
    if (period != run) continue;
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);
    
    for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
      int ch = map_name_ch[std::get<5>(*it1)];
      int obsch = get_obsch_name(std::get<5>(*it1));
      int covch = get_covch_name(std::get<5>(*it1));
      int flag_lee = std::get<7>(map_ch_hist[ch]);
      TString histoname = std::get<0>(*it1);
      //TH1F *htemp = map_histoname_hist[histoname];
      //
      map_histoname_infos[histoname] = std::make_tuple(ncount, covch, flag_lee, input_filename);
      map_no_histoname[ncount] = histoname;
      ncount ++;

      // std::cout << histoname << " " << obsch << " " << covch << " " << flag_lee << std::endl;
    }
  }

  std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > > map_passed_events; // last one is variable name ...
  std::map<TString, double> map_filename_pot;
  std::vector<int> max_lengths;
  std::vector<int> max_sup_lengths;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    //int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    if (period != run) continue;

    //map_all_events[input_filename];
    std::pair<std::vector<int>, std::vector<int>> lengths_pair = get_events_weights_xs(input_filename, map_passed_events, map_filename_pot, map_histoname_infos);
    std::vector<int> lengths = lengths_pair.first;
    std::vector<int> sup_lengths = lengths_pair.second;

    if (lengths.size() > max_lengths.size()) max_lengths.resize(lengths.size());
    for (size_t i = 0; i != lengths.size();i++){
      if (lengths.at(i) > max_lengths.at(i)) max_lengths.at(i) = lengths.at(i);
    }
    
    if (sup_lengths.size() > max_sup_lengths.size()) max_sup_lengths.resize(sup_lengths.size());
    for (size_t i = 0; i != sup_lengths.size();i++){
      if (sup_lengths.at(i) > max_sup_lengths.at(i)) max_sup_lengths.at(i) = sup_lengths.at(i);
    }
    
    //std::cout << input_filename << " " << lengths.size() << std::endl;
  // for (size_t j = 0; j!=max_lengths.size(); j++){
  //   int nsize = max_lengths.at(j);
  //   int sup_nsize = max_sup_lengths.at(j);
  //   std::cout << "knob# " << j << " input filename: " << input_filename << " nsize: " << nsize << " sup_nsize: " << sup_nsize << std::endl;
  // }

  }

  double data_pot = 5e19;
  const int rows = cov_xs_mat->GetNcols();
  float x[rows];
  (*cov_xs_mat).Zero();

  int acc_no = 0;
  // build covariance matrix ...
  
  for (size_t j = 0; j!=max_lengths.size(); j++){ // j: index of knobs
    int nsize = max_lengths.at(j);
    int sup_nsize = max_sup_lengths.at(j);

    TMatrixD temp_mat(rows, rows);
    temp_mat.Zero(); 
    
    // std::cout << "nsize: " << nsize << " sup_nsize: " << sup_nsize << std::endl;
    for (int i=0;i!=nsize;i++){ // i: index of universes of a knob
      //      if (i>=10) continue;
      
      for (int k = 0; k!= rows;k++){ // k: index of kinematic variable bins
     	x[k] = 0;
      }
      // if (nsize==600) {
      //   std::cout << "[wg] knob: " << j << " universe: " << i << " with nsize=600" << std::endl;
      // }
      fill_xs_histograms(j, max_lengths.size(), acc_no, i, nsize,  map_passed_events, map_histoname_infos, map_no_histoname, map_histoname_hists);

      // merge histograms according to POTs ...
      for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
     	// std::cout << it->first << std::endl;
     	int covch = it->first;
	auto tmp_results  = map_covch_hists[covch];
	TH1F *hpred = std::get<0>(tmp_results);
	TH1F *hsigma = std::get<1>(tmp_results);
	TH1F *hsigmabar = std::get<2>(tmp_results);
	TH2F *hR = std::get<3>(tmp_results);
	int num = std::get<4>(tmp_results);
	hpred->Reset();
	if (num!=1){
	  hsigma->Reset();
	  hsigmabar->Reset();
	  hR->Reset();
	}
	
     	for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	  TH1F *htemp = (TH1F*)hpred->Clone("htemp");
	  htemp->Reset();
	  TH1F *htemp1 = 0;
	  TH1F *htemp2 = 0;
	  TH2F *htemp3 = 0;
	  if (num!=1){
	    htemp1 = (TH1F*)hsigma->Clone("htemp1");
	    htemp2 = (TH1F*)hsigmabar->Clone("htemp2");
	    htemp3 = (TH2F*)hR->Clone("htemp3");
	    htemp1->Reset();
	    htemp2->Reset();
	    htemp3->Reset();
	  }
	  std::map<int, double> temp_map_mc_acc_pot;
	  
     	  for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
     	    TString histoname = (*it2).first;
     	    TString input_filename = map_histogram_inputfile[histoname];
     	    auto it3 = map_inputfile_info.find(input_filename);
     	    int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
     	    int norm_period = std::get<6>(it3->second);
     	    double mc_pot = map_filename_pot[input_filename];
     	    //std::cout << mc_pot << std::endl;
     	    if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
     	      temp_map_mc_acc_pot[norm_period] = mc_pot;
     	    }else{
     	      temp_map_mc_acc_pot[norm_period] += mc_pot;
     	    }
     	  }
	  
     	  for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
     	    TString histoname = (*it2).first;
     	    TString input_filename = map_histogram_inputfile[histoname];
     	    auto it3 = map_inputfile_info.find(input_filename);
     	    int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
     	    int norm_period = std::get<6>(it3->second);
     	    data_pot = std::get<5>(map_inputfile_info[input_filename]);
     	    double ratio = data_pot/temp_map_mc_acc_pot[norm_period];
	    auto tmp_hists = map_histoname_hists[histoname];
	    TH1F *hmc = std::get<0>(tmp_hists);
	    TH1F *hmc1 = std::get<1>(tmp_hists);
	    TH1F *hmc2 = std::get<2>(tmp_hists);
	    TH2F *hmc3 = std::get<3>(tmp_hists);
	    
	    htemp->Add(hmc, ratio);
	    if (num !=1){
	      htemp1->Add(hmc1,ratio);
	      htemp2->Add(hmc2,ratio);
	      htemp3->Add(hmc3,ratio);
	    }
	    //	std::cout << covch << " " << histoname << " " << ratio << std::endl;
	  }
	  
     	  hpred->Add(htemp);
     	  delete htemp;
	  if (num != 1){
	    hsigma->Add(htemp1);
	    hsigmabar->Add(htemp2);
	    hR->Add(htemp3);
	    delete htemp1;
	    delete htemp2;
	    delete htemp3;
	  }
	}

  // // hack: add a factor onto the CV hpred (reco) and hsigmabar (true)
  // if (nsize==600 and (covch==1 or covch==2)) {
  //     std::cout << "hack CV, covch: " << covch << " hpred: " << hpred->GetNbinsX() << " hsigmabar: " << hsigmabar->GetNbinsX() << std::endl;
  //     // sanity check
  //     double reco1[26] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // FC
  //     double signal1[10] = {1,1,1,1,1,1,1,1,1,1};
  //     double reco2[26] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // PC
  //     double signal2[10] = {1,1,1,1,1,1,1,1,1,1};
  //     for (int i=0; i<=hpred->GetNbinsX(); i++) { // overflow bin
  //       if (covch==1) hpred->SetBinContent(i+1, reco1[i] * hpred->GetBinContent(i+1));
  //       else if (covch==2) hpred->SetBinContent(i+1, reco2[i] * hpred->GetBinContent(i+1));
  //     }
  //     for (int i=0; i<hsigmabar->GetNbinsX(); i++) {
  //       if (covch==1) hsigmabar->SetBinContent(i+1, signal1[i] * hsigmabar->GetBinContent(i+1));
  //       else if (covch==2) hsigmabar->SetBinContent(i+1, signal2[i] * hsigmabar->GetBinContent(i+1));
  //     }
  // } // end of hack

     	int start_bin = map_covch_startbin[covch];
     	for (int k=0;k!=hpred->GetNbinsX()+1;k++){
	  if (num == 1){
	    x[start_bin+k] = hpred->GetBinContent(k+1) ;
	  }else{
	    x[start_bin+k] = - hpred->GetBinContent(k+1);
	  }
     	  //	  std::cout << i << " " << x[start_bin+i] << std::endl;
     	}
	if (num!=1){

	  for (int k=0;k!=hsigma->GetNbinsX();k++){
	    //	    int bin = std::round(hsigma->GetBinCenter(k+1));
	    for (int j=0; j!=hpred->GetNbinsX()+1;j++){
	      x[start_bin+j] += hR->GetBinContent(j+1,k+1)/hsigma->GetBinContent(k+1)*hsigmabar->GetBinContent(k+1); // hR(j,k)/hsigma(k): response in new universe
	    }
	  }
	  
    //// hack: save CV for UBGenieAll
    // if (nsize==600 and i==0) { // only save CV once
    //   auto ofile = new TFile("wgu_cv.root", "update");
    //   hpred->Write();
    //   hsigmabar->Write();
    //   ofile->Close();
    // }
    //// hack: save reweighted for UBGenieAll
    //// change accordingly in fill_xs_histograms()
    // auto ofile = TFile::Open("wgu_all_genie.root", "update");
    // auto hpred1 = (TH1F*)hpred->Clone(Form("%s_%d", hpred->GetName(), i));
    // auto hsigma1 = (TH1F*)hsigma->Clone(Form("%s_%d", hsigma->GetName(), i));
    // hpred1->Write();
    // hsigma1->Write();
    // ofile->Close();
	}
	
      }
      
      // add covariance matrix ...
      for (size_t n = 0;n!=rows; n++){
     	for (size_t m =0; m!=rows;m++){
     	  temp_mat(n,m) += x[n] * x[m];
     	}
      }  
    } // i

    
    if (nsize==2){  // second check
      temp_mat *= 1./sup_nsize;
    }else{
      temp_mat *= 1./nsize;
    }

   
    
    // //std::cout << j << " " << temp_mat(26+26+5,26+26+26) << std::endl;
    
    (*cov_xs_mat) +=  temp_mat;
    
    acc_no += nsize;
  }
  
  
  
  
  // build CV ...
  for (int i=0;i!=rows;i++){
    (*vec_mean)(i) = 0;
  }
  fill_xs_histograms(map_passed_events, map_histoname_infos, map_no_histoname, map_histoname_hists);

  
  // merge histograms according to POTs ...
  for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int covch = it->first;
    auto tmp_results  = map_covch_hists[covch];
    TH1F *hpred = std::get<0>(tmp_results);
    TH1F *hsigma = std::get<1>(tmp_results);
    TH1F *hsigmabar = std::get<2>(tmp_results);
    TH2F *hR = std::get<3>(tmp_results);
    int num = std::get<4>(tmp_results);
    hpred->Reset();
    if (num!=1){
      hsigma->Reset();
      hsigmabar->Reset();
      hR->Reset();
    }
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      TH1F *htemp1 = 0;
      TH1F *htemp2 = 0;
      TH2F *htemp3 = 0;
      if (num!=1){
	htemp1 = (TH1F*)hsigma->Clone("htemp1");
	htemp2 = (TH1F*)hsigmabar->Clone("htemp2");
	htemp3 = (TH2F*)hR->Clone("htemp3");
	htemp1->Reset();
	htemp2->Reset();
	htemp3->Reset();
      }
      std::map<int, double> temp_map_mc_acc_pot;
      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
     	TString histoname = (*it2).first;
     	TString input_filename = map_histogram_inputfile[histoname];
     	auto it3 = map_inputfile_info.find(input_filename);
     	int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
     	int norm_period = std::get<6>(it3->second);
     	double mc_pot = map_filename_pot[input_filename];
     	//std::cout << mc_pot << std::endl;
     	if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
     	  temp_map_mc_acc_pot[norm_period] = mc_pot;
     	}else{
     	  temp_map_mc_acc_pot[norm_period] += mc_pot;
     	}
     	//std::cout << histoname << " " << input_filename << " " << mc_pot << " " << period << std::endl;
      }

      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
     	TString histoname = (*it2).first;
     	TString input_filename = map_histogram_inputfile[histoname];
     	auto it3 = map_inputfile_info.find(input_filename);
     	int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
     	int norm_period = std::get<6>(it3->second);
     	data_pot = std::get<5>(map_inputfile_info[input_filename]);
     	double ratio = data_pot/temp_map_mc_acc_pot[norm_period];
	auto tmp_hists = map_histoname_hists[histoname];
	TH1F *hmc = std::get<0>(tmp_hists);
	TH1F *hmc1 = std::get<1>(tmp_hists);
	TH1F *hmc2 = std::get<2>(tmp_hists);
	TH2F *hmc3 = std::get<3>(tmp_hists);

     	htemp->Add(hmc, ratio);
	if (num !=1){
	  htemp1->Add(hmc1,ratio);
	  htemp2->Add(hmc2,ratio);
	  htemp3->Add(hmc3,ratio);
	  //std::cout << hmc1->GetSum() << " " << ratio << std::endl;
	}
	
     	//	std::cout << covch << " " << histoname << " " << ratio << " " << data_pot << std::endl;
      }
      
      hpred->Add(htemp);
      delete htemp;
      if (num != 1){
	//	std::cout << num << " " << htemp1->GetSum() << std::endl;
	hsigma->Add(htemp1);
	hsigmabar->Add(htemp2);
	hR->Add(htemp3);
	delete htemp1;
	delete htemp2;
	delete htemp3;
      }
    }
    
    int start_bin = map_covch_startbin[covch];
    for (int k=0;k!=hpred->GetNbinsX()+1;k++){
      (*vec_mean)(start_bin+k) = hpred->GetBinContent(k+1) ;
      //	  std::cout << i << " " << x[start_bin+i] << std::endl;
    }
    
    if (num!=1){
      // vec_signal, mat_R
      for (int k=0;k!=hsigma->GetNbinsX();k++){
	int bin = std::round(hsigma->GetBinCenter(k+1));
	(*vec_signal)(k) = hsigma->GetBinContent(k+1)/map_xs_bin_constant[bin];
      }
      // mat_R
      // loop real signal bin ...
      for (int k=0;k!=hsigma->GetNbinsX();k++){
	int bin = std::round(hsigma->GetBinCenter(k+1));
	for (int j=0; j!=hpred->GetNbinsX()+1;j++){
	  (*mat_R)(start_bin+j,k) = hR->GetBinContent(j+1,k+1)/(hsigma->GetBinContent(k+1)/map_xs_bin_constant[bin]);
	}
      }
    }
    
  }
  {
    // add additional uncertainties ...
    // POT ...
    auto it = map_xs_bin_errs.begin();
    double pot_err = (it->second).first;
    double target_err = (it->second).second;
    for (int i=0; i!= (*cov_xs_mat).GetNrows(); i++){
      for (int j=0;j!=(*cov_xs_mat).GetNcols(); j++){
    	(*cov_xs_mat)(i,j) += (*vec_mean)(i) * (*vec_mean)(j) *pot_err * pot_err; 
      }
    }
    // Target Nucleons
    std::vector<int> bins;
    for (auto it = xs_signal_ch_names.begin(); it != xs_signal_ch_names.end(); it++){
      TString ch_name = *it;
      int ch = map_name_ch[ch_name];
      int covch = map_ch_covch[ch];
      int nbin = map_covch_nbin[covch];
      int start_bin = map_covch_startbin[covch];
      for (int i=start_bin; i!= start_bin+nbin;i++){
     	bins.push_back(i);
      }
    }

    for (size_t i=0;i!=bins.size();i++){
      for (size_t j=0;j!=bins.size();j++){
    	(*cov_xs_mat)(bins.at(i),bins.at(j)) += (*vec_mean)(bins.at(i)) * (*vec_mean)(bins.at(j)) *target_err * target_err; 
      }
    }
    
  }
  
  
}

void LEEana::CovMatrix::fill_pred_R_signal(int run, TMatrixD* mat_R, TVectorD* vec_signal,  std::map<int, double>& map_data_period_pot, std::map<TString, std::tuple<TH1F*, TH2F*, double> >& map_name_xs_hists){

  // merge histograms according to POTs ...
  for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
    int covch = it->first;
    TH1F *hsigma = 0;
    TH2F *hR = 0;
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      std::map<int, double> temp_map_mc_acc_pot;
      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
      	TString histoname = (*it2).first;
	auto it4 = map_name_xs_hists.find(histoname);
	if (it4 == map_name_xs_hists.end()) continue;

	TH1F *h1 = std::get<0>(it4->second);
	TH2F *h2 = std::get<1>(it4->second);
	if (hsigma == 0){
	  hsigma = (TH1F*)h1->Clone("hsigma");
	  hR = (TH2F*)h2->Clone("hR");
	  hsigma->Reset();
	  hR->Reset();
	}
	TString input_filename = map_histogram_inputfile[histoname];
      	auto it3 = map_inputfile_info.find(input_filename);
      	int period = std::get<1>(it3->second);

	if (run !=0 && run !=period) continue;
	// if (period != run) continue; // skip ...
      	//int norm_period = std::get<6>(it3->second);

	double mc_pot = std::get<2>(it4->second);
	
      	if (temp_map_mc_acc_pot.find(period) == temp_map_mc_acc_pot.end()){
      	  temp_map_mc_acc_pot[period] = mc_pot;
      	}else{
      	  temp_map_mc_acc_pot[period] += mc_pot;
	}
      	//std::cout << histoname << " " << input_filename << " " << mc_pot << " " << period << std::endl;
      }
      if (hsigma ==0 ) continue;
      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
      	TString histoname = (*it2).first;
	auto it4 = map_name_xs_hists.find(histoname);
	if (it4 == map_name_xs_hists.end()) continue;
	TH1F *h1 = std::get<0>(it4->second);
	TH2F *h2 = std::get<1>(it4->second);
	
      	TString input_filename = map_histogram_inputfile[histoname];
      	auto it3 = map_inputfile_info.find(input_filename);
      	int period = std::get<1>(it3->second);
	if (period != run && run != 0) continue; // skip ...
      	double data_pot = map_data_period_pot[period];
      	double ratio = data_pot/temp_map_mc_acc_pot[period];

	//	std::cout << h1->GetSum() << " " << ratio << std::endl;
      	hsigma->Add(h1, ratio);
	hR->Add(h2, ratio);
      	//	std::cout << covch << " " << histoname << " " << ratio << " " << data_pot << std::endl;
      }
      
      
    }
    int start_bin = map_covch_startbin[covch];

    if (hsigma==0) continue;
    for (int k=0;k!=hsigma->GetNbinsX();k++){
      int bin = std::round(hsigma->GetBinCenter(k+1));
      (*vec_signal)(k) = hsigma->GetBinContent(k+1)/map_xs_bin_constant[bin];
    }
    // mat_R
    // loop real signal bin ...
    for (int k=0;k!=hsigma->GetNbinsX();k++){
      int bin = std::round(hsigma->GetBinCenter(k+1));
      for (int j=0; j!=hR->GetNbinsX()+1;j++){
	(*mat_R)(start_bin+j,k) = hR->GetBinContent(j+1,k+1)/(hsigma->GetBinContent(k+1)/map_xs_bin_constant[bin]);
      }
    }
    
    delete hsigma;
    delete hR;
  }
    

}


void LEEana::CovMatrix::fill_xs_histograms(int num, int tot_num, int acc_no, int no, int tot_no, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists){
  // int jinput = num; // backup the input
  for (auto it = map_histoname_hists.begin(); it != map_histoname_hists.end(); it++){
    int num = std::get<4>(it->second);
    TH1F *h1 = std::get<0>(it->second);
    TH1F *h2 = std::get<1>(it->second);
    TH1F *h3 = std::get<2>(it->second);
    TH2F *h4 = std::get<3>(it->second);
    if (num == 1){
      h1->Reset();
    }else{
      h1->Reset();
      h2->Reset();
      h3->Reset();
      h4->Reset();
    }
  }

  // std::cout << "acc_no: " << acc_no << " no: " << no << std::endl;
  for (auto it = map_passed_events.begin(); it != map_passed_events.end(); it++){
    TString filename = it->first;
    // loop over events ...
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      float weight = std::get<0>(*it1);
      float weight_lee = std::get<1>(*it1);
      if (std::get<3>(*it1).size() != tot_num) std::cout << "Incorrect Match Sys No! " << std::endl;
      if (std::get<3>(*it1).at(num) != tot_no) std::cout << "Mismatch No of Universe! " << std::endl;
      float rel_weight_diff = std::get<2>(*it1).at(acc_no+no);
      for (auto it2 = std::get<4>(*it1).begin(); it2 != std::get<4>(*it1).end(); it2++){
   	int no = std::get<0>(*it2);
   	float val = std::get<1>(*it2);
	bool flag_pass = std::get<2>(*it2);
	int nsignal_bin = std::get<3>(*it2);

   	TString histoname = map_no_histoname[no];
	auto tmp_hists = map_histoname_hists[histoname];
	TH1F *h1 = std::get<0>(tmp_hists);
	TH1F *h2 = std::get<1>(tmp_hists);
	TH1F *h3 = std::get<2>(tmp_hists);
	TH2F *h4 = std::get<3>(tmp_hists);
	int num = std::get<4>(tmp_hists);
	int flag_lee = std::get<2>(map_histoname_infos[histoname]);

   	if (std::isnan(rel_weight_diff) || std::isinf(rel_weight_diff)) continue;
	// seems to have extremely small cv weight
	if (fabs(rel_weight_diff)>100) continue;

	if (num==1){
	  if (flag_lee){
	    if (flag_pass) h1->Fill(val, rel_weight_diff * weight * weight_lee);
	  }else{
	    if (flag_pass) h1->Fill(val,  rel_weight_diff * weight); // rel = (genie-cv)/cv
	  }
	}else{
	  if (nsignal_bin != -1){
	    if (flag_lee){
	      if (flag_pass) h1->Fill(val, weight * weight_lee); // CV as the central one ...

	      h2->Fill(nsignal_bin, (1+rel_weight_diff) * weight*weight_lee); // signal
	      h3->Fill(nsignal_bin, weight*weight_lee); // nominal ...
	      
	      if (flag_pass) h4->Fill(val, nsignal_bin, (1+rel_weight_diff) * weight*weight_lee);
	    }else{
	      if (flag_pass) h1->Fill(val, weight); // CV as the central one
        // if (flag_pass) h1->Fill(val, (1+rel_weight_diff) * weight); // hack: alternative CV as the central one
	      
	      h2->Fill(nsignal_bin, (1+rel_weight_diff) * weight); // signal
	      h3->Fill(nsignal_bin, weight); // nominal ...
	      
	      if (flag_pass) h4->Fill(val, nsignal_bin, (1+rel_weight_diff) * weight);
	    }
	  }else{
	    std::cout << "Something wrong: cut/channel mismatch !" << std::endl;
	  }
	} // else	
      }
    }
  }
  
}


void LEEana::CovMatrix::fill_xs_histograms(std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, std::tuple<TH1F*, TH1F*, TH1F*, TH2F*, int> >& map_histoname_hists){
  for (auto it = map_histoname_hists.begin(); it != map_histoname_hists.end(); it++){
    int num = std::get<4>(it->second);
    TH1F *h1 = std::get<0>(it->second);
    TH1F *h2 = std::get<1>(it->second);
    TH1F *h3 = std::get<2>(it->second);
    TH2F *h4 = std::get<3>(it->second);
    if (num == 1){
      h1->Reset();
    }else{
      h1->Reset();
      h2->Reset();
      h3->Reset();
      h4->Reset();
    }
  }
  
  for (auto it = map_passed_events.begin(); it != map_passed_events.end(); it++){
    TString filename = it->first;
    // loop over events ...
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      float weight = std::get<0>(*it1);
      float weight_lee = std::get<1>(*it1);

      for (auto it2 = std::get<4>(*it1).begin(); it2 != std::get<4>(*it1).end(); it2++){
	int no = std::get<0>(*it2);
	float val = std::get<1>(*it2);
	bool flag_pass = std::get<2>(*it2);
	int nsignal_bin = std::get<3>(*it2);

	
	TString histoname = map_no_histoname[no];
	auto tmp_hists = map_histoname_hists[histoname];
	TH1F *h1 = std::get<0>(tmp_hists);
	TH1F *h2 = std::get<1>(tmp_hists);
	TH1F *h3 = std::get<2>(tmp_hists);
	TH2F *h4 = std::get<3>(tmp_hists);
	int num = std::get<4>(tmp_hists);
	int flag_lee = std::get<2>(map_histoname_infos[histoname]);

	//	if (no==0) std::cout << "Xin: " << " " << flag_pass << " " << nsignal_bin << " " << weight << " " << std::endl;
	
	if (num==1){
	  if (flag_lee){
	    if (flag_pass) h1->Fill(val, weight * weight_lee);
	  }else{
	    if (flag_pass) h1->Fill(val,  weight);
	  }
	}else{
	  if (nsignal_bin != -1){
	    if (flag_lee){
	      if (flag_pass) h1->Fill(val, weight * weight_lee);

	      h2->Fill(nsignal_bin, weight*weight_lee); 
	      h3->Fill(nsignal_bin, weight*weight_lee); // nominal ...
	      
	      if (flag_pass) h4->Fill(val, nsignal_bin, weight*weight_lee);
	    }else{
	      if (flag_pass) h1->Fill(val, weight);
	      
	      h2->Fill(nsignal_bin, weight);
	      h3->Fill(nsignal_bin, weight); // nominal ...
	      
	      if (flag_pass) h4->Fill(val, nsignal_bin, weight);
	    }
	  }else{
	    std::cout << "Something wrong: cut/channel mismatch !" << std::endl;
	  }
	} // else
	
      }
    }
  }

  // for (auto it = map_histoname_hists.begin(); it != map_histoname_hists.end(); it++){
  //   int num = std::get<4>(it->second);
  //   TH1F *h1 = std::get<0>(it->second);
  //   TH1F *h2 = std::get<1>(it->second);
  //   TH1F *h3 = std::get<2>(it->second);
  //   TH2F *h4 = std::get<3>(it->second);
  //   //    if (num==4)
  //   //  std::cout << h1->GetSum() << " " << h2->GetSum() << " " << h3->GetSum() << " " << h4->GetSum() << std::endl;
     
  // }


  
}

std::pair<std::vector<int>, std::vector<int> > LEEana::CovMatrix::get_events_weights_xs(TString input_filename, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > > >& map_passed_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){
  TFile *file = new TFile(input_filename);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
  
  
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars, tagger,2 );
  set_tree_address(T_eval, eval);
  set_tree_address(T_PFeval, pfeval);
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  // total POT calculations ...
  map_filename_pot[input_filename] = total_pot;

  // fill histogram ...
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_0",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_1",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_2",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_3",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_4",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_5",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_6",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_7",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_8",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_9",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_10",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_11",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_12",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_13",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_14",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_15",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_16",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_17",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_18",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_19",1);
  T_BDTvars->SetBranchStatus("mip_energy",1);
  T_BDTvars->SetBranchStatus("mip_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_tracks",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_showers",1);
  T_BDTvars->SetBranchStatus("gap_n_bad",1);
  T_BDTvars->SetBranchStatus("spt_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);

  if (tagger.flag_nc_gamma_bdt){
    T_BDTvars->SetBranchStatus("nc_delta_score", 1);
    T_BDTvars->SetBranchStatus("nc_pio_score", 1);
  }
  
  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);
  T_eval->SetBranchStatus("match_energy",1);
  T_eval->SetBranchStatus("match_isFC",1);
  T_eval->SetBranchStatus("match_found",1);
  if (T_eval->GetBranch("match_found_asInt")) T_eval->SetBranchStatus("match_found_asInt",1); 
  T_eval->SetBranchStatus("stm_eventtype",1);
  T_eval->SetBranchStatus("stm_lowenergy",1);
  T_eval->SetBranchStatus("stm_LM",1);
  T_eval->SetBranchStatus("stm_TGM",1);
  T_eval->SetBranchStatus("stm_STM",1);
  T_eval->SetBranchStatus("stm_FullDead",1);
  T_eval->SetBranchStatus("stm_clusterlength",1);
  
  T_eval->SetBranchStatus("weight_spline",1);
  T_eval->SetBranchStatus("weight_cv",1);
  T_eval->SetBranchStatus("weight_lee",1);
  T_eval->SetBranchStatus("weight_change",1);
  // MC enable truth information ...
  T_eval->SetBranchStatus("truth_isCC",1);
  T_eval->SetBranchStatus("truth_nuPdg",1);
  T_eval->SetBranchStatus("truth_vtxInside",1);
  T_eval->SetBranchStatus("truth_nuEnergy",1);
  T_eval->SetBranchStatus("truth_energyInside",1);
  T_eval->SetBranchStatus("truth_vtxX",1);
  T_eval->SetBranchStatus("truth_vtxY",1);
  T_eval->SetBranchStatus("truth_vtxZ",1);
  T_eval->SetBranchStatus("match_completeness_energy",1);
 
  
  T_KINEvars->SetBranchStatus("*",0);
  T_KINEvars->SetBranchStatus("kine_reco_Enu",1);
  T_KINEvars->SetBranchStatus("kine_energy_particle",1);
  T_KINEvars->SetBranchStatus("kine_particle_type",1);
  T_KINEvars->SetBranchStatus("kine_energy_info",1);
  T_KINEvars->SetBranchStatus("kine_energy_included",1);
  T_KINEvars->SetBranchStatus("kine_reco_add_energy",1);
  T_KINEvars->SetBranchStatus("kine_pio_mass",1);
  T_KINEvars->SetBranchStatus("kine_pio_flag",1);
  T_KINEvars->SetBranchStatus("kine_pio_vtx_dis",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_angle",1);
  if (T_KINEvars->GetBranch("vlne_v4_numu_full_primaryE")) {
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_totalE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_totalE",1);
  }

  T_PFeval->SetBranchStatus("*",0);
  T_PFeval->SetBranchStatus("reco_nuvtxX",1);
  T_PFeval->SetBranchStatus("reco_nuvtxY",1);
  T_PFeval->SetBranchStatus("reco_nuvtxZ",1);
  T_PFeval->SetBranchStatus("reco_showervtxX",1);
  T_PFeval->SetBranchStatus("reco_showervtxY",1);
  T_PFeval->SetBranchStatus("reco_showervtxZ",1);
  T_PFeval->SetBranchStatus("reco_muonMomentum",1);
  T_PFeval->SetBranchStatus("reco_showerKE",1);
  T_PFeval->SetBranchStatus("nuvtx_diff",1);
  T_PFeval->SetBranchStatus("showervtx_diff",1);
  T_PFeval->SetBranchStatus("muonvtx_diff",1);
  T_PFeval->SetBranchStatus("truth_muonMomentum",1);
  T_PFeval->SetBranchStatus("truth_pio_energy_1",1);
  T_PFeval->SetBranchStatus("truth_pio_energy_2",1);
  T_PFeval->SetBranchStatus("truth_pio_angle",1);
  T_PFeval->SetBranchStatus("truth_NprimPio",1);
  if (pfeval.flag_NCDelta){
    
    T_PFeval->SetBranchStatus("truth_NCDelta",1);
    T_PFeval->SetBranchStatus("truth_NprimPio",1);
  }
  if (pfeval.flag_recoprotonMomentum){
    T_PFeval->SetBranchStatus("reco_protonMomentum",1);
  }
  if (pfeval.flag_showerMomentum){
    T_PFeval->SetBranchStatus("reco_showerMomentum",1);
    T_PFeval->SetBranchStatus("reco_Nproton",1);
    T_PFeval->SetBranchStatus("truth_showerMomentum",1);
    T_PFeval->SetBranchStatus("truth_nuScatType",1);
  }
  if (T_PFeval->GetBranch("truth_startMomentum")){
    T_PFeval->SetBranchStatus("truth_Ntrack",1); 
    T_PFeval->SetBranchStatus("truth_pdg",1); 
    T_PFeval->SetBranchStatus("truth_mother",1); 
    T_PFeval->SetBranchStatus("truth_startMomentum",1); 
  }


  WeightInfo weight;
  TTree *T_weight = (TTree*)file->Get("wcpselection/T_weight");
  weight.file_type = new std::string();
  weight.expskin_FluxUnisim= new std::vector<float>;
  weight.horncurrent_FluxUnisim= new std::vector<float>;
  weight.kminus_PrimaryHadronNormalization= new std::vector<float>;
  weight.kplus_PrimaryHadronFeynmanScaling= new std::vector<float>;
  weight.kzero_PrimaryHadronSanfordWang= new std::vector<float>;
  weight.nucleoninexsec_FluxUnisim= new std::vector<float>;
  weight.nucleonqexsec_FluxUnisim= new std::vector<float>;
  weight.nucleontotxsec_FluxUnisim= new std::vector<float>;
  weight.piminus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  weight.pioninexsec_FluxUnisim= new std::vector<float>;
  weight.pionqexsec_FluxUnisim= new std::vector<float>;
  weight.piontotxsec_FluxUnisim= new std::vector<float>;
  weight.piplus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  
  weight.All_UBGenie= new std::vector<float>;
  weight.AxFFCCQEshape_UBGenie= new std::vector<float>;
  weight.DecayAngMEC_UBGenie= new std::vector<float>;
  weight.NormCCCOH_UBGenie= new std::vector<float>;
  weight.NormNCCOH_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_Reduced_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_UBGenie= new std::vector<float>;
  weight.RootinoFix_UBGenie= new std::vector<float>;
  weight.ThetaDelta2NRad_UBGenie= new std::vector<float>;
  weight.Theta_Delta2Npi_UBGenie= new std::vector<float>;
  weight.TunedCentralValue_UBGenie= new std::vector<float>;
  weight.VecFFCCQEshape_UBGenie= new std::vector<float>;
  weight.XSecShape_CCMEC_UBGenie= new std::vector<float>;
  weight.splines_general_Spline= new std::vector<float>;
  weight.xsr_scc_Fa3_SCC= new std::vector<float>;
  weight.xsr_scc_Fv3_SCC= new std::vector<float>;

  weight.reinteractions_piminus_Geant4 = new std::vector<float>;
  weight.reinteractions_piplus_Geant4 = new std::vector<float>;
  weight.reinteractions_proton_Geant4 = new std::vector<float>; 
  
  TString option;
  if (rw_type == 1){
    option = "reweight";
  }else if (rw_type == 2){
    option = "reweight_cor";
  }else if (T_weight->GetBranch("expskin_FluxUnisim")){
    option = "expskin_FluxUnisim";
  }else if (T_weight->GetBranch("horncurrent_FluxUnisim")){
    option = "horncurrent_FluxUnisim";
  }else if (T_weight->GetBranch("kminus_PrimaryHadronNormalization")){
    option = "kminus_PrimaryHadronNormalization";
  }else if (T_weight->GetBranch("kplus_PrimaryHadronFeynmanScaling")){
    option = "kplus_PrimaryHadronFeynmanScaling";
  }else if (T_weight->GetBranch("kzero_PrimaryHadronSanfordWang")){
    option = "kzero_PrimaryHadronSanfordWang";
  }else if (T_weight->GetBranch("nucleoninexsec_FluxUnisim")){
    option = "nucleoninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleonqexsec_FluxUnisim")){
    option = "nucleonqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleontotxsec_FluxUnisim")){
    option = "nucleontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piminus_PrimaryHadronSWCentralSplineVariation")){
    option = "piminus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("pioninexsec_FluxUnisim")){
    option = "pioninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("pionqexsec_FluxUnisim")){
    option = "pionqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("piontotxsec_FluxUnisim")){
    option = "piontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piplus_PrimaryHadronSWCentralSplineVariation")){
    option = "piplus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("All_UBGenie")){
    option = "UBGenieFluxSmallUni";
  }else if (T_weight->GetBranch("reinteractions_piminus_Geant4")){
    option = "reinteractions_piminus_Geant4";
  }else if (T_weight->GetBranch("reinteractions_piplus_Geant4")){
    option = "reinteractions_piplus_Geant4";
  }else if (T_weight->GetBranch("reinteractions_proton_Geant4")){
    option = "reinteractions_proton_Geant4";
  }
  
  set_tree_address(T_weight, weight, option);
  //std::cout << T_eval->GetEntries() << " " << T_weight->GetEntries() << " " << option << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > >& set_events = map_passed_events[input_filename];


  std::vector<int> max_lengths;
  std::vector<int> sup_lengths;

  int temp_sum = 0;
  
  for (size_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    T_weight->GetEntry(i);

    std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::tuple<int, float, bool, int> > > event_info;
    std::get<0>(event_info) = eval.weight_cv * eval.weight_spline;
    std::get<1>(event_info) = leeweight(eval.truth_nuEnergy);
     
    double reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info());//automatically 1 if reweighting is not applied
    std::get<0>(event_info) *= reweight;
 
     for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);

      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);

      //std::cout << no << " " << histoname << std::endl;
      
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);

      auto it3 = disabled_ch_names.find(ch_name);
      if (it3 != disabled_ch_names.end()) continue;
      
      float val = get_kine_var(kine, eval, pfeval, tagger, false, var_name);
      int flag_passall = (get_cut_pass(ch_name, add_cut, false, eval, pfeval, tagger, kine));
      bool flag_pass = flag_passall > 0;
      int signal_bin = -1;
      if (xs_signal_ch_names.find(ch_name) != xs_signal_ch_names.end()){
	signal_bin = get_xs_signal_no(cut_file, map_cut_xs_bin, eval, pfeval, tagger, kine);
      }

      //  std::cout << flag_pass << " " << signal_bin << " " << no << std::endl;
      if (flag_pass || signal_bin !=-1) {
        if(flag_passall >= 0)
          std::get<4>(event_info).insert(std::make_tuple(no, val, flag_pass, signal_bin));
	//	if (no == 0) std::cout << "Xin: " << " " << flag_pass << " " << signal_bin << " " << eval.weight_cv * eval.weight_spline << " " <<eval.run << " " << eval.subrun << " " << eval.event << std::endl;
      }
      //if (flag_pass || (signal_bin !=-1 && is_preselection(eval))) std::get<4>(event_info).insert(std::make_tuple(no, val, flag_pass, signal_bin));
     }
    
    if (std::get<4>(event_info).size()>0){
      if (option == "expskin_FluxUnisim"){
	std::get<2>(event_info).resize(weight.expskin_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.expskin_FluxUnisim->size());
	for (size_t j=0;j!=weight.expskin_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.expskin_FluxUnisim->at(j) - 1.0; // relative ...
	}
      }else if (option == "horncurrent_FluxUnisim"){
	std::get<2>(event_info).resize(weight.horncurrent_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.horncurrent_FluxUnisim->size());
	for (size_t j=0;j!= weight.horncurrent_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.horncurrent_FluxUnisim->at(j) - 1.0; // relative ...
	}
      }else if (option == "kminus_PrimaryHadronNormalization"){
	std::get<2>(event_info).resize(weight.kminus_PrimaryHadronNormalization->size());
	std::get<3>(event_info).push_back(weight.kminus_PrimaryHadronNormalization->size());
	for (size_t j=0;j!= weight.kminus_PrimaryHadronNormalization->size(); j++){
	  std::get<2>(event_info).at(j) = weight.kminus_PrimaryHadronNormalization->at(j) - 1.0; 
	}
      }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
	std::get<2>(event_info).resize(weight.kplus_PrimaryHadronFeynmanScaling->size());
	std::get<3>(event_info).push_back(weight.kplus_PrimaryHadronFeynmanScaling->size());
	for (size_t j=0;j!=weight.kplus_PrimaryHadronFeynmanScaling->size();j++){
	  std::get<2>(event_info).at(j) = weight.kplus_PrimaryHadronFeynmanScaling->at(j) - 1.0;
	}
      }else if (option == "kzero_PrimaryHadronSanfordWang"){
	std::get<2>(event_info).resize(weight.kzero_PrimaryHadronSanfordWang->size());
	std::get<3>(event_info).push_back(weight.kzero_PrimaryHadronSanfordWang->size());
	for (size_t j=0;j!=weight.kzero_PrimaryHadronSanfordWang->size();j++){
	  std::get<2>(event_info).at(j) = weight.kzero_PrimaryHadronSanfordWang->at(j) - 1.0;
	}
      }else if (option == "nucleoninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleoninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleoninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.nucleoninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.nucleoninexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "nucleonqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleonqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleonqexsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleonqexsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleonqexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "nucleontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleontotxsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleontotxsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleontotxsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piminus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piminus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
      }else if (option == "pioninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pioninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pioninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pioninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pioninexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "pionqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pionqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pionqexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pionqexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pionqexsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.piontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.piontotxsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.piontotxsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.piontotxsec_FluxUnisim->at(j) - 1.0;
	}
      }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piplus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piplus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
      }else if (option == "reinteractions_piminus_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_piminus_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_piminus_Geant4->size());
	for (size_t j=0; j!= weight.reinteractions_piminus_Geant4->size(); j++){
	  std::get<2>(event_info).at(j) = weight.reinteractions_piminus_Geant4->at(j) - 1.0;
	}
      }else if (option == "reinteractions_piplus_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_piplus_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_piplus_Geant4->size());
	for (size_t j=0; j!= weight.reinteractions_piplus_Geant4->size(); j++){
	  std::get<2>(event_info).at(j) = weight.reinteractions_piplus_Geant4->at(j) - 1.0;
	}
      }else if (option == "reinteractions_proton_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_proton_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_proton_Geant4->size());
	for (size_t j=0; j!= weight.reinteractions_proton_Geant4->size(); j++){
	  std::get<2>(event_info).at(j) = weight.reinteractions_proton_Geant4->at(j) - 1.0;
	}
      }else if (option == "reweight"){
        std::get<2>(event_info).resize(1000);
        std::get<3>(event_info).push_back(1000);
        if(!(flag_reweight)) reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
        for (size_t j=0;j!=1000;j++){
          if(flag_reweight){
            if (weight.weight_cv>0 && reweight!=1){
              gRandom->SetSeed(j*reweight*77777);
              double rand = gRandom->Gaus(reweight,abs(1-reweight));
              std::get<2>(event_info).at(j) = (rand-reweight)/reweight;
            }else std::get<2>(event_info).at(j) = 0;
          }else{
            gRandom->SetSeed(j*reweight*77777);
            double rand = gRandom->Gaus(1,abs(1-reweight));
            std::get<2>(event_info).at(j) = rand-1;
          }
        }
      }else if (option == "reweight_cor"){
        std::get<2>(event_info).resize(1);
        std::get<3>(event_info).push_back(1);
        if(flag_reweight){
          if (weight.weight_cv>0 && reweight!=1){
            std::get<2>(event_info).at(0) = (1-reweight)/reweight;
          }else{
            std::get<2>(event_info).at(0) = 0;
          }
        }else{
           reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
           std::get<2>(event_info).at(0) = reweight-1;
        }
      }else if (option == "UBGenieFluxSmallUni"){
	int acc_no = 0;
	std::get<2>(event_info).resize(weight.All_UBGenie->size());
	std::get<3>(event_info).push_back(weight.All_UBGenie->size());
	for (size_t j=0;j!=weight.All_UBGenie->size();j++){
	  if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.All_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.All_UBGenie->size();

        if(rw_type==3){
          std::get<2>(event_info).resize(acc_no+1000);
          std::get<3>(event_info).push_back(1000);
          if(!(flag_reweight)) reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
          for (size_t j=0;j!=1000;j++){
            if(flag_reweight){
              if (weight.weight_cv>0 && reweight!=1){
                gRandom->SetSeed(j*reweight*77777);
                double rand = gRandom->Gaus(reweight,abs(1-reweight));
                std::get<2>(event_info).at(acc_no+j) = (rand-reweight)/reweight;
              }else std::get<2>(event_info).at(acc_no+j) = 0;
            }else{
              gRandom->SetSeed(j*reweight*77777);
              double rand = gRandom->Gaus(1,abs(1-reweight));
              std::get<2>(event_info).at(acc_no+j) = rand-1;
            }
          }
          acc_no+=1000;
          std::get<2>(event_info).resize(acc_no+1);
          std::get<3>(event_info).push_back(1);
          if(flag_reweight){
            if (weight.weight_cv>0 && reweight!=1){
              std::get<2>(event_info).at(acc_no) = (1-reweight)/reweight;
            }else{
              std::get<2>(event_info).at(acc_no) = 0;
            }
          }else{
             reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
             std::get<2>(event_info).at(acc_no) = reweight-1;
          }
          acc_no++; 
        }

	std::get<2>(event_info).resize(acc_no + weight.AxFFCCQEshape_UBGenie->size());
	std::get<3>(event_info).push_back(weight.AxFFCCQEshape_UBGenie->size());
	for (size_t j=0; j!= weight.AxFFCCQEshape_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.AxFFCCQEshape_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.AxFFCCQEshape_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.DecayAngMEC_UBGenie->size());
	std::get<3>(event_info).push_back(weight.DecayAngMEC_UBGenie->size());
	for (size_t j=0; j!= weight.DecayAngMEC_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.DecayAngMEC_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.DecayAngMEC_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.NormCCCOH_UBGenie->size());
	std::get<3>(event_info).push_back(weight.NormCCCOH_UBGenie->size());
	for (size_t j=0; j!= weight.NormCCCOH_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.NormCCCOH_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.NormCCCOH_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.NormNCCOH_UBGenie->size());
	std::get<3>(event_info).push_back(weight.NormNCCOH_UBGenie->size());
	for (size_t j=0; j!= weight.NormNCCOH_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.NormNCCOH_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.NormNCCOH_UBGenie->size();

	if (!weight.flag_sep_28){
	  std::get<2>(event_info).resize(acc_no + weight.RPA_CCQE_Reduced_UBGenie->size());
	  std::get<3>(event_info).push_back(weight.RPA_CCQE_Reduced_UBGenie->size());
	  for (size_t j=0; j!= weight.RPA_CCQE_Reduced_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	      std::get<2>(event_info).at(acc_no+j) = (weight.RPA_CCQE_Reduced_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	    }else{
	      std::get<2>(event_info).at(acc_no+j) = 0;
	    }
	  }
	  acc_no += weight.RPA_CCQE_Reduced_UBGenie->size();
	}
	
	std::get<2>(event_info).resize(acc_no + weight.RPA_CCQE_UBGenie->size());
	std::get<3>(event_info).push_back(weight.RPA_CCQE_UBGenie->size());
	for (size_t j=0; j!= weight.RPA_CCQE_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.RPA_CCQE_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.RPA_CCQE_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.ThetaDelta2NRad_UBGenie->size());
	std::get<3>(event_info).push_back(weight.ThetaDelta2NRad_UBGenie->size());
	for (size_t j=0; j!= weight.ThetaDelta2NRad_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.ThetaDelta2NRad_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.ThetaDelta2NRad_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.Theta_Delta2Npi_UBGenie->size());
	std::get<3>(event_info).push_back(weight.Theta_Delta2Npi_UBGenie->size());
	for (size_t j=0; j!= weight.Theta_Delta2Npi_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.Theta_Delta2Npi_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.Theta_Delta2Npi_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.VecFFCCQEshape_UBGenie->size());
	std::get<3>(event_info).push_back(weight.VecFFCCQEshape_UBGenie->size());
	for (size_t j=0; j!= weight.VecFFCCQEshape_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.VecFFCCQEshape_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.VecFFCCQEshape_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.XSecShape_CCMEC_UBGenie->size());
	std::get<3>(event_info).push_back(weight.XSecShape_CCMEC_UBGenie->size());
	for (size_t j=0; j!= weight.XSecShape_CCMEC_UBGenie->size(); j++){
	    if (weight.weight_cv>0){
	    std::get<2>(event_info).at(acc_no+j) = (weight.XSecShape_CCMEC_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.XSecShape_CCMEC_UBGenie->size();

	
	std::get<2>(event_info).resize(acc_no + weight.xsr_scc_Fa3_SCC->size());
	std::get<3>(event_info).push_back(weight.xsr_scc_Fa3_SCC->size());
	for (size_t j=0; j!= weight.xsr_scc_Fa3_SCC->size(); j++){
	    if (weight.weight_cv>0){
	      std::get<2>(event_info).at(acc_no+j) = (weight.xsr_scc_Fa3_SCC->at(j) -1.0);//- weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.xsr_scc_Fa3_SCC->size();

	std::get<2>(event_info).resize(acc_no + weight.xsr_scc_Fv3_SCC->size());
	std::get<3>(event_info).push_back(weight.xsr_scc_Fv3_SCC->size());
	for (size_t j=0; j!= weight.xsr_scc_Fv3_SCC->size(); j++){
	    if (weight.weight_cv>0){
	      std::get<2>(event_info).at(acc_no+j) = (weight.xsr_scc_Fv3_SCC->at(j) - 1.0);//weight.weight_cv)/weight.weight_cv;
	  }else{
	    std::get<2>(event_info).at(acc_no+j) = 0;
	  }
	}
	acc_no += weight.xsr_scc_Fv3_SCC->size();
	
      }

      if (max_lengths.size() < std::get<3>(event_info).size()) max_lengths.resize(std::get<3>(event_info).size());
      for(size_t j = 0; j!= std::get<3>(event_info).size(); j++){
	if (max_lengths.at(j) < std::get<3>(event_info).at(j)) max_lengths.at(j) = std::get<3>(event_info).at(j);
      // if (max_length < std::get<3>(event_info).size()) max_length = std::get<3>(event_info).size();
      }
      
      set_events.insert(event_info);
      //  temp_sum ++;
    }
    
  }

  // std::cout << temp_sum << " " << set_events.size() << std::endl;



  if (option == "expskin_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "horncurrent_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "kminus_PrimaryHadronNormalization"){
    sup_lengths.push_back(1000);
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
    sup_lengths.push_back(1000);
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
    sup_lengths.push_back(1000);
  }else if (option == "nucleoninexsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "nucleonqexsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "nucleontotxsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
    sup_lengths.push_back(1000);
  }else if (option == "pioninexsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "pionqexsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "piontotxsec_FluxUnisim"){
    sup_lengths.push_back(1000);
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
    sup_lengths.push_back(1000);
  }else if (option == "reinteractions_piminus_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reinteractions_piplus_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reinteractions_proton_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reweight"){
    sup_lengths.push_back(1000);
  }else if (option == "reweight_cor"){
    sup_lengths.push_back(1);
  }else if (option == "UBGenieFluxSmallUni"){
    sup_lengths.push_back(600); // all_ubgenie
    if(rw_type==3){
      sup_lengths.push_back(1000); //reweighting syst
      sup_lengths.push_back(1); //cor reweighting syst
    } 
    sup_lengths.push_back(1);   // AxFFCCQEshape_UBGenie-
    sup_lengths.push_back(1);   // DecayAngMEC_UBGenie
    sup_lengths.push_back(1); // NormCCCOH_UBGenie
    sup_lengths.push_back(1); // NormNCCOH_UBGenie
    if (!weight.flag_sep_28){
      sup_lengths.push_back(1); // RPA_CCQE_Reduced_UBGenie
    }
    sup_lengths.push_back(2); // RPA_CCQE_UBGenie
    sup_lengths.push_back(1); // ThetaDelta2NRad_UBGenie
    sup_lengths.push_back(1); // Theta_Delta2Npi_UBGenie
    sup_lengths.push_back(1); // VecFFCCQEshape_UBGenie
    sup_lengths.push_back(1); // XSecShape_CCMEC_UBGenie
    sup_lengths.push_back(10); // xsr_scc_Fa3_SCC
    sup_lengths.push_back(10); // xsr_scc_Fv3_SCC

  }



  
  delete file;
  return std::make_pair(max_lengths, sup_lengths);
}


