void LEEana::CovMatrix::gen_pred_stat_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean, TMatrixD* cov_mat){
  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ; 
  std::map<int, TString> map_no_histoname;
  
  
  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);

    if (filetype == 5 || filetype == 15)   continue; // if this is data ...

    
    if (period != run && run!=0) continue;
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

     for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
      int ch = map_name_ch[std::get<5>(*it1)];
      int obsch = get_obsch_name(std::get<5>(*it1));
      int covch = get_covch_name(std::get<5>(*it1));
      int flag_lee = std::get<7>(map_ch_hist[ch]);
      TString histoname = std::get<0>(*it1);
      TH1F *htemp = map_histoname_hist[histoname];
      //
      map_histoname_infos[histoname] = std::make_tuple(ncount, covch, flag_lee, input_filename);
      map_no_histoname[ncount] = histoname;
      ncount ++;

      //std::cout << histoname << obsch << " " << covch << " " << flag_lee << std::endl;
    }
    
  }


  // results ... filename --> re --> variable, weight, lee weight, 
  std::map<TString, std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > > map_all_events;
  std::map<TString, double> map_filename_pot;
  std::map<int, double> map_data_period_pot;

  std::map<TString, std::set<std::pair<TString, double> > > map_filename_set_histoname_ratio;
  
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    
    if (filetype == 5 || filetype == 15) {
      map_data_period_pot[period] = get_ext_pot(input_filename);
      //      std::cout << period << " " << map_data_period_pot[period] << std::endl;
      continue; // skip data ...
    }
    if (period != run && run !=0) continue;

    get_pred_events_info(input_filename, map_all_events, map_filename_pot, map_histoname_infos);      
  }

 
  fill_pred_stat_histograms(map_all_events, map_histoname_infos, map_no_histoname, map_histoname_hist);
  
  // merge histograms according to POTs ...
  for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int covch = it->first;
    TH1F *hpred = map_covch_hist[covch];
    hpred->Reset();
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      std::map<int, double> temp_map_mc_acc_pot;
       
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
    	TString histoname = (*it2).first;
    	TString input_filename = map_histogram_inputfile[histoname];
    	auto it3 = map_inputfile_info.find(input_filename);
    	int period = std::get<1>(it3->second);
    	if (period != run && run != 0) continue;
       
    	double mc_pot = map_filename_pot[input_filename];
    	//std::cout << mc_pot << std::endl;
    	if (temp_map_mc_acc_pot.find(period) == temp_map_mc_acc_pot.end()){
    	  temp_map_mc_acc_pot[period] = mc_pot;
    	}else{
    	  temp_map_mc_acc_pot[period] += mc_pot;
    	}
      }
       
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
    	TString histoname = (*it2).first;
    	TString input_filename = map_histogram_inputfile[histoname];
    	auto it3 = map_inputfile_info.find(input_filename);
    	int period = std::get<1>(it3->second);
    	if (period != run && run !=0) continue; // skip ...
    	double data_pot = map_data_period_pot[period];
    	double ratio = data_pot/temp_map_mc_acc_pot[period];

	map_filename_set_histoname_ratio[input_filename].insert(std::make_pair(histoname, ratio));
	
    	TH1F *hmc = map_histoname_hist[histoname];
    	htemp->Add(hmc, ratio);
      }
      
      hpred->Add(htemp);
      delete htemp;
    }
     
    int start_bin = map_covch_startbin[covch];
    for (int i=0;i!=hpred->GetNbinsX()+1;i++){
      (*vec_mean)[start_bin+i] = hpred->GetBinContent(i+1);

      // std::cout << start_bin+i << " " << (*vec_mean_diff)(start_bin+i) << " " <<  hpred->GetBinContent(i+1) << std::endl;;
      // std::cout << x[start_bin+i] << std::endl;
    }
  }


  const int rows = cov_mat->GetNcols();
  float x[rows];
  (*cov_mat).Zero();

  std::map<TString, TH1D*> map_filename_histo;
  for (auto it = map_all_events.begin(); it != map_all_events.end(); it++){
    TString filename = it->first;
    int nsize = it->second.size();
    TH1D* htemp = new TH1D(filename, filename, nsize, 0.5, nsize+0.5);
    for (size_t i=0;i!=nsize;i++){
      htemp->SetBinContent(i+1, std::get<2>(it->second.at(i)) );
    }
    map_filename_histo[filename] = htemp;
  } 

  int nround = 5000;
  for (int qx = 0; qx != nround; qx++){
    if (qx % 500 ==0) std::cout << qx << std::endl;
    
    for (int i=0;i!=rows;i++){
      x[i] = 0;
    }
        
    // fill the histogram with variation ...
    fill_pred_stat_histograms(map_filename_histo, map_all_events, map_histoname_infos, map_no_histoname, map_histoname_hist, map_filename_set_histoname_ratio);
    
     // merge histograms according to POTs ...
    for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
      //std::cout << it->first << std::endl;
      int covch = it->first;
      TH1F *hpred = map_covch_hist[covch];
      hpred->Reset();
      
      for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
	TH1F *htemp = (TH1F*)hpred->Clone("htemp");
	htemp->Reset();
	
	
	for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	  TString histoname = (*it2).first;
	  TString input_filename = map_histogram_inputfile[histoname];
	  auto it3 = map_inputfile_info.find(input_filename);
	  int period = std::get<1>(it3->second);
	  if (period != run && run !=0) continue; // skip ...
	  	  
	  TH1F *hmc = map_histoname_hist[histoname];
	  htemp->Add(hmc, 1);
	}
	
	hpred->Add(htemp);
	delete htemp;
      }
      
      int start_bin = map_covch_startbin[covch];
      for (int i=0;i!=hpred->GetNbinsX()+1;i++){
	x[start_bin+i]= hpred->GetBinContent(i+1);
	//	std::cout << (*vec_mean)(start_bin+i) << " " << x[start_bin+i] << std::endl;
      }
    }
    
    
    // add covariance matrix ...
    for (size_t n = 0;n!=rows; n++){
      for (size_t m =0; m!=rows;m++){
	(*cov_mat)(n,m) += (x[n]-(*vec_mean)[n]) * (x[m] - (*vec_mean)[m]) * 1./nround;
      }
    }  
  }
  

  // clean up the memory ...
  for (auto it = map_filename_histo.begin(); it != map_filename_histo.end(); it++){
    delete it->second;
  }
  


  // update the final  ones ...
   fill_pred_stat_histograms(map_all_events, map_histoname_infos, map_no_histoname, map_histoname_hist);
  
  // merge histograms according to POTs ...
  for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int covch = it->first;
    TH1F *hpred = map_covch_hist[covch];
    hpred->Reset();
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      std::map<int, double> temp_map_mc_acc_pot;
       
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
    	TString histoname = (*it2).first;
    	TString input_filename = map_histogram_inputfile[histoname];
    	auto it3 = map_inputfile_info.find(input_filename);
    	int period = std::get<1>(it3->second);
    	if (period != run && run != 0) continue;
       
    	double mc_pot = map_filename_pot[input_filename];
    	//std::cout << mc_pot << std::endl;
    	if (temp_map_mc_acc_pot.find(period) == temp_map_mc_acc_pot.end()){
    	  temp_map_mc_acc_pot[period] = mc_pot;
    	}else{
    	  temp_map_mc_acc_pot[period] += mc_pot;
    	}
      }
       
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
    	TString histoname = (*it2).first;
    	TString input_filename = map_histogram_inputfile[histoname];
    	auto it3 = map_inputfile_info.find(input_filename);
    	int period = std::get<1>(it3->second);
    	if (period != run && run !=0) continue; // skip ...
    	double data_pot = map_data_period_pot[period];
    	double ratio = data_pot/temp_map_mc_acc_pot[period];
	 
    	TH1F *hmc = map_histoname_hist[histoname];
    	htemp->Add(hmc, ratio);
      }
      
      hpred->Add(htemp);
      delete htemp;
    }
     
    //    int start_bin = map_covch_startbin[covch];
    // for (int i=0;i!=hpred->GetNbinsX()+1;i++){
    //  (*vec_mean)[start_bin+i] = hpred->GetBinContent(i+1);

      // std::cout << start_bin+i << " " << (*vec_mean_diff)(start_bin+i) << " " <<  hpred->GetBinContent(i+1) << std::endl;;
      // std::cout << x[start_bin+i] << std::endl;
    //}
  }
  
}

void LEEana::CovMatrix::fill_pred_stat_histograms(std::map<TString, TH1D*> map_filename_histo, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist, std::map<TString, std::set<std::pair<TString, double> > >& map_filename_histoname_ratio){

  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
    it->second->Reset();
  }
  
  // loop over files
  for (auto it = map_all_events.begin(); it!=map_all_events.end(); it++){
    TString filename = it->first;
    TH1D *hweight = map_filename_histo[filename];
    
    auto it2 = map_filename_histoname_ratio.find(filename);
    if (it2 == map_filename_histoname_ratio.end()) continue;

    double max_sum = 0;
    std::map<TString, double> map_histoname_sum;

    std::map<double, double> map_ratio_sum;
    
    for (auto it3 = it2->second.begin(); it3 != it2->second.end(); it3++){
      double ratio = it3->second;
      TString temp_histoname = it3->first;
      
      auto it4 = map_ratio_sum.find(ratio);
      if (it4 == map_ratio_sum.end()){
	double sum = gRandom->Poisson(hweight->GetSum() * ratio);
	if (sum > max_sum) max_sum = sum;
	map_histoname_sum[temp_histoname] = sum;
	map_ratio_sum[ratio] = sum;
      }else{
	double sum = it4->second;
	map_histoname_sum[temp_histoname] = sum;
      }
    }
    
    
    for (size_t i=0;i<max_sum;i++){
      int global_index = hweight->FindBin(hweight->GetRandom())-1;
      double weight_lee = std::get<3>(it->second.at(global_index));
      
       for (auto it1 = std::get<4>(it->second.at(global_index)).begin(); it1 != std::get<4>(it->second.at(global_index)).end(); it1++){
	 int no = std::get<0>(*it1);
	 double val_cv = std::get<1>(*it1);
	 bool flag_cv = std::get<2>(*it1);

	 TString histoname = map_no_histoname[no];
	 TH1F *htemp = map_histoname_hist[histoname];
	 int flag_lee = std::get<2>(map_histoname_infos[histoname]);
	 
	 auto it4 = map_histoname_sum.find(histoname);
	 if (it4 == map_histoname_sum.end()) std::cout << "something wrong! " << std::endl;
	 double sum = it4->second;
	 
	 if (i>=sum) continue; // go over ...
	 
	  if (flag_cv){
	    if (flag_lee){
	      htemp->Fill(val_cv,weight_lee);
	    }else{
	      htemp->Fill(val_cv,1);
	    }
	  }
       }
      
    }
  }
  
}


void LEEana::CovMatrix::fill_pred_stat_histograms(std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist){
  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
    it->second->Reset();
  }
  
  // fill central value ...
  
  // loop over files
  for (auto it = map_all_events.begin(); it!=map_all_events.end(); it++){
    // loop over events ...
    //std::cout << it->first << " " << it->second.size() << std::endl;
   
	 
    for (size_t i=0;i!=it->second.size(); i++){
       double weight = std::get<2>(it->second.at(i));
       double weight_lee = std::get<3>(it->second.at(i));
      for (auto it1 = std::get<4>(it->second.at(i)).begin(); it1 != std::get<4>(it->second.at(i)).end(); it1++){
	int no = std::get<0>(*it1);
	double val_cv = std::get<1>(*it1);
	bool flag_cv = std::get<2>(*it1);

	TString histoname = map_no_histoname[no];
	TH1F *htemp = map_histoname_hist[histoname];
	 int flag_lee = std::get<2>(map_histoname_infos[histoname]);
	if (flag_cv){
	  if (flag_lee){
	    htemp->Fill(val_cv, weight * weight_lee);
	  }else{
	    htemp->Fill(val_cv, weight);
	  }
	}
      }
    }
  }
}

void LEEana::CovMatrix::get_pred_events_info(TString input_filename, std::map<TString, std::vector< std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > >&map_all_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){
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

  bool flag_data = true;
  if (T_eval->GetBranch("weight_cv")) flag_data = false;
  
  set_tree_address(T_BDTvars, tagger, 2);
  if (flag_data){
    set_tree_address(T_eval, eval,2); 
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval); 
    set_tree_address(T_PFeval, pfeval);
  }

  //  std::cout << flag_data << " " << input_filename << std::endl;
  
  set_tree_address(T_pot, pot); 
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;

  map_filename_pot[input_filename] = total_pot;
  //std::cout << input_filename << " " << total_pot << " " << std::endl;


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
  
  if (!flag_data){
    T_eval->SetBranchStatus("weight_spline",1);
    T_eval->SetBranchStatus("weight_cv",1);
    T_eval->SetBranchStatus("weight_lee",1);
    T_eval->SetBranchStatus("weight_change",1);
    // MC enable truth information ...
    T_eval->SetBranchStatus("truth_isCC",1);
    T_eval->SetBranchStatus("truth_nuPdg",1);
    T_eval->SetBranchStatus("truth_vtxInside",1);
    T_eval->SetBranchStatus("truth_nuEnergy",1);
    T_eval->SetBranchStatus("truth_vtxX",1);
    T_eval->SetBranchStatus("truth_vtxY",1);
    T_eval->SetBranchStatus("truth_vtxZ",1);
    // Xs related
    T_eval->SetBranchStatus("match_completeness_energy",1);
    T_eval->SetBranchStatus("truth_energyInside",1);
  }
  
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
  if (!flag_data){
    T_PFeval->SetBranchStatus("nuvtx_diff",1);
    T_PFeval->SetBranchStatus("showervtx_diff",1);
    T_PFeval->SetBranchStatus("muonvtx_diff",1);
    T_PFeval->SetBranchStatus("truth_muonMomentum",1);
  }
  if (pfeval.flag_NCDelta){
    
    if (!flag_data){
      T_PFeval->SetBranchStatus("truth_NCDelta",1);
      T_PFeval->SetBranchStatus("truth_NprimPio",1);
    }
  }
  if (pfeval.flag_recoprotonMomentum){
    T_PFeval->SetBranchStatus("reco_protonMomentum",1);
  }
  
  if (pfeval.flag_showerMomentum){
    T_PFeval->SetBranchStatus("reco_showerMomentum",1);
    T_PFeval->SetBranchStatus("reco_Nproton",1);
    if (!flag_data){
      T_PFeval->SetBranchStatus("truth_showerMomentum",1);
      T_PFeval->SetBranchStatus("truth_nuScatType",1);
      // oscillation formula ...
      T_PFeval->SetBranchStatus("truth_nu_momentum",1);
      T_PFeval->SetBranchStatus("neutrino_type",1);
      T_PFeval->SetBranchStatus("mcflux_ntype",1);
      T_PFeval->SetBranchStatus("mcflux_dk2gen",1);
      T_PFeval->SetBranchStatus("mcflux_gen2vtx",1);
      T_PFeval->SetBranchStatus("mcflux_ndecay",1);
    }
  }
  if (T_PFeval->GetBranch("truth_pdg")){
     T_PFeval->SetBranchStatus("truth_Ntrack",1); 
     T_PFeval->SetBranchStatus("truth_pdg",1); 
     T_PFeval->SetBranchStatus("truth_mother",1); 
     T_PFeval->SetBranchStatus("truth_startMomentum",1); 
  }

  std::vector<std::tuple<int, int, double, double, std::set<std::tuple<int, double, bool> > > > vec_events;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  vec_events.resize(T_eval->GetEntries());

  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);

    std::get<0>(vec_events.at(i)) = eval.run;
    std::get<1>(vec_events.at(i)) = eval.event;

    if (!flag_data){
      std::get<2>(vec_events.at(i)) = eval.weight_cv * eval.weight_spline;
      // hack for now ...
      std::get<3>(vec_events.at(i)) = leeweight(eval.truth_nuEnergy);
    }else{
      std::get<2>(vec_events.at(i)) = 1;
      std::get<3>(vec_events.at(i)) = 0;
    }

    double osc_weight = 1.0;
    bool flag_updated = false;
    
    for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      
      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);
      
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);

      double val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name);
      bool flag_pass = (get_cut_pass(ch_name, add_cut, flag_data, eval, pfeval, tagger, kine) > 0);

      if (flag_pass) std::get<4>(vec_events.at(i)).insert(std::make_tuple(no, val, flag_pass));
      
      if (flag_osc && is_osc_channel(ch_name) && (!flag_updated)){
	osc_weight = get_osc_weight(eval, pfeval);
	flag_updated = true;
      }
    }
    if (!flag_data){
      std::get<2>(vec_events.at(i)) *= osc_weight;
      std::get<2>(vec_events.at(i)) *= get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info());
    }  
}

  map_all_events[input_filename] = vec_events;

  delete file;
  
}
