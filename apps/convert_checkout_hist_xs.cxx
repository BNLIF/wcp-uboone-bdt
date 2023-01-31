#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

#include "WCPLEEANA/cuts.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{

  TString input_filename = argv[1];
  TString out_filename =  argv[2];

  
  bool flag_data = true;

  TFile *file = new TFile(input_filename,"READ");
  
  CovMatrix cov;
  cov.add_xs_config();
  
  std::cout << "Xs mode: " << std::endl;

  cov.print_rw(cov.get_rw_info());
  
  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");
 
  if (T_eval->GetBranch("weight_cv")) flag_data = false;

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

#include "init.txt"

  set_tree_address(T_BDTvars, tagger,2 );
  if (flag_data){
    set_tree_address(T_eval, eval,2);
    set_tree_address(T_PFeval, pfeval,2);
  }else{
    set_tree_address(T_eval, eval);
    set_tree_address(T_PFeval, pfeval);
  }
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  double ext_pot = cov.get_ext_pot(input_filename);
  if (ext_pot != 0) total_pot = ext_pot;
  
  std::cout << "Total POT: " << total_pot << " external POT: " << ext_pot << std::endl;

  
  // prepare histograms ...
  // declare histograms ...
  TH1F *htemp = 0; TH1F *htemp1 = 0; TH2F *htemp2 = 0; int num = 0;
  std::map<TString, std::tuple<TH1F*, TH1F*, TH2F*, int> > map_histoname_hists;
  
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
  
  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
  std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "CV:" << std::endl;

  for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    
    //    std::cout << std::get<0>( *it)  << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    num = 1;
   
    bool flag_spec = cov.is_xs_chname(ch_name) ;
    if (flag_spec){
      int nbins1 = cov.get_xs_nsignals();
      TString temp_histoname = histoname + "_signal";
      htemp1 = new TH1F(temp_histoname, temp_histoname,nbins1, 0.5,nbins1+0.5);
      temp_histoname = histoname + "_R";
      htemp2 = new TH2F(temp_histoname, temp_histoname,nbin, llimit, hlimit,nbins1, 0.5,nbins1+0.5);
      num = 3;
    }
    map_histoname_hists[histoname] = std::make_tuple(htemp, htemp1, htemp2, num);
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
  std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
  //  std::cout << "Error2: " << std::endl;
  for (auto it = histo_infos_err2.begin(); it != histo_infos_err2.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
        
    map_histoname_hists[histoname] = std::make_tuple(htemp, (TH1F*)(0), (TH2F*)(0), 1);
  }
  //  std::cout << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
  std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));
  //std::cout << "Cross: " << std::endl;
  for (auto it = histo_infos_cros.begin(); it != histo_infos_cros.end(); it++){
    TString histoname = std::get<0>(*it);
    Int_t nbin = std::get<1>(*it);
    float llimit = std::get<2>(*it);
    float hlimit = std::get<3>(*it);
    TString var_name = std::get<4>(*it);
    TString ch_name = std::get<5>(*it);
    TString add_cut = std::get<6>(*it);
    TString weight = std::get<7>(*it);
    //std::cout << std::get<0>( *it) << " " << std::get<1>(*it) << " " << std::get<4>(*it) << " " << std::get<5>(*it) << " " << std::get<6>(*it) << " " << std::get<7>(*it) << std::endl;
    htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
    map_histoname_hists[histoname] = std::make_tuple(htemp, (TH1F*)(0), (TH2F*)(0),1);
  }
  //  std::cout << std::endl;
  
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

  // T_eval->SetBranchStatus("run",1);
  // T_eval->SetBranchStatus("subrun",1);
  // T_eval->SetBranchStatus("event",1);
  
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
    T_eval->SetBranchStatus("truth_energyInside",1);
    T_eval->SetBranchStatus("truth_vtxX",1);
    T_eval->SetBranchStatus("truth_vtxY",1);
    T_eval->SetBranchStatus("truth_vtxZ",1);
    T_eval->SetBranchStatus("match_completeness_energy",1);
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
      T_PFeval->SetBranchStatus("truth_nuIntType",1);
      T_PFeval->SetBranchStatus("truth_muonMomentum",1);
      T_PFeval->SetBranchStatus("truth_pio_energy_1",1);
      T_PFeval->SetBranchStatus("truth_pio_energy_2",1);
      T_PFeval->SetBranchStatus("truth_pio_angle",1);
      T_PFeval->SetBranchStatus("truth_NprimPio",1);
      if(T_PFeval->GetBranch("truth_mother")){//prevents throwing an error for the non _PF files
        T_PFeval->SetBranchStatus("truth_Ntrack",1);
        T_PFeval->SetBranchStatus("truth_pdg",1); 
        T_PFeval->SetBranchStatus("truth_mother",1); 
        T_PFeval->SetBranchStatus("truth_startMomentum",1);
      }
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
    }
  }

  std::cout << "Total entries: " << T_eval->GetEntries() << std::endl;


  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    //    if (!is_preselection(eval)) continue;

    for (auto it = all_histo_infos.begin(); it != all_histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      Int_t nbin = std::get<1>(*it);
      float llimit = std::get<2>(*it);
      float hlimit = std::get<3>(*it);
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);
      TString weight = std::get<7>(*it);

      // get kinematics variable ...
      double val = get_kine_var(kine, eval, pfeval, tagger, flag_data, var_name);
      // get pass or not
      int flag_passall = get_cut_pass(ch_name, add_cut, flag_data, eval, pfeval, tagger, kine);
      bool flag_pass = flag_passall > 0;
      int signal_bin = -1;
      if (cov.is_xs_chname(ch_name))
	signal_bin = get_xs_signal_no(cov.get_cut_file(), cov.get_map_cut_xs_bin(), eval, pfeval, tagger, kine);

      if (!((signal_bin != -1) || flag_pass)) continue;
      // get weight ...
      double weight_val = get_weight(weight, eval, pfeval, kine, tagger, cov.get_rw_info(), flag_data);

      //  if (ch_name == "numuCC_signal_Enu_FC_overlay" && weight == "cv_spline") std::cout << "Xin: " << " " << flag_pass << " " << signal_bin << " " << weight_val << " " <<eval.run << " " << eval.subrun << " " << eval.event << std::endl;
      
      
      //     htemp = map_histoname_hist[histoname];
      auto tmp_hists = map_histoname_hists[histoname];
      TH1F *h1 = std::get<0>(tmp_hists);
      TH1F *h2 = std::get<1>(tmp_hists);
      TH2F *h3 = std::get<2>(tmp_hists);
      int num = std::get<3>(tmp_hists);
      if (num==1){
       	if (flag_pass) h1->Fill(val,  weight_val);
      }else{
       	if (signal_bin != -1){
       	  if (flag_pass) h1->Fill(val, weight_val);
       	  if (flag_passall >= 0) h2->Fill(signal_bin, weight_val);
       	  if (flag_pass) h3->Fill(val, signal_bin, weight_val);
       	}else{
          std::cout << num << " " << signal_bin << " " << ch_name << " " << flag_pass << " " << flag_passall << std::endl;
       	  std::cout << "[convt-hist-xs] Something wrong: cut/channel mismatch !" << std::endl;
       	}
      } // else
 


    }
  }
  

  // save histograms ...
  TFile *file1 = new TFile(out_filename,"RECREATE");
  file1->cd();
  TTree *T = new TTree("T","T");
  T->SetDirectory(file1);
  T->Branch("pot",&total_pot,"pot/D");
  T->Fill();
  
  for (auto it = map_histoname_hists.begin(); it!= map_histoname_hists.end(); it++){
  //   //std::cout<<"DEBUG: "<<it->first<<" "<<it->second->GetName()<<" "<<it->second->GetSum()<<"\n";
    TH1F *h1 = std::get<0>(it->second);
    TH1F *h2 = std::get<1>(it->second);
    TH2F *h3 = std::get<2>(it->second);
    int num = std::get<3>(it->second);
    if (num==1){
      h1->SetDirectory(file1);
    }else{
      h1->SetDirectory(file1);
      h2->SetDirectory(file1);
      h3->SetDirectory(file1);
    }
    //   it->second->SetDirectory(file1);
  }
  
  file1->Write();
  file1->Close();

  
     
  
  
  return 0;
}
