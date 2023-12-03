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

  bool flag_data = true;
  
  TFile *file = new TFile(input_filename,"READ");
  
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
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);
  
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

  std::cout << "Total entries: " << T_eval->GetEntries() << std::endl;

   T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);

  if (tagger.flag_nc_gamma_bdt){
    T_BDTvars->SetBranchStatus("nc_delta_score", 1);
    T_BDTvars->SetBranchStatus("nc_pio_score", 1);
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
      T_PFeval->SetBranchStatus("mcflux_ntype",1);
      T_PFeval->SetBranchStatus("mcflux_ptype",1);
      T_PFeval->SetBranchStatus("mcflux_tptype",1);
      T_PFeval->SetBranchStatus("mcflux_dk2gen",1);
      T_PFeval->SetBranchStatus("mcflux_gen2vtx",1);
      T_PFeval->SetBranchStatus("mcflux_ndecay",1);
      
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
    }
  }
  
  
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    
    if (!is_preselection(eval)) continue;

    double Enu = get_reco_Enu_corr(kine, flag_data);
      
    // if (tagger.nue_score>7.0 && eval.match_isFC==1 && Enu >= 800 && Enu <=1500){
    //   std::cout << eval.run << " \t" << eval.subrun << " \t" << eval.event << " \t" <<  Enu << "  \t" << pfeval.reco_nuvtxX << "  \t" << pfeval.reco_nuvtxY << "  \t" << pfeval.reco_nuvtxZ << std::endl;
    // }

    if (tagger.nue_score > 4.0){
      std::cout << eval.run << " \t" << eval.subrun << " \t" << eval.event << " \t" <<  Enu << "  \t" << pfeval.reco_nuvtxX << "  \t" << pfeval.reco_nuvtxY << "  \t" << pfeval.reco_nuvtxZ << " " << tagger.nue_score << " " << eval.match_isFC << std::endl;
    }
    
  }

  return 0;
}
