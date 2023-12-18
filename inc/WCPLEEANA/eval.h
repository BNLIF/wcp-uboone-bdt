#ifndef UBOONE_LEE_EVAL
#define UBOONE_LEE_EVAL

namespace LEEana{
  struct EvalInfo{
    bool is_match_found_int;
    bool is_file_type;
    
    Int_t run;
    Int_t subrun;
    Int_t event;
    Bool_t flash_found;
    Int_t flash_found_asInt;
    Float_t flash_time;
    
    std::string *file_type;
    
    Float_t flash_measPe;
    Float_t flash_predPe;
    Bool_t match_found;
    Int_t match_found_asInt;
    UInt_t match_type;
    Bool_t match_isFC;
    
    Bool_t match_isTgm;
    Bool_t match_notFC_FV;
    Bool_t match_notFC_SP;
    Bool_t match_notFC_DC;
    Float_t match_chargeU;
    
    Float_t match_chargeV;
    Float_t match_chargeY;
    Float_t match_energyY;
    Bool_t light_mismatch;
    Float_t match_charge;
    
    Float_t match_energy;
    Int_t stm_eventtype;
    Int_t stm_lowenergy;
    Int_t stm_LM;
    Int_t stm_TGM;
    
    Int_t stm_STM;
    Int_t stm_FullDead;
    Float_t stm_clusterlength;
    Float_t truth_nuEnergy;
    Float_t truth_energyInside;
    
    Float_t truth_electronInside;
    Int_t truth_nuPdg;
    Bool_t truth_isCC;
    Bool_t truth_isEligible;
    Bool_t truth_NCisEligible;
    
    Bool_t truth_isFC;
    Bool_t truth_vtxInside;
    Float_t truth_vtxX;
    Float_t truth_vtxY;
    Float_t truth_vtxZ;
    
    Float_t truth_nuTime;
    Float_t match_completeness;
    Float_t match_completeness_energy;
    Float_t match_purity;
    Float_t match_purity_xz;
    
    Float_t match_purity_xy;
    Float_t weight_spline;
    Float_t weight_cv;
    Float_t weight_flugg;
    Float_t weight_ppfx;
    Float_t weight_lee;
    
    Bool_t weight_change;


    // PeLEE
    Bool_t flag_pl;
    Float_t pl_shr_energy_tot_cali;
    Float_t pl_shr_tkfit_dedx_Y;
    Float_t pl_reco_nu_vtx_sce_x;
    Float_t pl_reco_nu_vtx_sce_y;
    Float_t pl_reco_nu_vtx_sce_z;
    Double_t pl_reco_e;
    Int_t pl_flag;

    // gLEE
    Bool_t flag_gl;
    Int_t gl_sel_type, gl_file_type;
    Float_t gl_true_Enu, gl_true_Elep, gl_reco_Eshower, gl_simple_pot_weight, gl_rem_orig_wc_pot_weight, gl_new_pot_weight, gl_overlap_weight,gl_wc_total_overlapped_weight;
};

 void set_tree_address(TTree *tree0, EvalInfo& eval_info, int flag = 1);
 void put_tree_address(TTree *tree0, EvalInfo& eval_info, int flag = 1);
}


 
void LEEana::set_tree_address(TTree *tree0, EvalInfo& eval_info, int flag){
  eval_info.is_match_found_int = false;
  eval_info.is_file_type = false;
  eval_info.flag_pl = false;
  eval_info.flag_gl = false;

  if (tree0->GetBranch("pl_flag")){
    eval_info.flag_pl = true;

    tree0->SetBranchAddress("pl_shr_energy_tot_cali",&eval_info.pl_shr_energy_tot_cali);
    tree0->SetBranchAddress("pl_shr_tkfit_dedx_Y",&eval_info.pl_shr_tkfit_dedx_Y);
    tree0->SetBranchAddress("pl_reco_nu_vtx_sce_x",&eval_info.pl_reco_nu_vtx_sce_x);
    tree0->SetBranchAddress("pl_reco_nu_vtx_sce_y",&eval_info.pl_reco_nu_vtx_sce_y);
    tree0->SetBranchAddress("pl_reco_nu_vtx_sce_z",&eval_info.pl_reco_nu_vtx_sce_z);
    tree0->SetBranchAddress("pl_reco_e",&eval_info.pl_reco_e);
    tree0->SetBranchAddress("pl_flag",&eval_info.pl_flag);
  }

  if (tree0->GetBranch("gl_sel_type")){
    eval_info.flag_gl = true;
    tree0->SetBranchAddress("gl_sel_type",&eval_info.gl_sel_type);
    tree0->SetBranchAddress("gl_file_type",&eval_info.gl_file_type);
    tree0->SetBranchAddress("gl_true_Enu",&eval_info.gl_true_Enu);
    tree0->SetBranchAddress("gl_true_Elep",&eval_info.gl_true_Elep);
    tree0->SetBranchAddress("gl_reco_Eshower",&eval_info.gl_reco_Eshower);
    tree0->SetBranchAddress("gl_simple_pot_weight",&eval_info.gl_simple_pot_weight);
    tree0->SetBranchAddress("gl_rem_orig_wc_pot_weight",&eval_info.gl_rem_orig_wc_pot_weight);
    tree0->SetBranchAddress("gl_new_pot_weight",&eval_info.gl_new_pot_weight);
    tree0->SetBranchAddress("gl_overlap_weight",&eval_info.gl_overlap_weight);
    tree0->SetBranchAddress("gl_wc_total_overlapped_weight",&eval_info.gl_wc_total_overlapped_weight);
  }
  
  
  tree0->SetBranchAddress("run", &eval_info.run);
  tree0->SetBranchAddress("subrun", &eval_info.subrun);
  tree0->SetBranchAddress("event", &eval_info.event);
  tree0->SetBranchAddress("flash_found", &eval_info.flash_found);
  if (tree0->GetBranch("flash_found_asInt"))   {
    tree0->SetBranchAddress("flash_found_asInt", &eval_info.flash_found_asInt);
    eval_info.is_match_found_int = true;
  }
  if (tree0->GetBranch("file_type")){
    tree0->SetBranchAddress("file_type",&eval_info.file_type);
    eval_info.is_file_type = true;
  }
  tree0->SetBranchAddress("flash_time", &eval_info.flash_time);

  tree0->SetBranchAddress("flash_measPe", &eval_info.flash_measPe);
  tree0->SetBranchAddress("flash_predPe", &eval_info.flash_predPe);
  tree0->SetBranchAddress("match_found", &eval_info.match_found);
  if (tree0->GetBranch("match_found_asInt"))   tree0->SetBranchAddress("match_found_asInt", &eval_info.match_found_asInt); 
  tree0->SetBranchAddress("match_type", &eval_info.match_type);
  tree0->SetBranchAddress("match_isFC", &eval_info.match_isFC);

  tree0->SetBranchAddress("match_isTgm", &eval_info.match_isTgm);
  tree0->SetBranchAddress("match_notFC_FV", &eval_info.match_notFC_FV);
  tree0->SetBranchAddress("match_notFC_SP", &eval_info.match_notFC_SP);
  tree0->SetBranchAddress("match_notFC_DC", &eval_info.match_notFC_DC);
  tree0->SetBranchAddress("match_chargeU", &eval_info.match_chargeU);
  
  tree0->SetBranchAddress("match_chargeV", &eval_info.match_chargeV);
  tree0->SetBranchAddress("match_chargeY", &eval_info.match_chargeY);
  tree0->SetBranchAddress("match_energyY", &eval_info.match_energyY);
  tree0->SetBranchAddress("light_mismatch", &eval_info.light_mismatch);
  tree0->SetBranchAddress("match_charge", &eval_info.match_charge);

  tree0->SetBranchAddress("match_energy", &eval_info.match_energy);
  tree0->SetBranchAddress("stm_eventtype", &eval_info.stm_eventtype);
  tree0->SetBranchAddress("stm_lowenergy", &eval_info.stm_lowenergy);
  tree0->SetBranchAddress("stm_LM", &eval_info.stm_LM);
  tree0->SetBranchAddress("stm_TGM", &eval_info.stm_TGM);

  tree0->SetBranchAddress("stm_STM", &eval_info.stm_STM);
  tree0->SetBranchAddress("stm_FullDead", &eval_info.stm_FullDead);
  tree0->SetBranchAddress("stm_clusterlength", &eval_info.stm_clusterlength);

  if (flag==1){
    tree0->SetBranchAddress("truth_nuEnergy", &eval_info.truth_nuEnergy);
    tree0->SetBranchAddress("truth_energyInside", &eval_info.truth_energyInside);
    
    tree0->SetBranchAddress("truth_electronInside", &eval_info.truth_electronInside);
    tree0->SetBranchAddress("truth_nuPdg", &eval_info.truth_nuPdg);
    tree0->SetBranchAddress("truth_isCC", &eval_info.truth_isCC);
    tree0->SetBranchAddress("truth_isEligible", &eval_info.truth_isEligible);
    tree0->SetBranchAddress("truth_NCisEligible", &eval_info.truth_NCisEligible);
    
    tree0->SetBranchAddress("truth_isFC", &eval_info.truth_isFC);
    tree0->SetBranchAddress("truth_vtxInside", &eval_info.truth_vtxInside);
    tree0->SetBranchAddress("truth_vtxX", &eval_info.truth_vtxX);
    tree0->SetBranchAddress("truth_vtxY", &eval_info.truth_vtxY);
    tree0->SetBranchAddress("truth_vtxZ", &eval_info.truth_vtxZ);
    
    tree0->SetBranchAddress("truth_nuTime", &eval_info.truth_nuTime);
    tree0->SetBranchAddress("match_completeness", &eval_info.match_completeness);
    tree0->SetBranchAddress("match_completeness_energy", &eval_info.match_completeness_energy);
    tree0->SetBranchAddress("match_purity", &eval_info.match_purity);
    tree0->SetBranchAddress("match_purity_xz", &eval_info.match_purity_xz);
    
    tree0->SetBranchAddress("match_purity_xy", &eval_info.match_purity_xy);
    tree0->SetBranchAddress("weight_spline", &eval_info.weight_spline);
    tree0->SetBranchAddress("weight_cv", &eval_info.weight_cv);
    tree0->SetBranchAddress("weight_flugg", &eval_info.weight_flugg);
    tree0->SetBranchAddress("weight_ppfx", &eval_info.weight_ppfx);
    tree0->SetBranchAddress("weight_lee", &eval_info.weight_lee);
    
    if (tree0->GetBranch("weight_change")) tree0->SetBranchAddress("weight_change",&eval_info.weight_change);
  }
}

void LEEana::put_tree_address(TTree *tree0, EvalInfo& eval_info, int flag){
  tree0->Branch("run", &eval_info.run,"data/I");
  tree0->Branch("subrun", &eval_info.subrun,"data/I");
  tree0->Branch("event", &eval_info.event,"data/I");
  tree0->Branch("flash_found", &eval_info.flash_found,"data/O");

  tree0->Branch("flash_time", &eval_info.flash_time,"data/F");

  if (eval_info.is_match_found_int){
    tree0->Branch("flash_found_asInt", &eval_info.flash_found_asInt, "data/I");
    tree0->Branch("match_found_asInt", &eval_info.match_found_asInt, "data/I"); 
  }
  if (eval_info.is_file_type){
    tree0->Branch("file_type",&eval_info.file_type);
  }
  
  tree0->Branch("flash_measPe", &eval_info.flash_measPe,"data/F");
  tree0->Branch("flash_predPe", &eval_info.flash_predPe,"data/F");
  tree0->Branch("match_found", &eval_info.match_found,"data/O");
  tree0->Branch("match_type", &eval_info.match_type,"data/i");
  tree0->Branch("match_isFC", &eval_info.match_isFC,"data/O");

  tree0->Branch("match_isTgm", &eval_info.match_isTgm,"data/O");
  tree0->Branch("match_notFC_FV", &eval_info.match_notFC_FV,"data/O");
  tree0->Branch("match_notFC_SP", &eval_info.match_notFC_SP,"data/O");
  tree0->Branch("match_notFC_DC", &eval_info.match_notFC_DC,"data/O");
  tree0->Branch("match_chargeU", &eval_info.match_chargeU,"data/F");
  
  tree0->Branch("match_chargeV", &eval_info.match_chargeV,"data/F");
  tree0->Branch("match_chargeY", &eval_info.match_chargeY,"data/F");
  tree0->Branch("match_energyY", &eval_info.match_energyY,"data/F");
  tree0->Branch("light_mismatch", &eval_info.light_mismatch,"data/O");
  tree0->Branch("match_charge", &eval_info.match_charge,"data/F");

  tree0->Branch("match_energy", &eval_info.match_energy,"data/F");
  tree0->Branch("stm_eventtype", &eval_info.stm_eventtype,"data/I");
  tree0->Branch("stm_lowenergy", &eval_info.stm_lowenergy,"data/I");
  tree0->Branch("stm_LM", &eval_info.stm_LM,"data/I");
  tree0->Branch("stm_TGM", &eval_info.stm_TGM,"data/I");

  tree0->Branch("stm_STM", &eval_info.stm_STM,"data/I");
  tree0->Branch("stm_FullDead", &eval_info.stm_FullDead,"data/I");
  tree0->Branch("stm_clusterlength", &eval_info.stm_clusterlength,"data/F");
  if (flag==1){
    tree0->Branch("truth_nuEnergy", &eval_info.truth_nuEnergy,"data/F");
    tree0->Branch("truth_energyInside", &eval_info.truth_energyInside,"data/F");
    
    tree0->Branch("truth_electronInside", &eval_info.truth_electronInside,"data/F");
    tree0->Branch("truth_nuPdg", &eval_info.truth_nuPdg,"data/I");
    tree0->Branch("truth_isCC", &eval_info.truth_isCC,"data/O");
    tree0->Branch("truth_isEligible", &eval_info.truth_isEligible,"data/O");
    tree0->Branch("truth_NCisEligible", &eval_info.truth_NCisEligible,"data/O");
    
    tree0->Branch("truth_isFC", &eval_info.truth_isFC,"data/O");
    tree0->Branch("truth_vtxInside", &eval_info.truth_vtxInside,"data/O");
    tree0->Branch("truth_vtxX", &eval_info.truth_vtxX,"data/F");
    tree0->Branch("truth_vtxY", &eval_info.truth_vtxY,"data/F");
    tree0->Branch("truth_vtxZ", &eval_info.truth_vtxZ,"data/F");
    
    tree0->Branch("truth_nuTime", &eval_info.truth_nuTime,"data/F");
    tree0->Branch("match_completeness", &eval_info.match_completeness,"data/F");
    tree0->Branch("match_completeness_energy", &eval_info.match_completeness_energy,"data/F");
    tree0->Branch("match_purity", &eval_info.match_purity,"data/F");
    tree0->Branch("match_purity_xz", &eval_info.match_purity_xz,"data/F");
    
    tree0->Branch("match_purity_xy", &eval_info.match_purity_xy,"data/F");
    tree0->Branch("weight_spline", &eval_info.weight_spline,"data/F");
    tree0->Branch("weight_cv", &eval_info.weight_cv,"data/F");
    tree0->Branch("weight_flugg", &eval_info.weight_flugg,"data/F");
    tree0->Branch("weight_ppfx", &eval_info.weight_ppfx,"data/F");
    tree0->Branch("weight_lee", &eval_info.weight_lee,"data/F");
    
    tree0->Branch("weight_change",&eval_info.weight_change,"data/O");
  }

  if (eval_info.flag_pl){

    tree0->Branch("pl_shr_energy_tot_cali",&eval_info.pl_shr_energy_tot_cali,"data/F");
    tree0->Branch("pl_shr_tkfit_dedx_Y",&eval_info.pl_shr_tkfit_dedx_Y,"data/F");
    tree0->Branch("pl_reco_nu_vtx_sce_x",&eval_info.pl_reco_nu_vtx_sce_x,"data/F");
    tree0->Branch("pl_reco_nu_vtx_sce_y",&eval_info.pl_reco_nu_vtx_sce_y,"data/F");
    tree0->Branch("pl_reco_nu_vtx_sce_z",&eval_info.pl_reco_nu_vtx_sce_z,"data/F");
    tree0->Branch("pl_reco_e",&eval_info.pl_reco_e,"data/D");
    tree0->Branch("pl_flag",&eval_info.pl_flag,"data/I");
  }

  if (eval_info.flag_gl){
    tree0->Branch("gl_sel_type",&eval_info.gl_sel_type,"gl_sel_type/I");
    tree0->Branch("gl_file_type",&eval_info.gl_file_type,"gl_file_type/I");
    tree0->Branch("gl_true_Enu",&eval_info.gl_true_Enu,"gl_true_Enu/F");
    tree0->Branch("gl_true_Elep",&eval_info.gl_true_Elep,"gl_true_Elep/F");
    tree0->Branch("gl_reco_Eshower",&eval_info.gl_reco_Eshower,"gl_reco_Eshower/F");
    tree0->Branch("gl_simple_pot_weight",&eval_info.gl_simple_pot_weight,"gl_simple_pot_weight/F");
    tree0->Branch("gl_rem_orig_wc_pot_weight",&eval_info.gl_rem_orig_wc_pot_weight,"gl_rem_orig_wc_pot_weight/F");
    tree0->Branch("gl_new_pot_weight",&eval_info.gl_new_pot_weight,"gl_new_pot_weight/F");
    tree0->Branch("gl_overlap_weight",&eval_info.gl_overlap_weight,"gl_overlap_weight/F");
    tree0->Branch("gl_wc_total_overlapped_weight",&eval_info.gl_wc_total_overlapped_weight,"gl_wc_total_overlapped_weight/F");
  }
  
}


#endif
