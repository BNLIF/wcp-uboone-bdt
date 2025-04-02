#ifndef UBOONE_LEE_SPACE
#define UBOONE_LEE_SPACE

namespace LEEana{
struct SpaceInfo{
  Int_t run;
  Int_t subrun;
  Int_t event;
  std::vector<double> *Trec_spacepoints_x;
  std::vector<double> *Trec_spacepoints_y;
  std::vector<double> *Trec_spacepoints_z;
  std::vector<double> *Trec_spacepoints_q;
  std::vector<double> *Trec_spacepoints_cluster_id;
  std::vector<double> *Trec_spacepoints_real_cluster_id;
  std::vector<double> *Trec_spacepoints_sub_cluster_id;

  std::vector<double> *Treccharge_spacepoints_x;
  std::vector<double> *Treccharge_spacepoints_y;
  std::vector<double> *Treccharge_spacepoints_z;
  std::vector<double> *Treccharge_spacepoints_q;
  std::vector<double> *Treccharge_spacepoints_cluster_id;
  std::vector<double> *Treccharge_spacepoints_real_cluster_id;
  std::vector<double> *Treccharge_spacepoints_sub_cluster_id;

  std::vector<double> *Trecchargeblob_spacepoints_x;
  std::vector<double> *Trecchargeblob_spacepoints_y;
  std::vector<double> *Trecchargeblob_spacepoints_z;
  std::vector<double> *Trecchargeblob_spacepoints_q;
  std::vector<double> *Trecchargeblob_spacepoints_cluster_id;
  std::vector<double> *Trecchargeblob_spacepoints_real_cluster_id;
  std::vector<double> *Trecchargeblob_spacepoints_sub_cluster_id;
};

 void set_tree_address(TTree *tree0, SpaceInfo& space_info, bool flag_runinfo = true);
 void put_tree_address(TTree *tree0, SpaceInfo& space_info);
 void clear_space_info(SpaceInfo& space_info);
 void init_pointers(SpaceInfo& space_info);
 void del_pointers(SpaceInfo& space_info);

}

void LEEana::init_pointers(SpaceInfo& space_info) {
  space_info.Trec_spacepoints_x = new std::vector<double>;
  space_info.Trec_spacepoints_y = new std::vector<double>;
  space_info.Trec_spacepoints_z = new std::vector<double>;
  space_info.Trec_spacepoints_q = new std::vector<double>;
  space_info.Trec_spacepoints_cluster_id = new std::vector<double>;
  space_info.Trec_spacepoints_real_cluster_id = new std::vector<double>;
  space_info.Trec_spacepoints_sub_cluster_id = new std::vector<double>;

  space_info.Treccharge_spacepoints_x = new std::vector<double>;
  space_info.Treccharge_spacepoints_y = new std::vector<double>;
  space_info.Treccharge_spacepoints_z = new std::vector<double>;
  space_info.Treccharge_spacepoints_q = new std::vector<double>;
  space_info.Treccharge_spacepoints_cluster_id = new std::vector<double>;
  space_info.Treccharge_spacepoints_real_cluster_id = new std::vector<double>;
  space_info.Treccharge_spacepoints_sub_cluster_id = new std::vector<double>;

  space_info.Trecchargeblob_spacepoints_x = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_y = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_z = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_q = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_cluster_id = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_real_cluster_id = new std::vector<double>;
  space_info.Trecchargeblob_spacepoints_sub_cluster_id = new std::vector<double>;
}

void LEEana::del_pointers(SpaceInfo& space_info) {
  delete space_info.Trec_spacepoints_x;
  delete space_info.Trec_spacepoints_y;
  delete space_info.Trec_spacepoints_z;
  delete space_info.Trec_spacepoints_q;
  delete space_info.Trec_spacepoints_cluster_id;
  delete space_info.Trec_spacepoints_real_cluster_id;
  delete space_info.Trec_spacepoints_sub_cluster_id;

  delete space_info.Treccharge_spacepoints_x;
  delete space_info.Treccharge_spacepoints_y;
  delete space_info.Treccharge_spacepoints_z;
  delete space_info.Treccharge_spacepoints_q;
  delete space_info.Treccharge_spacepoints_cluster_id;
  delete space_info.Treccharge_spacepoints_real_cluster_id;
  delete space_info.Treccharge_spacepoints_sub_cluster_id;

  delete space_info.Trecchargeblob_spacepoints_x;
  delete space_info.Trecchargeblob_spacepoints_y;
  delete space_info.Trecchargeblob_spacepoints_z;
  delete space_info.Trecchargeblob_spacepoints_q;
  delete space_info.Trecchargeblob_spacepoints_cluster_id;
  delete space_info.Trecchargeblob_spacepoints_real_cluster_id;
  delete space_info.Trecchargeblob_spacepoints_sub_cluster_id;
}

void LEEana::clear_space_info(SpaceInfo& space_info) {
  space_info.Trec_spacepoints_x->clear();
  space_info.Trec_spacepoints_y->clear();
  space_info.Trec_spacepoints_z->clear();
  space_info.Trec_spacepoints_q->clear();
  space_info.Trec_spacepoints_cluster_id->clear();
  space_info.Trec_spacepoints_real_cluster_id->clear();
  space_info.Trec_spacepoints_sub_cluster_id->clear();

  space_info.Treccharge_spacepoints_x->clear();
  space_info.Treccharge_spacepoints_y->clear();
  space_info.Treccharge_spacepoints_z->clear();
  space_info.Treccharge_spacepoints_q->clear();
  space_info.Treccharge_spacepoints_cluster_id->clear();
  space_info.Treccharge_spacepoints_real_cluster_id->clear();
  space_info.Treccharge_spacepoints_sub_cluster_id->clear();

  space_info.Trecchargeblob_spacepoints_x->clear();
  space_info.Trecchargeblob_spacepoints_y->clear();
  space_info.Trecchargeblob_spacepoints_z->clear();
  space_info.Trecchargeblob_spacepoints_q->clear();
  space_info.Trecchargeblob_spacepoints_cluster_id->clear();
  space_info.Trecchargeblob_spacepoints_real_cluster_id->clear();
  space_info.Trecchargeblob_spacepoints_sub_cluster_id->clear();
}

void LEEana::set_tree_address(TTree *tree0, SpaceInfo& space_info, bool flag_runinfo) {
  init_pointers(space_info);
  if(flag_runinfo){
    tree0->SetBranchAddress("run", &space_info.run);
    tree0->SetBranchAddress("subrun", &space_info.subrun);
    tree0->SetBranchAddress("event", &space_info.event);
  }
  tree0->SetBranchAddress("Trec_spacepoints_x", &space_info.Trec_spacepoints_x);
  tree0->SetBranchAddress("Trec_spacepoints_y", &space_info.Trec_spacepoints_y);
  tree0->SetBranchAddress("Trec_spacepoints_z", &space_info.Trec_spacepoints_z);
  tree0->SetBranchAddress("Trec_spacepoints_q", &space_info.Trec_spacepoints_q);
  tree0->SetBranchAddress("Trec_spacepoints_cluster_id", &space_info.Trec_spacepoints_cluster_id);
  tree0->SetBranchAddress("Trec_spacepoints_real_cluster_id", &space_info.Trec_spacepoints_real_cluster_id);
  tree0->SetBranchAddress("Trec_spacepoints_sub_cluster_id", &space_info.Trec_spacepoints_sub_cluster_id);
  tree0->SetBranchAddress("Treccharge_spacepoints_x", &space_info.Treccharge_spacepoints_x);
  tree0->SetBranchAddress("Treccharge_spacepoints_y", &space_info.Treccharge_spacepoints_y);
  tree0->SetBranchAddress("Treccharge_spacepoints_z", &space_info.Treccharge_spacepoints_z);
  tree0->SetBranchAddress("Treccharge_spacepoints_q", &space_info.Treccharge_spacepoints_q);
  tree0->SetBranchAddress("Treccharge_spacepoints_cluster_id", &space_info.Treccharge_spacepoints_cluster_id);
  tree0->SetBranchAddress("Treccharge_spacepoints_real_cluster_id", &space_info.Treccharge_spacepoints_real_cluster_id);
  tree0->SetBranchAddress("Treccharge_spacepoints_sub_cluster_id", &space_info.Treccharge_spacepoints_sub_cluster_id);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_x", &space_info.Trecchargeblob_spacepoints_x);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_y", &space_info.Trecchargeblob_spacepoints_y);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_z", &space_info.Trecchargeblob_spacepoints_z);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_q", &space_info.Trecchargeblob_spacepoints_q);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_cluster_id", &space_info.Trecchargeblob_spacepoints_cluster_id);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_real_cluster_id", &space_info.Trecchargeblob_spacepoints_real_cluster_id);
  tree0->SetBranchAddress("Trecchargeblob_spacepoints_sub_cluster_id", &space_info.Trecchargeblob_spacepoints_sub_cluster_id);
}

void LEEana::put_tree_address(TTree *tree0, SpaceInfo& space_info) {
  // add new branches for run, subrun, event
  tree0->Branch("run", &space_info.run, "run/I");
  tree0->Branch("subrun", &space_info.subrun, "subrun/I");
  tree0->Branch("event", &space_info.event, "event/I");
  tree0->Branch("Trec_spacepoints_x", &space_info.Trec_spacepoints_x);
  tree0->Branch("Trec_spacepoints_y", &space_info.Trec_spacepoints_y);
  tree0->Branch("Trec_spacepoints_z", &space_info.Trec_spacepoints_z);
  tree0->Branch("Trec_spacepoints_q", &space_info.Trec_spacepoints_q);
  tree0->Branch("Trec_spacepoints_cluster_id", &space_info.Trec_spacepoints_cluster_id);
  tree0->Branch("Trec_spacepoints_real_cluster_id", &space_info.Trec_spacepoints_real_cluster_id);
  tree0->Branch("Trec_spacepoints_sub_cluster_id", &space_info.Trec_spacepoints_sub_cluster_id);
  tree0->Branch("Treccharge_spacepoints_x", &space_info.Treccharge_spacepoints_x);
  tree0->Branch("Treccharge_spacepoints_y", &space_info.Treccharge_spacepoints_y);
  tree0->Branch("Treccharge_spacepoints_z", &space_info.Treccharge_spacepoints_z);
  tree0->Branch("Treccharge_spacepoints_q", &space_info.Treccharge_spacepoints_q);
  tree0->Branch("Treccharge_spacepoints_cluster_id", &space_info.Treccharge_spacepoints_cluster_id);
  tree0->Branch("Treccharge_spacepoints_real_cluster_id", &space_info.Treccharge_spacepoints_real_cluster_id);
  tree0->Branch("Treccharge_spacepoints_sub_cluster_id", &space_info.Treccharge_spacepoints_sub_cluster_id);
  tree0->Branch("Trecchargeblob_spacepoints_x", &space_info.Trecchargeblob_spacepoints_x);
  tree0->Branch("Trecchargeblob_spacepoints_y", &space_info.Trecchargeblob_spacepoints_y);
  tree0->Branch("Trecchargeblob_spacepoints_z", &space_info.Trecchargeblob_spacepoints_z);
  tree0->Branch("Trecchargeblob_spacepoints_q", &space_info.Trecchargeblob_spacepoints_q);
  tree0->Branch("Trecchargeblob_spacepoints_cluster_id", &space_info.Trecchargeblob_spacepoints_cluster_id);
  tree0->Branch("Trecchargeblob_spacepoints_real_cluster_id", &space_info.Trecchargeblob_spacepoints_real_cluster_id);
  tree0->Branch("Trecchargeblob_spacepoints_sub_cluster_id", &space_info.Trecchargeblob_spacepoints_sub_cluster_id);
}

#endif
