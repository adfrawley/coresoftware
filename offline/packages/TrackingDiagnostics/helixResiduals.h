#ifndef HELIXRESIDUALS_H
#define HELIXRESIDUALS_H

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>

#include <fun4all/SubsysReco.h>

class TFile;
class TTree;
class TH2D;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

class helixResiduals : public SubsysReco
{
 public:
  helixResiduals(const std::string& name = "helixResiduals");
  ~helixResiduals() {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* /*topNode*/) override;
  int End(PHCompositeNode* /*topNode*/) override;
  void set_output_file(const std::string& outputfile) { filepath = outputfile; }

 private:
  int getNodes(PHCompositeNode* topNode);

  void fill_residuals(TrackSeed* seed, int seed_id, bool isTpc);
  void getTrackletClusterList(TrackSeed *tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec);
  unsigned int getSensor(TrkrDefs::hitsetkey hitsetkey);

  ActsGeometry* tGeometry = nullptr;

  //  TNtuple* ntp_residuals = nullptr;

  SvtxTrackMap* _tracks = nullptr;
  TrkrClusterContainer* _clusters = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;

  TFile* fout = nullptr;

  std::string filepath = "";

  float clusx,clusy,clusz,pcax,pcay,pcaz,pt,px,py,pz;
  int seed_id, isTpc,nclus;
  int crossing;
  unsigned int layer;
  TrkrDefs::cluskey cluskey;
  TrkrDefs::hitsetkey hitsetkey;

  TTree *tree_residuals;
  TH2D *hpar[57][6];

};

#endif  // HELIXRESIDUALS_H
