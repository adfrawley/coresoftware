#include "helixResiduals.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>

#include <TFile.h>
#include <TNtuple.h>

helixResiduals::helixResiduals(const std::string &name)
  : SubsysReco(name)
{
}

int helixResiduals::InitRun(PHCompositeNode *topNode)
{
  const char *cfilepath = filepath.c_str();
  fout = new TFile(cfilepath, "recreate");
  ntp_residuals = new TNtuple("ntp_residuals", "Seed Residuals", "seed_id:layer:nclus:dphi:dx:dy:dz:pcax:pcay:pcaz:x:y:z:pt:px:py:pz:crossing:isSilicon:isTpc");

  getNodes(topNode);

  return 0;
}

int helixResiduals::End(PHCompositeNode * /**topNode*/)
{
  fout->cd();
  ntp_residuals->Write();
  fout->Close();

  return 0;
}

int helixResiduals::getNodes(PHCompositeNode *topNode)
{
  _tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_tracks)
  {
    std::cerr << PHWHERE << "No SvtxTrackMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _tpc_seeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_seeds)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TpcTrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _si_seeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_si_seeds)
  {
    std::cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clusters)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cerr << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int helixResiduals::process_event(PHCompositeNode * /*topNode*/)
{
  for (auto tpcseed_iter = _tpc_seeds->begin(); tpcseed_iter != _tpc_seeds->end(); ++tpcseed_iter)
  {
    int id = _tpc_seeds->index(tpcseed_iter);
    TrackSeed *tpcseed = _tpc_seeds->get(id);
    if (!tpcseed)
    {
      continue;
    }
    fill_residuals(tpcseed, id, true);
  }
  for (auto siseed_iter = _si_seeds->begin(); siseed_iter != _si_seeds->end(); ++siseed_iter)
  {
    int id = _si_seeds->index(siseed_iter);
    TrackSeed *siseed = _si_seeds->get(id);
    if (!siseed)
    {
      continue;
    }
    fill_residuals(siseed, id, false);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void helixResiduals::fill_residuals(TrackSeed *seed, int seed_id, bool isTpc)
{
  if (seed->size_cluster_keys() == 0)
  {
    return;
  }

  std::vector<Acts::Vector3> clusterPositions;
  std::vector<TrkrDefs::cluskey> clusterKeys;
  getTrackletClusterList(seed, clusterKeys); 
  auto nclus = clusterKeys.size();
  // A silicon seed with 3 clusters is not useful
  // the helical fit will be close to perfect if there are only 3 clusters
  if(nclus < 5) return;

 TrackFitUtils::getTrackletClusters(tGeometry, _clusters, clusterPositions, clusterKeys);   
  std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
  if(fitparams.size() == 0) 
    {
      std::cout << "Fit failed " << std::endl;
      return;
    }

  float pt = seed->get_pt();
  float px = seed->get_px(_clusters, tGeometry);
  float py = seed->get_py(_clusters, tGeometry);
  float pz = seed->get_pz();
  unsigned int crossing = seed->get_crossing();

  for (size_t i = 0; i < clusterPositions.size(); i++)
  {
    unsigned int layer = TrkrDefs::getLayer(clusterKeys[i]);
    Acts::Vector3 position = clusterPositions[i];
   Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitparams, position);
    float cluster_phi = atan2(position(1), position(0));
    float pca_phi = atan2(pca(1), pca(0));
    float dphi = cluster_phi - pca_phi;
    if (dphi > M_PI)
    {
      dphi = 2 * M_PI - dphi;
    }
    if (dphi < -M_PI)
    {
      dphi = 2 * M_PI + dphi;
    }
    float dx = position(0) - pca(0);
    float dy = position(1) - pca(1);
    float dz = position(2) - pca(2);
    float pcax = pca(0);
    float pcay = pca(1);
    float pcaz = pca(2);
    float x = position(0);
    float y = position(1);
    float z = position(2);

    if(Verbosity() > 0)
      {
	std::cout << "   layer " << layer << " nclus " << nclus << " position(microns) " << position(0)*10000 << "  " << position(1)*10000 << "  " << position(2)*10000 << " pca " << pca(0)*10000 << "  " << pca(1)*10000 << "  " << pca(2)*10000 << " dx " << dx*10000 << "  dy " << dy*1000 << " dz " << dz*10000 << std::endl;
      }

    float data[20] = {(float)seed_id, (float)layer, (float)nclus, dphi, dx, dy, dz, pcax, pcay, pcaz, x, y, z, pt, px, py, pz, (float)crossing, (float)!isTpc, (float)isTpc};

    ntp_residuals->Fill(data);
  }
}

void helixResiduals::getTrackletClusterList(TrackSeed *tracklet, std::vector<TrkrDefs::cluskey>& cluskey_vec)
{
  for (auto clusIter = tracklet->begin_cluster_keys();
       clusIter != tracklet->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = _clusters->findCluster(key);
      if(!cluster)
	{
	  std::cout << "Failed to get cluster with key " << key << std::endl;
	  continue;
	}	  
      
      /// Make a safety check for clusters that couldn't be attached to a surface
      auto surf = tGeometry->maps().getSurface(key, cluster);
      if(!surf)  
	{
	  std::cout << "Failed to find surface for cluster key " << key << std::endl; 
	  continue; 
	}

      // drop some bad layers in the TPC completely
      unsigned int layer = TrkrDefs::getLayer(key);
      if(layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39) {continue;}

      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}
