#include "helixResiduals.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>
#include <micromegas/MicromegasDefs.h>
#include <trackbase/TrackFitUtils.h>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

helixResiduals::helixResiduals(const std::string &name)
  : SubsysReco(name)
{
}

int helixResiduals::InitRun(PHCompositeNode *topNode)
{
  const char *cfilepath = filepath.c_str();
  fout = new TFile(cfilepath, "recreate");
  tree_residuals = new TTree("tree_residuals", "Seed Residuals");
  tree_residuals->Branch("cluskey",&cluskey,"cluskey/l");
   tree_residuals->Branch("hitsetkey",&hitsetkey,"hitsetkey/i");
  tree_residuals->Branch("seed_id",&seed_id,"seed_id/I");
  tree_residuals->Branch("crossing",&crossing,"crossing/I");
  tree_residuals->Branch("isTpc",&isTpc,"isTpc/I");
  tree_residuals->Branch("layer",&layer,"layer/i");
  tree_residuals->Branch("nclus",&nclus,"nclus/I");
  tree_residuals->Branch("clusx",&clusx,"clusx/F");
  tree_residuals->Branch("clusy",&clusy,"clusy/F");
  tree_residuals->Branch("clusz",&clusz,"clusz/F");
  tree_residuals->Branch("pcax",&pcax,"pcax/F");
  tree_residuals->Branch("pcay",&pcay,"pcay/F");
  tree_residuals->Branch("pcaz",&pcaz,"pcaz/F");
  tree_residuals->Branch("pt",&pt,"pt/F");
  tree_residuals->Branch("px",&px,"px/F");
  tree_residuals->Branch("py",&py,"py/F");
  tree_residuals->Branch("pz",&pz,"pz/F");

  for(int ilayer=0;ilayer<57;++ilayer)
    {
      for(int ipar = 0; ipar < 6; ++ipar)
	{
	  double range = 0.5;  // mm
	  double range_angles = 0.03;  // rad
	  if( ilayer > 6 && ilayer < 56 )
	    {
	      if(ipar == 5)
		range = 3.0;
	      else 
		range = 2.0;
	    }
	  if( ilayer > 2 && ilayer < 7 )
	    {
              if(ipar == 5)
                range = 20.0;
            }
	  if(ipar < 3) range = range_angles;
	  
	  char name[500];
	  char title[500];
	  sprintf(name,"hpar_%i_%i", ilayer, ipar);
	  if(ilayer < 3)  sprintf(title,"MVTX parameter %i", ipar);
	  else if (ilayer > 2 && ilayer < 7) sprintf(title,"INTT parameter %i", ipar);
	  else if (ilayer > 6 && ilayer < 55) sprintf(title,"TPC parameter %i", ipar);
	  else  sprintf(title,"MMS parameters %i", ipar);
	  
	  hpar[ilayer][ipar] = new TH2D(name, title, 600, 0, 200, 2000, -range, +range);  // sensor number, parameter range
	  
	  hpar[ilayer][ipar]->GetXaxis()->SetNdivisions(504);
	  hpar[ilayer][ipar]->GetXaxis()->SetLabelSize(0.05);
	  hpar[ilayer][ipar]->GetXaxis()->SetTitleSize(0.05);
	  if(ipar < 3)
	    hpar[ilayer][ipar]->GetXaxis()->SetTitle("sensor");
	  else
	    hpar[ilayer][ipar]->GetXaxis()->SetTitle("sensor");	    	    
	}
    }
 
  getNodes(topNode);

  return 0;
}

int helixResiduals::End(PHCompositeNode * /**topNode*/)
{
  fout->cd();
  tree_residuals->Write();
  for(int lyr = 0; lyr < 57; ++lyr)
    for(int ipar=0; ipar<6; ++ipar)
      {
	hpar[lyr][ipar]->Write();
      }

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

void helixResiduals::fill_residuals(TrackSeed *seed, int id, bool tpc)
{
  if (seed->size_cluster_keys() == 0)
  {
    return;
  }

  std::vector<Acts::Vector3> clusterPositions;
  std::vector<TrkrDefs::cluskey> clusterKeys;
  getTrackletClusterList(seed, clusterKeys); 
  nclus = clusterKeys.size();
  // A silicon seed with 3 or 4 clusters is not useful, the helical fit will be close to perfect if there are only 3 clusters
  if(!tpc && nclus < 5) return;
  if(tpc && nclus < 20) return;

 TrackFitUtils::getTrackletClusters(tGeometry, _clusters, clusterPositions, clusterKeys);   
  std::vector<float> fitparams = TrackFitUtils::fitClusters(clusterPositions, clusterKeys);
  if(fitparams.size() == 0) 
    {
      std::cout << "Fit failed " << std::endl;
      return;
    }

  pt = seed->get_pt();
  px = seed->get_px(_clusters, tGeometry);
  py = seed->get_py(_clusters, tGeometry);
  pz = seed->get_pz();
  crossing = seed->get_crossing();
  seed_id = id;

  for (size_t i = 0; i < clusterPositions.size(); i++)
  {
    layer = TrkrDefs::getLayer(clusterKeys[i]);
    cluskey = clusterKeys[i];
    hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(clusterKeys[i]);
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

    pcax = pca(0);
    pcay = pca(1);
    pcaz = pca(2);
    clusx = position(0);
    clusy = position(1);
    clusz = position(2);
    isTpc = tpc;

    if(Verbosity() > 0)
      {
	float dx = position(0) - pca(0);
	float dy = position(1) - pca(1);
	float dz = position(2) - pca(2);
	std::cout << "   layer " << layer << " nclus " << nclus << " position(microns) " << position(0)*10000 << "  " << position(1)*10000 << "  " << position(2)*10000 << " pca " << pca(0)*10000 << "  " << pca(1)*10000 << "  " << pca(2)*10000 << " dx " << dx*10000 << "  dy " << dy*1000 << " dz " << dz*10000 << std::endl;
      }

    tree_residuals->Fill();

    // make a 2D histogram suitable for plotting by "detailed_alignment_parameter_plots.C" 

    // convert to mm for plotting
    clusx *= 10; clusy *= 10; clusz *= 10;
    pcax *= 10; pcay *= 10; pcaz *= 10;

    unsigned int sensor = getSensor(hitsetkey);
    float pars[6] = {0,0,0,0,0,0};
    pars[3] = clusx-pcax;
    pars[4] = clusy-pcay;
    pars[5] = clusz-pcaz;
    
    for(int ipar=0;ipar<6;++ipar)
      {
	hpar[layer][ipar]->Fill(sensor, pars[ipar]);
      }
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
      unsigned int lyr = TrkrDefs::getLayer(key);
      if(lyr == 7 || lyr == 22 || lyr == 23 || lyr == 38 || lyr == 39) {continue;}

      cluskey_vec.push_back(key);
      
    } // end loop over clusters for this track 
}


unsigned int helixResiduals::getSensor(TrkrDefs::hitsetkey hskey)
{
  unsigned int sensor = 0;
  unsigned int trkrid = TrkrDefs::getTrkrId(hskey);

  if(trkrid == TrkrDefs::mvtxId)
    {
      unsigned int staveid = MvtxDefs::getStaveId(hskey);
      unsigned int chipid = MvtxDefs::getChipId(hskey);
      sensor = staveid*9 + chipid;
    }
  else if(trkrid == TrkrDefs::inttId)
    {
      unsigned int ladderzid = InttDefs::getLadderZId(hskey);
      unsigned int ladderphiid = InttDefs::getLadderPhiId(hskey);
      sensor = ladderphiid*4 + ladderzid;      
    }
  else if(trkrid == TrkrDefs::tpcId)
    {
      unsigned int sector = TpcDefs::getSectorId(hskey);
      unsigned int side = TpcDefs::getSide(hskey);
      sensor = side * 12 + sector; 
    }
  else
    {
      // trkrid == TrkrDefs::micromegasId
      unsigned int tile = MicromegasDefs::getTileId(hskey);
      sensor = tile;
    }

  return sensor;

}
