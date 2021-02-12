#include "TpcClustering.h"   

/// Tracking includes

#include <trackbase/TrkrHit.h>            // for TrkrHit
#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TMatrixFfwd.h>    // for TMatrixF
#include <TMatrixT.h>       // for TMatrixT, ope...
#include <TMatrixTUtils.h>  // for TMatrixTRow

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
TpcClustering::TpcClustering(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
TpcClustering::~TpcClustering()
{

}

//____________________________________________________________________________..
int TpcClustering::InitRun(PHCompositeNode *topNode)
{
  return Setup(topNode);
}

int TpcClustering::process_event(PHCompositeNode *topNode)
{
  return Process(topNode);
}

int TpcClustering::End(PHCompositeNode *topNode)
{
  return End();
}

//____________________________________________________________________________..

int TpcClustering::Setup(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return  Fun4AllReturnCodes::EVENT_OK;
}


/*
int TpcClustering::Process(PHCompositeNode *topNode)
{

  if(Verbosity() > 0) 
    std::cout << std::endl << "original size of hit map: " << _hit_map->size() << std::endl;  

  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange = _hit_map->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    // Every hitset corresponds to a single layer
    TrkrHitSet *hitset = hitsetitr->second;
    int layer = TrkrDefs::getLayer(hitsetitr->first);
    if (Verbosity() > 2)
      if (layer == print_layer)
      {
	cout << endl << PHWHERE << hitsetitr->first
	     << " layer " << (int) TrkrDefs::getLayer(hitsetitr->first)
	     << " side " << (int) TpcDefs::getSide(hitsetitr->first)
	     << " sector " << (int) TpcDefs::getSectorId(hitsetitr->first)
	     << " nhits " << (int) hitset->size() 
	     << endl;
      }

    // We do not cluster across sector boundaries, or across the central membrane.
    // So we want to do the clustering locally within a module
    // we have a single hitset, get the info that identifies the sector
    int sector = TpcDefs::getSectorId(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);

    // we will need the geometry object for this layer to get the number of Z and phi bins in the layer

    PHG4CylinderCellGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    int NPhiBins = layergeom->get_phibins();
    NPhiBinsMin = 0;
    NPhiBinsMax = NPhiBins;

    int NZBins = layergeom->get_zbins();
    if (side == 0)
    {
      NZBinsMin = 0;
      NZBinsMax = NZBins / 2;
    }
    else
    {
      NZBinsMin = NZBins / 2 + 1;
      NZBinsMax = NZBins;
    }

    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<double>> adcval(NPhiBins, std::vector<double>(NZBins, 0));






  return Fun4AllReturnCodes::EVENT_OK;
}
*/

int  TpcClustering::GetNodes(PHCompositeNode* topNode)
{
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
void TpcClustering::rotate_error(double erphi, double ez, double clusphi, double error[][3])
{
 TMatrixF ROT(3, 3);
  TMatrixF ERR(3,3);
 for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      {
	ROT[i][j] = 0;
	ERR[i][j] = 0;
      }

  ROT[0][0] = cos(clusphi);
  ROT[0][1] = -sin(clusphi);
  ROT[1][0] = sin(clusphi);
  ROT[1][1] = cos(clusphi);
  ROT[2][2] = 1.0;

  ERR[1][1] = erphi*erphi; 
  ERR[2][2] = ez*ez;
 
  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF COVAR_ERR(3, 3);
  COVAR_ERR = ROT * ERR * ROT_T;

  for(int i=0;i<3;++i)
    for(int j=0;j<3;++j)
      error[i][j] = COVAR_ERR[i][j];

}
*/
