// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		 TpcClustering.h
 *  \brief		 Base class for TPC clustering 
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TPCCLUSTERING_H
#define TPCCLUSTERING_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrClusterHitAssoc;

class TpcClustering : public SubsysReco
{
 public:

  TpcClustering(const std::string &name = "TpcClustering");

  virtual ~TpcClustering();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

protected:

  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// override in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process(PHCompositeNode *topNode) = 0;

  /// Called in SubsysReco::End
  virtual int End() = 0;

  TrkrHitSetContainer *m_hits = nullptr;
  TrkrClusterContainer *m_clusterlist = nullptr;
  TrkrClusterHitAssoc *m_clusterhitassoc = nullptr;


 private:

  int GetNodes(PHCompositeNode* topNode);


};

#endif // TPCCLUSTERING_H
