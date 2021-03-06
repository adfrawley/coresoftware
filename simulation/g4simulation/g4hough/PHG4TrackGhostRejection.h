#ifndef __PHG4TRACKGHOSTREJECTION_H__
#define __PHG4TRACKGHOSTREJECTION_H__

//===========================================================
/// \file PHG4TrackGhostRejection.h
/// \brief Quick N-shared hit rejection routine
/// \author Mike McCumber
//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>

// standard includes
#include <vector>
#include <map>

// forward declarations
class PHCompositeNode;
class SvtxTrackMap;

class PHG4TrackCandidate
{
  
 public:

  int trackid;
  short nhits;
  std::vector<unsigned int> hitids;
  float chisq;
  bool keep;

};

/// \class PHG4TrackGhostRejection
///
/// \brief Quick N-shared hit rejection routine
///
/// This module runs after the pattern recognition to remove
/// track candidates with a user defined overlap. The steps are:
/// (1) Sort the hits on each track by index
/// (2) Fill a multi-map between a track and all of the overlapping tracks
/// (3) Keep only the best track from the set of overlapping tracks
///
class PHG4TrackGhostRejection : public SubsysReco
{

 public:
 
  PHG4TrackGhostRejection(int nlayers, const std::string &name = "PHG4TrackGhostRejection");
  virtual ~PHG4TrackGhostRejection() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
  void set_max_shared_hits(unsigned int nhits) { _max_shared_hits = nhits; }
  unsigned int get_max_shared_hits() { return _max_shared_hits; }

  void set_layer_enabled(int layer, bool enabled) {_layer_enabled[layer] = enabled;}
  bool get_layer_enabled(int layer) {return _layer_enabled[layer];}

 private:

  SvtxTrackMap *_g4tracks;
  std::vector< PHG4TrackCandidate > _candidates;

  unsigned int _nlayers;
  unsigned int _max_shared_hits;
  std::vector<bool> _layer_enabled;

  std::multimap< unsigned int, unsigned int > _overlapping;

};

#endif // __PHG4TRACKGHOSTREJECTION_H__
