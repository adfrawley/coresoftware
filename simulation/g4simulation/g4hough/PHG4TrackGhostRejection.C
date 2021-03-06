#include "PHG4TrackGhostRejection.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <SvtxTrackMap.h>
#include <SvtxTrack.h>

// standard includes
#include <iostream>
#include <algorithm>

using namespace std;

bool hit_sort(unsigned int i, unsigned int j) { return (i < j);}

PHG4TrackGhostRejection::PHG4TrackGhostRejection(int nlayers, const string &name) :
SubsysReco(name)
{
  verbosity = 0;
  _nlayers = nlayers;
  _max_shared_hits = _nlayers;
  _layer_enabled.assign(_nlayers,true);
  _overlapping.clear();
  _candidates.clear();
}

int PHG4TrackGhostRejection::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::InitRun(PHCompositeNode *topNode)
{
  if (verbosity >= 0) {
    cout << "================== PHG4TrackGhostRejection::InitRun() =====================" << endl;
    cout << " CVS Version: $Id: PHG4TrackGhostRejection.C,v 1.12 2015/04/21 23:47:10 pinkenbu Exp $" << endl;
    cout << " Maximum allowed shared hits: " << _max_shared_hits << endl;
    for (unsigned int i=0;i<_layer_enabled.size();++i) {
      cout << " Enabled for hits in layer #" << i << ": " << boolalpha << _layer_enabled[i] << noboolalpha << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::process_event(PHCompositeNode *topNode)
{
  if(verbosity > 0) cout << "PHG4TrackGhostRejection::process_event -- entered" << endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  // Pull the reconstructed track information off the node tree...
  _g4tracks = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if(!_g4tracks) 
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if (verbosity > 1) {
    _g4tracks->identify();
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = &iter->second;
      track->identify();
    }
  }

  //----------------------------
  // Sort the hits on each track
  //----------------------------

  _candidates.clear();
  
  for (SvtxTrackMap::Iter iter = _g4tracks->begin();
       iter != _g4tracks->end();
       ++iter) {

    SvtxTrack* track = &iter->second;
  
    PHG4TrackCandidate combo;

    combo.trackid = track->getTrackID();
    combo.nhits = track->getNhits();

    for (unsigned int j = 0; j < _nlayers; ++j) {
      if (!_layer_enabled[j]) continue;
      
      if (track->hasCluster(j)) {
	combo.hitids.push_back(track->getClusterID(j));
      }
    }

    if (track->getNDF() != 0) {
      combo.chisq = track->getChisq()/track->getNDF();
    }

    combo.keep = true;

    // sort the hits by index
    stable_sort(combo.hitids.begin(),combo.hitids.end(),hit_sort);

    _candidates.push_back(combo);
  }

  //---------------------
  // Fill the overlap map
  //---------------------

  _overlapping.clear();
  for(unsigned int i = 0; i < _candidates.size(); i++)
    {
      for(unsigned int j = i+1; j < _candidates.size(); j++)
  	{ 
	  // determine the maximum length of the anticipated storage
	  unsigned maxhits = _candidates[i].hitids.size();
	  if(_candidates[j].hitids.size() > maxhits)
	    {
	      maxhits = _candidates[j].hitids.size();
	    }

	  // create the difference storage
	  std::vector<unsigned int> diff;
	  diff.assign(maxhits,0);

	  // run the difference algorithm
	  std::vector<unsigned int>::iterator it = diff.begin();
	  it = std::set_difference(_candidates[i].hitids.begin(),_candidates[i].hitids.end(),
				   _candidates[j].hitids.begin(),_candidates[j].hitids.end(),
				   diff.begin());

	  // calculate the overlap
	  unsigned int overlap = maxhits - int(it - diff.begin());

	  // insert an overlapping pair into the map
	  if(overlap > _max_shared_hits) _overlapping.insert(std::make_pair(i,j));	  
   	}
    }

  //----------------------
  // Flag the ghost tracks
  //----------------------

  std::multimap<unsigned,unsigned int>::iterator iter;
  for (iter = _overlapping.begin(); iter != _overlapping.end(); iter++) {

    unsigned int key = iter->first;
    unsigned int value = iter->second;

    if (_candidates[key].nhits > _candidates[value].nhits) {
      // prefer longer track
      _candidates[value].keep = false;
    } else if (_candidates[key].nhits < _candidates[value].nhits) {
      // prefer longer track
      _candidates[key].nhits = false;
    } else {
      // choose between equal length tracks by chisq/dof
      if (_candidates[key].chisq < _candidates[value].chisq) {
	_candidates[value].keep = false;
      } else {
	_candidates[key].keep = false;
      }
    }
  }

  //------------------------
  // Remove the ghost tracks
  //------------------------

  int initial_size = _g4tracks->size();

  // loop over container and delete!
  for (unsigned int i = 0; i < _candidates.size(); i++) {
    if (!_candidates[i].keep) {
      // look for the track to delete
      if (_g4tracks->find(_candidates[i].trackid) != _g4tracks->end()) {
	_g4tracks->erase(_candidates[i].trackid);
      }
    }
  }

  if (verbosity > 1) {
    _g4tracks->identify();
    for (SvtxTrackMap::Iter iter = _g4tracks->begin();
	 iter != _g4tracks->end();
	 ++iter) {
      SvtxTrack *track = &iter->second;
      track->identify();
    }
  }

  if(verbosity > 0)
    cout << "PHG4TrackGhostRejection - rejected and removed " 
         << initial_size - _g4tracks->size()
         << " tracks" << endl;;
  
  if(verbosity > 0) cout << "PHG4TrackGhostRejection::process_event -- exited" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TrackGhostRejection::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

