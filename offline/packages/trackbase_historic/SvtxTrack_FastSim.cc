/*
 * SvtxTrack_FastSim.C
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#include "SvtxTrack_FastSim.h"

#include "SvtxTrack.h"  // for SvtxTrack::ConstClusterIter, SvtxTrack

#include <climits>
#include <map>          // for _Rb_tree_const_iterator
#include <ostream>      // for operator<<, basic_ostream, basic_ostream<>::_...

using namespace std;

SvtxTrack_FastSim::SvtxTrack_FastSim()
  : _truth_track_id(UINT_MAX)
  , _nmeas(0)
{
}

SvtxTrack_FastSim::~SvtxTrack_FastSim()
{
}

void SvtxTrack_FastSim::identify(std::ostream& os) const
{
  os << "SvtxTrack_FastSim Object ";
  os << "truth_track_id:" << get_truth_track_id() << endl;
  os << "id: " << get_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << endl;


  os << "clusters: ";
  for (SvtxTrack::ConstClusterKeyIter iter = begin_cluster_keys();
       iter != end_cluster_keys();
       ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;
      os << cluster_key << " ";
    }
  os << endl;

  return;
}

int SvtxTrack_FastSim::isValid() const
{
  return 1;
}
