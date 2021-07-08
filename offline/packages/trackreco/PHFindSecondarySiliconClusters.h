// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHFINDSECONDARYSILICONCLUSTERS_H
#define PHFINDSECONDARYSILICONCLUSTERS_H

#include <trackreco/PHTrackPropagating.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class TrkrCluster;
class TF1;


class PHFindSecondarySiliconClusters : public PHTrackPropagating
{
 public:

  PHFindSecondarySiliconClusters(const std::string &name = "PHFindSecondarySiliconClusters");

  ~PHFindSecondarySiliconClusters() override;

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  void  line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b);
  void  line_fit(std::vector<std::pair<double,double>> points, double &a, double &b);
  void CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0);
  //  std::vector<double> GetCircleClusterResiduals(std::vector<std::pair<double,double>> points, double R, double X0, double Y0);
  // std::vector<double> GetLineClusterResiduals(std::vector<std::pair<double,double>> points, double A, double B);
  std::vector<TrkrCluster*> getTrackTpcClusters(SvtxTrack *_tracklet);
  std::vector<TrkrCluster*> getTrackSiliconClusters(SvtxTrack *_tracklet);
  std::vector<TrkrCluster*> findMatchedSiliconClusters(std::vector<TrkrCluster*> all_si_clusters, double R, double X0, double Y0, double A, double B, SvtxTrack*  track);
  std::vector<TrkrCluster*> getAllSiliconClusters();
  void findRoot(const double R, const double X0, const double Y0, double& x, double& y);				 

  std::string _track_map_name_silicon;
  
  SvtxTrack *_track{nullptr};

  unsigned int _min_tpc_layer = 7;
  unsigned int _max_tpc_layer = 54;

  double _z_proj= 0;


  double _xy_residual_cut = 0.15;
  double _z_residual_cut = 0.15;

};

#endif // PHFINDSECONDARYSILICONCLUSTERS_H
