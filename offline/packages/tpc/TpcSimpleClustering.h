#ifndef TPC_TPCSIMPLECLUSTERING_H
#define TPC_TPCSIMPLECLUSTERING_H

#include "TpcClustering.h"  // base class

#include <fun4all/SubsysReco.h>
#include <map> 
#include <vector>
#include <string>

class PHCompositeNode;
class TrkrHitSet;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4CylinderCellGeom;

typedef std::pair<int, int> iphiz;
typedef std::pair<double, iphiz> ihit;

class TpcSimpleClustering : public TpcClustering
{
 public:
  TpcSimpleClustering(const std::string &name = "TpcSimpleClustering");
  virtual ~TpcSimpleClustering(){}

  void set_sector_fiducial_cut(const double cut){SectorFiducialCut = cut; }
  void set_search_bins(const int bins){NSearch = bins;}

 protected:

  virtual int Setup(PHCompositeNode *topNode);
  virtual int Process(PHCompositeNode *topNode);
  virtual int End();


 private:
  bool is_local_maximum(int phi, int z, std::vector<std::vector<double>> &adcval);
  bool is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom);
  void get_cluster(int phibin, int zbin, std::vector<std::vector<double>> &adcval, std::vector<ihit> &ihit_list);
  void find_z_range(int phibin, int zbin, std::vector<std::vector<double>> &adcval, int& zdown, int& zup);
  void find_phi_range(int phibin, int zbin, std::vector<std::vector<double>> &adcval, int& phidown, int& phiup);
  void remove_hit(double adc, int phibin, int zbin, std::multimap<double, ihit> &all_hit_map, std::vector<std::vector<double>> &adcval );
  void remove_hits(std::vector<ihit> &ihit_list, std::multimap<double, ihit> &all_hit_map, std::vector<std::vector<double>> &adcval);
  void calc_cluster_parameter(std::vector<ihit> &ihit_list, int iclus, PHG4CylinderCellGeom *layergeom, TrkrHitSet *hitset);
  void print_cluster(std::vector<ihit> &ihit_list);

  double zz_shaping_correction;
  double pedestal;

  double SectorFiducialCut;
  int NSearch;

  int NPhiBinsMax;
  int NPhiBinsMin;
  int NZBinsMax;
  int NZBinsMin;

};

#endif
