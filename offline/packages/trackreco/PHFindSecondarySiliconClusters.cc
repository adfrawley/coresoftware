#include "PHFindSecondarySiliconClusters.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>
using namespace std;

//____________________________________________________________________________..
PHFindSecondarySiliconClusters::PHFindSecondarySiliconClusters(const std::string &name):
 PHTrackPropagating(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHFindSecondarySiliconClusters::PHFindSecondarySiliconClusters(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHFindSecondarySiliconClusters::~PHFindSecondarySiliconClusters()
{

}

//____________________________________________________________________________..
int PHFindSecondarySiliconClusters::Setup(PHCompositeNode *topNode)
{
  std::cout << PHWHERE << " Parameters:  _xy_residual_cut " << _xy_residual_cut
	    << " _z_residual_cut " << _z_residual_cut 
	    << std::endl; 

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHFindSecondarySiliconClusters::Process()
{
  // _track_map contains the tracks
  // We want to identify and associate the TPC only tracks with silicon clusters

  // get a list of all clusteres in the silicon layers
  std::vector<TrkrCluster*> all_si_clusters = getAllSiliconClusters();

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size()  << endl;

  // loop over the TPC track seeds
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      _track = phtrk_iter->second;

      // get the silicon clusters for this track, skip if there are some
      std::vector<TrkrCluster*> si_clusters = getTrackSiliconClusters(_track);
      if(si_clusters.size() > 0) continue;

      // This is a track with TPC clusters only      

      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing track: " << phtrk_iter->first
	    << ": nhits: " << _track-> size_cluster_keys()
	    << ": pT: " << _track->get_pt()
	    << ": phi: " << _track->get_phi()
	    << ": eta: " << _track->get_eta()
	    << endl;
	}

      // get the tpc cluster positions in z and r

      // Get the TPC clusters for this track
      std::vector<TrkrCluster*> clusters = getTrackTpcClusters(_track);

      // count TPC layers for this track
      std::set<unsigned int> layers;
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  unsigned int layer = TrkrDefs::getLayer(clusters[i]->getClusKey());
	  layers.insert(layer);
	}
      unsigned int nlayers = layers.size();
      if(Verbosity() > 2) std::cout << "    TPC layers this track: " << nlayers << std::endl;

      if(clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc tracklet, not enough TPC clusters " << std::endl; 
	  continue;  // skip to the next TPC tracklet
	}

      //================
      // Get track Z intercept
      //================

      std::vector<std::pair<double, double>> points;      
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  double z = clusters[i]->getZ();
	  double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));
	  
	  points.push_back(make_pair(r,z));
	}
            
      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit(points, A, B);
      _z_proj = B;
      if(Verbosity() > 2) std::cout << " First fitted line has _z_proj " << _z_proj << std::endl;

      //===================
      // Get tracklet PCA in (x,y)
      //===================

      // make circle fit to TPC clusters
      std::vector<std::pair<double, double>> cpoints;
      for (unsigned int i=0; i<clusters.size(); ++i)
	{
	  cpoints.push_back(make_pair(clusters[i]->getX(), clusters[i]->getY()));
	}
      double R, X0, Y0;
      CircleFitByTaubin(cpoints, R, X0, Y0);
      if(Verbosity() > 2) 
      std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // we now have a TPC track with no silicon match, and a circle fit to the clusters
      // see if any silicon clusters are along the projection
      std::vector<TrkrCluster*> matched_si_clusters = findMatchedSiliconClusters(all_si_clusters, R, X0, Y0, A, B, _track);
      if(matched_si_clusters.size() == 0) continue;
      // add to the track and the association map
      for(auto cluster : matched_si_clusters)
	{
	  if(Verbosity() > 1) std::cout << "  inserting cluster " << cluster->getClusKey() << " into track " << _track->get_id() << std::endl;
	  _track->insert_cluster_key(cluster->getClusKey());
	  _assoc_container->SetClusterTrackAssoc(cluster->getClusKey(), _track->get_id());
	  clusters.push_back(cluster);

	  double z = cluster->getZ();
	  double r = sqrt(pow(cluster->getX(),2) + pow(cluster->getY(), 2));	  
	  points.push_back(make_pair(r,z));
	  cpoints.push_back(make_pair(cluster->getX(), cluster->getY()));
	}

      // finally, redo the helical fit with all clusters and update the track PCA
      //====================================================         
      // get the straight line representing the z trajectory in the form of z vs radius
      line_fit(points, A, B);
      _z_proj = B;
      if(Verbosity() > 2) std::cout << " Final fitted line has _z_proj " << _z_proj << std::endl;

      // make circle fit to all clusters
      CircleFitByTaubin(cpoints, R, X0, Y0);
      if(Verbosity() > 2) 
	std::cout << " Final fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // Get PCA from circle fit
      double pcax = NAN;
      double pcay = NAN;
      findRoot(R, X0, Y0, pcax, pcay);
      if(std::isnan(pcax))
	{
	  std::cout << PHWHERE << " bad PCA found " << std::endl;
	  pcax = 9999.0;
	  pcay = 9999.0;
	  continue;
	}

      // update track position
      _track->set_x(pcax);
      _track->set_y(pcay);
      _track->set_z(_z_proj);
      
    }  // end loop over tracks

  if (Verbosity() > 0)
    cout << "PHFindSecondarySiliconClusters::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHFindSecondarySiliconClusters::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHFindSecondarySiliconClusters::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  return Fun4AllReturnCodes::EVENT_OK;
}


void  PHFindSecondarySiliconClusters::line_fit(std::vector<std::pair<double,double>> points, double &a, double &b)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  
   double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
   for (unsigned int i=0; i<points.size(); ++i)
    {
      double r = points[i].first;
      double z = points[i].second;

      xsum=xsum+r;                        //calculate sigma(xi)
      ysum=ysum+z;                        //calculate sigma(yi)
      x2sum=x2sum+pow(r,2);                //calculate sigma(x^2i)
      xysum=xysum+r*z;                    //calculate sigma(xi*yi)
    }
   a=(points.size()*xysum-xsum*ysum)/(points.size()*x2sum-xsum*xsum);            //calculate slope
   b=(x2sum*ysum-xsum*xysum)/(x2sum*points.size()-xsum*xsum);            //calculate intercept

   if(Verbosity() > 10)
     {
       for (unsigned int i=0;i<points.size(); ++i)
	 {
	   double r = points[i].first;
	   double z_fit = a * r + b;                    //to calculate z(fitted) at given r points
	   std::cout << " r " << r << " z " << points[i].second << " z_fit " << z_fit << std::endl; 
	 } 
     }

    return;
}   

void  PHFindSecondarySiliconClusters::line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b)
{
  std::vector<std::pair<double,double>> points;
  
   for (unsigned int i=0; i<clusters.size(); ++i)
     {
       double z = clusters[i]->getZ();
       double r = sqrt(pow(clusters[i]->getX(),2) + pow(clusters[i]->getY(), 2));

       points.push_back(make_pair(r,z));
     }

   line_fit(points, a, b);

    return;
}

void PHFindSecondarySiliconClusters::CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0)
/*  
      Circle fit to a given set of data points (in 2D)
      This is an algebraic fit, due to Taubin, based on the journal article
      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With 
                  Applications To Edge And Range Image Segmentation",
                  IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
*/
{
  int iter,IterMAX=99;
  
  double Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
  double A0,A1,A2,A22,A3,A33;
  double x,y;
  double DET,Xcenter,Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;
  for(unsigned int i = 0; i < points.size(); ++i)
    {
      meanX += points[i].first;
      meanY += points[i].second;
      weight++;
    }
  meanX /= weight;
  meanY /= weight;

  //     computing moments 
  
  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;
  
  for (unsigned int i=0; i<points.size(); i++)
    {
      double Xi = points[i].first - meanX;   //  centered x-coordinates
      double Yi = points[i].second - meanY;   //  centered y-coordinates
      double Zi = Xi*Xi + Yi*Yi;
      
      Mxy += Xi*Yi;
      Mxx += Xi*Xi;
      Myy += Yi*Yi;
      Mxz += Xi*Zi;
      Myz += Yi*Zi;
      Mzz += Zi*Zi;
    }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;
  
  //  computing coefficients of the characteristic polynomial
  
  Mz = Mxx + Myy;
  Cov_xy = Mxx*Myy - Mxy*Mxy;
  Var_z = Mzz - Mz*Mz;
  A3 = 4*Mz;
  A2 = -3*Mz*Mz - Mzz;
  A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;
  
  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0  
  //    (it is guaranteed to converge to the right root)
  
  for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x*(A22 + A33*x);
      double xnew = x - y/Dy;
      if ((xnew == x)||(!isfinite(xnew))) break;
      double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  assembling the output
  
  X0 = Xcenter + meanX;
  Y0 = Ycenter + meanY;
  R = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);
}

std::vector<TrkrCluster*> PHFindSecondarySiliconClusters::getTrackTpcClusters(SvtxTrack *_tracklet_tpc)
{
  std::vector<TrkrCluster*> clusters;
  
  for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
       key_iter != _tracklet_tpc->end_cluster_keys();
       ++key_iter)
    {
      TrkrDefs::cluskey cluster_key = *key_iter;
      unsigned int layer = TrkrDefs::getLayer(cluster_key);
      
      if(layer < _min_tpc_layer) continue;
      if(layer > _max_tpc_layer) continue;
      
      // get the cluster
      TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);
      
      //tpc_clusters_map.insert(std::make_pair(layer, tpc_clus));
      clusters.push_back(tpc_clus);
      
      if(Verbosity() > 5) 
	std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		  << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " clusters.size() " << clusters.size() << std::endl;
    }
  return clusters;
}

std::vector<TrkrCluster*> PHFindSecondarySiliconClusters::getTrackSiliconClusters(SvtxTrack *_tracklet_tpc)
{
  std::vector<TrkrCluster*> clusters;
  
  for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_tpc->begin_cluster_keys();
       key_iter != _tracklet_tpc->end_cluster_keys();
       ++key_iter)
    {
      TrkrDefs::cluskey cluster_key = *key_iter;
      unsigned int layer = TrkrDefs::getLayer(cluster_key);
      
      if(layer >= _min_tpc_layer) continue;
      
      // get the cluster
      TrkrCluster *si_clus =  _cluster_map->findCluster(cluster_key);
      
      //tpc_clusters_map.insert(std::make_pair(layer, tpc_clus));
      clusters.push_back(si_clus);
      
      if(Verbosity() > 5) 
	std::cout << "  Solicon cluster in layer " << layer << " with position " << si_clus->getX() 
		  << "  " << si_clus->getY() << "  " << si_clus->getZ() << " clusters.size() " << clusters.size() << std::endl;
    }
  return clusters;
}

void PHFindSecondarySiliconClusters::findRoot(const double R, const double X0,
				    const double Y0, double& x,
				    double& y)
{
  /**
   * We need to determine the closest point on the circle to the origin
   * since we can't assume that the track originates from the origin
   * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
   * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the 
   * origin and some (currently, unknown) point on the circle x,y.
   * 
   * Solving the circle eqn for x and substituting into d gives an eqn for
   * y. Taking the derivative and setting equal to 0 gives the following 
   * two solutions. We take the smaller solution as the correct one, as 
   * usually one solution is wildly incorrect (e.g. 1000 cm)
   */
  
  double miny = (sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double miny2 = (-sqrt(pow(X0, 2) * pow(R, 2) * pow(Y0, 2) + pow(R, 2) 
		      * pow(Y0,4)) + pow(X0,2) * Y0 + pow(Y0, 3)) 
    / (pow(X0, 2) + pow(Y0, 2));

  double minx = sqrt(pow(R, 2) - pow(miny - Y0, 2)) + X0;
  double minx2 = -sqrt(pow(R, 2) - pow(miny2 - Y0, 2)) + X0;
  
  if(Verbosity() > 3)
    std::cout << "minx1 and x2 : " << minx << ", " << minx2 << std::endl
	      << "miny1 and y2 : " << miny << ", " << miny2 << std::endl;

  /// Figure out which of the two roots is actually closer to the origin
  if(fabs(minx) < fabs(minx2))
    x = minx;
  else
    x = minx2;

  if(fabs(miny) < fabs(miny2))
    y = miny;
  else
    y = miny2;
  
  if(Verbosity() > 3)
    {
      std::cout << "Minimum x and y positions " << x << ",  " 
		<< y << std::endl;
    }

}

std::vector<TrkrCluster*> PHFindSecondarySiliconClusters::findMatchedSiliconClusters(std::vector<TrkrCluster*> all_si_clusters, double R, double X0, double Y0, double A, double B, SvtxTrack* track)
{
  std::vector<TrkrCluster*> matched_clusters;

  for(auto cluster : all_si_clusters)
    {
      TrkrDefs::cluskey cluskey = cluster->getClusKey();
      unsigned int layer = TrkrDefs::getLayer(cluskey);

      // xy residual
      double x = cluster->getX();
      double y = cluster->getY();
      // To ensure clusters are on same side as track?
      double clusphi = atan2(y - track->get_y(), x - track->get_x());
      if(fabs(clusphi - track->get_phi()) > M_PI / 2.0) continue;
      // The shortest distance of a point from a circle is along the radial line from the circle center to the point
      double xydca = sqrt( (x - X0)*(x-X0) + (y-Y0)*(y-Y0) )  -  R;  
      if(fabs(xydca) > _xy_residual_cut) continue;

      // line residual
      double z = cluster->getZ();
      double r = sqrt(x*x+y*y);
      double a = -A;
      double b = 1.0;
      double c = -B;
      double zdca = sqrt(pow(a*r+b*z+c, 2)) / sqrt(a*a+b*b);
      
      if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::mvtxId)
	if(fabs(zdca) > _z_residual_cut) continue;
            
      // if we are still here, this cluster passed the cuts
      matched_clusters.push_back(cluster);

      if(Verbosity() > 1)
	std::cout << "   found cluster in layer " << layer  << " with radius " << r << " xydca " << xydca << " zdca " << zdca << " x " << x << " y " << y  << " clusphi " << clusphi << " track phi " << track->get_phi() << std::endl;
    }

  return matched_clusters;  

}

std::vector<TrkrCluster*> PHFindSecondarySiliconClusters::getAllSiliconClusters()
{
  std::vector<TrkrCluster*> clusters;
  
  auto mvtx_hitsetrange = _hitsets->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (auto hitsetitr = mvtx_hitsetrange.first;
       hitsetitr != mvtx_hitsetrange.second;
       ++hitsetitr)
    {
      auto range = _cluster_map->getClusters(hitsetitr->first);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  const auto cluster = clusIter->second;
	  clusters.push_back(cluster);	  
	}
    }

  auto intt_hitsetrange = _hitsets->getHitSets(TrkrDefs::TrkrId::inttId);
  for (auto hitsetitr = intt_hitsetrange.first;
       hitsetitr != intt_hitsetrange.second;
       ++hitsetitr)
    {
      auto range = _cluster_map->getClusters(hitsetitr->first);
      for( auto clusIter = range.first; clusIter != range.second; ++clusIter )
	{
	  const auto cluster = clusIter->second;
	  clusters.push_back(cluster);	  
	}
    }
  
  return clusters;
}
