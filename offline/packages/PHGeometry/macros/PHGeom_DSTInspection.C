// $Id: $

/*!
 * \file PHGeom_DSTInspection.C
 * \brief Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <phgeom/PHGeomUtility.h>
#include <phgeom/PHGeomIOTGeo.h>
#include <phgeom/PHGeomTGeo.h>
#include <phool/recoConsts.h>

#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TGLClip.h>
#include <TGLUtil.h>
#include <TGLViewer.h>
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TROOT.h>

#include <cassert>
#include <string>
#include <vector>

using namespace std;

R__LOAD_LIBRARY(libphgeom.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libfun4all.so)

// This function edit and builds TGeoManager from the DST's Geometry IO node.
// It need to be called before any module using TGeo, i.e. before modules that calls PHGeomUtility::GetTGeoManager(se->topNode())
// i.e. after PHG4Reco and before any Kalman fitter calls

void AddActsTpcSurfaces(  TGeoManager *geoManager, TGeoVolume *tpc_gas_vol, int verbosity = 1)
{
  // TPC surface subdivisions
  double m_minSurfZ = 0.0;
  double m_maxSurfZ = 105.5;
  unsigned int m_nSurfZ = 11;
  unsigned int m_nSurfPhi = 10;
  double m_surfStepPhi;
  double m_surfStepZ;
  double m_moduleStepPhi;
  double m_modulePhiStart;

  // these don't change, we are building the tpc this way!
  const unsigned int m_nTpcLayers = 48;
  const unsigned int m_nTpcModulesPerLayer = 12;
  const double m_minRadius[3] = {30.0, 40.0, 60.0};
  const double m_maxRadius[3] = {40.0, 60.0, 77.0};
  double layer_thickness_sector[3] = {0};
  double layer_radius[48]= {0};
  double layer_thickness[48] = {0};

  // These are arbitrary subdivisions, and may change
  m_surfStepZ = (m_maxSurfZ - m_minSurfZ) / (double) m_nSurfZ;
  m_moduleStepPhi = 2.0 * M_PI / 12.0;
  m_modulePhiStart = -M_PI;
  m_surfStepPhi = 2.0 * M_PI / (double) (m_nSurfPhi * m_nTpcModulesPerLayer);

  // prevent boxes from touching when placed
  double half_width_clearance_thick = 0.4999;
  double half_width_clearance_phi = 0.4999;
  double half_width_clearance_z = 0.4999;
 
  for(unsigned int isector = 0; isector < 3; ++isector)
    {
      layer_thickness_sector[isector] = (m_maxRadius[isector] - m_minRadius[isector]) / 16.0;

      for(unsigned int ilayer =0; ilayer < 16; ++ilayer)
	{
	  layer_radius[isector*16 + ilayer] = m_minRadius[isector] + layer_thickness_sector[isector]*(double) ilayer + layer_thickness_sector[isector] / 2.0;
	  layer_thickness[isector*16 + ilayer] = layer_thickness_sector[isector];
	}
    }
  
  //  TGeoTube *tpc_gas_shape = dynamic_cast<TGeoTube *>(tpc_gas_vol->GetShape());
  // assert(tpc_gas_shape);
  
  TGeoMedium *tpc_gas_medium = tpc_gas_vol->GetMedium();
  assert(tpc_gas_medium);

  TGeoVolume *tpc_gas_measurement_vol[48];
  double tan_half_phi = tan(m_surfStepPhi / 2.0);
  for(unsigned int ilayer = 0; ilayer < m_nTpcLayers; ++ilayer)
    {
      // make a box for this layer
      char bname[500];
      sprintf(bname,"tpc_gas_measurement_%i",ilayer);

      // Because we use a box, not a section of a cylinder, we need this to prevent overlaps
      // set the nominal r*phi dimension of the box so they just touch at the inner edge when placed 
      double box_r_phi = 2.0 * tan_half_phi * (layer_radius[ilayer] - layer_thickness[ilayer] / 2.0);

      tpc_gas_measurement_vol[ilayer] = geoManager->MakeBox(bname, tpc_gas_medium, 
							    layer_thickness[ilayer]*half_width_clearance_thick, 
							    box_r_phi*half_width_clearance_phi, 
							    m_surfStepZ*half_width_clearance_z);

      tpc_gas_measurement_vol[ilayer]->SetLineColor(kBlack);
      tpc_gas_measurement_vol[ilayer]->SetFillColor(kYellow);
      tpc_gas_measurement_vol[ilayer]->SetVisibility(kTRUE);
      cout << "Made box for layer " << ilayer << " with dx " << layer_thickness[ilayer] << " dy " 
	   << box_r_phi << " ref arc " << m_surfStepPhi*layer_radius[ilayer] << " dz " << m_surfStepZ << endl;

      tpc_gas_measurement_vol[ilayer]->Print();
    }

  int copy = 0;	      
  for (unsigned int iz = 0; iz < m_nSurfZ; ++iz)
    {
      // The (half) tpc gas volume is 105.5 cm long and is symmetric around (x,y,z) = (0,0,0) in its frame
      double z_center = -105.5/2.0 + m_surfStepZ / 2.0 + (double) iz * m_surfStepZ;
      
      for (unsigned int imod = 0; imod < m_nTpcModulesPerLayer; ++imod)
	{
	  for (unsigned int iphi = 0; iphi < m_nSurfPhi; ++iphi)
	    {

	      double min_phi = m_modulePhiStart + (double) imod * m_moduleStepPhi + (double) iphi * m_surfStepPhi;
	      double phi_center = min_phi + m_surfStepPhi / 2.0;
	      double phi_center_degrees = phi_center * 180.0 / 3.14159;
	      
	      for (unsigned int ilayer = 0; ilayer < m_nTpcLayers; ++ilayer)
		{
		  copy++;
		  
		  // place copies of the gas volume to fill up the layer
		  
		  double x_center = layer_radius[ilayer] * cos(phi_center);
		  double y_center = layer_radius[ilayer] * sin(phi_center);
		  
		  char rot_name[500];
		  sprintf(rot_name,"tpc_gas_rotation_%i", copy);
		  TGeoCombiTrans *tpc_gas_measurement_location = new TGeoCombiTrans(x_center, y_center, z_center,
										    new TGeoRotation(rot_name,phi_center_degrees, 0, 0));
		  
		  tpc_gas_vol->AddNode(tpc_gas_measurement_vol[ilayer], copy, tpc_gas_measurement_location);
		  
		  if(verbosity && ilayer == 30) 
		    {
		      cout << " Made copy " << copy << " iz " << iz << " imod " << imod << " ilayer " << ilayer << " iphi " << iphi << endl;
		      cout << "    x_center " << x_center << " y_center " << y_center << " z_center " << z_center << " phi_center_degrees " << phi_center_degrees << endl;
		    }
		}
	    }
	}
      
    }
}


void EditTPCGeometry(PHCompositeNode *topNode, const int verbosity = 1)
{
  // clean up existing TGeoManager before building a new one
  PHGeomTGeo *geom_node =  PHGeomUtility::GetGeomTGeoNode(topNode, true);
  assert(geom_node);
  if (geom_node)
  {
    if (geom_node->isValid())
    {
      geom_node->Reset();
    }
  }

  // read back geometry again from the DST persistent IO node
  PHGeomIOTGeo *dst_geom_io = PHGeomUtility::GetGeomIOTGeoNode(topNode, false);
  if (not dst_geom_io)
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomTGeo node RUN/GEOMETRY due to missing PHGeomIOTGeo node at RUN/GEOMETRY_IO"
         << endl;
    exit(1) ;
  }
  if (not dst_geom_io->isValid())
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - failed to update PHGeomTGeo node RUN/GEOMETRY due to invalid PHGeomIOTGeo node at RUN/GEOMETRY_IO"
         << endl;
    exit(1) ;
  }

  // build new TGeoManager
  TGeoManager *geoManager = dst_geom_io->ConstructTGeoManager();
  geom_node->SetGeometry(geoManager);
  assert(geoManager);

  TGeoVolume *World_vol = geoManager->GetTopVolume();
  TGeoNode *tpc_envelope_node = nullptr;
  TGeoNode *tpc_gas_north_node = nullptr;

  // find tpc north gas volume at path of World*/tpc_envelope*
  if (verbosity)
  {
    cout << "EditTPCGeometry - searching under volume: ";
    World_vol->Print();
  }
  for (int i = 0; i < World_vol->GetNdaughters(); i++)
  {
    TString node_name = World_vol->GetNode(i)->GetName();

    if (verbosity)
      cout << "EditTPCGeometry - searching node " << node_name << endl;

    if (node_name.BeginsWith("tpc_envelope"))
    {
      if (verbosity)
        cout << "EditTPCGeometry - found! " << endl;

      tpc_envelope_node = World_vol->GetNode(i);
      break;
    }
  }
  assert(tpc_envelope_node);

  // find tpc north gas volume at path of World*/tpc_envelope*/tpc_gas_north*
  TGeoVolume *tpc_envelope_vol = tpc_envelope_node->GetVolume();
  assert(tpc_envelope_vol);
  if (verbosity)
  {
    cout << "EditTPCGeometry - searching under volume: ";
    tpc_envelope_vol->Print();
  }
  for (int i = 0; i < tpc_envelope_vol->GetNdaughters(); i++)
  {
    TString node_name = tpc_envelope_vol->GetNode(i)->GetName();

    if (verbosity)
      cout << "EditTPCGeometry - searching node " << node_name << endl;

    if (node_name.BeginsWith("tpc_gas_north"))
    {
      if (verbosity)
        cout << "EditTPCGeometry - found! " << endl;

      tpc_gas_north_node = tpc_envelope_vol->GetNode(i);
      break;
    }
  }
  assert(tpc_gas_north_node);
  TGeoVolume *tpc_gas_north_vol = tpc_gas_north_node->GetVolume();
  assert(tpc_gas_north_vol);

  if (verbosity)
  {
    cout << "EditTPCGeometry - gas volume: ";
    tpc_gas_north_vol->Print();
  }

  // north = side 0?
  AddActsTpcSurfaces(geoManager, tpc_gas_north_vol, verbosity);

  // Closing geometry implies checking the geometry validity,
  // fixing shapes with negative parameters (run-time shapes)building the cache manager, voxelizing all volumes,
  // counting the total number of physical nodes and registering the manager class to the browser.
  geoManager->CloseGeometry();

  // save the edited geometry to DST persistent IO node for downstream DST files
  PHGeomUtility::UpdateIONode(topNode);
}


//! Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
//! Based on abhisek's display macro
void PHGeom_DSTInspection(string DST_file_name = "/phenix/u/jinhuang/links/sPHENIX_work/G4sPHENIX.root",
                          bool do_clip = true)
{
  TEveManager::Create();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 12345);

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->fileopen(DST_file_name);
  se->registerInputManager(hitsin);

  // run one event as example
  se->run(1);


  EditTPCGeometry(se->topNode());


//  PHGeomUtility::GetTGeoManager(se->topNode());
  assert(gGeoManager);

  if (!gROOT->GetListOfGeometries()->FindObject(gGeoManager))
    gROOT->GetListOfGeometries()->Add(gGeoManager);
  if (!gROOT->GetListOfBrowsables()->FindObject(gGeoManager))
    gROOT->GetListOfBrowsables()->Add(gGeoManager);
  //  gGeoManager->UpdateElements();

  TGeoNode *current = gGeoManager->GetCurrentNode();
  //Alternate drawing
  //current->GetVolume()->Draw("ogl");
  //Print the list of daughters
  //current->PrintCandidates();
  for (int igeom = 0; igeom < current->GetNdaughters(); igeom++)
  {
    // hide black holes
    TGeoNode *geo_node = (TGeoNode *) current->GetNodes()->UncheckedAt(igeom);
    if (string(geo_node->GetName()).find("BH_") != string::npos)
      geo_node->GetVolume()->SetVisibility(false);
  }
  
  TEveGeoTopNode *eve_node = new TEveGeoTopNode(gGeoManager, current);
  // eve_node->SetVisLevel(6);
  gEve->AddGlobalElement(eve_node);
  gEve->FullRedraw3D(kTRUE);

  // EClipType not exported to CINT (see TGLUtil.h):
  // 0 - no clip, 1 - clip plane, 2 - clip box
  TGLViewer *v = gEve->GetDefaultGLViewer();
//  if (do_clip)
//  {
//    v->GetClipSet()->SetClipType(TGLClip::kClipPlane);
//  }z
  //  v->ColorSet().Background().SetColor(kMagenta + 4);
  v->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);

  v->CurrentCamera().RotateRad(3.14/2, 0.);

  v->DoDraw();
  

  PHGeomUtility::ExportGeomtry(se->topNode(),"TGeo_export.root");
  cout <<"Done export Geometry to TGeo_export.root"<<endl;
}
