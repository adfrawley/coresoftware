#include "PHG4TruthSubsystem.h"
#include "PHG4TruthEventAction.h"
#include "PHG4TruthSteppingAction.h"
#include "PHG4TruthTrackingAction.h"


#include "PHG4TruthInfoContainer.h"

#include "PHG4VtxPointv1.h"
#include "PHG4Particlev2.h"

#include "PHG4InEvent.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4SystemOfUnits.hh>

using namespace std;

//_______________________________________________________________________
PHG4TruthSubsystem::PHG4TruthSubsystem( const string &name ):
  PHG4Subsystem( name ),
  eventAction_( 0 ),
  steppingAction_( 0 ),
  trackingAction_( 0 )
{}

//_______________________________________________________________________
int PHG4TruthSubsystem::InitRun( PHCompositeNode* topNode )
{

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create truth information container
  PHG4TruthInfoContainer* truthInfoList =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );
  if ( !truthInfoList )
    {
      dstNode->addNode( new PHIODataNode<PHObject>( truthInfoList = new PHG4TruthInfoContainer(), "G4TruthInfo", "PHObject" ));
    }

  // event action
  eventAction_ = new PHG4TruthEventAction();

  // create stepping action
  //steppingAction_ = new PHG4TruthSteppingAction( eventAction_ );

  // create tracking action
  trackingAction_ = new PHG4TruthTrackingAction( eventAction_ );

  return 0;

}

//_______________________________________________________________________
int
PHG4TruthSubsystem::process_event( PHCompositeNode* topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if ( eventAction_ )
    {
      eventAction_->SetInterfacePointers( topNode );
    }
  else
    {
      cout << PHWHERE << " No EventAction registered" << endl;
      exit(1);
    }

  if ( trackingAction_ )
    {
      trackingAction_->SetInterfacePointers( topNode );
    }
  else
    {
      cout << PHWHERE << " No TrackingAction registered" << endl;
      exit(1);
    }

  // this is called before G4 is kicked into gear. So we can fill in the information of the
  // G4 input particles from the InEvent Node.
  PHG4InEvent *inEvent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!inEvent) // if this node doesn't exist, they get the particles from somewhere else (not good)
    {
      cout << "Could not locate PHG4INEVENT node, where do you get your Geant 4 input from???" << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  PHG4TruthInfoContainer* truthInfoList =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );

//     cout << "truthInfoList identify" << endl;
//     truthInfoList->identify();
//     cout << "truthInfoList identify done" << endl;
//     cout << "truthInfoList maxkey: " << truthInfoList->maxindex() << endl;
//     cout << "truthInfoList minkey: " << truthInfoList->minindex() << endl;
  trackingAction_->TrackIdOffset(truthInfoList->maxtrkindex());
  eventAction_->TrackIdOffset(truthInfoList->maxtrkindex());
  map<int, PHG4VtxPoint *>::const_iterator vtxiter;
  multimap<int, PHG4Particle *>::const_iterator particle_iter;
  std::pair< std::map<int, PHG4VtxPoint *>::const_iterator, std::map<int, PHG4VtxPoint *>::const_iterator > vtxbegin_end = inEvent->GetVertices();
  for (vtxiter = vtxbegin_end.first; vtxiter != vtxbegin_end.second; vtxiter++)
    {
      PHG4VtxPoint *vtx = new PHG4VtxPointv1(vtxiter->second);
      int my_vtx_id = truthInfoList->AddPrimaryVertex(vtx);
	   
      pair<multimap<int, PHG4Particle *>::const_iterator, multimap<int, PHG4Particle *>::const_iterator > particlebegin_end = inEvent->GetParticles(vtxiter->first);
      for (particle_iter = particlebegin_end.first; particle_iter != particlebegin_end.second; particle_iter++)
	{
	  PHG4Particle *particle = new PHG4Particlev2(particle_iter->second);
	  particle->set_vtx_id(my_vtx_id);
	  G4ParticleTable* particleTable( G4ParticleTable::GetParticleTable() );
	  G4ParticleDefinition* g4particle( particleTable->FindParticle( particle->get_pid() ) );
	  if ( !g4particle && (particle->get_name()).find("geantino") == string::npos ) // keep geantinos
	    {
	      std::cout << PHWHERE << ": unable to find particle properties for "
			<< particle->get_name() << " with PDG id "
			<< particle->get_pid() << std::endl;
	    }
	  double mass = 0;
	  if (g4particle)
	    {
	      mass = ( g4particle->GetPDGMass() / GeV );
	    }
	  double energy = sqrt( pow(mass, 2) + pow(particle->get_px(), 2) + pow(particle->get_py(), 2) + pow(particle->get_pz(), 2) );
	  particle->set_e( energy );

	  PHG4TruthInfoContainer::ConstIterator piter = truthInfoList->AddPrimaryParticle(particle);
	  particle->set_track_id(piter->first); // update track id with returned id from insert
	}

    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4TruthSubsystem::ResetEvent(PHCompositeNode *topNode)
{
  trackingAction_->ResetEvent(topNode);
  eventAction_->ResetEvent(topNode);
  return 0;
}

//_______________________________________________________________________
PHG4EventAction* PHG4TruthSubsystem::GetEventAction( void ) const
{ return eventAction_; }


//_______________________________________________________________________
PHG4SteppingAction* PHG4TruthSubsystem::GetSteppingAction( void ) const
{ return steppingAction_; }

PHG4TrackingAction* 
PHG4TruthSubsystem::GetTrackingAction( void ) const
{ 
  return trackingAction_; 
}
