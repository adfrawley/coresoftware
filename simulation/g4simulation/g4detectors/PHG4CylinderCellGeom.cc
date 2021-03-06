#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellDefs.h"

#include <cmath>
#include <cstdlib>

ClassImp(PHG4CylinderCellGeom)

using namespace std;

PHG4CylinderCellGeom::PHG4CylinderCellGeom():
  layer(-9999),
  binning(0),
  radius(NAN),
  nzbins(-1),
  zmin(NAN),
  zstep(NAN),
  nphibins(-1),
  phimin(-M_PI),
  phistep(NAN),
  thickness(NAN)
{
  return;
}

void
PHG4CylinderCellGeom::set_zbins(const int i)
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  nzbins = i;
}

void
PHG4CylinderCellGeom::set_zmin(const double z)
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  zmin = z;
}

int
PHG4CylinderCellGeom::get_zbins() const
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  return nzbins;
}

double
PHG4CylinderCellGeom::get_zmin() const
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  return zmin;
}

double
PHG4CylinderCellGeom::get_zstep() const
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  return zstep;
}

void
PHG4CylinderCellGeom::set_zstep(const double z)
{
  check_binning_method(phg4cylindercelldefs::sizebinning);
  zstep = z;
}

int
PHG4CylinderCellGeom::get_phibins() const
{
  check_binning_method_phi("PHG4CylinderCellGeom::get_phibins");
return  nphibins;
}

double
PHG4CylinderCellGeom::get_phistep() const
{
  check_binning_method_phi("PHG4CylinderCellGeom::get_phistep");
  return  phistep;
}

double
PHG4CylinderCellGeom::get_phimin() const
{
  check_binning_method_phi("PHG4CylinderCellGeom::get_phimin");
  return  phimin;
}

void
PHG4CylinderCellGeom::set_phibins(const int i)
{
  check_binning_method_phi("PHG4CylinderCellGeom::set_phibins");
  nphibins = i;
}

void
PHG4CylinderCellGeom::set_phistep(const double r)
{
  check_binning_method_phi("PHG4CylinderCellGeom::set_phistep");
  phistep = r;
}

void
PHG4CylinderCellGeom::set_phimin(const double r)
{
  check_binning_method_phi("PHG4CylinderCellGeom::set_phimin");
  phimin = r;
}

int
PHG4CylinderCellGeom::get_etabins() const
{
  check_binning_method_eta("PHG4CylinderCellGeom::get_etabins");
  return nzbins;
}

double
PHG4CylinderCellGeom::get_etastep() const
{
  check_binning_method_eta("PHG4CylinderCellGeom::get_etastep");
  return zstep;
}
double
PHG4CylinderCellGeom::get_etamin() const
{
  check_binning_method_eta("PHG4CylinderCellGeom::get_etamin");
  return zmin;
}

void
PHG4CylinderCellGeom::set_etamin(const double z)
{
  check_binning_method_eta("PHG4CylinderCellGeom::set_etamin");
  zmin = z;
}

void
PHG4CylinderCellGeom::set_etastep(const double z)
{
  check_binning_method_eta("PHG4CylinderCellGeom::set_etastep");
  zstep = z;
}

void
PHG4CylinderCellGeom::set_etabins(const int i)
{
  check_binning_method_eta("PHG4CylinderCellGeom::set_etabins");
  nzbins = i;
}

void
PHG4CylinderCellGeom::identify(std::ostream& os) const
{
  os << "layer: " << layer
     << ", radius: " << radius
     << ", thickness: " << thickness;
  switch (binning)
    {
    case phg4cylindercelldefs::sizebinning:
      os << ", zbins: " << nzbins
	 << ", zmin: " << zmin
	 << ", zstepsize: " << zstep;
      break;
    case phg4cylindercelldefs::etaphibinning:
      os << ", etabins: " << nzbins
	 << ", etamin: " << zmin
	 << ", etastepsize: " << zstep;
      break;
    case phg4cylindercelldefs::etaslatbinning:
      os << ", etabins: " << nzbins
	 << ", etamin: " << zmin
	 << ", etastepsize: " << zstep;
      break;
    default:
      os << "no valid binning method: " << binning << endl;
      return;
      break;
    }
  os << ", phimin: " << phimin
     << ", phibins: " << nphibins
     << ", phistep: " << phistep
     << endl;
  return;
}

pair<double, double>
PHG4CylinderCellGeom::get_zbounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method(phg4cylindercelldefs::sizebinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return make_pair(zlow, zhigh);
}

pair<double, double>
PHG4CylinderCellGeom::get_etabounds(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method_eta("PHG4CylinderCellGeom::get_etabounds");
//  check_binning_method(phg4cylindercelldefs::etaphibinning);
  double zlow = zmin + ibin * zstep;
  double zhigh = zlow + zstep;
  return make_pair(zlow, zhigh);
}


pair<double, double>
PHG4CylinderCellGeom::get_phibounds(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  double philow = phimin + ibin * phistep;
  double phihigh = philow + phistep;
  return make_pair(philow, phihigh);
}

int
PHG4CylinderCellGeom::get_zbin(const double z) const
{
  if (z < zmin || z > (zmin+nzbins*zstep))
  {
    //    cout << "Asking for bin for z outside of z range: " << z << endl;
    return -1;
  }
  
  check_binning_method(phg4cylindercelldefs::sizebinning);
  return floor( (z-zmin)/zstep );
}

int
PHG4CylinderCellGeom::get_etabin(const double eta) const
{
  if (eta < zmin || eta > (zmin+nzbins*zstep))
  {
    //    cout << "Asking for bin for eta outside of eta range: " << eta << endl;
    return -1;
  }
  check_binning_method_eta();
  return floor( (eta-zmin)/zstep );
}

int
PHG4CylinderCellGeom::get_phibin(const double phi) const
{
  double norm_phi = phi;
  if(phi < phimin || phi > (phimin+nphibins*phistep))
  {
    int nwraparound = -floor((phi-phimin) * 0.5 / M_PI);
    norm_phi += 2*M_PI*nwraparound;
  }
  check_binning_method_phi();
  return floor( (norm_phi-phimin)/phistep );
}

double
PHG4CylinderCellGeom::get_zcenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in z: " << ibin << endl;
      exit(1);
    }
  check_binning_method(phg4cylindercelldefs::sizebinning);
  return zmin + (ibin + 0.5)*zstep;
}

double
PHG4CylinderCellGeom::get_etacenter(const int ibin) const
{
  if (ibin < 0 || ibin > nzbins)
    {
      cout << "Asking for invalid bin in eta: " << ibin << endl;
      cout << "minbin: 0, maxbin " << nzbins << endl; 
      exit(1);
    }
  check_binning_method_eta();
  return zmin + (ibin + 0.5)*zstep;
}

double
PHG4CylinderCellGeom::get_phicenter(const int ibin) const
{
  if (ibin < 0 || ibin > nphibins)
    {
      cout << "Asking for invalid bin in phi: " << ibin << endl;
      exit(1);
    }

  check_binning_method_phi();
  return (phimin + (ibin + 0.5)*phistep);
}

string
PHG4CylinderCellGeom::methodname(const int i) const
{
  switch (i)
    {
    case phg4cylindercelldefs::sizebinning:
      return "Bins in cm";
      break;
    case phg4cylindercelldefs::etaphibinning:
      return "Eta/Phi bins";
      break;
    case phg4cylindercelldefs::etaslatbinning:
      return "Eta/numslat bins";
      break;
    default:
      break;
    }
  return "Unknown";
}

void
PHG4CylinderCellGeom::check_binning_method(const int i) const
{
  if (binning != i)
    {
      cout << "different binning method used " << methodname(binning)
           << ", not : " << methodname(i)
           << endl;
      exit(1);
    }
  return;
}

void
PHG4CylinderCellGeom::check_binning_method_eta(const std::string & src) const
{
  if (binning != phg4cylindercelldefs::etaphibinning && 
      binning != phg4cylindercelldefs::etaslatbinning)
    {
      if (src.size())
        cout << src<<" : ";

      cout << "different binning method used " << methodname(binning)
           << ", not : " << methodname(phg4cylindercelldefs::etaphibinning)
	   << " or " << methodname(phg4cylindercelldefs::etaslatbinning)
           << endl;
      exit(1);
    }
  return;
}

void
PHG4CylinderCellGeom::check_binning_method_phi(const std::string & src) const
{
  if (binning != phg4cylindercelldefs::etaphibinning && 
      binning != phg4cylindercelldefs::sizebinning &&
      binning != phg4cylindercelldefs::etaslatbinning)
    {
      if (src.size())
        cout << src<<" : ";

      cout << "different binning method used " << methodname(binning)
           << ", not : " << methodname(phg4cylindercelldefs::etaphibinning)
	   << " or " << methodname(phg4cylindercelldefs::sizebinning)
           << endl;
      exit(1);
    }
  return;
}
