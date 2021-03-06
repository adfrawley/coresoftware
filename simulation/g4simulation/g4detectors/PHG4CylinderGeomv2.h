#ifndef PHG4CylinderGeomv2_H__
#define PHG4CylinderGeomv2_H__

#include "PHG4CylinderGeomv1.h"

class PHG4CylinderGeomv2: public PHG4CylinderGeomv1
{
 public:
  PHG4CylinderGeomv2();
  PHG4CylinderGeomv2(const double r, const double zmi, const double zma, const double thickn, const int n_scint):
    PHG4CylinderGeomv1(r,zmi,zma,thickn),
    nscint(n_scint)
      {}

  virtual ~PHG4CylinderGeomv2() {}

  void identify(std::ostream& os = std::cout) const;
  void set_nscint(const int i) {nscint = i;}
  int get_nscint() const {return nscint;}

 protected:
  int nscint;

  ClassDef(PHG4CylinderGeomv2,1)
};

#endif
