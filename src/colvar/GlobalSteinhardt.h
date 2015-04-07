/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_colvar_GlobalSteinhardt_h
#define __PLUMED_colvar_GlobalSteinhardt_h
#include "Colvar.h"
#include "tools/SwitchingFunction.h"

namespace PLMD{

class NeighborList;

namespace colvar{

class GlobalSteinhardt : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  bool invalidateList;
  bool firsttime;
  unsigned tmom;
  std::vector<double> coeff_poly;
  std::vector<double> normaliz;  
  double deriv_poly( const unsigned& m, const double& val, double& df );
protected:
  SwitchingFunction switchingFunction;
  void setAngularMomentum( const unsigned& ang ) {tmom = ang;}
  void setPolyCoeff (const std::vector<double > & coeff) {coeff_poly = coeff;}
  void setNormalization (const std::vector<double > & normaliz_) {normaliz = normaliz_;}
public:
  GlobalSteinhardt(const ActionOptions&);
  ~GlobalSteinhardt();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};

}
}
#endif
