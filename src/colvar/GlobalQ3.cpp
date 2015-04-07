/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "GlobalSteinhardt.h"
#include "core/ActionRegister.h"
#include <string>

using namespace std;

namespace PLMD {
namespace colvar {

class GlobalQ3 : public GlobalSteinhardt {  
public:
  static void registerKeywords( Keywords& keys );
  GlobalQ3( const ActionOptions& ao );
};

  PLUMED_REGISTER_ACTION(GlobalQ3,"GlobalQ3")

  void GlobalQ3::registerKeywords( Keywords& keys ){
    GlobalSteinhardt::registerKeywords( keys );
  }

  GlobalQ3::GlobalQ3(const ActionOptions& ao ):
      Action(ao),
      GlobalSteinhardt(ao)
  {
    setAngularMomentum(3);

    std::vector<double> normaliz, coeff_poly;

    normaliz.resize( 4 );
    normaliz[0] = sqrt( ( 7.0*6.0 ) / (4.0*pi*6.0) );
    normaliz[1] = -sqrt( ( 7.0*2.0 ) / (4.0*pi*24.0) );
    normaliz[2] = sqrt( ( 7.0*1.0) / (4.0*pi*120.0) );
    normaliz[3] = -sqrt( ( 7.0*1.0) / (4.0*pi*720.0) );
  
    coeff_poly.resize( 4 ); 
    coeff_poly[0]=0.0; 
    coeff_poly[1]=-1.5;
    coeff_poly[2]=0.0; 
    coeff_poly[3]=2.5;

    setPolyCoeff (coeff_poly);
    setNormalization (normaliz);
  
    checkRead();
   
    log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
  }
}
}

