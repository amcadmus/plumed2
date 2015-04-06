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

class GlobalQ4 : public GlobalSteinhardt {  
public:
  static void registerKeywords( Keywords& keys );
  GlobalQ4( const ActionOptions& ao );
};

  PLUMED_REGISTER_ACTION(GlobalQ4,"GlobalQ4")

  void GlobalQ4::registerKeywords( Keywords& keys ){
    GlobalSteinhardt::registerKeywords( keys );
    keys.addFlag("IMAG",false,"Compute the imaginary part of the Steinhardt parameter ");
    keys.add("compulsory","MVALUE","0","The value m in angular momentum ");
  }

  GlobalQ4::GlobalQ4(const ActionOptions& ao ):
      Action(ao),
      GlobalSteinhardt(ao)
  {
    setAngularMomentum(4);

    std::vector<double> normaliz, coeff_poly;

    normaliz.resize( 5 );
    normaliz[0] = sqrt( ( 9.0*24.0 ) / (4.0*pi*24.0) );
    normaliz[1] = -sqrt( ( 9.0*6.0 ) / (4.0*pi*120.0) );
    normaliz[2] = sqrt( ( 9.0*2.0) / (4.0*pi*720.0) );
    normaliz[3] = -sqrt( ( 9.0*1) / (4.0*pi*5040.0) );
    normaliz[4] = sqrt( (9.0*1) / (4.0*pi*40320.0) );

    coeff_poly.resize( 5 ); 
    coeff_poly[0]=0.375; coeff_poly[1]=0.0;
    coeff_poly[2]=-3.75; coeff_poly[3]=0.0;
    coeff_poly[4]=4.375;

    setPolyCoeff (coeff_poly);
    setNormalization (normaliz);

    parse("MVALUE", mvalue);

    parseFlag("IMAG",doImag);
  
    checkRead();
   
    log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
  }
}
}

