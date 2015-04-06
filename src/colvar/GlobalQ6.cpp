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

class GlobalQ6 : public GlobalSteinhardt {  
public:
  static void registerKeywords( Keywords& keys );
  GlobalQ6( const ActionOptions& ao );
};

  PLUMED_REGISTER_ACTION(GlobalQ6,"GlobalQ6")

  void GlobalQ6::registerKeywords( Keywords& keys ){
    GlobalSteinhardt::registerKeywords( keys );
    keys.addFlag("IMAG",false,"Compute the imaginary part of the Steinhardt parameter ");
    keys.add("compulsory","MVALUE","0","The value m in angular momentum ");
    keys.add("compulsory","NN","6","The n parameter of the switching function ");
    keys.add("compulsory","MM","12","The m parameter of the switching function ");
    keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
    keys.add("compulsory","R_0","The r_0 parameter of the switching function");
    keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
	     "The following provides information on the \\ref switchingfunction that are available. " 
	     "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 
  }

  GlobalQ6::GlobalQ6(const ActionOptions& ao ):
      Action(ao),
      GlobalSteinhardt(ao)
  {
    setAngularMomentum(6);

    std::vector<double> normaliz, coeff_poly;

    normaliz.resize( 7 );
    normaliz[0] = sqrt( ( 13.0*720.0 ) / (4.0*pi*720.0) );
    normaliz[1] = -sqrt( ( 13.0*120.0 ) / (4.0*pi*5040) );
    normaliz[2] = sqrt( ( 13.0*24) / (4.0*pi*40320) );
    normaliz[3] = -sqrt( ( 13.0*6) / (4.0*pi*362880) );
    normaliz[4] = sqrt( (13.0*2) / (4.0*pi*3628800) );
    normaliz[5] = -sqrt( (13.0*1) / (4.0*pi*39916800) );
    normaliz[6] = sqrt( (13.0*1) / (4.0*pi*479001600) );

    coeff_poly.resize( 7 ); 
    coeff_poly[0]=-0.3125; coeff_poly[1]=0.0;
    coeff_poly[2]=6.5625; coeff_poly[3]=0.0;
    coeff_poly[4]=-19.6875; coeff_poly[5]=0.0;
    coeff_poly[6]=14.4375; 

    setPolyCoeff (coeff_poly);
    setNormalization (normaliz);

    string sw,errors;
    parse("SWITCH",sw);
    if(sw.length()>0){
      switchingFunction.set(sw,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
    } else {
      int nn=6;
      int mm=12;
      double d0=0.0;
      double r0=0.0;
      parse("R_0",r0);
      if(r0<=0.0) error("R_0 should be explicitly specified and positive");
      parse("D_0",d0);
      parse("NN",nn);
      parse("MM",mm);
      switchingFunction.set(nn,mm,r0,d0);
    }
    parse("MVALUE", mvalue);

    parseFlag("IMAG",doImag);
  
    checkRead();
   
    log<<"  contacts are counted with cutoff "<<switchingFunction.description()<<"\n";
  }
}
}

