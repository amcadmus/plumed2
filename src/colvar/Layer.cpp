/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
// #include "Layer.h"
#include "Colvar.h"
#include "ActionRegister.h"
// #include "tools/NeighborList.h"
// #include "tools/Communicator.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

class Layer : public Colvar {
  bool serial;
  double nn;
  unsigned dir;
  bool comp_shift;
  double shift;
public:
  Layer(const ActionOptions&);
  ~Layer();
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Layer,"LAYER")

void Layer::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("COMP_SHIFT",false,"Compute the shift of the cos function automatically");
  keys.add("atoms","ATOMS","The group of atoms that you are calculating the Layer parameter");
  keys.add("compulsory","NUMB_LAYER","The number of layer in the system");
  keys.add("compulsory","DIRECTION","The direction along which the layer forms");
  keys.add("optional","SHIFT", "The shift angle of the cos function");
}

  Layer::Layer(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      serial(false),
      nn(0.),
      dir(0),
      comp_shift (false),
      shift (0)
  {

    parseFlag("SERIAL",serial);

    std::vector<AtomNumber> atoms;
    parseAtomList("ATOMS",atoms);

    parse("NUMB_LAYER",nn);
    string direction;
    parse("DIRECTION",direction);
    if (direction == "X"){
      dir = 0;
    }
    else if (direction == "Y"){
      dir = 1;
    }
    else if (direction == "Z"){
      dir = 2;
    }
    else {
      error ("unrecongnized direction input!");
    }
    string boolstring;
    parse("COMP_SHIFT",boolstring);
    if (boolstring == "TRUE" || boolstring == "true"){
      comp_shift = true;
    }
    else if (boolstring == "FALSE" || boolstring == "false") {
      comp_shift = false;
    }
    else {
      error (" unrecongnized bool input at COMP_SHIFT!");
    }
    if (!comp_shift) {
      parse("SHIFT",shift);
    }
    else {
      shift = 0.;
    }

    checkRead();

    addValueWithDerivatives(); setNotPeriodic();  
    requestAtoms(atoms);
 
    log.printf("  consider groups of %u atoms\n",static_cast<unsigned>(atoms.size()));
    log.printf("  first group:\n");
    for(unsigned int i=0;i<atoms.size();++i){
      if ( (i+1) % 25 == 0 ) log.printf("  \n");
      log.printf("  %d", atoms[i].serial());
    }
    log.printf("  \n");
    if (comp_shift){
      log.printf("  auto shift comp: yes\n");
    }
    else {
      log.printf("  auto shift comp: no, using provided shift value of %f\n", shift);
    }
  }

  Layer::~Layer(){
  }
  
// calculator
  void Layer::calculate()
  {

    // unsigned stride=comm.Get_size();
    // unsigned rank=comm.Get_rank();
    // if(serial){
    //   stride=1;
    //   rank=0;
    // }else{
    //   stride=comm.Get_size();
    //   rank=comm.Get_rank();
    // }

    double ninv = 1./double(getNumberOfAtoms());
    double Li = 1./(getBox()[dir][dir]);

    double valueA = 0;
    double valueB = 0;
    double twoPiNLi = 2. * M_PI * nn * Li;
    if (comp_shift){
      for (unsigned ii = 0; ii < getNumberOfAtoms(); ++ii){
	Vector posi = getPosition(ii) ;
	valueA += cos (twoPiNLi * posi[dir]);
	valueB += sin (twoPiNLi * posi[dir]);
      }
      if (fabs(valueA) < 1e-12){
	shift = M_PI * 0.5;
      }
      else {
	shift = atan (-valueB / valueA);
      }
    }
 
    double value0=0.;
    double value1=0.;
    // double check = 0;
    vector<Vector> deriv(getNumberOfAtoms());
    // vector<Vector> deriv0(getNumberOfAtoms());
    // vector<Vector> deriv1(getNumberOfAtoms());
    Tensor virial;
    // cout << "used shift is " << shift <<endl;

    for(unsigned int i = 0; i < getNumberOfAtoms(); i++) {
      Vector posi = getPosition(i) ;
      value0 += ninv * cos(twoPiNLi * posi[dir] + shift);
      value1 += ninv * cos(twoPiNLi * posi[dir] + shift + M_PI);
    }
    double value = 0.;
    bool used0;
    if (fabs(value0) > fabs(value1)){
      used0 = true;  
      value =  value0;
    }
    else {
      used0 = false;
      value = value1;
    }
    for(unsigned int i = 0; i < getNumberOfAtoms(); i++) {
      Vector posi = getPosition(i) ;
      if (used0){
	deriv[i][dir] += value * -ninv * twoPiNLi * sin(twoPiNLi * posi[dir] + shift);
      }
      else {
	deriv[i][dir] += value * -ninv * twoPiNLi * sin(twoPiNLi * posi[dir] + shift + M_PI);
      }
    }
    value = 0.5 * value * value;
    
    setValue (value);
    for(unsigned i = 0; i < deriv.size(); ++ i) {
      setAtomsDerivatives (i, deriv[i]);
    }

    
    // for(unsigned int i = 0; i < getNumberOfAtoms(); i++) {
    //   Vector posi = getPosition(i) ;
    //   // check += sin(twoPiNLi * posi[dir] + shift);
    //   value0 += ninv * (0.5 * cos(twoPiNLi * posi[dir] + shift) + .5);
    //   deriv0[i][dir] = deriv0[i][dir] + -ninv * 0.5 * twoPiNLi * sin(twoPiNLi * posi[dir] + shift);
    //   if (comp_shift){
    // 	value1 += ninv * (0.5 * cos(twoPiNLi * posi[dir] + shift + M_PI) + .5);
    // 	deriv1[i][dir] = deriv1[i][dir] + -ninv * 0.5 * twoPiNLi * sin(twoPiNLi * posi[dir] + shift + M_PI);
    //   }
    // }

    // if (comp_shift){
    //   if (value0 > value1){
    // 	setValue           (value0);
    // 	for(unsigned i = 0; i < deriv0.size(); ++ i) {
    // 	  setAtomsDerivatives(i,deriv0[i]);
    // 	}
    //   }
    //   else {
    // 	setValue           (value1);
    // 	for(unsigned i = 0; i < deriv1.size(); ++ i) {
    // 	  setAtomsDerivatives(i,deriv1[i]);
    // 	}
    //   }
    // }
    // else {
    //   setValue           (value0);
    //   for(unsigned i = 0; i < deriv0.size(); ++ i) {
    // 	setAtomsDerivatives(i,deriv0[i]);
    //   }      
    // }
  }
}
}
