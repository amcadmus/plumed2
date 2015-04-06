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
#include "GlobalSteinhardt.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"

#include <string>
#include <complex>
#include <iostream>

using namespace std;

namespace PLMD{
  namespace colvar{

    void GlobalSteinhardt::registerKeywords( Keywords& keys ){
      Colvar::registerKeywords(keys);
      keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
      keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
      keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
      keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
      keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
      keys.add("atoms","GROUPA","First list of atoms");
      keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
    }

    GlobalSteinhardt::GlobalSteinhardt(const ActionOptions&ao):
	PLUMED_COLVAR_INIT(ao),
	pbc(true),
	serial(false),
	invalidateList(true),
	firsttime(true)
    {

      parseFlag("SERIAL",serial);

      vector<AtomNumber> ga_lista,gb_lista;
      parseAtomList("GROUPA",ga_lista);
      parseAtomList("GROUPB",gb_lista);

      bool nopbc=!pbc;
      parseFlag("NOPBC",nopbc);
      pbc=!nopbc;

// pair stuff
      bool dopair=false;
      parseFlag("PAIR",dopair);

// neighbor list stuff
      bool doneigh=false;
      double nl_cut=0.0;
      int nl_st=0;
      parseFlag("NLIST",doneigh);
      if(doneigh){
	parse("NL_CUTOFF",nl_cut);
	if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
	parse("NL_STRIDE",nl_st);
	if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
      }
  
      addValueWithDerivatives(); setNotPeriodic();
      if(gb_lista.size()>0){
	if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
	else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
      } else {
	if(doneigh)  nl= new NeighborList(ga_lista,pbc,getPbc(),nl_cut,nl_st);
	else         nl= new NeighborList(ga_lista,pbc,getPbc());
      }
  
      requestAtoms(nl->getFullAtomList());
 
      log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
      log.printf("  first group:\n");
      for(unsigned int i=0;i<ga_lista.size();++i){
	if ( (i+1) % 25 == 0 ) log.printf("  \n");
	log.printf("  %d", ga_lista[i].serial());
      }
      log.printf("  \n  second group:\n");
      for(unsigned int i=0;i<gb_lista.size();++i){
	if ( (i+1) % 25 == 0 ) log.printf("  \n");
	log.printf("  %d", gb_lista[i].serial());
      }
      log.printf("  \n");
      if(pbc) log.printf("  using periodic boundary conditions\n");
      else    log.printf("  without periodic boundary conditions\n");
      if(dopair) log.printf("  with PAIR option\n");
      if(doneigh){
	log.printf("  using neighbor lists with\n");
	log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
      }

    }

    GlobalSteinhardt::~GlobalSteinhardt(){
      delete nl;
    }

    void GlobalSteinhardt::prepare(){
      if(nl->getStride()>0){
	if(firsttime || (getStep()%nl->getStride()==0)){
	  requestAtoms(nl->getFullAtomList());
	  invalidateList=true;
	  firsttime=false;
	}else{
	  requestAtoms(nl->getReducedAtomList());
	  invalidateList=false;
	  if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
	}
	if(getExchangeStep()) firsttime=true;
      }
    }

    inline double GlobalSteinhardt::
    deriv_poly( const unsigned& m, const double& val, double& df ){
      double fact=1.0;
      for(unsigned j=1;j<=m;++j) fact=fact*j;
      double res=coeff_poly[m]*fact;

      double pow=1.0, xi=val, dxi=1.0; df=0.0;
      for(unsigned i=m+1;i<=tmom;++i){
	double fact=1.0;
	for(unsigned j=i-m+1;j<=i;++j) fact=fact*j;
	res=res+coeff_poly[i]*fact*xi;
	df = df + pow*coeff_poly[i]*fact*dxi;
	xi=xi*val; dxi=dxi*val; pow+=1.0;
      }
      df = df*normaliz[m];
      return normaliz[m]*res;
    }

// calculator
    void GlobalSteinhardt::calculate()
    {

      double ncoord=0.;
      Tensor virial;
      vector<Vector> deriv(getNumberOfAtoms());
      double steinhardtPrefactor = 1./(6.*getNumberOfAtoms());
      // double steinhardtPrefactor = 1.;
      cout << steinhardtPrefactor << endl;

      if(nl->getStride()>0 && invalidateList){
	nl->update(getPositions());
      }

      unsigned stride=comm.Get_size();
      unsigned rank=comm.Get_rank();
      if(serial){
	stride=1;
	rank=0;
      }else{
	stride=comm.Get_size();
	rank=comm.Get_rank();
      }

      for(unsigned int i=rank;i<nl->size();i+=stride) {                   // sum over close pairs
 
	Vector distance;
	unsigned i0=nl->getClosePair(i).first;
	unsigned i1=nl->getClosePair(i).second;

	if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;

	if(pbc){
	  distance=pbcDistance(getPosition(i0),getPosition(i1));
	} else {
	  distance=delta(getPosition(i0),getPosition(i1));
	}

	// calculates the Steinhardt parameter at value of m
	{
	  double dfunc, dpoly_ass, md, tq6, itq6, real_z, imag_z; 
	  Vector dz, myrealvec, myimagvec, real_dz, imag_dz;
	  // The square root of -1
	  std::complex<double> ii( 0.0, 1.0 ), dp_x, dp_y, dp_z;
	  double sw, poly_ass, dlen;
	  std::complex<double> powered;
	  dlen=distance.modulo(); 
	  sw = switchingFunction.calculate( dlen, dfunc );
	  if( sw>=1e-12 ){
	    double dlen3 = dlen*dlen*dlen;
	    // Derivatives of z/r wrt x, y, z
	    dz = -( distance[2] / dlen3 )*distance;
	    dz[2] += (1.0 / dlen);
	    if (mvalue == 0){
	      // Do stuff for m=0
	      poly_ass=deriv_poly( 0, distance[2]/dlen, dpoly_ass );
	      // Derivative wrt to the vector connecting the two atoms
	      myrealvec = (+sw)*dpoly_ass*dz + poly_ass*(+dfunc)*distance;
	      // Accumulate the derivatives
	      deriv[i0] = deriv[i0] + (-myrealvec);
	      deriv[i1] = deriv[i1] + (myrealvec);
	      virial = virial + Tensor( -myrealvec,distance );
	      // And store the vector function
	      ncoord += sw*poly_ass;	    
	    }
	    else {
	      bool posiFlag = true;
	      int m = mvalue;
	      if (mvalue < 0) {
		m = -mvalue;
		posiFlag = false;
	      }
	      std::complex<double> com1( distance[0]/dlen ,distance[1]/dlen );
	      // Calculate Legendre Polynomial
	      poly_ass=deriv_poly( m, distance[2]/dlen, dpoly_ass );
	      // Calculate powe of complex number
	      powered=pow(com1,m-1); md=static_cast<double>(m);
	      // Real and imaginary parts of z
	      real_z = real(com1*powered); imag_z = imag(com1*powered );
 
	      // Calculate steinhardt parameter
	      tq6=poly_ass*real_z;   // Real part of steinhardt parameter
	      itq6=poly_ass*imag_z;  // Imaginary part of steinhardt parameter

	      // Derivatives wrt ( x/r + iy )^m
	      dp_x = md*powered*( (1.0/dlen)-(distance[0]*distance[0])/dlen3-ii*(distance[0]*distance[1])/dlen3 );
	      dp_y = md*powered*( ii*(1.0/dlen)-(distance[0]*distance[1])/dlen3-ii*(distance[1]*distance[1])/dlen3 );
	      dp_z = md*powered*( -(distance[0]*distance[2])/dlen3-ii*(distance[1]*distance[2])/dlen3 );

	      // Derivatives of real and imaginary parts of above
	      real_dz[0] = real( dp_x ); real_dz[1] = real( dp_y ); real_dz[2] = real( dp_z );
	      imag_dz[0] = imag( dp_x ); imag_dz[1] = imag( dp_y ); imag_dz[2] = imag( dp_z );  

	      // Complete derivative of steinhardt parameter
	      myrealvec = (+sw)*dpoly_ass*real_z*dz + (+dfunc)*distance*tq6 + (+sw)*poly_ass*real_dz; 
	      myimagvec = (+sw)*dpoly_ass*imag_z*dz + (+dfunc)*distance*itq6 + (+sw)*poly_ass*imag_dz;

	      // Real part
	      if (posiFlag){
		ncoord += sw*tq6 ;
		deriv[i0] = deriv[i0] + (-myrealvec);
		deriv[i1] = deriv[i1] + (myrealvec);
		virial = virial + Tensor( -myrealvec,distance );
		// addComponent( tmom+m, sw*tq6 );
		// addAtomsDerivative( tmom+m, 0, -myrealvec );
		// addAtomsDerivative( tmom+m, i, myrealvec );
		// addBoxDerivatives( tmom+m, Tensor( -myrealvec,distance ) );
	      }
	      else {
		double pref=1.;
		if (m % 2 != 0) pref = -1.;
		ncoord += pref*sw*tq6 ;
		deriv[i0] = deriv[i0] + (-pref * myrealvec);
		deriv[i1] = deriv[i1] + (pref * myrealvec);
		virial = virial + pref*Tensor( -myrealvec,distance );
		// Store -m part of vector
		// -m part of vector is just +m part multiplied by (-1.0)**m and multiplied by complex
		// conjugate of Legendre polynomial
		// Real part
		// addComponent( tmom-m, pref*sw*tq6 );
		// addAtomsDerivative( tmom-m, 0, -pref*myrealvec );
		// addAtomsDerivative( tmom-m, i, pref*myrealvec );
		// addBoxDerivatives( tmom-m, pref*Tensor( -myrealvec,distance ) );
	      }
	    }// end m != 0
	  }// end sw larger than tolerance
	}// end Steinhardt
      }// end neighbor loop

      if(!serial){
	comm.Sum(ncoord);
	if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
	comm.Sum(virial);
      }

      for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i, steinhardtPrefactor * deriv[i]);
      setValue           (ncoord * steinhardtPrefactor);
      setBoxDerivatives  (virial * steinhardtPrefactor);
    }
  }
}


