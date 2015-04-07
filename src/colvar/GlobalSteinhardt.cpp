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
      keys.add("compulsory","NN","6","The n parameter of the switching function ");
      keys.add("compulsory","MM","12","The m parameter of the switching function ");
      keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
      keys.add("compulsory","R_0","The r_0 parameter of the switching function");
      keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
	       "The following provides information on the \\ref switchingfunction that are available. " 
	       "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 
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

      vector<double > qvalue_real (tmom + 1, 0.);
      vector<double > qvalue_imag (tmom + 1, 0.);
      vector<Tensor > virial_real (tmom + 1);
      vector<Tensor > virial_imag (tmom + 1);
      vector<vector<Vector > > deriv_real (tmom + 1);
      vector<vector<Vector > > deriv_imag (tmom + 1);
      for (unsigned ii = 0; ii < tmom + 1; ++ ii){
	deriv_real[ii].resize (getNumberOfAtoms());
	deriv_imag[ii].resize (getNumberOfAtoms());
      }
      double ncoord_value = 0;
      Tensor ncoord_virial;
      vector<Vector> ncoord_deriv(getNumberOfAtoms());

      double steinhardtPrefactor = 1. * sqrt(4*M_PI / (2. * tmom + 1.));

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

	{
	  double dfunc, dpoly_ass, md, tq6, itq6, real_z, imag_z; 
	  // The square root of -1
	  std::complex<double> ii( 0.0, 1.0 ), dp_x, dp_y, dp_z;
	  double sw, poly_ass, dlen;
	  std::complex<double> powered;
	  dlen=distance.modulo(); 
	  sw = switchingFunction.calculate( dlen, dfunc );

	  if( sw>=1e-12 ){
	    ncoord_value += sw;
	    ncoord_deriv[i0] = ncoord_deriv[i0] + (-dfunc)*distance ;
	    ncoord_deriv[i1] = ncoord_deriv[i1] + dfunc*distance ;
	    ncoord_virial = ncoord_virial + (-dfunc)*Tensor(distance,distance);

	    double dlen3 = dlen*dlen*dlen;
	    // Derivatives of z/r wrt x, y, z
	    Vector dz (-( distance[2] / dlen3 )*distance);
	    dz[2] += (1.0 / dlen);

	    // Do stuff for m=0
	    poly_ass=deriv_poly( 0, distance[2]/dlen, dpoly_ass );
	    // Derivative wrt to the vector connecting the two atoms
	    Vector myrealvec = (+sw)*dpoly_ass*dz + poly_ass*(+dfunc)*distance;
	    // And store the vector function
	    qvalue_real[0] += sw*poly_ass;
	    // printf ("atom %d with %d, \t distance %f %f %f, \t value_acc %f\n", i0, i1, distance[0], distance[1], distance[2], sw*poly_ass);
	    // Accumulate the derivatives
	    deriv_real[0][i0] = deriv_real[0][i0] + (-myrealvec);
	    deriv_real[0][i1] = deriv_real[0][i1] + (myrealvec);
	    virial_real[0] = virial_real[0] + Tensor( -myrealvec,distance );

	    for (unsigned m = 1; m <= tmom; ++m) {

	      std::complex<double> com1( distance[0]/dlen ,distance[1]/dlen );
	      // Calculate Legendre Polynomial
	      poly_ass=deriv_poly( m, distance[2]/dlen, dpoly_ass );
	      // Calculate powe of complex number
	      powered=pow(com1,m-1); md=static_cast<double>(m);

	      // Derivatives wrt ( x/r + iy )^m
	      dp_x = md*powered*( (1.0/dlen)-(distance[0]*distance[0])/dlen3-ii*(distance[0]*distance[1])/dlen3 );
	      dp_y = md*powered*( ii*(1.0/dlen)-(distance[0]*distance[1])/dlen3-ii*(distance[1]*distance[1])/dlen3 );
	      dp_z = md*powered*( -(distance[0]*distance[2])/dlen3-ii*(distance[1]*distance[2])/dlen3 );
	      // Derivatives of real and imaginary parts of above
	      Vector real_dz, imag_dz;
	      real_dz[0] = real( dp_x ); real_dz[1] = real( dp_y ); real_dz[2] = real( dp_z );
	      imag_dz[0] = imag( dp_x ); imag_dz[1] = imag( dp_y ); imag_dz[2] = imag( dp_z );

	      // Real part
	      // Real and imaginary parts of z
	      real_z = real(com1*powered); 
	      // Calculate steinhardt parameter
	      tq6=poly_ass*real_z;   // Real part of steinhardt parameter
	      // Complete derivative of steinhardt parameter
	      Vector myrealvec ( (+sw)*dpoly_ass*real_z*dz + (+dfunc)*distance*tq6 + (+sw)*poly_ass*real_dz ); 
	      qvalue_real[m] += sw*tq6 ;
	      deriv_real[m][i0] = deriv_real[m][i0] +  (-myrealvec);
	      deriv_real[m][i1] = deriv_real[m][i1] +  (myrealvec);
	      virial_real[m] = virial_real[m] + Tensor( -myrealvec,distance );

	      imag_z = imag(com1*powered );
	      itq6=poly_ass*imag_z;  // Imaginary part of steinhardt parameter
	      Vector myimagvec ( (+sw)*dpoly_ass*imag_z*dz + (+dfunc)*distance*itq6 + (+sw)*poly_ass*imag_dz );
	      qvalue_imag[m] += sw*itq6;
	      deriv_imag[m][i0] = deriv_imag[m][i0] +  (-myimagvec);
	      deriv_imag[m][i1] = deriv_imag[m][i1] +  (myimagvec);
	      virial_imag[m] = virial_imag[m] + Tensor( -myimagvec,distance ) ;
	    }
	  }// end sw larger than tolerance
	}// end Steinhardt
      }// end neighbor loop

      if(!serial){
	for (unsigned m = 0; m <= tmom; ++m){
	  comm.Sum(qvalue_real[m]);
	  if(!deriv_real[m].empty()) comm.Sum(&deriv_real[m][0][0],3*deriv_real[m].size());
	  comm.Sum(virial_real[m]);
	  comm.Sum(qvalue_imag[m]);
	  if(!deriv_imag[m].empty()) comm.Sum(&deriv_imag[m][0][0],3*deriv_imag[m].size());
	  comm.Sum(virial_imag[m]);
	}
	comm.Sum(ncoord_value);
	if(!ncoord_deriv.empty()) comm.Sum(&ncoord_deriv[0][0],3*ncoord_deriv.size());
	comm.Sum(ncoord_virial);
      }
      
      // for(unsigned i=0;i<deriv.size();++i) {
      // 	setAtomsDerivatives(i, steinhardtPrefactor * deriv[i] );
      // }
      // setValue           (steinhardtPrefactor * qvalue  );
      // setBoxDerivatives  (steinhardtPrefactor * (virial) );

      double qvalue=0.;
      Tensor virial;
      vector<Vector> deriv(getNumberOfAtoms());
      qvalue = qvalue_real[0] * qvalue_real[0];
      for (unsigned m = 1; m <= tmom; ++m){
	qvalue += 2. * (qvalue_real[m] * qvalue_real[m] + qvalue_imag[m] * qvalue_imag[m]);
      }
      qvalue = sqrt(qvalue);

      for (unsigned ii = 0; ii < deriv.size(); ++ii){
	deriv[ii] += qvalue_real[0] * deriv_real[0][ii];
	for (unsigned m = 1; m <= tmom; ++ m){
	  deriv[ii] += 2. * (qvalue_real[m] * deriv_real[m][ii] + qvalue_imag[m] * deriv_imag[m][ii]);
	}
	deriv[ii] = deriv[ii] * (1./qvalue);
      }
      
      virial += qvalue_real[0] * virial_real[0];
      for (unsigned m = 1; m <= tmom; ++ m){
	virial += 2. * (qvalue_real[m] * virial_real[m] + qvalue_imag[m] * virial_imag[m]);
      }
      virial = virial * (1./qvalue);      

      double ncoord_valuei = 1./(ncoord_value);
      double ncoord_valuei2 = ncoord_valuei * ncoord_valuei;
      for(unsigned i=0;i<deriv.size();++i) {
      	setAtomsDerivatives(i,
      			    steinhardtPrefactor * ( deriv[i] * ncoord_value - qvalue * ncoord_deriv[i] ) * ncoord_valuei2 );
      }
      setValue           (steinhardtPrefactor * qvalue * ncoord_valuei );
      setBoxDerivatives  (steinhardtPrefactor * (virial * ncoord_value - qvalue * ncoord_virial) * ncoord_valuei2 );
    }
  }
}


