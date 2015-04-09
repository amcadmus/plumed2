/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "Bias.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/Random.h"


using namespace std;


namespace PLMD{
  namespace bias{

    class AFED : public Bias{
      double			kb ;
      int			seed;
      double			dt;
      std::vector<double>	temp;
      std::vector<double>	tau;
      std::vector<double>	kappa;
      std::vector<double>	gamma;
      std::vector<double>	sigma;
      std::vector<double>	mass;

      std::vector<double>	oldaa;
      std::vector<double>	oldaav;
      std::vector<double>	oldf;
      std::vector<double>	work;

      Random			rand;

      void   stepA (const double & mydt,
		    std::vector<double > & myaa,
		    std::vector<double > & myaav,
		    std::vector<double > & myf) ;
      double stepB (const double & mydt,
		    std::vector<double > & myaa,
		    std::vector<double > & myaav,
		    std::vector<double > & myf) ;
      void   stepO (const double & mydt,
		    std::vector<double > & myaa,
		    std::vector<double > & myaav,
		    std::vector<double > & myf) ;
  public:
      AFED(const ActionOptions&);
      void calculate();
      static void registerKeywords( Keywords& keys );
    };

    PLUMED_REGISTER_ACTION(AFED,"AFED")

    void AFED::registerKeywords( Keywords& keys ){
      Bias::registerKeywords(keys);
      keys.use("ARG");
      keys.add ("compulsory", "SEED", "0", "The random seed ");
      keys.add ("compulsory", "DT", "0.001", "The time step ");
      keys.add ("compulsory", "TEMP", "set the target temperature for CVs ");
      keys.add ("compulsory", "TAU",  "set time-scale of temperature control for CVs ");
      keys.add ("compulsory", "KAPPA",  "restraint constant ");
      componentsAreNotOptional(keys);
      keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
      keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
      keys.addOutputComponent("_cntr","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
			      "these quantities will named with  the arguments of the bias followed by "
			      "the character string _cntr. These quantities give the instantaneous position "
			      "of the center of the harmonic potential.");
      keys.addOutputComponent("_cntrv","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
			      "these quantities will named with  the arguments of the bias followed by "
			      "the character string _cntrv. These quantities give the instantaneous velocity "
			      "of the center of the harmonic potential.");
      keys.addOutputComponent("_work","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
			      "These quantities will named with the arguments of the bias followed by "
			      "the character string _work. These quantities tell the user how much work has "
			      "been done by the potential in dragging the system along the various colvar axis.");
    }

    AFED::AFED(const ActionOptions&ao):
	PLUMED_BIAS_INIT(ao),
	kb(0.008314),
	temp(getNumberOfArguments()),
	tau(getNumberOfArguments()),
	kappa(getNumberOfArguments()),
	gamma(getNumberOfArguments()),
	sigma(getNumberOfArguments()),
	mass(getNumberOfArguments()),
	oldaa(getNumberOfArguments(), 0.),
	oldaav(getNumberOfArguments(), 0.),
	oldf(getNumberOfArguments(), 0.),
	work(getNumberOfArguments(), 0.)
    {
      parse ("DT",dt);
      parseVector ("TEMP",temp);
      parseVector ("TAU",tau);
      parseVector ("KAPPA",kappa);
      parse ("SEED", seed);
      
      unsigned rank = 0;
      // unsigned rank=comm.Get_rank();
      rand.setSeed (seed + rank);
      
      for (unsigned ii = 0; ii < getNumberOfArguments(); ++ ii){
	mass[ii] = tau[ii] * tau[ii] * kb * temp[ii];	// from Yu et.cl. JCP 2014	
	gamma[ii] = 1./tau[ii];
	sigma[ii] = sqrt (2. * gamma[ii] * (kb * temp[ii]));
      }

      log.printf ("  rand seed is %d\n", seed);
      if (getNumberOfArguments() != temp.size()){
	error ("the size of TEMP should match the number of arguments");
      }
      if (getNumberOfArguments() != tau.size()){
	error ("the size of TAU should match the number of arguments");
      }
      if (getNumberOfArguments() != kappa.size()){
	error ("the size of KAPPA should match the number of arguments");
      }
      log.printf ("  using temperature: ");
      for (unsigned ii = 0; ii < temp.size(); ++ii){
	log.printf ("%f ", temp[ii]);
      }
      log.printf ("\n");
      log.printf ("  using tau: ");
      for (unsigned ii = 0; ii < tau.size(); ++ii){
	log.printf ("%f ", tau[ii]);
      }
      log.printf ("\n");
      log.printf ("  using kappa: ");
      for (unsigned ii = 0; ii < kappa.size(); ++ii){
	log.printf ("%f ", kappa[ii]);
      }
      log.printf ("\n");      
      log.printf ("  using mass: ");
      for (unsigned ii = 0; ii < kappa.size(); ++ii){
	log.printf ("%f ", mass[ii]);
      }
      log.printf ("\n");      

      checkRead();

      addComponent("bias"); componentIsNotPeriodic("bias");
      addComponent("force2"); componentIsNotPeriodic("force2");

      // add the centers of the restraint as additional components that can be retrieved (useful for debug)

      std::string comp;
      for(unsigned i=0;i< getNumberOfArguments() ;i++){
	comp=getPntrToArgument(i)->getName()+"_cntr"; // each spring has its own center 
        addComponent(comp); componentIsNotPeriodic(comp);
	comp=getPntrToArgument(i)->getName()+"_work"; // each spring has its own work
        addComponent(comp); componentIsNotPeriodic(comp);
	comp=getPntrToArgument(i)->getName()+"_cntrv"; // each spring has its own velocity 
        addComponent(comp); componentIsNotPeriodic(comp);
      }

      log<<"  Bibliography ";
      log<<cite("Grubmuller, Heymann, and Tavan, Science 271, 997 (1996)")<<"\n";

    }

    double AFED::stepB (const double & mydt,
		      std::vector<double > & myaa,
		      std::vector<double > & myaav,
		      std::vector<double > & myf) 
    {
      unsigned narg=getNumberOfArguments();
      double ener = 0.;
      for (unsigned ii = 0; ii < narg; ++ ii ){
	const double cv=difference(ii, myaa[ii], getArgument(ii)); // this gives: getArgument(i) - aa[i]
	myf[ii] = - kappa[ii] * cv;	// this is the f to the system, 
	myaav[ii] += -myf[ii] * mydt;	// to evlove the cv, minus the f.
	ener += 0.5 * kappa[ii] * cv * cv;
      }
      return ener;
    }

    void AFED::stepA (const double & mydt,
		      std::vector<double > & myaa,
		      std::vector<double > & myaav,
		      std::vector<double > & myf) 
    {
      unsigned narg=getNumberOfArguments();
      for (unsigned ii = 0; ii < narg; ++ ii ){
	myaa[ii] += myaav[ii] * dt / mass[ii];
      }
    }

    void AFED::stepO (const double & mydt,
		      std::vector<double > & myaa,
		      std::vector<double > & myaav,
		      std::vector<double > & myf) 
    {
      unsigned narg=getNumberOfArguments();
      for (unsigned ii = 0; ii < narg; ++ ii ){
	double expgammadt = exp (-gamma[ii] * dt) ;
	double rand_gaussian = rand.Gaussian();
	// std::cout << "rank : " << comm.Get_rank()  << " rand " << rand_gaussian << std::endl;
	// std::cout << std::fflush;
	myaav[ii] = 
	    expgammadt * myaav[ii] + 
	    sigma[ii] / (sqrt(2. * gamma[ii])) * sqrt(1. - expgammadt * expgammadt) * sqrt(mass[ii]) * rand_gaussian;
      }   
    }

    void AFED::calculate(){

      double ene=0.0;
      double totf2=0.0;
      unsigned narg=getNumberOfArguments();

      std::vector<double > aa(oldaa), aav(oldaav), f(oldf);

      long int now=getStep();
      if (now == 0){
	for(unsigned i=0;i<narg;++i){
	  aa[i] = getArgument(i);
	  aav[i] = 0.;
	}
      }

      stepB (0.5 * dt, aa, aav, f);
      stepA (0.5 * dt, aa, aav, f);
      stepO (1.0 * dt, aa, aav, f);
      stepA (0.5 * dt, aa, aav, f);
      ene = stepB (0.5 * dt, aa, aav, f);

      for(unsigned i=0;i<narg;++i){
	// std::cout << "write to : " << getPntrToArgument(i)->getName() << std::endl;
	getPntrToComponent(getPntrToArgument(i)->getName()+"_cntr")->set(aa[i]); 
	getPntrToComponent(getPntrToArgument(i)->getName()+"_cntrv")->set(aav[i]/mass[ii]); 
	if(oldaa.size()==aa.size() && oldf.size()==f.size()) {
	  work[i] += 0.5 * (oldf[i] + f[i]) * (aa[i] - oldaa[i]) ;
	}
	getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(work[i]); 
	setOutputForce(i,f[i]);
	totf2+=f[i]*f[i];
      };
      oldf=f;
      oldaa=aa;
      oldaav = aav;
      getPntrToComponent("bias")->set(ene);
      getPntrToComponent("force2")->set(totf2);
    }

  }
}


