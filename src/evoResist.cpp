#include <EvoResistance.h>

//' Probability of resistance through de novo mutation.
//'
//' p_mu calculates the probability of resistance from de novo mutations. This is the probability of any mutation being resistant in t<tfix.
//'
//' @param P Parameter object. Create this using the "parameters()" constructor.
//' @param T Trajectory object. Create this using the "trajectory()" constructor.
//' @param F Fixation_probability object. Create this using the "fixation_probability()" constructor.
//' @export
// [[Rcpp::export]]
double p_mu(parameters& P, trajectory& T, fixation_probability& F)
  //Eqn. 13, Probability of resistance from de novo mutations. This is the probability of any mutation being resistant in t<tfix
{
  double p_cum = 1;
  int t=1;
  
  while(1)
    {
	//Below is Eqn. 5. This defines the rate of de novo mutations after the driver is introduced. 
	//2 alleles
	//P.N = population size
	//P.mu = population mutation rate, given as a parameter to the eqn.
	//(1-T.x(t))*)1-T.x(t)) = fraction of homozygous wild type alleles. Notice, this assumes few resistance alleles.
	//T.x(t)*(1-T.x(t))*(1-c) = fraction of heterozygotes times the chance that there is no drive (no conversion, but not due
	//mistakes, failed conversion, or broekn driver allele. Simply the chance that drive doesn't happen.)
      double m = 2*P.N*P.mu*( (1-T.x(t))*(1-T.x(t))+T.x(t)*(1-T.x(t))*(1-P.c) );
        //Below is the calculation in Eqn. 13, probability of mutation times probability of fixation
      double p = exp(-m*F.pi(t));
      if(p>1-P.sigma) { break; }
      else { p_cum = p_cum*p; }	//sum in Eqn. 13
      t++;
    }
  return 1-p_cum;	//Final solution to Eqn. 13
}

//' Probability of resistance through NHEJ.
//'
//' p_delta calculates the probability of resistance from NHEJ during driver homing. This is the probability of any mutation being resistant in t<tfix.
//'
//' @param P Parameter object. Create this using the "parameters()" constructor.
//' @param T Trajectory object. Create this using the "trajectory()" constructor.
//' @param F Fixation_probability object. Create this using the "fixation_probability()" constructor.
//' @export
// [[Rcpp::export]]
double p_delta(parameters& P, trajectory& T, fixation_probability& F)
  //Eqn. 15, probability of resistance from NHEJ, from t0 to t<tfix
{
  double p_cum = 1;
  int t=1;
  
  while(1)
    {
	//Eqn. 6. This defines the rate of NHEJ mutations after the driver is introduced
	//2 = alleles
	//P.N = Population size
	//P.delta = Fraction of cleavage events resulting in restant allels through NHEJ
	//P.c = probability of drive attempting to convert WT allele to driver allele. This is the only time NHEJ will be active
	//T.x(t)*(1-T.x(t)) = fraction of heterozygotes in the population. This is the only time conversion can occur. 
      double m = 2*P.N*P.delta*P.c*T.x(t)*(1-T.x(t));
        //Actual calculation in Eqn. 15
      double p = exp(-m*F.pi(t));
      if(p>1-P.sigma) { break; }
      else { p_cum = p_cum*p; } //Sum in Eqn. 15
      t++;
    }
  return 1-p_cum;	//Final solution to Eqn. 15
}

//' Probability of resistance through standing genetic variation.
//'
//' p_sgv calculates the probability of resistance existing as standing genetic variation. This is the probability of resistance becoming fixed during t<tfix.
//'
//' @param P Parameter object. Create this using the "parameters()" constructor.
//' @param F Fixation_probability object. Create this using the "fixation_probability()" constructor.
//' @export
// [[Rcpp::export]]
double p_sgv(parameters& P, fixation_probability& F)
  //Eqn. 12, probability of resistance from standing genetic variation. These are derived from Hermisson and Pennings 2005, Eqn. 8
{
  if(P.sro>0) // codominant, the resistant/WT heterozygote has a fitness cost greater than 0 but less than the resistant homozygote
    {
	// 4*P.Ne*P.mu -> theta, twice the population-level mutation rate
	// 2*P.Ne*F.pi(0) = establishment probability of mutation present at t0
	// 4*P.Ne*P.sro = fitness cost of heterozygote
      return 1-exp(-4*P.Ne*P.mu*log(1+2*P.Ne*F.pi(0)/(4*P.Ne*P.sro+1)));
    }
  else // recessive, resistant/WT heterozygote has no fitness cost, resistant homozygote does have one.
    {
	//sqrt(2*P.Ne*P.srr) = discussed in A11 of Hermissson and Pennings 2005.
      return 1-exp(-4*P.Ne*P.mu*log(1+2*P.Ne*F.pi(0)/(sqrt(2*P.Ne*P.srr)+1)));
    }
}

