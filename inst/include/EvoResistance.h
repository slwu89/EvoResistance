// -*-c++-*-
#ifndef _EVORESISTANCE_H_
#define _EVORESISTANCE_H_


#include <RcppCommon.h>

// parameters class
class parameters
{
public:
  parameters(double x0, double c, double delta, double mu, double sdo, double sdr, double sro, double sdd, double srr, double Ne, double N, double sigma): x0(x0), c(c), delta(delta), mu(mu), sdo(sdo), sdr(sdr), sro(sro), sdd(sdd), srr(srr), Ne(Ne), N(N), sigma(sigma) {}
  double x0;	//Introduction frequency of driver allele
  double c;	//Cleavage Rate in drive/wild-type heterozygotes
  double delta;	//Fraction of cleavage events resulting in resistance allele by NHEJ
  double mu;	//Rate at which wild-type alleles mutation into resistance alleles
  double sdo;	//Fitness cost of driver/wild-type heterozygotes
  double sdr;	//Fitness cost of driver/resistance heterozygotes
  double sro;	//Fitness cost of resistance/wild-type heterozygotes
  double sdd;	//Fitness cost of driver homozygotes
  double srr;	//Fitness cost of resistance homozygotes
  double Ne;	//Variance effective population size
  double N;	//Census population size
  double sigma;	//Precision for simulations
};

//thoughts
//mu does not account for mutations in driver allele that block or modify homing, ie
//	mutations to PAM sequence (abolish homing?), or mutations in gRNA(miss-target events)

// w function
inline double w(double x, parameters& P)
{
  // Average population fitness
  // This ignores NHEJ, and assumes that resistance alleles are low in frequency, such that no homozygotes exist and the frequency
  // of driver allels is 1-(wild type alleles). 
  return 2*x*(1-x)*((1-P.c)*(1-P.sdo)+P.c*(1-P.sdd))+x*x*(1-P.sdd)+(1-x)*(1-x);
}


class trajectory
{
private:
	std::vector<double> X;
public:
  int t_fix;  

  trajectory(parameters& P)
  {
    double x = P.x0;	//starting frequency of driver allele

    t_fix = 0;
    while(x<(1-1/(2*P.N)))	//This is a thing. Time to fixation? I can't remember what this is, but it's defined as something
      {
	t_fix++;
	X.push_back(x);
	double x_next = (x*(1-x)*((1-P.c)*(1-P.sdo)+2*P.c*(1-P.sdd)) +x*x*(1-P.sdd))/w(x,P);	//eqn. 9, frequency of driver allele
	x = x_next;
	if(t_fix>1e4) { Rcpp::Rcout << "NA NA NA NA"<< std::endl; Rcpp::stop("you did a bad."); }
      }
  }
  
  double x(int t)
  {
    if(t>=t_fix) { return 1; } 
    else         { return X[t]; }
  }
};


inline double s(double x, parameters& P)
  //Effective selection coefficient, assuming resistance allele is rare and only exists in heterozygotes.
  //This is a simplification of eqn. 2
{
  return ((1-x)*(1-P.sro)+x*(1-P.sdr))/w(x,P) - 1;	//eqn. 8
}


inline double g(int t, int t_end, parameters& P, trajectory& T)
  //Preparing for eqn. 10. This is the integral of the time-depended effective selection coefficient.
{
  double G = 0.0;
  for(int i=t; i<t_end; i++)
    {
      G += s(T.x(i),P);	//Integral from t0 to tfix of the time dependend effective selection coefficient
    }
  return G;
}


inline double f(int t, int t_end, parameters& P, trajectory& T)
  //Preparing for eqn. 10. This is the integral of a lot of shit. Basically, a continuation of Uecker and Hermisson (2011),
  //it is the approximation for a time period t0 to tend<tfix. Derived in Appendix A of Unckless and Messer.
{
  double F = 0.0;
  for(int i=t; i<t_end; i++)
    {
      double S = s(T.x(i),P);	//selection coefficient at time t, determined through the fraction of driver allele at that time.
      if(S==0.0) { return -1; }	//problems if zero, so define at -1 as a cop-out
      F += exp(-g(t,i,P,T)) * ( 1 - exp(-S) )/S;	//integral of stuff. Basically, this makes eqn. 10 look prettier. 
    }
  return F;
}


class fixation_probability
  //This class calculates the probability of things occuring. The first method is is the fixation probability, eqn 10.
  //p_mu is the probability of a resistance occuring from a de-novo mutation
  //p_delta is the probability of resistance from NHEJ
  //p_sgv is the probability of resistance from standing genetic variation
{
private:
  std::vector<double> p;
  int t_fix;
public:

  fixation_probability(parameters& P, trajectory& T)
	//Eqn. 10. The for loop calculates fixation probability for t<tfix. The final if/else calculates it for t=tfix
	//This part fills a vector of probabilities to be called later. 
  {
    t_fix = T.t_fix;
    
    for(int t=0; t<t_fix; t++)
      {
	double F = f(t,t_fix,P,T); //See above
	double G = g(t,t_fix,P,T); //See above
	double S = s(1.0,P);	   //See above

	if(F==-1 || S==0.0) { p.push_back( 0.0 ); }	//Take care of problem when se(t)=0
	else      { p.push_back( 2.0/( 1.0 + (P.N/P.Ne)*(F+exp(-G)/S) ) ); }	//eqn. 10 for t<tfix. This is derived in appendix A
      }

    double S = s(1.0,P);	//Determine the selection coefficient at fixation

    if(S==0) { p.push_back( 0.0 ); }	//if zero, equation fails, insert 0
    else     { p.push_back( 2.0/( 1.0 + (P.N/P.Ne)*(1.0/S) ) ); }	//Eqn. 10 for t>=tfix. Notice, this is not time-varying.
  }

  double pi(int t)
	//This returns the fixation probability when asked.
  {
    if(t>=t_fix) { return p.back(); }	//If t>tfix, pi(t) isn't time dependent, so just pull the last cell from the vector.
    else         { return p[t]; }	//If t<tfix, pull the appropriate time. 
  }
};

#include <EvoResistance/RcppR6_pre.hpp>
#include <EvoResistance/RcppR6_post.hpp>

#endif
