// -*-c++-*-
#ifndef _EVORESISTANCE_H_
#define _EVORESISTANCE_H_


#include <RcppCommon.h>

// parameters class
class parameters{
public:
  parameters(double x0, double c, double delta, double mu, double sdo, double sdr, double sdd, double srr, double Ne, double N, double sigma): x0(x0), c(c), delta(delta), mu(mu), sdo(sdo), sdr(sdr), sdd(sdd), srr(srr), Ne(Ne), N(N), sigma(sigma) {}
  double x0;
  double c;
  double delta;
  double mu;
  double sdo;
  double sdr;
  double sro;
  double sdd;
  double srr;
  double Ne;
  double N;
  double sigma;
};


// w function
inline double w(double x, parameters& P){
  return 2*x*(1-x)*((1-P.c)*(1-P.sdo)+P.c*(1-P.sdd))+x*x*(1-P.sdd)+(1-x)*(1-x);
}

// trajectory class
class trajectory
{
private:
  std::vector<double> X;
public:
  int t_fix;

  trajectory(parameters& P)
  {
    double x = P.x0;

    t_fix = 0;
    while(x<(1-1/(2*P.N)))
      {
	t_fix++;
	X.push_back(x);
	double x_next = (x*(1-x)*((1-P.c)*(1-P.sdo)+2*P.c*(1-P.sdd)) +x*x*(1-P.sdd))/w(x,P);
	x = x_next;
	if(t_fix>1e4) { Rcpp::Rcout << "NA NA NA NA"<< std::endl; Rcpp::stop("you did a bad"); }
      }
  }

  double x(int t)
  {
    if(t>=t_fix) { return 1; }
    else         { return X[t]; }
  }
};

// s function
inline double s(double x, parameters& P)
{
  return ((1-x)*(1-P.sro)+x*(1-P.sdr))/w(x,P) - 1;
}

// g function
inline double g(int t, int t_end, parameters& P, trajectory& T)
{
  double G = 0.0;
  for(int i=t; i<t_end; i++)
    {
      G += s(T.x(i),P);
    }
  return G;
}

// f function
inline double f(int t, int t_end, parameters& P, trajectory& T)
{
  double F = 0.0;
  for(int i=t; i<t_end; i++)
    {
      double S = s(T.x(i),P);
      if(S==0.0) { return -1; }
      F += exp(-g(t,i,P,T)) * ( 1 - exp(-S) )/S;
    }
  return F;
}

// fixation_probability class
class fixation_probability
{
private:
  std::vector<double> p;
  int t_fix;
public:

  fixation_probability(parameters& P, trajectory& T)
  {
    t_fix = T.t_fix;

    for(int t=0; t<t_fix; t++)
      {
	double F = f(t,t_fix,P,T);
	double G = g(t,t_fix,P,T);
	double S = s(1.0,P);

	if(F==-1 || S==0.0) { p.push_back( 0.0 ); }
	else      { p.push_back( 2.0/( 1.0 + (P.N/P.Ne)*(F+exp(-G)/S) ) ); }
      }

    double S = s(1.0,P);

    if(S==0) { p.push_back( 0.0 ); }
    else     { p.push_back( 2.0/( 1.0 + (P.N/P.Ne)*(1.0/S) ) ); }
  }

  double pi(int t)
  {
    if(t>=t_fix) { return p.back(); }
    else         { return p[t]; }
  }
};

#include <EvoResistance/RcppR6_pre.hpp>
#include <EvoResistance/RcppR6_post.hpp>

#endif
