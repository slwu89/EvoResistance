#include <EvoResistance.h>

//' p_mu
//'
//' p_mu
//'
//' @param P parameters
//' @param T trajectory
//' @param F fixation_probability
//' @export
// [[Rcpp::export]]
double p_mu(parameters& P, trajectory& T, fixation_probability& F)
{
  double p_cum = 1;
  int t=1;

  while(1)
    {
      double m = 2*P.N*P.mu*( (1-T.x(t))*(1-T.x(t))+T.x(t)*(1-T.x(t))*(1-P.c) );
      double p = exp(-m*F.pi(t));
      if(p>1-P.sigma) { break; }
      else { p_cum = p_cum*p; }
      t++;
    }
  return 1-p_cum;
}

//' p_delta
//'
//' p_delta
//'
//' @param P parameters
//' @param T trajectory
//' @param F fixation_probability
//' @export
// [[Rcpp::export]]
double p_delta(parameters& P, trajectory& T, fixation_probability& F)
{
  double p_cum = 1;
  int t=1;

  while(1)
    {
      double m = 2*P.N*P.delta*P.c*T.x(t)*(1-T.x(t));
      double p = exp(-m*F.pi(t));
      if(p>1-P.sigma) { break; }
      else { p_cum = p_cum*p; }
      t++;
    }
  return 1-p_cum;
}

//' p_sgv
//'
//' p_sgv
//'
//' @param P parameters
//' @param F fixation_probability
//' @export
// [[Rcpp::export]]
double p_sgv(parameters& P, fixation_probability& F)
{
  if(P.sro>0) // codominant
    {
      return 1-exp(-4*P.Ne*P.mu*log(1+2*P.Ne*F.pi(0)/(4*P.Ne*P.sro+1)));
    }
  else // recessive
    {
      return 1-exp(-4*P.Ne*P.mu*log(1+2*P.Ne*F.pi(0)/(sqrt(2*P.Ne*P.srr)+1)));
    }
}
