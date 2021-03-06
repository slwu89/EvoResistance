% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppR6.R
\name{parameters}
\alias{parameters}
\title{Parameters object for EvoResistance simulations.}
\usage{
parameters(x0, c, delta, mu, sdo, sdr, sro, sdd, srr, Ne, N, sigma)
}
\arguments{
\item{x0:}{Introduction frequency of driver allele.}

\item{c:}{Cleavage rate in drive/wild-type heterozygotes.}

\item{delta:}{Fraction of cleavage events resulting in a resistant allele by NHEJ.}

\item{mu:}{Mutation rate of wild-type alleles, resulting in resistant allele.}

\item{sdo:}{Fitness cost of driver/wild-type heterozygote.}

\item{sdr:}{Fitness cost of driver/resistant heterozygote.}

\item{sro:}{Fitness cost of resistant/wild-type heterozygote.}

\item{sdd:}{Fitness cost of driver homozote.}

\item{srr:}{Fitness cost of resistant homozygote.}

\item{Ne:}{Variance effective population size.}

\item{N:}{Census population size.}

\item{sigma:}{Precision for simulations.}
}
\description{
This creates an object "Parameters" that contains all of the parameters for the simulation. It is the input for all further classes and functions used to calulate the allele frequences and probability of mutation from different sources.
}
