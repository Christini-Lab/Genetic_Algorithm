/*
  Christini lab genetic algorithm
  Modified from IlliGAL GA, reference below
*/

/*
 Single & Multi-Objective Real-Coded Genetic Algorithms Code
 Author: Kumara Sastry
 Illinois Genetic Algorithms Laboratory (IlliGAL)
 Deparment of General Engineering
 University of Illinois at Urbana-Champaign
 104 S. Mathews Ave, Urbana, IL 61801
*/

#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include "globalSetup.hpp"
#include "random.hpp"
#include <iostream>

extern GlobalSetup *globalSetup;
extern Random myRandom;

class Chromosome
{
  private:
    double *genes;
  public:
  Chromosome();
  Chromosome(const Chromosome &chrom);
  ~Chromosome();
  inline void setValue(const int &locus, const double &value) {
    genes[locus] = value; }
  inline double &operator[](const int &locus) const { return genes[locus]; }
  void mutateMinMax(const int *freezeMask);
  void mutateNormal(const int *freezeMask);
  void mutatePolynomial(const int *freezeMask);
  void mutatePolynomial(const int index);
  Chromosome &operator= (const Chromosome &sourceChrom);
  friend std::ostream &operator<< (std::ostream &out, Chromosome &chrom);
};

#endif
