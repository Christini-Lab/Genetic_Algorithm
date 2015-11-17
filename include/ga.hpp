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

#ifndef _GA_H
#define _GA_H

#include <iostream>
#include <fstream>

#include "globalSetup.hpp"
#include "population.hpp"
#include "random.hpp"

extern GlobalSetup *globalSetup;
extern Random myRandom;

class GA
{
  int genID;
  long noOfGlobalEvals;
  long noOfLocalEvals;
  Population *population;
  bool stoppingCriteria(void);
  int successiveNoChangeInBestFitness;
  int successiveNoChangeInBestObjective;
  int successiveNoChangeInNoOfFronts;
  int successiveNoChangeInAvgObjective;
  int successiveNoChangeInAvgFitness;
  int successiveNoChangeInFitnessVar;
 public:
  GA(void);
  ~GA(void);

  void reinitialize(void);
  bool generate();
  bool nsgaGenerate();
  Individual getBestGuy() {return *(population->bestInd);}
  void freeze(int locus) {population->freeze(locus);}
  void flood(int locus) {population->flood(locus);}
  void gaOutput(std::ostream &out) {out << *population;}

  Population *getPopulation(void) {return population;}

  inline int getGenID (void) const {return genID;}
  inline long getNoOfGlobalEvals(void) const {return noOfGlobalEvals;}
  inline long getNoOfLocalEvals(void) const {return noOfLocalEvals;}
  inline int getSuccessiveNoChangeInBestFitness(void) const {
    return successiveNoChangeInBestFitness;
  }
  inline int getSuccessiveNoChangeInBestObjective(void) const {
    return successiveNoChangeInBestObjective;
  }
  inline int getSuccessiveNoChangeInNoOfFronts(void) const {
    return successiveNoChangeInNoOfFronts;
  }
  inline int getSuccessiveNoChangeInAvgObjective(void) const {
    return successiveNoChangeInAvgObjective;
  }
  inline int getSuccessiveNoChangeInAvgFitness(void) const {
    return successiveNoChangeInAvgFitness;
  }
  inline int getSuccessiveNoChangeInFitnessVar(void) const {
    return successiveNoChangeInFitnessVar;
  }

  inline void setGenID (int iValue) {genID = iValue;}
  inline void setNoOfGlobalEvals(long lValue) {noOfGlobalEvals = lValue;}
  inline void setNoOfLocalEvals(long lValue) {noOfLocalEvals = lValue;}
  inline void setSuccessiveNoChangeInBestFitness(int iValue) {
    successiveNoChangeInBestFitness = iValue;
  }
  inline void setSuccessiveNoChangeInBestObjective(int iValue) {
    successiveNoChangeInBestObjective = iValue;
  }
  inline void setSuccessiveNoChangeInNoOfFronts(int iValue) {
    successiveNoChangeInNoOfFronts = iValue;
  }
  inline void setSuccessiveNoChangeInAvgObjective(int iValue) {
    successiveNoChangeInAvgObjective = iValue;
  }
  inline void setSuccessiveNoChangeInAvgFitness(int iValue) {
    successiveNoChangeInAvgFitness = iValue;
  }
  inline void setSuccessiveNoChangeInFitnessVar(int iValue) {
    successiveNoChangeInFitnessVar = iValue;
  }
};

#endif
