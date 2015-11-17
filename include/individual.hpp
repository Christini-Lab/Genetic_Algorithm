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

#ifndef _INDIVIDUAL_H
#define _INDIVIDUAL_H

#include "globalSetup.hpp"
#include "random.hpp"
#include "chromosome.hpp"
#include <iostream>
#include <stdlib.h>

extern GlobalSetup *globalSetup;

class Individual {
protected:
  Chromosome chrom;
  double *fitness;
  double *objFunction;
  double penalty;
  int noOfViolations;
  double *violation;
public:
  Individual(void);
  Individual(const Individual &sourceInd );
  Individual(const Individual *sourceInd );
  virtual ~Individual(void);

  void mutate(int *freezeMask);
  inline void mutate(int index) {chrom.mutatePolynomial(index);}
  inline double *getViolation() { return violation; }
  inline double getViolation(int ii) const {return violation[ii];}
  inline void setViolation(int index, double newViolation) {
    violation[index] = newViolation;
  }
  inline double getPenalty() const { return penalty; }
  inline int getNoOfViolations() const { return noOfViolations; }
  inline void setPenalty(double dValue) { penalty = dValue; }
  inline void setNoOfViolations(int iValue) { noOfViolations = iValue; }
  inline double getObjective() const { return *objFunction; }
  inline double getFitness() const {return *fitness; }
  inline void setObjective(double newObjective) {
    if (globalSetup->gaType==SGA) *objFunction = newObjective; else exit(0);
  }
  inline void setFitness(double newFitness) {
    if (globalSetup->gaType==SGA) *fitness = newFitness; else exit(0);
  }
  inline void setValue(int index, double value) { chrom.setValue(index, value);}
  inline void freeze(int index, double value) { chrom.setValue(index, value); }
  inline void flood(int index) {
    double value;
    value = myRandom.boundedRandom(
        globalSetup->variableRanges[index][0],
        globalSetup->variableRanges[index][1]);
    if(globalSetup->variableTypes[index] == Integer)
      value = int(value+0.5);
    chrom.setValue(index,value);
  }
  void evaluateFitness(void);
  void loadIndividual(double *varValues, double *objValues,
                      double *constViolValues, double penaltyValue);
  friend bool constrTournComp (const Individual *guy1, const Individual *guy2,
                               int compareWhat);
  friend bool constrPenaltyComp(const Individual *guy1, const Individual *guy2);
  friend bool individualComp(const Individual *guy1, const Individual *guy2);
  friend bool isBetter(const Individual *guy1, const Individual *guy2);
  friend std::ostream &operator<< (std::ostream &out, Individual &x);
  inline double &operator[] (const int index) {return chrom[index];}

  Individual &operator= (const Individual &sourceInd);
};

class NsgaIndividual: public Individual {

protected:
  int rank;
  double crowdingDistance;

public:
  NsgaIndividual(void);
  NsgaIndividual(const Individual &sourceInd);
  NsgaIndividual(const Individual *sourceInd);
  NsgaIndividual(const NsgaIndividual &sourceInd);
  NsgaIndividual(const NsgaIndividual *sourceInd);
  ~NsgaIndividual(void);

  // Other member functions
  friend int dominates(const NsgaIndividual *guy1, const NsgaIndividual *guy2);
  friend bool crowdingComp(const NsgaIndividual *guy1,
                           const NsgaIndividual *guy2);
  inline void setFitness(int index, double newFitness) {
    //    cout << index << " " << newFitness << endl;
    fitness[index] = newFitness;}
  inline void setObjective(int index, double newObjective) {
    //    cout << index << " " << newObjective << endl;
    objFunction[index] = newObjective;}
  inline double getFitness(int ii) const {return fitness[ii];}
  inline double getObjective(int ii) const {return objFunction[ii];}
  inline double getViolation(int ii) const {return violation[ii];}
  inline int getRank(void) const { return rank; }
  inline double getCrowdingDistance (void) const { return crowdingDistance; }
  inline void setRank(int frontRank) { rank = frontRank; }
  inline void setCrowdingDistance(double distance) {
    crowdingDistance = distance;
  }

  // Different type of operator=.
  NsgaIndividual & operator= (const NsgaIndividual &sourceInd);
  NsgaIndividual & operator= (const Individual &sourceInd);

};

#endif
