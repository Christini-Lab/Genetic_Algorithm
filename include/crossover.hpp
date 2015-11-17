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

#ifndef _CROSSOVER_H
#define _CROSSOVER_H

#include "individual.hpp"
#include "globalSetup.hpp"
#include "random.hpp"
#include <iostream>
#include <math.h>

class Crossover {
protected:
  void deterministicCrowding(Individual *parent1, Individual *parent2,
                             Individual *child1, Individual *child2);

public:
  virtual void crossover(Individual *parent1, Individual *parent2)=0;
};

class OneTwoPointCrossover : public Crossover {
private:
  int noOfPoints;
public:
  OneTwoPointCrossover(int numPoints);
  void crossover (Individual *parent1, Individual *parent2);
};

class UniformCrossover: public Crossover {
private:
  double genewiseProbability;
public:
  UniformCrossover(void);
  UniformCrossover(double geneProb);
  void crossover (Individual *parent1, Individual *parent2);
};

class SimulatedBinaryCrossover: public Crossover {
private:
  double genewiseProbability;
public:
  SimulatedBinaryCrossover(void);
  SimulatedBinaryCrossover(double geneProb);
  void crossover (Individual *parent1, Individual *parent2);
};

#endif
