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

#ifndef _POPULATION_H
#define _POPULATION_H
class GA;                 // Predeclaration for making friends
class Selection;
class Crossover;
class Population;

#include <assert.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "random.hpp"
#include "globalSetup.hpp"
#include "individual.hpp"
#include "selection.hpp"
#include "crossover.hpp"
#include "localsearch.hpp"

extern GlobalSetup *globalSetup;
extern Random myRandom;

class Population
{
 protected:
  Individual **guys;              // The chromosomes
  Individual **newGuys;           // Individuals after crossover and mutation
  int        *mpool;              // mating pool
  Individual *bestInd;            // individual in the population

  Selection *selection;           // Pointer to selection function
  Crossover *crossover;           // Pointer to xover function
  LocalSearch *localSearch;
  int *freezeMask;

  int	     noOfFeasible;  	  // number of feasible individuals
  double     *bestobj;              // maximum objective
  double     *worstobj;              // minimum objective
  double     *avgobj;              // average objective
  double     *maxfit;              // maximum fitness
  double     *minfit;              // minimum fitness
  double     *avgfit;              // average fitness
  double     *varfit;
  double     *bestFitChange;
  double     *avgFitChange;
  double     *fitVarChange;
  double     *bestObjChange;
  double     *avgObjChange;

 public:

  // Big 3
  Population();
  virtual ~Population();

  // Declaration of friendship!
  friend class GA;

  // Interface functions
  inline int    getNoOfFeasible(void) const { return noOfFeasible; }
  inline double getMaxObj(void) const { return *bestobj; }
  inline double getMinObj(void) const { return *worstobj; }
  inline double getAvgObj(void) const { return *avgobj; }
  inline double getMaxFit(void) const { return *maxfit; }
  inline double getMinFit(void) const { return *minfit; }
  inline double getAvgFit(void) const { return *avgfit; }
  inline double getFitVar(void) const { return *varfit; }
  inline double getAvgObjChange(void) const { return *avgObjChange; }
  inline double getBestObjChange(void) const { return *bestObjChange; }
  inline double getFitVarChange(void) const { return *fitVarChange; }
  inline double getAvgFitChange(void) const { return *avgFitChange; }
  inline double getBestFitChange(void) const { return *bestFitChange; }
  
  inline void setNoOfFeasible(int iValue) { noOfFeasible = iValue; }
  inline void setMaxObj(double dValue) { *bestobj = dValue; }
  inline void setMinObj(double dValue) { *worstobj = dValue; }
  inline void setAvgObj(double dValue) { *avgobj = dValue; }
  inline void setMaxFit(double dValue) { *maxfit = dValue; }
  inline void setMinFit(double dValue) { *minfit = dValue; }
  inline void setAvgFit(double dValue) { *avgfit = dValue; }
  inline void setFitVar(double dValue) { *varfit = dValue; }
  inline void setAvgObjChange(double dValue) { *avgObjChange = dValue; }
  inline void setBestObjChange(double dValue) { *bestObjChange = dValue; }
  inline void setFitVarChange(double dValue) { *fitVarChange = dValue; }
  inline void setAvgFitChange(double dValue) { *avgFitChange = dValue; }
  inline void setBestFitChange(double dValue) { *bestFitChange = dValue; }

  inline double getFitness(const int index) const {
    return guys[index]->getFitness();
  }
  Individual *getBestIndividual(void) {
    return bestInd;
  }

  inline int*   getMPool(void) { return mpool; }
  inline int*	getFreezeMask(void) { return freezeMask; }

  //Common GA operations called on each and every individual
  int doEvaluate(void);
  void loadPopulationFromFile(void);
  void doSelect(void);
  void doCrossover(void);
  void doMutate(void);
  int doLocalSearch(void);
  void replacePopulation(void);
  void mapObjectiveToFitness(void);
  void computeObjStatistics(void);
  void computeFitnessStatistics(void);
  void scaleFitness(void);
  void rankingQuickSort(int *output, int left, int right);
  void swap(int& ii, int& jj);

  //Freeze and flood operations
  void	freeze(int locus);
  void	flood( int locus);

  // Some operations
  Individual *operator[]( int index ) const {
    return guys[index];
  }

  friend std::ostream &operator<< ( std::ostream &out, const Population &pop );

  //Some GA operations defined in Population
  void shareFitness(void);
  void doRTS(void);

};

class NsgaPopulation : public Population {

 protected:
  NsgaIndividual **combinedGuys;
  int *numIndsFront;
  int **paretoFront;
  int numFronts;
  int numFrontChange;
 public:

  // The Big 3
  NsgaPopulation(void);
  NsgaPopulation(const NsgaPopulation &sourcePop);
  ~NsgaPopulation(void);
  double getCrowdingDistance(int index) {
    return combinedGuys[index]->getCrowdingDistance();
  }
  // I think friendship is not derived - so I am redeclaring
  friend class GA;
  void doNonDominatedSort(int whichGuys);
  void computeCrowdingDistance(int whichGuys);
  void quickSort(NsgaIndividual **theGuys, int *output, int left, int right,
                 int objID);
  void regQSort(double *crowdingDistance, int *output, int left, int right);
  void swap(int &ii, int &jj);
  void computeObjStatistics(int whichGuys);
  void computeFitnessStatistics(int whichGuys);
  void mapObjectiveToFitness(int whichGuys);

  inline int getNoOfIndsBestFront(void) const { return numIndsFront[0]; }

  inline int* getNoOfIndsFront(void) { return numIndsFront; }
  inline int** getParetoFront(void) { return paretoFront; }
  inline NsgaIndividual** getCombinedGuys(void) { return combinedGuys; }
  
  inline int getNoOfFronts(void) const { return numFronts; }
  inline int getNoOfFrontChange(void) const { return numFrontChange; }
  inline double getMaxObj(int index) const { return bestobj[index]; }
  inline double getMinObj(int index) const { return worstobj[index]; }
  inline double getAvgObj(int index) const { return avgobj[index]; }
  inline double getMaxFit(int index) const { return maxfit[index]; }
  inline double getMinFit(int index) const { return minfit[index]; }
  inline double getAvgFit(int index) const { return avgfit[index]; }
  inline double getFitVar(int index) const { return varfit[index]; }
  
  inline void setNoOfFronts(int iValue) { numFronts = iValue; }
  inline void setNoOfFrontChange(int iValue) { numFrontChange = iValue; }
  inline void setMaxObj(int index, double dValue) { bestobj[index] = dValue; }
  inline void setMinObj(int index, double dValue) { worstobj[index] = dValue; }
  inline void setAvgObj(int index, double dValue) { avgobj[index] = dValue; }
  inline void setMaxFit(int index, double dValue) { maxfit[index] = dValue; }
  inline void setMinFit(int index, double dValue) { minfit[index] = dValue; }
  inline void setAvgFit(int index, double dValue) { avgfit[index] = dValue; }
  inline void setFitVar(int index, double dValue) { varfit[index] = dValue; }
};
#endif
