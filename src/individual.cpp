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

#include "../include/individual.hpp"
// Started Modification 02/14/2001:21:19:59 - PVP

extern GlobalSetup *globalSetup;
extern void globalEvaluate(double *values, double *objFunction,
			   double *violation, double *penalty,
			   int *noOfViolations);

Individual::Individual(void) {

  // This is the default constructor. Here, the properties
  // of the Individual will be set appropriately. It should be
  // noted that the default constructor of a parent class is
  // automatically called when a child class object is constructed.
  // The reason why I am saying this is that we have double *fitness
  // etc and we would have to initialize it in the NsgaIndividual
  // class constructor to an array of fitness, etc. Here, we should
  // take care that we do not allocate it again. Therefore, we have
  // to check what type of GA this is, before we allocate memory.

  // The things to be initialized:
  //    1.   Chromosome chrom; - just call default constructor - will be
  //         initialized to a random chromosome based on bounds
  //    2.   double *fitness; - initialize to new doubl if GA is SGA
  //    3.   double *objFunction; -          -do-
  //    4.   double penalty;  - initialize to ZERO.
  //    5.   double *violation, totalViolation; - Set violation to NULL if
  //         unconstrained.
  //    6.   int *freezeMask; - initialize to int[globalSetup->noOfVariables]
  //         and set everything to 0.

  // Do the evaluation after everything is done.

  // Chromosomes default constructor is called automatically.

  if (!(globalSetup->finalNoOfConstraints))
    violation = NULL;
  else
    violation = new double[globalSetup->finalNoOfConstraints];
  noOfViolations = 0;
  penalty = 0.0;
  if (globalSetup->gaType==SGA) {
    fitness = new double(0.00);
    objFunction = new double(0.00);
    // In order to parallelize fitness evaluation, evaluate fitness must now be
    // explictly called, see Population class default constructor
    // if(!(globalSetup->loadPopulation)) evaluateFitness(); // If Nsga, it's
    // constructor will call it's evaluate.
  }
}

Individual::Individual(const Individual &sourceInd ) {
  // Call operator= of chromosomes and just copy all other values after
  // appropriate allocations
  int ii;
  chrom = sourceInd.chrom;
  if (globalSetup->gaType==SGA) {
    fitness = new double(*(sourceInd.fitness));
    objFunction = new double(*(sourceInd.objFunction));
  }
  if (globalSetup->finalNoOfConstraints) {
    noOfViolations = sourceInd.noOfViolations;
    violation = new double[globalSetup->finalNoOfConstraints];
    for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++)
      violation[ii] = sourceInd.violation[ii];
    penalty = sourceInd.penalty;
  }
  else
    violation = NULL;

}

Individual::Individual(const Individual *sourceInd ) {
  int ii;
  chrom = sourceInd->chrom;
  if (globalSetup->gaType==SGA) {
    fitness = new double(*(sourceInd->fitness));
    objFunction = new double(*(sourceInd->objFunction));
  }
  if (globalSetup->finalNoOfConstraints) {
    noOfViolations = sourceInd->noOfViolations;
    violation = new double[globalSetup->finalNoOfConstraints];
    for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++)
      violation[ii] = sourceInd->violation[ii];
    penalty = sourceInd->penalty;
  }
  else
    violation = NULL;

}

Individual & Individual::operator= (const Individual &sourceInd) {
  // Difference between copy constructor and this is
  // that this guy has already been allocated some memory
  // for the dynamically allocated stuff and so it has to be freed
  // and copied or copied without reallocation.

  int ii;
  chrom = sourceInd.chrom;
  if (globalSetup->gaType==SGA) {
    *fitness = (*(sourceInd.fitness));
    *objFunction = (*(sourceInd.objFunction));
  }
  if(globalSetup->finalNoOfConstraints) {
    noOfViolations = sourceInd.noOfViolations;
    for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++)
      violation[ii] = sourceInd.violation[ii];
    penalty = sourceInd.penalty;
  }
  return *this;
}

void Individual::loadIndividual(double *varValues, double *objValues,
                                double *constViolValues, double penaltyValue) {
  int ii;
  for(ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
    chrom[ii] = varValues[ii];
  for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++)
    objFunction[ii] = objValues[ii];
  if(globalSetup->finalNoOfConstraints) {
    noOfViolations = 0;
    for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++)
      violation[ii] = constViolValues[ii];
    if(violation[ii] >= 1.0e-5) {
      noOfViolations++;
    }
    penalty = penaltyValue;
  }
}

void Individual::evaluateFitness(void) {
  double *values;
  int ii;

  values = new double[globalSetup->noOfDecisionVariables];

  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
    values[ii] = chrom[ii];
  globalEvaluate(values, objFunction, violation,
		 &penalty, &noOfViolations);
  delete [] values;
}

void Individual::mutate(int *freezeMask) {
  switch(globalSetup->mutationType) {
    case Selective:
      chrom.mutateMinMax(freezeMask);
      break;
    case Genewise:
      chrom.mutateNormal(freezeMask);
      break;
    case Polynomial:
      chrom.mutatePolynomial(freezeMask);
      break;
    default: exit(0);
  }
}


Individual::~Individual(void) {
  if (globalSetup->gaType==SGA) {
    delete fitness;
    delete objFunction;
  }
  if (globalSetup->finalNoOfConstraints)
    delete [] violation;
}

// Some stuff for NsgaIndividual

NsgaIndividual::NsgaIndividual(void) {
  // At this point, chrom is ready, violation is ready and freezeMask is there
  // even though it is not necessary. We have to initialize fitness and
  // objFunction and do the evaluation

  fitness = new double[globalSetup->finalNoOfObjectives];
  objFunction = new double[globalSetup->finalNoOfObjectives];
  rank = 0;
  crowdingDistance=0.0;
  if(!(globalSetup->loadPopulation)) evaluateFitness();
}

NsgaIndividual::NsgaIndividual(const Individual &sourceInd):
    Individual(sourceInd) {
  int ii;
  NsgaIndividual *newSource = (NsgaIndividual *) &sourceInd;
  fitness = new double[globalSetup->finalNoOfObjectives];
  objFunction = new double[globalSetup->finalNoOfObjectives];

  for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = newSource->fitness[ii];
    objFunction[ii] = newSource->objFunction[ii];
  }
  rank = newSource->rank;
  crowdingDistance = newSource->crowdingDistance;
}

NsgaIndividual::NsgaIndividual(const Individual *sourceInd):
    Individual(*sourceInd){
  int ii;

  NsgaIndividual *newSource = (NsgaIndividual *) sourceInd;
  fitness = new double[globalSetup->finalNoOfObjectives];
  objFunction = new double[globalSetup->finalNoOfObjectives];
  for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = newSource->fitness[ii];
    objFunction[ii] = newSource->objFunction[ii];
  }
  rank = newSource->rank;
  crowdingDistance = newSource->crowdingDistance;
}
NsgaIndividual::NsgaIndividual(const NsgaIndividual &sourceInd):
    Individual(sourceInd) {
  int ii;
  fitness = new double[globalSetup->finalNoOfObjectives];
  objFunction = new double[globalSetup->finalNoOfObjectives];
  for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = sourceInd.fitness[ii];
    objFunction[ii] = sourceInd.objFunction[ii];
  }
  rank = sourceInd.rank;
  crowdingDistance = sourceInd.crowdingDistance;
}

NsgaIndividual::NsgaIndividual(const NsgaIndividual *sourceInd):
    Individual(*sourceInd){


  int ii;
  fitness = new double[globalSetup->finalNoOfObjectives];
  objFunction = new double[globalSetup->finalNoOfObjectives];
  for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = sourceInd->fitness[ii];
    objFunction[ii] = sourceInd->objFunction[ii];
  }
  rank = sourceInd->rank;
  crowdingDistance = sourceInd->crowdingDistance;
}


NsgaIndividual::~NsgaIndividual(void) {
  // Here Individual's destructor is called automatically.
  // No need to call Individual::~Individual explicitly.
  delete [] fitness;
  delete [] objFunction;
}

NsgaIndividual & NsgaIndividual::operator= (const NsgaIndividual &sourceInd) {

  int ii;
  *((Individual *)this) = sourceInd;
  for (ii=0; ii<globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = sourceInd.fitness[ii];
    objFunction[ii] = sourceInd.objFunction[ii];
  }
  rank = sourceInd.rank;
  crowdingDistance = sourceInd.crowdingDistance;
  return *this;
}

NsgaIndividual & NsgaIndividual::operator= (const Individual &sourceInd){

  int ii;
  Individual::operator=(sourceInd);
  NsgaIndividual *newSource = (NsgaIndividual *)&sourceInd;
  for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    fitness[ii] = newSource->fitness[ii];
    objFunction[ii] = newSource->objFunction[ii];
  }
  rank = newSource->rank;
  crowdingDistance = newSource->crowdingDistance;
  return *this;
}

std::ostream &operator << (std::ostream &out, Individual &x) {
  out << x.chrom;
  out << "Objective Function value = " << *(x.objFunction) << std::endl;
  for (int i=0; i<globalSetup->finalNoOfConstraints; i++)
    out << "Constraint Violation #" << i << " " << x.violation[i] << std::endl;
  return out;
}

/***
    Crowding comparison operator: Determines if NsgaIndvidual 1 is better than
    NsgaIndividual 2. If the rank of individual 1 is less than that of
    individual 2, then return 1. If the rank of individual 1 and 2 are same and
    the crowding distance of individual 1 is greater than that of individual 2
    then return 1. For all other cases return 0.
*/
bool crowdingComp (const NsgaIndividual *guy1, const NsgaIndividual *guy2) {
  if((guy1->rank < guy2->rank)||
     ((guy1->rank == guy2->rank)&&
      (guy1->crowdingDistance > guy2->crowdingDistance))) return 1;
  else return 0;
}

/***
    Constraint tournament comparison operator: Determines if the individual 1 is
    better than individual 2. The flag compareWhat selects whether the objective
    alue or the fitness is to be compared. If objective value is to be compared
    and we are required to maximise the objective function then the fitness is
    set to be equal to the objective value. On the other hand if the objective
    function is to be minimized, then the fitness is set to negative of
    objective function. If both the individuals are feasible or if they have the
    same amount of infeasibility (quantified by total constraint violation),
    then return 1 if the fitness of individual 1 is greater than that of
    individual 2. On the other hand if individual 1 is feasible and individual 2
    is infeasible then return 1. If both individuals are infeasible but
    individual 1 has lower infeasibility (lower constraint violation value) then
    return 1. For all other cases return 0
*/
bool constrTournComp (const Individual *guy1, const Individual *guy2,
                      int compareWhat)
{
  double guy1Fitness, guy2Fitness;

  switch(compareWhat) {
    case FITNESS:
      guy1Fitness = *(guy1->fitness);
      guy2Fitness = *(guy2->fitness);
      break;
    case OBJECTIVE:
      if(*(globalSetup->typeOfOptimizations) == Maximization) {
        guy1Fitness = *(guy1->objFunction);
        guy2Fitness = *(guy2->objFunction);
      }
      else {
        guy1Fitness = -(*(guy1->objFunction));
        guy2Fitness = -(*(guy2->objFunction));
      }
      break;
    default: exit(0);
  }
  if((guy1->penalty <= ZERO) && (guy2->penalty > ZERO))
    return 1;
  else if(((guy1->penalty <= ZERO)&&(guy2->penalty <= ZERO))||
	  (fabs(guy1->penalty-guy2->penalty) <= ZERO)) {
    if(guy1Fitness > guy2Fitness) return 1;
  }
  else
    if(guy1->penalty < guy2->penalty) return 1;
  return 0;
}

/**
   Domination determinator. Determines if NsgaIndvidual 1 dominates
   NsgaIndvidual 2.
   If the constraint handling method is either penalty method or if there are no
   constraints, then
   -- If every fitness value of individual 1 is greater than or equal to the
      fitness value of individual 2 then return 1.
   -- If every fitness value of individual 2 is greater than or equal to the
      fitness value of individual 2 then return -1.
   -- For all other cases return 0.
      If the constraint handling method is through tournament selection, then
   -- If both individuals are feasible or have the same amount of infeasibility,
      then if every fitness value of individual 1 is greater than or equal to
      the fitness value of individual 2 then return 1.
   -- If indiviudal 1 is feasible and individual 2 is infeasible return 1.
   -- If both individuals are infeasible, but the infeasiblity of individual 1
      is lesser than individual 2 then return 1
   -- If both individuals are feasible or have the same amount of infeasibility,
      then if every fitness value of individual 2 is greater than or equal to
      the fitness value of individual 1 then return -1.
   -- If indiviudal 2 is feasible and individual 1 is infeasible then return -1
   -- If both individuals are infeasible, but the infeasiblity of individual 2
      is lesser than individual 1 then return -1.
   -- For all other cases return 0.

*/
int dominates (const NsgaIndividual *guy1, const NsgaIndividual *guy2) {
  int ii, idomj = 0, jdomi = 0;

  switch(globalSetup->constraintMethod) {
    case NoConstraints:
      for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
        if(guy1->fitness[ii] > guy2->fitness[ii]) idomj = 1;
        else if(guy1->fitness[ii] < guy2->fitness[ii]) jdomi = 1;
      }
      break;
    case Penalty:
      for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
        if(guy1->fitness[ii] > guy2->fitness[ii]) idomj = 1;
        else if(guy1->fitness[ii] < guy2->fitness[ii]) jdomi = 1;
      }
      break;
    case Tournament:
      if((guy1->penalty <= ZERO)&&(guy2->penalty <= ZERO)) {
        for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
          if(guy1->fitness[ii] > guy2->fitness[ii]) idomj = 1;
          else if(guy1->fitness[ii] < guy2->fitness[ii]) jdomi = 1;
        }
      }
      else if(((guy1->penalty <= ZERO)&&(guy2->penalty <= ZERO))||
              (fabs(guy1->penalty-guy2->penalty) <= ZERO)) {
        for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
          if(guy1->fitness[ii] > guy2->fitness[ii]) idomj = 1;
          else if(guy1->fitness[ii] < guy2->fitness[ii]) jdomi = 1;
        }
      }
      else {
        if(guy1->penalty < guy2->penalty) idomj = 1;
        else if(guy2->penalty < guy1->penalty) jdomi = 1;
      }
      break;
    default: exit(0);
  }
  if((idomj == 1)&&(jdomi == 1)) return INONDOMJ;
  else if((idomj == 1)&&(jdomi == 0)) return IDOMJ;
  else if((idomj == 0)&&(jdomi == 1)) return JDOMI;
  else return INONDOMJ;
}


// Penaly constraint comparison operator: Checks if indiviudal 1 is better than
// individual 2 if the constraint handling is through penalty method.
bool constrPenaltyComp(const Individual *guy1, const Individual *guy2) {
  if(*(globalSetup->typeOfOptimizations) == Maximization) {
    if(*(guy1->objFunction) - guy1->penalty >
       *(guy2->objFunction) - guy2->penalty) return 1;
  }
  else {
    if(*(guy1->objFunction) + guy1->penalty <
       *(guy2->objFunction) + guy2->penalty) return 1;
  }
  return 0;
}


// Individual comparison operator: Determines if individual 1 is better than
// individual 2 when no constraints are present.
bool individualComp(const Individual *guy1, const Individual *guy2) {
  if(*(globalSetup->typeOfOptimizations) == Maximization) {
    if(*(guy1->objFunction) > *(guy2->objFunction)) return 1;
  }
  else {
    if(*(guy1->objFunction) < *(guy2->objFunction)) return 1;
  }
  return 0;
}

bool isBetter(const Individual *guy1, const Individual *guy2) {
  bool status;

  switch(globalSetup->constraintMethod) {
    case NoConstraints:
      status = individualComp(guy1, guy2);
      break;
    case Penalty:
      status = constrPenaltyComp(guy1, guy2);
      break;
    case Tournament:
      status = constrTournComp(guy1, guy2, OBJECTIVE);
      break;
    default:
      exit(0);
  }
  return status;
}
