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

#include "../include/selection.hpp"

// Constructor for the base selection class.
// This has to be called explicitly from any
// derived class constructor since the default
// will be overridden

Selection::Selection (const Population *Pop) {
  pop=Pop;
}

// Now for the derived class cosntructors

TournamentSelection::TournamentSelection (
    const int size, const Population *pop) :
    Selection(pop) , tournamentSize(size) {
}

TournamentSelectionWithReplacement::TournamentSelectionWithReplacement (
    const int size, const Population *pop) :
    Selection(pop) , tournamentSize(size) {
}

StochasticUniversalSelection::StochasticUniversalSelection(
    const Population *pop) :
    Selection(pop) {
}

RouletteWheelSelection::RouletteWheelSelection (const Population *pop) :
    Selection(pop) {
}

TruncationSelection::TruncationSelection(const int selP, const Population *pop) :
    Selection(pop), selectionPressure(selP) {
}

/***
    Tournament Selection Without Replacement: Randomly select S individuals
    without replacement and copy the best individual to the mating pool. Repeat
    this process till N individuals are selected.
    Implementation:
    1. Shuffle the Individuals of the current population
    2. Select S individuals without replacement
    3. Copy the index of the best individual among those S individuals
    to the mating pool
    4. Go to step 2 till all the individuals are chosen
    5. If the number of individuals in mating pool is less than N
    then go to step 1.

    Note: Tournament selection without replacement should be preferred over
    tournament selection with replacement.
*/
void TournamentSelection::select(int *matingPool) {
  int ii, jj, kk, ll = 0;
  int p1, p2, *randomArray, *winner;

  randomArray = new int[globalSetup->populationSize];
  winner = new int[globalSetup->populationSize];

  for(ii = 0; ii < globalSetup->populationSize; ii++)
    randomArray[ii] = ii;

  for(ii = 0; ii < tournamentSize; ii++) {
    myRandom.shuffleArray(randomArray, globalSetup->populationSize);

    for(jj = 0; jj < globalSetup->populationSize; jj += tournamentSize) {
      p1 = randomArray[jj];

      for(kk = 1; kk < tournamentSize; kk++) {
        p2 = randomArray[jj+kk];

        if(betterIndividual((*pop)[p1],(*pop)[p2]))
          winner[ll] = p1;
        else
          winner[ll] = p2;

        p1 = winner[ll];
      }
      ++ll;
    }
  }

  for(ii = 0; ii < globalSetup->populationSize; ii++)
    matingPool[ii] = winner[ii];
  delete []randomArray;
  delete []winner;
}

/***
    Tournament Selection With Replacement: Randomly select S individuals with
    replacement and copy the best individual to the mating pool. Repeat this
    process till N individuals are selected.
    Implementation:
    1. Shuffle the Individuals of the current population
    2. Select S individuals with replacement
    3. Copy the index of the best individual among those S individuals
    to the mating pool
    4. Go to step 1 till the mating pool has N individuals.

    Note: Tournament selection without replacement should be preferred over
    tournament selection with replacement.
*/
void TournamentSelectionWithReplacement::select(int *matingPool) {
  int ii, jj, p1, p2;
  int *winner;

  winner = new int[globalSetup->populationSize];

  for(ii = 0; ii < globalSetup->populationSize; ii++) {
    p1 = myRandom.boundedIntegerRandom(0,globalSetup->populationSize);
    for(jj = 1; jj < tournamentSize; jj++) {
      do{ p2 = myRandom.boundedIntegerRandom(0,globalSetup->populationSize); }
      while(p2 == p1);
      if(betterIndividual((*pop)[p1],(*pop)[p2])) winner[ii] = p1;
      else winner[ii] = p2;
      p1 = winner[ii];
    }
  }
  for(ii = 0; ii < globalSetup->populationSize; ii++)
    matingPool[ii] = winner[ii];
  delete []winner;
}

/***
    Truncation Selection: Sort the individuals in the population according to
    the fitness (or some performace criteria) and allocate S copies to the top
    N/S individuals.
*/
void TruncationSelection::select(int *matingPool) {
  int ii, jj, kk, *randomArray, *winner;

  randomArray = new int[globalSetup->populationSize];
  winner = new int[globalSetup->populationSize];

  for(ii = 0; ii < globalSetup->populationSize; ii++) randomArray[ii] = ii;
  selectionQuickSort(pop, randomArray, 0, globalSetup->populationSize);
  for(ii = 0, kk = 0; kk < globalSetup->populationSize; ii++)
    for(jj = 0; jj < selectionPressure; jj++)
      winner[kk++] = randomArray[globalSetup->populationSize-1-ii];
  for(ii = 0; ii < globalSetup->populationSize; ii++)
    matingPool[ii] = winner[ii];

  delete []randomArray;
  delete []winner;
}

/***
    Roulette Wheel Selection: Each individual is assigned a slice proportional
    to its fitness. The wheel is spun N times, where N is the number of
    individuals in the population.

    Note: Ordinal selection schemes like tournament or trunction selection
    schemes should always be preferred over roulette wheel selection
    Note: Stochastic universal selection should be preferred over roulette wheel
    selection.
    Note: Roulette wheel selection should always be used with some scaling
    method.
*/
void RouletteWheelSelection::select(int *matingPool) {
  int ii, jj, *winner;
  double rndNo, *sumFit;

  winner = new int[globalSetup->populationSize];
  sumFit = new double[globalSetup->populationSize];

  sumFit[0] = pop->getFitness(0);
  for(ii = 1; ii < globalSetup->populationSize; ii++)
    sumFit[ii] = sumFit[ii-1] + pop->getFitness(ii);
  for(ii = 0; ii < globalSetup->populationSize; ii++) {
    rndNo = myRandom.boundedRandom(0,sumFit[globalSetup->populationSize-1]);
    for(jj = 0; jj < globalSetup->populationSize; jj++) {
      if(sumFit[jj] >= rndNo) { winner[ii] = jj; break; }
    }
  }
  for(ii = 0; ii < globalSetup->populationSize; ii++)  matingPool[ii] =
                                                           winner[ii];
  delete []winner;
  delete []sumFit;
}

/***
    Stochasting Universal Selection: It is a method for selecting a population
    according to some given probability in a way that minimizes chance
    fluctuations associated with Roulette Wheel Selection. Instead of spinning
    the roulette wheel N times to select N parents, SUS spins the wheel once,
    with N equally spaced pointers, which are used to select N individuals.
    Under this method, each individual is guaranteed to reporduce at least
    floor(Expected number of copies of the individual) and at most
    ceil(Expected number of copies of the individual). SUS also has the
    advantage of being easier and quicker to implement that roulette wheel
    sampling. In practice it should always be preferred over roulette wheel
    selection.
    Note: This should almost always be used with some scaling method.
    Note: Currently this works only with SGA.

    Reference: Baker, J.E. (1985). "Adaptive Selection Methods for Genetic
    Algorithms", In Grefenstte, J. (Ed.), Proceedings of the International
    Conference on Genetic Algorithms and Their Applications (pp. 101-111).
    Hillsdale, NJ:Lawrence Erlbaum Associates (TCGA No. 00460).
*/
void StochasticUniversalSelection::select(int *matingPool) {
  int ii, jj = 0, *winner;
  double rndNo, sum, avgRank = 0.0;

  winner = new int[globalSetup->populationSize];

  rndNo = myRandom.random01();
  for(sum = 0.0, ii = 0; ii < globalSetup->populationSize; ii++)
    for(sum += (pop->getFitness(ii))/(pop->getAvgFit()); sum > rndNo; rndNo++)
      winner[jj++] = ii;
  for(ii = 0; ii < globalSetup->populationSize; ii++)
    matingPool[ii] = winner[ii];

  delete []winner;
}

/***
    Returns 1 if individual 1 is better than individual 2 depending on
    constraint handling method and the GA type, else returns 0.
    Methods Invoked: constrTournComp(Individual*, Individual*, compareWhat)
    crowdingCom(NsgaIndvidual*, NsgaIndividual*);
    Details: If the GA type is SGA then
    * If the constraint handling method is either penalty or if
    there are no constraints then return 1 if individual 1
    has higher fitness compared to individual 2, otherwise
    return 0.
    * If the constraint handling method is tournament selection,
    then return 1 if the function compTournComp returns 1, else
    return 0.
    If the GA type is NSGA then, return 1 if crowdingComp returns 1,
    otherwise return 0.
*/
int Selection::betterIndividual(Individual *guy1, Individual *guy2) {
  if(globalSetup->gaType == SGA) {
    switch(globalSetup->constraintMethod) {
      case NoConstraints: case Penalty:
        if(guy1->getFitness() > guy2->getFitness()) return 1;
        else return 0;
        break;
      case Tournament:
        if(constrTournComp(guy1,guy2, FITNESS)) return 1;
        else return 0;
        break;
      default: exit(0);
    };
  }
  else if(globalSetup->gaType == NSGA) {
    if(crowdingComp((NsgaIndividual*)guy1,(NsgaIndividual*)guy2)) return 1;
    else return 0;
  }
  return 0;
}

///This is a quick sort routine used for truncation selection
void Selection::selectionQuickSort(const Population *pop, int *output,
                                   int left, int right) {
  int ii, xx;
  Individual *target;

  if(right > left) {
    (target) = (*pop)[output[right-1]];
    ii = left-1;
    xx = right-1;
    for(;;) {
      while(!betterIndividual((*pop)[output[++ii]],target))
        if(ii >= right-1) break;
      if(ii >= xx) break;
      while(betterIndividual((*pop)[output[--xx]],target))
        if(xx <= 0) break;
      if(ii >= xx) break;
      swap(output[ii],output[xx]);
    }
    swap(output[ii],output[right-1]);
    selectionQuickSort(pop, output, left, ii);
    selectionQuickSort(pop, output, ii+1, right);
  }
}

///This is a swap routine used by the quick sort subroutine
void Selection::swap(int& ii, int& jj) {
  int temp;
  temp = jj; jj = ii; ii = temp;
}
