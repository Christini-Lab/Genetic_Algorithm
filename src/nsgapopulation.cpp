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

#include "../include/population.hpp"

NsgaPopulation::NsgaPopulation(void) {
  int ii, combinedPopSize = 2*(globalSetup->populationSize);
  numFronts = 0;
  numFrontChange = 0;
  combinedGuys = new NsgaIndividual*[combinedPopSize];
  numIndsFront = new int[combinedPopSize];
  paretoFront = new int*[combinedPopSize];
  for(ii = 0; ii < combinedPopSize; ii++) {
    combinedGuys[ii] = new NsgaIndividual;
    paretoFront[ii] = new int[combinedPopSize];
  }
  bestobj = new double[globalSetup->finalNoOfObjectives];
  worstobj = new double[globalSetup->finalNoOfObjectives];
  avgobj = new double[globalSetup->finalNoOfObjectives];
  maxfit = new double[globalSetup->finalNoOfObjectives];
  minfit = new double[globalSetup->finalNoOfObjectives];
  avgfit = new double[globalSetup->finalNoOfObjectives];
  varfit = new double[globalSetup->finalNoOfObjectives];
  computeObjStatistics(GUYS);
  mapObjectiveToFitness(GUYS);
  computeFitnessStatistics(GUYS);
}

NsgaPopulation::~NsgaPopulation(void) {
  int ii, combinedPopSize = 2*(globalSetup->populationSize);
  for(ii = 0; ii < combinedPopSize; ii++) {
    delete combinedGuys[ii];
    delete paretoFront[ii];
  }
  delete []combinedGuys;
  delete []numIndsFront;
  delete []paretoFront;
  delete []bestobj;
  delete []worstobj;
  delete []avgobj;
  delete []minfit;
  delete []maxfit;
  delete []avgfit;
  delete []varfit;
}


/***
    Peforms a nondominated sort of individuals to create a set of pareto fronts.
    Details:
    *  The flag whichGuys refer to either the current population (for the
    initial generation) or the combined population.
    *  domCount[i]: (domination count) Number of the individuals that dominate
    individual i
    *  numIndDom[i]: Number of individuals that individual i dominates
    *  indDomByMe[i]: Array of indices of individuals that individual i
    dominates
    *  paretoFront[i]: Array of indices of individuals that are members of p
    areto front i.
    *  numIndsFront[i]: Number of individuals in the ith front.
    *  For every unique pair of individuals (i,j)
    --Check if i dominates j, If yes then add j to indDomByMe[i], increment
    numIndDom[i] by 1 and increment domCount[j] by 1.
    --On the other hand if j dominates i, then add i to indDomByMe[j],
    increment numIndDom[j] by 1 and increment domCount[i] by 1.
    * For every individual i
    --If domCount[i] equals 0, then add individual i to paretoFront[0],
    increment numIndsFront[0], and set the rank of individual i to 0
    * Set frontID = 0. While number of individuals in the pareto Front is
    greater than zero, do
    --For every individual i in the pareto front frontID
    *  For every individual j that individual i dominates
    --Decrement the domination count of j by i
    --If the domination count of j is equal to zero, then add j to the
    paretoFron[frontID+1], increment numIndsFront[frontID+1], and set
    the rank of individual j to frontID+1.
    --Increment the frontID by one
    * Set total number of pareto fronts to frontID.
    */
void NsgaPopulation::doNonDominatedSort(int whichGuys) {
  int ii, jj, popSize, fid, idomj, nif;
  int *domCount, **indDomByMe, *numIndDom;
  int oldNoOfFronts;
  NsgaIndividual **theGuys;

  oldNoOfFronts = numFronts;

  if(whichGuys == GUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)guys;
  }
  else if(whichGuys == COMBINEDGUYS) {
    popSize = 2*(globalSetup->populationSize);
    theGuys = combinedGuys;
  }
  domCount = new int[popSize];
  numIndDom = new int[popSize];
  indDomByMe = new int*[popSize];
  for(ii = 0; ii < popSize; ii++) {
    indDomByMe[ii] = new int[popSize];
    domCount[ii] = numIndDom[ii] = 0;
  }
  for(ii = 0; ii < popSize-1; ii++) {
    for(jj = ii+1; jj < popSize; jj++) {
      idomj = dominates(theGuys[ii],theGuys[jj]);
      if(idomj == IDOMJ) {
	indDomByMe[ii][numIndDom[ii]] = jj;
	numIndDom[ii] += 1;
	domCount[jj] += 1;
      }
      else if(idomj == JDOMI) {
	indDomByMe[jj][numIndDom[jj]] = ii;
	numIndDom[jj] += 1;
	domCount[ii] += 1;
      }
    }
  }
  nif = 0;

  for(ii = 0; ii < popSize; ii++) {
    if(domCount[ii] == 0) {
      paretoFront[0][nif++] = ii;
      theGuys[ii]->setRank(0);
    }
  }
  numIndsFront[0] = nif;

  fid = 0;
  while(numIndsFront[fid] != 0) {
    nif = 0;
    for(ii = 0; ii < numIndsFront[fid]; ii++) {
      for(jj = 0; jj < numIndDom[paretoFront[fid][ii]]; jj++) {
	domCount[indDomByMe[paretoFront[fid][ii]][jj]] -= 1;
	if(domCount[indDomByMe[paretoFront[fid][ii]][jj]] == 0) {
	  paretoFront[fid+1][nif++] = indDomByMe[paretoFront[fid][ii]][jj];
	  theGuys[indDomByMe[paretoFront[fid][ii]][jj]]->setRank(fid+1);
	}
      }
    }
    fid += 1;
    numIndsFront[fid] = nif;
  }
  numFronts = fid;
  numFrontChange = abs(oldNoOfFronts - numFronts);
  if(whichGuys == GUYS)
    for(ii = 0; ii < popSize; ii++) guys[ii] = theGuys[ii];
  else if(whichGuys == COMBINEDGUYS)
    for(ii = 0; ii < popSize; ii++) combinedGuys[ii] = theGuys[ii];

  for(ii = 0; ii < popSize; ii++)
    delete indDomByMe[ii];
  delete []domCount;
  delete []numIndDom;
  delete []indDomByMe;
}

/***
    Computes the crowding distance of each individual.
    Details:
    * The flag whichGuys refer to either the current population (for the initial
    generation) or the combined population.
    * For every front i
    -- For every objective j
    * Sort the individuals in front i in ascending order of objective j
    * Set the crowding distance of the first and the last individual to infinity
    * For every individual k in the front i except the first and the last
    individuals
    --Normalize the fitness j by dividing the fitness j of individual k-1  and
    individual k+1 by maximum jth fitness.
    --Add absolute value of (Normalized jth fitness of individual k+1 -
    Normalized jth fitness of individual k-1) to the crowding distance of kth
    individual.
*/
void NsgaPopulation::computeCrowdingDistance(int whichGuys) {
  int ii, jj, kk, firstInd, lastInd;
  int indId1, indId2, indId, popSize;
  int *sortListByObj;
  double normFit1, normFit2, *crowdDist;
  NsgaIndividual **theGuys;

  if(whichGuys == GUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)guys;
  }
  else if(whichGuys == COMBINEDGUYS) {
    popSize = 2*(globalSetup->populationSize);
    theGuys = combinedGuys;
  }

  crowdDist = new double[popSize];
  for(ii = 0; ii < popSize; ii++) crowdDist[ii] = 0.0;

  for(ii = 0; ii < numFronts; ii++) {
    sortListByObj = new int[numIndsFront[ii]];
    for(jj = 0; jj < globalSetup->finalNoOfObjectives; jj++) {
      for(kk = 0; kk < numIndsFront[ii]; kk++)
	sortListByObj[kk] = paretoFront[ii][kk];
      quickSort(theGuys, sortListByObj, 0, numIndsFront[ii],jj);
      firstInd = sortListByObj[0];
      lastInd = sortListByObj[numIndsFront[ii]-1];
      crowdDist[firstInd] = crowdDist[lastInd] = INFTY;
      for(kk = 1; kk < numIndsFront[ii]-1; kk++) {
	indId = sortListByObj[kk];
	indId1 = sortListByObj[kk+1];
	indId2 = sortListByObj[kk-1];
	normFit1 = theGuys[indId1]->getFitness(jj)/(1.0+maxfit[jj]);
	normFit2 = theGuys[indId2]->getFitness(jj)/(1.0+maxfit[jj]);
	crowdDist[indId] += fabs(normFit1 - normFit2);
      }
    }
    delete []sortListByObj;
  }

  if(whichGuys == GUYS)
    for(ii = 0; ii < popSize; ii++)
      ((NsgaIndividual**)guys)[ii]->setCrowdingDistance(crowdDist[ii]);
  else if (whichGuys == COMBINEDGUYS)
    for(ii = 0; ii < popSize; ii++)
      combinedGuys[ii]->setCrowdingDistance(crowdDist[ii]);
  delete [] crowdDist;
}

/// Quick sort routine used for crowding distance computation. Sorts the
// individuals in a given pareto front in ascending order of objId th objective
// value.
void NsgaPopulation::quickSort(NsgaIndividual **theGuys, int *output,
			       int left, int right, int objId) {
  int ii, xx;
  double target;

  if(right > left) {
    target = theGuys[output[right-1]]->getFitness(objId);
    ii = left-1;
    xx = right-1;
    for(;;) {
      while(theGuys[output[++ii]]->getFitness(objId) < target)
	if(ii >= right-1) break;
      if(ii >= xx) break;
      while(theGuys[output[--xx]]->getFitness(objId) > target)
	if(xx <= 0) break;
      if(ii >= xx) break;
      swap(output[ii],output[xx]);
    }
    swap(output[ii],output[right-1]);
    quickSort(theGuys, output, left, ii, objId);
    quickSort(theGuys, output, ii+1, right, objId);
  }
}

///Swap used by the sort function
void NsgaPopulation::swap(int& ii, int& jj) {
  int temp;
  temp = jj; jj = ii; ii = temp;
}

///Quick sort function to sort individuals in a pareto front according to the
// crowding distance.
void NsgaPopulation::regQSort(double *crowdDist, int *output,
                              int left, int right) {
  int ii, xx;
  double target;
  if(right > left) {
    target = crowdDist[output[right-1]];
    ii = left-1;
    xx = right-1;
    for(;;) {
      while(crowdDist[output[++ii]] < target)
	if(ii >= right-1) break;
      if(ii >= xx) break;
      while(crowdDist[output[--xx]] > target)
	if(xx <= 0) break;
      if(ii >= xx) break;
      swap(output[ii],output[xx]);
    }
    swap(output[ii],output[right-1]);
    regQSort(crowdDist, output, left, ii);
    regQSort(crowdDist, output, ii+1, right);
  }
}

/*============================================================
** Function Name: NsgaPopulation::computeObjStatistics(int whichGuys)
** The flag whichGuys refer to either the current population (for the
** initial generation) or the combined population.
** Function Task: this method generates statistics according to
** each objective function values. It finds out the average objective
** value for the population and the best and worst objective
** values for each objective function. This depends on whether the
** type of optimization in the problem of each objective function
** is minimization or maximization.
** Output:None
** Functions called: Individual::getObjective(int jj)
**========================================================*/
void NsgaPopulation::computeObjStatistics(int whichGuys) {
  int ii, jj, popSize;
  double oneOverPopulationSize;
  NsgaIndividual **theGuys;

  if(whichGuys == GUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)guys;
  }
  else if(whichGuys == NEWGUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)newGuys;
  }
  else if(whichGuys == COMBINEDGUYS) {
    popSize = 2*(globalSetup->populationSize);
    theGuys = combinedGuys;
  }

  oneOverPopulationSize = 1.0/(1.0*popSize);

  for(jj = 0; jj < globalSetup->finalNoOfObjectives; jj++) {
    bestobj[jj] = worstobj[jj] = theGuys[0]->getObjective(jj);
    avgobj[jj] = oneOverPopulationSize*(theGuys[0]->getObjective(jj));
    for(ii = 1; ii < popSize; ii++) {
      avgobj[jj] += oneOverPopulationSize*(theGuys[ii]->getObjective(jj));
      if(globalSetup->typeOfOptimizations[jj] == Maximization) {
	if(theGuys[ii]->getObjective(jj) > bestobj[jj])
          bestobj[jj] = theGuys[ii]->getObjective(jj);
	if(theGuys[ii]->getObjective(jj) < worstobj[jj])
          worstobj[jj] = theGuys[ii]->getObjective(jj);
      }
      if(globalSetup->typeOfOptimizations[jj]==Minimization) {
	if(theGuys[ii]->getObjective(jj) < bestobj[jj])
          bestobj[jj] = theGuys[ii]->getObjective(jj);
	if(theGuys[ii]->getObjective(jj) > worstobj[jj])
          worstobj[jj] = theGuys[ii]->getObjective(jj);
      }
    }
  }
}

/*============================================================
** Function Name: NsgaPopulation::computeFitnessStatistics(int whichGuys)
** The flag whichGuys refer to either the current population
** (for the initial generation) or the combined population.
** Function Task: this method generates statistics according to
** every fitness function values. It finds out the average fitness
** value for the population and the maximum and minimum fitness
** values as well as the fitness variance for each fitness function.
** Output:None
** Functions called: Individual::getFitness(int jj)
**========================================================*/
void NsgaPopulation::computeFitnessStatistics(int whichGuys) {
  int ii, jj, popSize;
  double oneOverPopulationSize;
  NsgaIndividual **theGuys;

  if(whichGuys == GUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)guys;
  }
  else if(whichGuys == NEWGUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual**)newGuys;
  }
  else if(whichGuys == COMBINEDGUYS) {
    popSize = 2*(globalSetup->populationSize);
    theGuys = combinedGuys;
  }

  oneOverPopulationSize = 1.0/(1.0*popSize);

  for(jj = 0; jj < globalSetup->finalNoOfObjectives; jj++) {
    maxfit[jj] = minfit[jj] = theGuys[0]->getFitness(jj);
    avgfit[jj] = oneOverPopulationSize*(theGuys[0]->getFitness(jj));
    varfit[jj] = oneOverPopulationSize*(theGuys[0]->getFitness(jj)) *
        (theGuys[0]->getFitness(jj));

    for(ii = 1; ii < popSize; ii++) {
      avgfit[jj] += oneOverPopulationSize*(theGuys[ii]->getFitness(jj));
      varfit[jj] += oneOverPopulationSize*(theGuys[ii]->getFitness(jj)) *
          (theGuys[ii]->getFitness(jj));
      if(theGuys[ii]->getFitness(jj) > maxfit[jj])
        maxfit[jj] = theGuys[ii]->getFitness(jj);
      if(theGuys[ii]->getFitness(jj) < minfit[jj])
        minfit[jj] = theGuys[ii]->getFitness(jj);
    }
    varfit[jj] -= avgfit[jj]*avgfit[jj];
  }
}

/***
    Map each objective to corresponding fitness value depending on minimization
    or maximization. The flag whichGuys refer to either the current population
    (for the initial generation) or the combined population. If the jth
    objective value is a maximization, then the fitness is set equal to the
    objective value. On the other hand, if the problem is a minimization one,
    then set the fitness equal to the negative of the objective value.
*/
void NsgaPopulation::mapObjectiveToFitness(int whichGuys) {
  int ii,jj, popSize;
  NsgaIndividual **theGuys;

  if(whichGuys == GUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual **)guys;
  }
  else if(whichGuys == NEWGUYS) {
    popSize = globalSetup->populationSize;
    theGuys = (NsgaIndividual **)newGuys;
  }
  else if(whichGuys == COMBINEDGUYS) {
    popSize = 2*(globalSetup->populationSize);
    theGuys = combinedGuys;
  }

  for (jj=0; jj<globalSetup->finalNoOfObjectives; jj++) {
    if (globalSetup->typeOfOptimizations[jj] == Maximization)
      for(ii = 0; ii < popSize; ii++)
        theGuys[ii]->setFitness(jj, theGuys[ii]->getObjective(jj));
    else
      for (ii=0; ii< popSize; ii++)
	theGuys[ii]->setFitness(jj, -theGuys[ii]->getObjective(jj));
  }
}
