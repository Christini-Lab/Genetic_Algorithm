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


#include "../include/chromosome.hpp"

extern GlobalSetup *globalSetup;
extern Random myRandom;

// The Big 3

Chromosome::Chromosome(void) {
  // Initialize to random values within the limits of each variable
  // I am going to follow the same naming conventions as used previously
  // for random number generators although we will be using different
  // functions.
  int ii;

  genes = new double[globalSetup->noOfDecisionVariables];
  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
    genes[ii] =
        myRandom.boundedRandom(globalSetup->variableRanges[ii][0],
                               globalSetup->variableRanges[ii][1]);
    if (globalSetup->variableTypes[ii]==Integer)
      genes[ii] = int(genes[ii]+0.5);
  }
}

Chromosome::Chromosome(const Chromosome &chrom) {
  int ii;
  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
    genes[ii] = chrom[ii];

}

Chromosome::~Chromosome(void) {
  delete []genes;
}

std::ostream &operator<< (std::ostream &out, Chromosome &chrom) {
  int ii;
  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
    out<<chrom[ii]<<std::endl;
  return out;
}

Chromosome &Chromosome::operator= (const Chromosome &sourceChrom)
{
  int ii;
  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
    genes[ii] = sourceChrom.genes[ii];
  return *this;
}


/*
 *==========================================================
 * Function Name: mutateMinMax
 * Function Task: Selective mutation. Randomly select one of
 *    the genes and replace it with a random variables sampled
 *    uniformly between the gene ranges if the gene is not
 *    frozen. Perform this process with a probability equal
 *    to mutation probability
 * Input parameters: freezeMask[noOfDecisionVariables]
 * Output: None
 * Functions Invoked: myRandom::boundedRandom(Min, Max)
 *=========================================================
 */

void Chromosome::mutateMinMax (const int *freezeMask) {
  int position;

  // select a gene position randomly
  position = myRandom.boundedIntegerRandom(0,globalSetup->noOfDecisionVariables);

  // If the chosen gene is not frozen (only for SGA)
  if(((globalSetup->gaType == SGA)&&(freezeMask[position] == OFF))||
     (globalSetup->gaType == NSGA)) {
    //Flip a biased coin (bias = mutation probability)
    if(myRandom.flip(globalSetup->mutationProbability)) {
      //Replace the selected gene with a uniform random variate within gene
      // ranges
      genes[position] = myRandom.boundedRandom(
          globalSetup->variableRanges[position][0],
          globalSetup->variableRanges[position][1]);
      //If the gene is an integer round it off
      if(globalSetup->variableTypes[position] == Integer)
	genes[position] = int(genes[position]+0.5);
    }
  }
}

/*
 *==========================================================
 * Function Name: mutateNormal
 * Function Task: Genewise mutation. To each of the gene
 *    add a value sampled from a normal distribution with
 *    zero mean and a user specified standard deviation.
 *    if the gene is not frozen. Perform this process with a
 *    probability equal to mutation probability
 * Input parameters: freezeMask[noOfDecisionVariables]
 *                   globalSetup::mutationParameters
 * Output: None
 * Functions Invoked: myRandom::normalRandom(standardDeviation)
 *=========================================================
 */
void Chromosome::mutateNormal(const int *freezeMask) {
  int ii;

  //Copy the standard deviations from global setup

  //For each decision variable
  for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {

    // If the GA type is SGA then check if freeze mask of the variable is off
    if(((globalSetup->gaType == SGA)&&(freezeMask[ii] == OFF))||
       (globalSetup->gaType == NSGA)) {
      //Flip a biased coin with mutation probability
      if(myRandom.flip(globalSetup->mutationProbability)) {
	// Added a random guassian noise to the variable value
	genes[ii] = genes[ii] +
            myRandom.normalRandom(((double *)globalSetup->
                                   mutationParameters)[ii]);
	// If the variable is an integer, round it to closest integer
	if(globalSetup->variableTypes[ii]==Integer)
	  genes[ii] = int(genes[ii]+0.5);
	//If the variables are beyond the range then perturb it only till the
        // range
	if(genes[ii] > globalSetup->variableRanges[ii][1])
          genes[ii] = globalSetup->variableRanges[ii][1];
	if(genes[ii] < globalSetup->variableRanges[ii][0])
          genes[ii] = globalSetup->variableRanges[ii][0];
      }
    }
  }
}

/*
 *==========================================================
 * Function Name: mutatePolynomial
 * Function Task: Polynomial mutation. This procedure is identical to simulated
 * binary crossover, except that the slope between the current gene value and
 * the gene range (either minimum or maximum whichever is closer to the gene
 * value).
 * Input parameters: freezeMask[noOfDecisionVariables]
 *                   globalSetup::mutationParameters
 * Output: None
 * Functions Invoked: myRandom::random01()
 *                    myRandom::flip(double probability)
 *                    myRandom::boundedRandom(Min, Max)
 * Reference: Deb, K., & Agarwal, R.B. (1995) "Simulated Binary Crossover for
 * Continuous Search Space", Complex Systems. 9. pp. 115-148. (TCGA No. 05896).
 *            Deb, K., & Kumar, A. (1995) "Real-Coded Genetic Algorithms with
 * Simulated Binary Crossover: Studies on Multimodal and Multiobjective
 * Problems", Complex Systems. 9. pp. 431-454. (TCGA No. 05897).
 *=========================================================
 */

void Chromosome::mutatePolynomial(const int *freezeMask) {
  int ii;
  double geneMin, geneMax, rndNo;
  double delta, value, deltaq;
  int etaM; //The order of the polynomial
  if(globalSetup->mutationParameters == NULL)
    etaM = 20;
  else
    etaM = ((int *)globalSetup->mutationParameters)[0];

  // For each variable do the following
  for(ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
    // If the GA type is SGA check if the freeze mask for the current variable
    // is off
    if(((globalSetup->gaType == SGA)&&(freezeMask[ii] == OFF))||
       (globalSetup->gaType == NSGA)) {
      // Flip a biased coin with mutation probability
      if(myRandom.flip(globalSetup->mutationProbability)) {
	geneMin = globalSetup->variableRanges[ii][0];
	geneMax = globalSetup->variableRanges[ii][1];
	// If the gene value is greater than the minimum
	if(genes[ii] > geneMin) {
	  // if the gene value is closer to the minimum, compute the slope with
          // the gene value and minimum as points

	  if((genes[ii]-geneMin) < (geneMax-genes[ii]))
	    delta = (genes[ii] - geneMin)/(geneMax-geneMin);
	  // if the gene value is closer to the maximum, compute the slope with
          // the gene value and maximum as points
	  else delta = (geneMax - genes[ii])/(geneMax-geneMin);
	  rndNo = myRandom.random01();
	  if(rndNo <= 0.5) {
	    value = 2.0*rndNo +
                (1.0 - 2.0*rndNo)*pow(1.0-delta,(double)(etaM+1));
	    deltaq = pow(value, 1.0/((double)(etaM+1))) - 1.0;
	  }
	  else {
	    value = 2.0*(1.0-rndNo) +
                2.0*(rndNo - 0.5)*pow(1.0-delta,(double)(etaM+1));
	    deltaq = 1.0 - pow(value, 1.0/((double)(etaM+1)));
	  }
	  genes[ii] = genes[ii] + deltaq*(geneMax - geneMin);
	}
	// If the gene value is less than the minimum, assign a value sampled
        // from a uniform distribution between the gene ranges.
	else genes[ii] = myRandom.boundedRandom(geneMin, geneMax);
	// If the gene is an integer, then round it off to the closest integer
	if(globalSetup->variableTypes[ii] == Integer)
          genes[ii] = int(genes[ii]+0.5);
	// If the gene value is not within the range then bounded either to its
        // minimum or maximum value
	if(genes[ii] > globalSetup->variableRanges[ii][1])
          genes[ii] = globalSetup->variableRanges[ii][1];
	if(genes[ii] < globalSetup->variableRanges[ii][0])
          genes[ii] = globalSetup->variableRanges[ii][0];
      }
    }
  }
}

/*
 *==========================================================
 * Function Name: mutatePolynomial
 * Function Task: Polynomial mutation. Mutates a given variable
 *    using polynomial mutation. This routine is used by local
 *    search method (simplex and complex) to create initial
 *    points near a given point.
 * Input parameters: int index
 * Input range: [0, noOfDecisionVariables)
 * Output: None
 * Functions Invoked: myRandom::random01()
 * Reference: Deb, K., & Agarwal, R.B. (1995) "Simulated Binary Crossover for
 * Continuous Search Space", Complex Systems. 9. pp. 115-148. (TCGA No. 05896).
 *            Deb, K., & Kumar, A. (1995) "Real-Coded Genetic Algorithms with
 * Simulated Binary Crossover: Studies on Multimodal and Multiobjective
 * Problems", Complex Systems. 9. pp. 431-454. (TCGA No. 05897).
 *=========================================================
 */

void Chromosome::mutatePolynomial(const int index) {
  double geneMin, geneMax, rndNo;
  double delta, value, deltaq;
  int etaM; //The order of the polynomial
  etaM = ((int *)globalSetup->mutationParameters)[0];

  geneMin = globalSetup->variableRanges[index][0];
  geneMax = globalSetup->variableRanges[index][1];
  if(genes[index] > geneMin) {
    if((genes[index]-geneMin) < (geneMax-genes[index]))
      delta = (genes[index] - geneMin)/(geneMax-geneMin);
    else delta = (geneMax - genes[index])/(geneMax-geneMin);
    rndNo = myRandom.random01();
    if(rndNo <= 0.5) {
      value = 2.0*rndNo + (1.0 - 2.0*rndNo)*pow(1.0-delta,(double)(etaM+1));
      deltaq = pow(value, 1.0/((double)(etaM+1))) - 1.0;
    }
    else {
      value = 2.0*(1.0-rndNo) +
          2.0*(rndNo - 0.5)*pow(1.0-delta,(double)(etaM+1));
      deltaq = 1.0 - pow(value, 1.0/((double)(etaM+1)));
    }
    genes[index] = genes[index] + deltaq*(geneMax - geneMin);
  }
  else genes[index] = myRandom.boundedRandom(geneMin, geneMax);
  if(globalSetup->variableTypes[index] == Integer)
    genes[index] = int(genes[index]+0.5);
  if((genes[index] > geneMax) || (genes[index] < geneMin))
    genes[index] = myRandom.boundedRandom(geneMin, geneMax);
}
