/*
  Christini lab genetic algorithm
  Modified from IlliGAL GA, reference below

  Most of file derived from original userDefinables.cpp
  Added Timestamp feature, which requires C++11 standard.
*/

/*
  Single & Multi-Objective Real-Coded Genetic Algorithms Code
  Author: Kumara Sastry
  Illinois Genetic Algorithms Laboratory (IlliGAL)
  Deparment of General Engineering
  University of Illinois at Urbana-Champaign
  104 S. Mathews Ave, Urbana, IL 61801
*/

#include "../include/ga.hpp"

#include <iostream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>

#define BLANK_STR " \t\n\r"

/*
  Uncomment to add a time stamp at the beginning of each data file. Requires
  the use of chrono library, which requires code to be compiled using the
  C++11 standard.
*/
// #define TIMESTAMP

#ifdef TIMESTAMP
#include <chrono>
#endif

extern GlobalSetup *globalSetup;
extern Random myRandom;

/*
  Read a non-comment, non-blank line from the specified file stream
  At EOF, it returns NULL.
  Otherwise, it returns the pointer to the first token.
  (The string just read in is altered by strtok().)
*/
static char* readOneLine(char *pcBuf, int iMaxSize, FILE *fStream) {
  char *pToken;

  do {
    pToken = NULL;

    *pcBuf = '\0';
    if (fgets(pcBuf, iMaxSize, fStream) == NULL) {
      std::cerr << "Error: unable to read line" << std::endl;
      exit(1);
    }

    if (feof(fStream))
      break;

    // get the first token
    pToken = strtok(pcBuf, BLANK_STR);

    // if there is no first token, it is a blank line.
    // if the first token starts with '#', it is a comment line.
  } while ((pToken == NULL) || (*pToken == '#'));

  return (pToken);
}

int run_GA(char *inputFileName) {
  int ii;
  const int ciBufSize = 1024;
  char *pToken, caBuf[ciBufSize];
  FILE *fInput, *fOutput;

  fInput = fopen(inputFileName, "r");
  if (fInput == NULL) {
    printf("Error! opening file %s\n", inputFileName);
    exit(1);
  }

  globalSetup = new GlobalSetup;

  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }

  //GA type
  if (strcmp("SGA", pToken) == 0) {
    globalSetup->gaType = SGA;
  }
  else if (strcmp("NSGA", pToken) == 0) {
    globalSetup->gaType = NSGA;
  }
  else {
    fclose(fInput);
    printf("Unknown parameter! It should be either SGA or NSGA\n");
    exit(1);
  }

  // decision variables

  // read the number of decision variables
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  // the number can't be less than or equal to zero
  if ((globalSetup->noOfDecisionVariables = atoi(pToken)) <= 0) {
    fclose(fInput);
    printf("Error! number of decision variables should be > 0\n");
    exit(1);
  }

  globalSetup->variableTypes =
      new VariableType[globalSetup->noOfDecisionVariables];
  globalSetup->variableRanges =
      new double*[(globalSetup->noOfDecisionVariables)];
  for(ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
    // read a line
    if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }
    // variable type
    if (strcmp("double", pToken) == 0) {
      globalSetup->variableTypes[ii] = Real;
    }
    else if (strcmp("int", pToken) == 0) {
      globalSetup->variableTypes[ii] = Integer;
    }
    else {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }

    // variable ranges
    // allocate memory
    globalSetup->variableRanges[ii] = new double[2];

    // lower bound
    if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }
    globalSetup->variableRanges[ii][0] = atof(pToken);

    // upper bound
    if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }
    globalSetup->variableRanges[ii][1] = atof(pToken);
    if(globalSetup->variableRanges[ii][1] <=
       globalSetup->variableRanges[ii][0]) {
      fclose(fInput);
      printf("Error! lower bound, %f, must be lower than the upper bound, %f\n",
             globalSetup->variableRanges[ii][0],
             globalSetup->variableRanges[ii][1]);
      exit(1);
    }
  }

  // objectives

  // read the number of objectives
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  // the number can't be less than or equal to zero
  if ((globalSetup->noOfRawObjectives = atoi(pToken)) <= 0) {
    fclose(fInput);
    printf("Error! number of objectives should be > 0\n");
    exit(1);
  }
  globalSetup->finalNoOfObjectives = globalSetup->noOfRawObjectives;
  globalSetup->noOfLinearObjectiveCombinations = 0;

  globalSetup->typeOfOptimizations =
      new OptimType[globalSetup->finalNoOfObjectives];
  for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
    // read a line
    if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }
    // optimization type
    if (strcmp("Min", pToken) == 0) {
      globalSetup->typeOfOptimizations[ii] = Minimization;
    }
    else if (strcmp("Max", pToken) == 0) {
      globalSetup->typeOfOptimizations[ii] = Maximization;
    }
    else {
      fclose(fInput);
      printf("Error! optimization type can either be Min or Max\n");
      exit(1);
    }
  }

  // constrained variables

  // read the number of constrained variables
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  // the number can't be less than zero
  if ((globalSetup->noOfRawConstraints = atoi(pToken)) < 0) {
    fclose(fInput);
    printf("Error! number of constraints should be >= 0\n");
    exit(1);
  }
  else if (globalSetup->noOfRawConstraints == 0) {
    if ((strlen(pToken) != 1) || (*pToken != '0')) {
      fclose(fInput);
      printf("Error! number of constraints should be >= 0\n");
      exit(1);
    }
  }
  globalSetup->noOfLinearConstraintCombinations = 0;
  globalSetup->finalNoOfConstraints = globalSetup->noOfRawConstraints;

  // penalty weights
  globalSetup->penaltyWeights = new double[globalSetup->finalNoOfConstraints];
  for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++) {
    // read a line
    if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
      fclose(fInput);
      printf("Error in the input file, please refer to the documentation\n");
      exit(1);
    }
    globalSetup->penaltyWeights[ii] = atof(pToken);
  }

  // general parameters
  // population size
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    // Use a default value of 30*ell*log(ell);
    printf("Using default population-sizing thumbrule: n = 30*ell*log(ell)\n");

    // Take care of small problem sizes where log(ell) < 1
    if(globalSetup->noOfDecisionVariables > 2)
      globalSetup->populationSize =
          (int)(30*(globalSetup->noOfDecisionVariables) *
                log((double)(globalSetup->noOfDecisionVariables)));
    else
      globalSetup->populationSize =
          (int)(30*(globalSetup->noOfDecisionVariables));

    //Round it to next nearest tenth number
    if((globalSetup->populationSize)%10)
      globalSetup->populationSize += (globalSetup->populationSize)%10;
    printf("The population size used is: %d\n", globalSetup->populationSize);
  }
  // the number can't be less than or equal to zero
  else if ((globalSetup->populationSize = atoi(pToken)) <= 0) {
    fclose(fInput);
    printf("The population size must be > 0\n");
    exit(1);
  }
  else if ((globalSetup->populationSize % 2) != 0) {
    // the number can't be an odd number
    fclose(fInput);
    printf("Error! population size must be an even number\n");
    exit(1);
  }
  // maximum generations
  // the number can't be less than or equal to zero
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    // Use a default value of 6*ell;
    printf("Using default convergence-time thumbrule: tc = 6*ell\n");
    globalSetup->maxGenerations = 6*(globalSetup->noOfDecisionVariables);
    if((globalSetup->maxGenerations)%10)
      globalSetup->maxGenerations += 10-(globalSetup->maxGenerations)%10;
    printf("The maximum number of generations set is: %d\n",
           globalSetup->maxGenerations);
  }
  else if ((globalSetup->maxGenerations = atoi(pToken)) <= 0) {
    fclose(fInput);
    printf("Error! maximum number of generations must be > 0\n");
    exit(1);
  }
  // replace proportion
  // the number should be in (0.0, 1.0]
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    // Setting a default value of 0.9
    printf("Using default replacement proportion of 0.9\n");
    globalSetup->replaceProportion = 0.9;
  }
  else if (((globalSetup->replaceProportion = atof(pToken)) <= 0.0) ||
           (globalSetup->replaceProportion > 1.0)) {
    fclose(fInput);
    printf("Error! proportion of parent population that should be replaced must"
           " be > 0 and <= 1\n");
    exit(1);
  }

  // niching (multimodal handling)
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    // Using NoNiching by default
    printf("No niching method is used by default\n");
    globalSetup->nichingType = NoNiching;
  }
  // niching type
  else if (strcmp("NoNiching", pToken) == 0) {
    globalSetup->nichingType = NoNiching;  globalSetup->nichingParameters = NULL;

  }
  else if (strcmp("Sharing", pToken) == 0) {
    globalSetup->nichingType = Sharing;
  }
  else if (strcmp("RTS", pToken) == 0) {
    globalSetup->nichingType = RTS;
  }
  else if (strcmp("DeterministicCrowding", pToken) == 0) {
    globalSetup->nichingType = DeterministicCrowding;
  }
  else {
    fclose(fInput);
    printf("Error! valid niching types are: NoNiching, Sharing, RTS, and "
           "DeterministicCrowding\n");
    exit(1);
  }
  // check niching type
  if ((globalSetup->gaType == NSGA) &&
      (globalSetup->nichingType != NoNiching)) {
    fclose(fInput);
    printf("Error! valid choice for niching types with NSGA is: NoNiching\n");
    exit(1);
  }
  // read niching parameters
  switch(globalSetup->nichingType) {
    case NoNiching:
    case DeterministicCrowding:
      // no extra parameters
      break;
    case Sharing:
      {
        globalSetup->nichingParameters = new double[2];
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default sharing radius of 4.24\n");
          ((double *)(globalSetup->nichingParameters))[0] = 4.24;
        }
        else
          ((double *)globalSetup->nichingParameters)[0] = atof(pToken);
        if (((double*)globalSetup->nichingParameters)[0] <= 0.0) {
          fclose(fInput);
          printf("Error! niching radius must be > 0\n");
          exit(1);
        }
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default scaling of 1 for fitness sharing\n");
          ((double*)(globalSetup->nichingParameters))[1] = 1.0;
        }
        else
          ((double*)(globalSetup->nichingParameters))[1] = atof(pToken);
      }
      break;
    case RTS:
      {
        globalSetup->nichingParameters = new int[1];

        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Setting the window size for RTS using default rule: w = "
                 "min(ell, n/20)\n");
          if(globalSetup->noOfDecisionVariables <
             (globalSetup->populationSize/20))
            ((int*)(globalSetup->nichingParameters))[0] =
                globalSetup->noOfDecisionVariables;
          else
            ((int*)(globalSetup->nichingParameters))[0] =
                (globalSetup->populationSize)/20;

          //Adjust for small problem sizes: w = n/20;
          if(globalSetup->noOfDecisionVariables < 20)
            ((int*)(globalSetup->nichingParameters))[0] =
                (globalSetup->populationSize)/20;

          //Check if the window size is greater than the population size
          if(((int*)globalSetup->nichingParameters)[0] >
             globalSetup->populationSize)
            ((int*)(globalSetup->nichingParameters))[0] =
                globalSetup->populationSize;

          printf("The window size used for RTR is: %d\n",
                 ((int*)globalSetup->nichingParameters)[0]);
        }
        else
          ((int*)(globalSetup->nichingParameters))[0] = atoi(pToken);

        // the window size should be in (0, populationSize]
        if ((((int*)globalSetup->nichingParameters)[0] <= 0) ||
            (((int*)globalSetup->nichingParameters)[0] >
             globalSetup->populationSize)) {
          fclose(fInput);
          printf("Error! window size for RTR should be > 0 and less than the "
                 "population size\n");
          exit(1);
        }
      }
      break;
    default:
      {
        fclose(fInput);
        printf("Error! valid choices for niching type are: NoNiching, Sharing, "
               "RTS, and DeterministicCrowding\n");
        exit(1);
      }
      break;
  }

  // selection
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Using tournament selection w/o replacement as a default selection "
           "method\n");
    globalSetup->selectionType = TournamentWOR;
  }
  // selection type
  else if (strcmp("TournamentWOR", pToken) == 0) {
    globalSetup->selectionType = TournamentWOR;
  }
  else if (strcmp("SUS", pToken) == 0) {
    globalSetup->selectionType = SUS;
  }
  else if (strcmp("Truncation", pToken) == 0) {
    globalSetup->selectionType = Truncation;
  }
  else if (strcmp("RouletteWheel", pToken) == 0) {
    globalSetup->selectionType = RouletteWheel;
  }
  else if (strcmp("TournamentWR", pToken) == 0) {
    globalSetup->selectionType = TournamentWR;
  }
  else {
    fclose(fInput);
    printf("Error! valid selection methods are: RouletteWheel, SUS, "
           "TournamentWOR, TournamentWR, and Truncation\n");
    exit(1);
  }

  // check selection type
  if ((globalSetup->gaType == NSGA) &&
      ((globalSetup->selectionType == SUS) ||
       (globalSetup->selectionType == RouletteWheel))) {
    fclose(fInput);
    printf("Error! with NSGA, valid selection methods are: TournamentWOR, "
           "TournamentWR, and Truncation\n");
    exit(1);
  }
  // read selection parameters
  switch(globalSetup->selectionType) {
    case TournamentWOR:
    case Truncation:
    case TournamentWR:
      {
        globalSetup->selectionParameters = new int[1];

        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default tournament size of 2\n");
          ((int*)globalSetup->selectionParameters)[0] = 2;
        }
        else
          ((int*)globalSetup->selectionParameters)[0] = atoi(pToken);
      }
      break;
    case SUS:
    case RouletteWheel:
      // no extra parameters
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid selection parameter\n");
        exit(1);
      }
      break;
  }

  // Crossover
  // crossover probability
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Using a default crossover probability of 0.9\n");
    globalSetup->xOverProbability = 0.9;
  }
  // the number should be [0.0, 1.0]
  else if (((globalSetup->xOverProbability = atof(pToken)) < 0.0) ||
           (globalSetup->xOverProbability > 1.0)) {
    fclose(fInput);
    printf("Error! crossover probability must be >= 0.0 and <= 1.0\n");
    exit(1);
  }
  // crossover type
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Using SBX as the default crossover method. Note that this might be "
           "an inappropriate choice if your variables are binary\n");
    globalSetup->xOverType = SBX;
  }
  else if (strcmp("OnePoint", pToken) == 0) {
    globalSetup->xOverType = OnePoint;
  }
  else if (strcmp("TwoPoint", pToken) == 0) {
    globalSetup->xOverType = TwoPoint;
  }
  else if (strcmp("Uniform", pToken) == 0) {
    globalSetup->xOverType = Uniform;
  }
  else if (strcmp("SBX", pToken) == 0) {
    globalSetup->xOverType = SBX;
  }
  else {
    fclose(fInput);
    printf("Error! valid crossover types are: OnePoint, TwoPoint, Uniform, and "
           "SBX\n");
    exit(1);
  }
  // read crossover parameters
  switch(globalSetup->xOverType) {
    case OnePoint:
    case TwoPoint:
      //no extra parameters
      break;
    case Uniform:
      {
        globalSetup->xOverParameters = new double[1];
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default genewise swap probability of 0.5\n");
          ((double*)globalSetup->xOverParameters)[0] = 0.5;
        }
        else
          ((double*)globalSetup->xOverParameters)[0] = atof(pToken);
        if((((double*)globalSetup->xOverParameters)[0] <= 0.0)||
           (((double*)globalSetup->xOverParameters)[0] >= 1.0)) {
          fclose(fInput);
          printf("Genewise probability must be > 0.0 and < 1.0\n");
          exit(1);
        }
      }
    case SBX:
      {
        globalSetup->xOverParameters = new double[2];

        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default genewise swap probability of 0.5\n");
          ((double*)globalSetup->xOverParameters)[0] = 0.5;
        }
        else
          ((double*)globalSetup->xOverParameters)[0] = atof(pToken);
        if((((double*)globalSetup->xOverParameters)[0] <= 0.0)||
           (((double*)globalSetup->xOverParameters)[0] >= 1.0)) {
          fclose(fInput);
          printf("Error! genewise probability must be > 0.0 and < 1.0\n");
          exit(1);
        }
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using default polynomial order for SBX: 10\n");
          ((double*)globalSetup->xOverParameters)[1] = 10;
        }
        else
          ((double*)globalSetup->xOverParameters)[1] = atof(pToken);
        if(((double*)globalSetup->xOverParameters)[1] < 0.0) {
          fclose(fInput);
          printf("Error! genewise probability must be >= 0.0\n");
          exit(1);
        }
      }
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid crossover parameter\n");
        exit(1);
      }
      break;
  }

  // mutation
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Using a default mutation probability of 0.1\n");
    globalSetup->mutationProbability = 0.1;
  }
  // mutation probability
  // the number should be [0.0, 1.0]
  else if (((globalSetup->mutationProbability = atof(pToken)) < 0.0)
           || (globalSetup->mutationProbability > 1.0)) {
    fclose(fInput);
    printf("Error! mutation probability must be >= 0.0 and <= 1.0.\n");
    exit(1);
  }
  // mutation type
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Using Polynomial as the default mutation method. Note that this "
           "might be an inappropriate choice if your variables are binary\n");
    globalSetup->mutationType = Polynomial;
  }
  else if (strcmp("Selective", pToken) == 0) {
    globalSetup->mutationType = Selective;
  }
  else if (strcmp("Genewise", pToken) == 0) {
    globalSetup->mutationType = Genewise;
  }
  else if (strcmp("Polynomial", pToken) == 0) {
    globalSetup->mutationType = Polynomial;
  }
  else {
    fclose(fInput);
    printf("Error! valid mutation types are: Selective, Genewise, and "
           "Polynomial\n");
    exit(1);
  }
  // read mutation parameters
  switch(globalSetup->mutationType) {
    case Selective:
      // no extra parameters
      break;
    case Polynomial:
      {
        globalSetup->mutationParameters = new int[1];
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using a default value for the polynomial probability: 20\n");
          ((int*)globalSetup->mutationParameters)[0] = 20;
        }
        else
          ((int*)globalSetup->mutationParameters)[0] = atoi(pToken);
        if(((int*)globalSetup->mutationParameters)[0] < 0) {
          fclose(fInput);
          printf("Error! polynomial order for polynomial mutation must be > "
                 "0\n");
          exit(1);
        }
      }
      break;
    case Genewise:
      {
        globalSetup->mutationParameters =
            new double[globalSetup->noOfDecisionVariables];

        for (ii = 0 ; ii < globalSetup->noOfDecisionVariables ; ii++) {
          if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
            printf("Using default std. deviation of 10 percent of the variable "
                   "range\n");
            ((double*)globalSetup->mutationParameters)[ii] = 0.1 *
                (globalSetup->variableRanges[ii][1] -
                 globalSetup->variableRanges[ii][0]);
          }
          else
            ((double*)globalSetup->mutationParameters)[ii] = atof(pToken);
          if(((double*)globalSetup->mutationParameters)[ii] <= 0.0) {
            fclose(fInput);
            printf("Error! standard deviation for gene %d for genewise mutation"
                   " must be > 0\n", ii);
            exit(1);
          }
        }
      }
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid mutation parameter\n");
        exit(1);
      }
      break;
  }

  // scaling method
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Scaling is not used by default\n");
    globalSetup->scalingMethod = NoScaling;
  }
  // scaling method
  else if (strcmp("NoScaling", pToken) == 0) {
    globalSetup->scalingMethod = NoScaling;
  }
  else if (strcmp("Ranking", pToken) == 0) {
    globalSetup->scalingMethod = Ranking;
  }
  else if (strcmp("SigmaScaling", pToken) == 0) {
    globalSetup->scalingMethod = SigmaScaling;
  }
  else {
    fclose(fInput);
    printf("Error! valid scaling methods are: NoScaling, Ranking, and "
           "SigmaScaling\n");
    exit(1);
  }
  // read scaling parameters
  switch(globalSetup->scalingMethod) {
    case NoScaling:
    case Ranking:
      // no extra parameters
      break;
    case SigmaScaling:
      {
        globalSetup->scalingParameters = new double[1];

        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          ((double*)globalSetup->scalingParameters)[0] = 1.0;
        }
        else
          ((double*)globalSetup->scalingParameters)[0] = atof(pToken);
        if(((double*)globalSetup->scalingParameters)[0] <= 0.0) {
          fclose(fInput);
          printf("Error! scaling parameter for SigmaScaling must be > 0\n");
          exit(1);
        }
      }
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid scaling parameter\n");
        exit(1);
      }
      break;
  }

  // constraint method
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    if(globalSetup->finalNoOfConstraints == 0) {
      printf("Using no constraint handling method by default\n");
      globalSetup->constraintMethod = NoConstraints;
    }
    else {
      printf("Using tournament selection as the default constraint handling "
             "method\n");
      globalSetup->constraintMethod = Tournament;
    }
  }
  // constraint method
  else if (strcmp("NoConstraints", pToken) == 0) {
    globalSetup->constraintMethod = NoConstraints;
  }
  else if (strcmp("Penalty", pToken) == 0) {
    globalSetup->constraintMethod = Penalty;
  }
  else if (strcmp("Tournament", pToken) == 0) {
    globalSetup->constraintMethod = Tournament;
  }
  else {
    fclose(fInput);
    printf("Error! valid constraint handling methods are: NoConstraint, Penalty"
           ", and Tournament\n");
    exit(1);
  }
  // check constraint method
  if ((globalSetup->gaType == NSGA) &&
      (globalSetup->constraintMethod == Penalty)) {
    fclose(fInput);
    printf("Error! penalty based constraint handling method cannot be used with"
           " NSGA\n");
    exit(1);
  }
  if ((globalSetup->finalNoOfConstraints == 0) &&
      (globalSetup->constraintMethod != NoConstraints)) {
    fclose(fInput);
    printf("Error! valid constraint-handling method when there are no "
           "constraints is NoConstraints\n");
    exit(1);
  }
  // read penalty function
  switch(globalSetup->constraintMethod) {
    case NoConstraints:
    case Tournament:
      // no extra parameters
      break;
    case Penalty:
      {
        // penalty function
        if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
          printf("Using linear penalty be default\n");
          globalSetup->penaltyFunction = Linear;
        }
        else if (strcmp("Linear", pToken) == 0) {
          globalSetup->penaltyFunction = Linear;
        }
        else if (strcmp("Quadratic", pToken) == 0) {
          globalSetup->penaltyFunction = Quadratic;
        }
        else {
          fclose(fInput);
          printf("Error! valid penalty function methods are: Linear and "
                 "Quadratic\n");
          exit(1);
        }
      }
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid constraint-handling method parameter\n");
        exit(1);
      }
      break;
  }

  // local search method
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0)
  {
    printf("Using no local search by default\n");
    globalSetup->localSearchMethod = NoLocalSearch;
  }
  // local search method
  else if (strcmp("NoLocalSearch", pToken) == 0) {
    globalSetup->localSearchMethod = NoLocalSearch;
  }
  else if (strcmp("SimplexSearch", pToken) == 0) {
    globalSetup->localSearchMethod = SimplexSearch;
  }
  else {
    fclose(fInput);
    printf("Error! valid local search methods are: NoLocalSearch and "
           "SimplexSearch\n");
    exit(1);
  }
  // check local search method
  if ((globalSetup->localSearchMethod != NoLocalSearch) &&
      (globalSetup->gaType == NSGA)) {
    fclose(fInput);
    printf("Error! cannot use local search with NSGA\n");
    exit(1);
  }
  switch(globalSetup->localSearchMethod) {
    case NoLocalSearch:
      // no extra parameters
      break;
    case SimplexSearch:
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Using default tolerance of 0.001\n");
        globalSetup->maxLocalTolerance = 1.0E-3;
      }
      else
        globalSetup->maxLocalTolerance = atof(pToken);
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Setting maximum local search evaluations to 20\n");
        globalSetup->maxLocalEvaluations = 20;
      }
      else
        globalSetup->maxLocalEvaluations = atoi(pToken);
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Using default local penalty parameter of 1.0\n");
        globalSetup->initialLocalPenaltyParameter = 1.0;
      }
      else
        globalSetup->initialLocalPenaltyParameter = atof(pToken);
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Using default local update parameter of 2.0\n");
        globalSetup->localUpdateParameter = 2.0;
      }
      else
        globalSetup->localUpdateParameter = atof(pToken);
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Using default Lamarckian probability of 0.15\n");
        globalSetup->lamarckianProbability = 0.15;
      }
      else
        globalSetup->lamarckianProbability = atof(pToken);
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Using default overall local search probability of 0.5\n");
        globalSetup->localSearchProbability = 0.5;
      }
      else
        globalSetup->localSearchProbability = atof(pToken);
      break;
    default:
      {
        fclose(fInput);
        printf("Error! invalid local search parameters\n");
        exit(1);
      }
      break;
  }

  // stopping criteria

  // noOfStoppingCriterias
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0)
  {
    printf("Using no extra stopping criterias by default.\n");
    globalSetup->noOfStoppingCriterias = 0;
  }
  // the number can't be less than zero
  else if ((globalSetup->noOfStoppingCriterias = atoi(pToken)) < 0) {
    fclose(fInput);
    printf("Error! number of stopping criterias must be > 0\n");
    exit(1);
  }
  else if (globalSetup->noOfStoppingCriterias == 0) {
    if ((strlen(pToken) != 1) || (*pToken != '0')) {
      fclose(fInput);
      printf("Error! number of stopping criterias must be a number!\n");
      exit(1);
    }
  }

  // allocate memory
  if(globalSetup->noOfStoppingCriterias > 0) {
    // number of generation window
    if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
      printf("The default window size used in stopping critiria is: 5.\n");
      globalSetup->genNumWindow = 5;
    }
    else
      globalSetup->genNumWindow = atoi(pToken);

    // the number should be greater than zero
    if (globalSetup->genNumWindow <= 0) {
      fclose(fInput);
      printf("Error! window size used in stopping critieria must be > 0\n");
      exit(1);
    }
    globalSetup->otherStoppingCriteria =
        new StoppingCriterias[globalSetup->noOfStoppingCriterias];
    globalSetup->stoppingParameter =
        new double[globalSetup->noOfStoppingCriterias];

    for (ii = 0 ; ii < globalSetup->noOfStoppingCriterias ; ii++) {
      // read a line
      if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
        fclose(fInput);
        printf("Error in the input file, Please see the documentation.\n");
        exit(1);
      }
      // stopping criterion type
      if (strcmp("NoOfEvaluations", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = NoOfEvaluations;
      }
      else if (strcmp("FitnessVariance", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = FitnessVariance;
      }
      else if (strcmp("AverageFitness", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = AverageFitness;
      }
      else if (strcmp("AverageObjective", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = AverageObjective;
      }
      else if (strcmp("ChangeInBestFitness", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInBestFitness;
      }
      else if (strcmp("ChangeInAvgFitness", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInAvgFitness;
      }
      else if (strcmp("ChangeInFitnessVar", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInFitnessVar;
      }
      else if (strcmp("ChangeInBestObjective", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInBestObjective;
      }
      else if (strcmp("ChangeInAvgObjective", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInAvgObjective;
      }
      else if (strcmp("NoOfFronts", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = NoOfFronts;
      }
      else if (strcmp("NoOfGuysInFirstFront", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = NoOfGuysInFirstFront;
      }
      else if (strcmp("ChangeInNoOfFronts", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = ChangeInNoOfFronts;
      }
      else if (strcmp("BestFitness", pToken) == 0) {
        globalSetup->otherStoppingCriteria[ii] = BestFitness;
      }
      else {
        fclose(fInput);
        printf("Error! invalid stopping criteria parameter\n");
        exit(1);
      }

      // check criteria
      switch(globalSetup->otherStoppingCriteria[ii]) {
        case NoOfEvaluations:
        case FitnessVariance:
          // 1 - nothing to check
          break;

        case AverageFitness:
        case AverageObjective:
        case ChangeInBestFitness:
        case ChangeInAvgFitness:
        case ChangeInFitnessVar:
        case ChangeInBestObjective:
        case ChangeInAvgObjective:
        case BestFitness:
          // 2 - SGA w/o MM only
          if ((globalSetup->gaType != SGA) ||
              (globalSetup->nichingType != NoNiching)) {
            fclose(fInput);
            printf("Cannot use the following stopping critier with NSGA or when"
                   " Niching is used:\n");
            printf("\t AverageFitness\n\t AverageObjective\n");
            printf("\t ChangeInBestFitness\n\t ChangeInAvgFitness\n");
            printf("\t ChangeInFitnessVar\n\t ChangeInBestObjective\n");
            printf("\t ChangeInAvgObjective\n\t BestFitness\n");
            printf("Error! invalid choice of stopping criteria");
            exit(1);
          }
          break;
        case NoOfFronts:
        case NoOfGuysInFirstFront:
        case ChangeInNoOfFronts:
          // 3 - NSGA only
          if (globalSetup->gaType != NSGA) {
            fclose(fInput);
            printf("Error! following stopping criteria can be used only with "
                   "NSGA:\n");
            printf("\t NoOfFronts\n\t NoOfGuysInFirstFront\n\t "
                   "ChangeInNoOfFronts\n");
            exit(1);
          }
          break;

        default:
          {
            fclose(fInput);
            printf("Error! invalid stopping criteria\n");
            exit(1);
          }
          break;
      }
      // parameter
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        fclose(fInput);
        printf("Error! invalid stopping criteria parameter\n");
        exit(1);
      }
      ((double*)globalSetup->stoppingParameter)[ii] = atof(pToken);
    }
  }
  // Initial population generation

  // Load population
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0)
  {
    printf("Using random initialization by default.\n");
    globalSetup->loadPopulation = false;
    globalSetup->evaluateAgain = true;
  }
  // Check if the load population is 0 or 1
  else {
    // load the population from a file
    if(atoi(pToken) == 1) {
      globalSetup->loadPopulation = true;

      globalSetup->populationFileName = new char[256];
      // Read the name of the file to load the population from
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        fclose(fInput);
        printf("Error! invalid file name to load the population from, see "
               "documentation for more details.\n");
        exit(1);
      }
      strcpy(globalSetup->populationFileName, pToken);

      // Evaluate the population or not
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        printf("Evaluating the initial population by default.\n");
        globalSetup->evaluateAgain = true;
      }
      else if(atoi(pToken) == 0) {
        globalSetup->evaluateAgain = false;
      }
      else
        globalSetup->evaluateAgain = true;
    }
    // Use random initialization
    else {
      globalSetup->loadPopulation = false;
      globalSetup->evaluateAgain = true;
    }
  }

  // Save population
  if ((pToken = readOneLine(caBuf, ciBufSize, fInput)) == NULL) {
    fclose(fInput);
    printf("Error in the input file, please refer to the documentation\n");
    exit(1);
  }
  if(strcmp("default", pToken) == 0) {
    printf("Saving the evaluated individuals to a file by default.\n");
    globalSetup->savePopulation = true;
    globalSetup->saveEvalSolutions = new char[256];
    printf("Enter the filename you want to save the population to.\n");
    fflush(stdout);
    if (scanf("%s", globalSetup->saveEvalSolutions) == EOF) {
      std::cerr << "Error: unable to read filename" << std::endl;
      exit(1);
    }
  }
  // Check if the save population is 0 or 1
  else {
    fflush(stdout);
    // Save the population to a file
    if(atoi(pToken) == 1) {
      globalSetup->savePopulation = true;
      // Filename to save the evaluated solutions to.
      globalSetup->saveEvalSolutions = new char[256];
      if ((pToken = strtok(NULL, BLANK_STR)) == NULL) {
        fclose(fInput);
        printf("Error in the input file, please refer to the"
               " documentation\n");
        exit(1);
      }
      // Add time to file name
#ifdef TIMESTAMP
      std::chrono::system_clock::time_point p =
          std::chrono::system_clock::now();
      std::time_t t = std::chrono::system_clock::to_time_t(p);
      char time[18]; // holds "2013-12-01_21.31_"
      strftime( time, 18, "%Y.%m-%d_%H.%M_", localtime( &t ) );
      strcpy(globalSetup->saveEvalSolutions, time);
      strcat(globalSetup->saveEvalSolutions, pToken);
#else
      strcpy(globalSetup->saveEvalSolutions, pToken);
#endif
    }
  }
  fclose(fInput);
  fflush(stdout);

  if(globalSetup->savePopulation) {
    fOutput = fopen(globalSetup->saveEvalSolutions, "w");
    for(ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
      fprintf(fOutput, "var#%d\t", ii);
    }
    for(ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
      fprintf(fOutput, "obj#%d\t", ii);
    }
    if(globalSetup->finalNoOfConstraints > 0) {
      for(ii = 0; ii < globalSetup->finalNoOfConstraints; ii++) {
        fprintf(fOutput, "const#%d\t", ii);
      }
      fprintf(fOutput, "penalty");
    }
    fprintf(fOutput, "\n");
    fclose(fOutput);
  }

  GA ga;
  if (globalSetup->gaType == SGA) {
    while(ga.generate());
  }
  else {
    while(ga.nsgaGenerate());
  }

  delete [](globalSetup->penaltyWeights);
  delete [](globalSetup->variableTypes);
  for(ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
    delete [](globalSetup->variableRanges[ii]);
  }
  delete [](globalSetup->variableRanges);
  delete [](globalSetup->typeOfOptimizations);
  if((globalSetup->selectionType == TournamentWR)||
     (globalSetup->selectionType == TournamentWOR)||
     (globalSetup->selectionType == Truncation)) {
    delete [](int*)(globalSetup->selectionParameters);
  }
  if((globalSetup->xOverType==Uniform)||(globalSetup->xOverType == SBX))
    delete [](double*)(globalSetup->xOverParameters);
  if(globalSetup->mutationType == Polynomial)
    delete [](int*)(globalSetup->mutationParameters);
  else if(globalSetup->mutationType == Genewise)
    delete [](double*)(globalSetup->mutationParameters);
  if(globalSetup->nichingType == RTS)
    delete [](int*)(globalSetup->nichingParameters);
  if(globalSetup->nichingType == Sharing)
    delete [](double*)(globalSetup->nichingParameters);
  if(globalSetup->scalingMethod == SigmaScaling)
    delete [](double*)(globalSetup->scalingParameters);
  if(globalSetup->noOfStoppingCriterias > 0) {
    delete [] globalSetup->otherStoppingCriteria;
    delete [](double*)(globalSetup->stoppingParameter);
  }
  if(globalSetup->loadPopulation)
    delete [](globalSetup->populationFileName);
  if(globalSetup->savePopulation)
    delete [](globalSetup->saveEvalSolutions);

  return 1;
}
