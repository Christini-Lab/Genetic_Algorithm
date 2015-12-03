/*
  Christini lab genetic algorithm
  Modified from IlliGAL GA, reference below

  Example main derived from original IlliGAL code.

  Original code housed main and objective function within the
  userDefinables.cpp. For easier readabiity and development, these functions
  have been separated.
*/

/*
  Single & Multi-Objective Real-Coded Genetic Algorithms Code
  Author: Kumara Sastry
  Illinois Genetic Algorithms Laboratory (IlliGAL)
  Deparment of General Engineering
  University of Illinois at Urbana-Champaign
  104 S. Mathews Ave, Urbana, IL 61801
*/

#include "../../include/main_ga.hpp"
#include "../../include/random.hpp"
#include "../../include/globalSetup.hpp"

#include <fstream>
#include <cstdlib>

GlobalSetup *globalSetup;
Random myRandom;

/* Objective function and constraints go here */
void globalEvaluate(double *x, double *objArray, double *constraintViolation,
                    double *penalty, int *noOfViolations) {
  int ii;
  FILE *outEvals;

  for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++)
    objArray[ii] = 0.0;

  if (globalSetup->gaType == SGA) {
    *objArray = (x[0]*x[0] + x[1] - 11.0)*(x[0]*x[0] + x[1] - 11.0) +
        (x[0] + x[1]*x[1] - 7.0)*(x[0] + x[1]*x[1] - 7.0);

    constraintViolation[0] = (x[0]-5.0)*(x[0]-5.0) + x[1]*x[1] - 26.0;
    if (constraintViolation[0] <= 0.0) constraintViolation[0] = 0.0;

    constraintViolation[1] = 4*x[1] + x[2] - 20.0;
    if (constraintViolation[1] <= 0.0) constraintViolation[1] = 0.0;
  }
  else {
    objArray[0] = -10*exp(-0.2*sqrt(x[0]*x[0] + x[1]*x[1])) -
        10*exp(-0.2*sqrt(x[1]*x[1] + x[2]*x[2]));
    objArray[1] = pow(fabs(x[0]),0.8) + pow(fabs(x[1]),0.8) +
        pow(fabs(x[2]),0.8) + 5.0*sin(x[0]*x[0]*x[0]) +
        5.0*sin(x[1]*x[1]*x[1]) + 5.0*sin(x[2]*x[2]*x[2]);
  }

#pragma omp ordered
  {
    if (globalSetup->savePopulation) {
      outEvals = fopen(globalSetup->saveEvalSolutions, "a");
      for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
        fprintf(outEvals, "%f\t", x[ii]);
      }
      for (ii = 0; ii < globalSetup->finalNoOfObjectives; ii++) {
        fprintf(outEvals, "%f\t", objArray[ii]);
      }
      if (globalSetup->finalNoOfConstraints > 0) {
        for (ii = 0; ii < globalSetup->finalNoOfConstraints; ii++) {
          fprintf(outEvals, "%f\t", constraintViolation[ii]);
        }
        fprintf(outEvals, "%f", *penalty);
      }
      fprintf(outEvals, "\n");
      fflush(outEvals);
      fclose(outEvals);
    }
  }
}

int main(int argc, char *argv[]) {
  // Requires one argument, GA input file
  if(argc != 2) {
    std::cout << "Error: GA input file required" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  run_GA(argv[1]);
  return 1;
}
