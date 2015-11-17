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

#include "../include/globalSetup.hpp"

GlobalSetup::~GlobalSetup() {
  int ii;
  delete []variableTypes;
  for (ii = 0; ii < noOfDecisionVariables; ii++)
    delete [] variableRanges[ii];
  delete [] variableRanges;
}
