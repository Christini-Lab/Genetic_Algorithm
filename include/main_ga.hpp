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

#ifndef MAINGA_HPP
#define MAINGA_HPP

#include <fstream>

static char* readOneLine(char *pcBuf, int iMaxSize, FILE *fStream);
void setup_GA(char *argv);
int run_GA(char *argv);
int run_GA();

#endif
