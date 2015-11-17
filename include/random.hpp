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

#ifndef _random_hpp
#define _random_hpp
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS1 1.2e-7
#define RNMX (1.0-EPS1)

class Random
{
 private:
  long seed;
  int iset;
  double gset;
 public:
  Random();
  Random(const long seedValue);
  ~Random();
  long getSeed() const { return seed; }
  void setSeed(long seedValue) {seed = seedValue;}
  double random01(void);
  bool flip(double probability);
  double boundedRandom(double minimum, double maximum);
  int boundedIntegerRandom(int minimum, int maximum);
  double normalRandom(double standardDeviation);
  void shuffleArray(int *theArray, int arrayLength);
  long makeSeed();
};

#endif
