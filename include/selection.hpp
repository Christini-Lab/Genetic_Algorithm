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

#ifndef _SELECTION_H
#define _SELECTION_H

#include "population.hpp"
#include "globalSetup.hpp"
#include "random.hpp"
#include <iostream>

class Selection {
 protected:
  const Population *pop;
  void selectionQuickSort(const Population *pop, int *randomArray, int left, int right);
  void swap(int &ii, int &jj);
 public:
  Selection(const Population *pop);
  virtual void select(int *matingPool)=0;
  int betterIndividual(Individual *guy1, Individual *guy2);
};

class TournamentSelection: public Selection {
 private:
  int tournamentSize;
 public:
  TournamentSelection(const int size, const Population *pop);
  void select(int *matingPool);
};

class TournamentSelectionWithReplacement: public Selection {
 private:
  int tournamentSize;
 public:
  TournamentSelectionWithReplacement(const int size, const Population *pop);
  void select(int *matingPool);
};

class StochasticUniversalSelection: public Selection {
 public:
  StochasticUniversalSelection(const Population *pop);
  void select(int *matingPool);
};

class RouletteWheelSelection: public Selection {
 public:
  RouletteWheelSelection(const Population *pop);
  void select(int *matingPool);

};

class TruncationSelection: public Selection {
 private:
  int selectionPressure;
 public:
  TruncationSelection(const int selP, const Population *pop);
  void select(int *matingPool);

};

#endif
