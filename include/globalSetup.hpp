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

#ifndef _GLOBALSETUP_H
#define _GLOBALSETUP_H

#include <iostream>
#include <fstream>

#define SWAP(A,B)    (A^=B^=A^=B)
#define ZERO         1.0E-10
#define INFTY        1.0E+20;
#define IDOMJ        1
#define JDOMI       -1
#define INONDOMJ     0
#define GUYS         0
#define NEWGUYS      1
#define COMBINEDGUYS 2
#define NO           0
#define YES          1
#define OFF          0
#define ON           1
#define OBJECTIVE    0
#define FITNESS      1
#define REFLECT -1.0
#define EXPAND   2.0
#define CONTRACT 0.5
#define RECONSTRUCT 0.0

enum OptimType {Minimization, Maximization};
enum GAType {SGA, NSGA};
enum VariableType {Integer, Real};
enum ConstraintMethod {NoConstraints, Penalty, Tournament};
enum SelectionType {TournamentWOR, SUS, Truncation, RouletteWheel,
                    TournamentWR};
enum XOverType {OnePoint, TwoPoint, Uniform, SBX};
enum MutationType {Selective, Genewise, Polynomial};
enum NichingType {NoNiching, Sharing, RTS, DeterministicCrowding};
enum LocalSearchMethod {NoLocalSearch, SimplexSearch};
enum PenaltyFunction{Linear, Quadratic};
enum ScalingMethod {NoScaling,Ranking, SigmaScaling};
enum StoppingCriterias {
    NoOfEvaluations,		// 1
    FitnessVariance,		// 1
    BestFitness,			// 2
    AverageFitness,		// 2
    AverageObjective,		// 2
    ChangeInBestFitness,	        // 2
    ChangeInAvgFitness,		// 2
    ChangeInFitnessVar,		// 2
    ChangeInBestObjective,	// 2
    ChangeInAvgObjective,	        // 2
    NoOfFronts,			// 3
    NoOfGuysInFirstFront,	        // 3
    ChangeInNoOfFronts		// 3
};


class GlobalSetup {
public:
    // Configurator specific information
    int noOfDecisionVariables;
    int noOfRawConstraints;
    int noOfRawObjectives;
    int noOfLinearObjectiveCombinations;
    int noOfLinearConstraintCombinations;
    int finalNoOfObjectives;
    int finalNoOfConstraints;
    double **linearObjectiveCombinationWeights;
    double **linearConstraintCombinationWeights;
    double *penaltyWeights;
    double maxLocalTolerance;
    int maxLocalEvaluations;
    double initialLocalPenaltyParameter;
    double localUpdateParameter;
    // GA Specific Information

    VariableType *variableTypes;
    double **variableRanges;

    OptimType *typeOfOptimizations;
    GAType gaType;

    int populationSize;
    int maxGenerations;

    double xOverProbability, mutationProbability;
    double lamarckianProbability, localSearchProbability;
    LocalSearchMethod localSearchMethod;
    SelectionType selectionType;
    XOverType xOverType;
    MutationType mutationType;
    NichingType nichingType;
    ConstraintMethod constraintMethod;
    PenaltyFunction penaltyFunction;

    ScalingMethod scalingMethod;
    void	*scalingParameters;

    void *selectionParameters, *xOverParameters;
    void *nichingParameters, *mutationParameters;
    int	genNumWindow;			// 5
    int	noOfStoppingCriterias;
    double *stoppingParameter;
    StoppingCriterias *otherStoppingCriteria;
    double replaceProportion;		// 0.9
    bool loadPopulation;
    bool evaluateAgain;
    char *populationFileName;
    bool savePopulation;
    char *saveEvalSolutions;

    ~GlobalSetup();
};


#endif
