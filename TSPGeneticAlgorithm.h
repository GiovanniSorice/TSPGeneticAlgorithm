//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#define TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#include "graph.h"
#include <unordered_map>
#include <vector>
#include <random>
#include "geneticAlgorithm.h"
template<typename TId, typename TValue>
class TSPGeneticAlgorithm : public geneticAlgorithm {
 private:
  graph<TId, TValue>* graph_;
  //! population of chromosome
  std::vector<std::vector<std::pair<TId, double>>> population;
  std::vector<std::vector<std::pair<TId, double>>> intermediatePopulation;
  //! random probability generator
  std::random_device rd;
  std::mt19937 gen{rd()};
  std::uniform_real_distribution<double> unif{0, 1};
  //! Chromosome current value
  std::vector<std::pair<size_t, double>> chromosomeEvals;
  //! Crossover probability
  double crossoverProbability;
  //! Elite percentage
  double elitePercentage;

  double randomProbabilityGenerator();
 public:
  void initializer(int totalPopulation) override;
  void evaluate() override;
  void selectionReproduction() override;
  void crossover() override;
  void mutation() override;
  void run() override;

};

template
class TSPGeneticAlgorithm<int, int>;
template
class TSPGeneticAlgorithm<int, double>;
template
class TSPGeneticAlgorithm<unsigned long long, double>;

#endif //TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
