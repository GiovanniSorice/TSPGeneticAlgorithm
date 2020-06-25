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
  graph<TId, TValue> *graph_;
  //! TODO: Creare la top list dei migliori cromosomi
  //! population of chromosome
  std::vector<std::vector<TId >> rankedPopulation;
  std::vector<std::vector<TId >> population;
  std::vector<std::vector<TId>> intermediatePopulation;
  //! random probability generator
  std::random_device rd;
  std::mt19937 gen{rd()};
  std::uniform_real_distribution<double> unif{0, 1};
  //! Chromosome current value
  std::vector<std::pair<size_t, double>> chromosomeEvals;
  //! Crossover probability
  double crossoverProbability;
  //! Crossover probability
  double mutationProbability;
  //! Elite percentage
  double elitePercentage;
  int seed;
  double randomProbabilityGenerator();
  void adjustPopulation();

 public:
  void initializer(size_t totalPopulation) override;
  void evaluate() override;
  void selectionReproduction() override;
  void crossover() override;
  void mutation() override;
  void run() override;
  void setGraph(graph<TId, TValue>* graph);
  void setRandomGraph(size_t nNodes, int seed=0);

};

template
class TSPGeneticAlgorithm<int, int>;
template
class TSPGeneticAlgorithm<int, double>;
template
class TSPGeneticAlgorithm<unsigned long long, double>;

#endif //TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
