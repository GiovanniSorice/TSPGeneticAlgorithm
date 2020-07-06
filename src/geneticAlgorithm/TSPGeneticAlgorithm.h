//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#define TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#include "../graph/graph.h"
#include <unordered_map>
#include <vector>
#include <random>
#include "geneticAlgorithm.h"
template<typename TId, typename TValue>
class TSPGeneticAlgorithm : public geneticAlgorithm {
 private:
  graph<TId, TValue> *graph_;
  //! population of chromosome
  std::vector<std::vector<TId >> population;
  std::vector<std::vector<TId >> rankedPopulation;
  std::vector<std::vector<TId>> intermediatePopulation;
  //! Chromosome current value
  std::vector<std::pair<size_t, double>> chromosomeEvals;
  //! random probability generator
  std::random_device rd;
  std::mt19937 gen{rd()};
  std::uniform_real_distribution<double> unif{0, 1};
  int seed;
  //! Crossover probability
  double crossoverProbability;
  //! Crossover probability
  double mutationProbability;
  int multiplier;
  size_t totalPopulation;
  double randomProbabilityGenerator();
  void adjustPopulation();
  void initializer();
  void evaluate();
  void selectionReproduction();
  void crossover();
  void mutation();

 public:
  explicit TSPGeneticAlgorithm();
  explicit TSPGeneticAlgorithm(int seed_);
  explicit TSPGeneticAlgorithm(int seed_, double crossoverProbability_, double mutationProbability_);
  void run(int iteration) override;
  void setGraph(graph<TId, TValue>* graph);
  void setRandomGraph(size_t nNodes);
  void SetCrossoverProbability(double crossover_probability);
  void SetMutationProbability(double mutation_probability);
  void SetMultiplier(int multiplier);
  void SetTotalPopulation(size_t total_population);

};

template
class TSPGeneticAlgorithm<int, int>;
template
class TSPGeneticAlgorithm<int, double>;
template
class TSPGeneticAlgorithm<unsigned long long, double>;

#endif //TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
