//
// Created by gs1010 on 26/06/20.
//

#ifndef TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMST_H_
#define TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMST_H_
#include "../graph/graph.h"
#include <unordered_map>
#include <vector>
#include <random>
#include "geneticAlgorithm.h"
template<typename TId, typename TValue>
class TSPGeneticAlgorithmST : public geneticAlgorithm {
 private:
  graph<TId, TValue> *graph_;
  //! TODO: Creare la top list dei migliori cromosomi
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
  int nWorker;
  double randomProbabilityGenerator();
  void adjustPopulation();

 public:
  explicit TSPGeneticAlgorithmST();
  explicit TSPGeneticAlgorithmST(int seed_);
  explicit TSPGeneticAlgorithmST(int seed_, double crossoverProbability_, double mutationProbability_);
  void initializer() override;
  void evaluate() override;
  void selectionReproduction() override;
  void crossover() override;
  void mutation() override;
  void run(int iteration) override;
  void setGraph(graph<TId, TValue> *graph);
  void setRandomGraph(size_t nNodes);
  void SetCrossoverProbability(double crossover_probability);
  void SetMutationProbability(double mutation_probability);
  void SetMultiplier(int multiplier);
  void SetTotalPopulation(size_t total_population);
  void SetNWorker(int n_worker);

};
template
class TSPGeneticAlgorithmST<int, int>;
template
class TSPGeneticAlgorithmST<int, double>;
template
class TSPGeneticAlgorithmST<unsigned long long, double>;

#endif //TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMST_H_
