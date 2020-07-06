//
// Created by gs1010 on 27/06/20.
//


#include "TSPGeneticAlgorithmFF.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <deque>
#include <thread>
#include "../graph/undirectedGraph.h"

// #define VALUES
#define TIME

template<typename TId, typename TValue>
double TSPGeneticAlgorithmFF<TId, TValue>::randomProbabilityGenerator() {
  return unif(gen);
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::run(int iteration) {
  std::vector<std::unique_ptr<ff::ff_node>> initializerWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    initializerWorkers.push_back(std::unique_ptr<Initializer>(new Initializer(graph_->getNodes(), gen)));
  auto emitterInitializer = std::unique_ptr<EmitterInitializer>(new EmitterInitializer(totalPopulation, nWorker));
  ff::ff_Farm<size_t> initializationFarm(std::move(initializerWorkers), std::move(emitterInitializer));


  std::vector<std::unique_ptr<ff::ff_node>> evaluateWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    evaluateWorkers.push_back(std::unique_ptr<Evaluate>(new Evaluate(graph_, gen)));
  auto cInitEEval = std::unique_ptr<collectorInitializerEmitterEvaluate>(new collectorInitializerEmitterEvaluate(totalPopulation,
                                                             iteration,
                                                             totalPopulation * crossoverProbability
                                                                 + totalPopulation * mutationProbability));
  auto sR = std::unique_ptr<selectionReproductionNode>(new selectionReproductionNode(totalPopulation, multiplier, gen));

  ff::ff_Farm evaluateFarm(std::move(evaluateWorkers), std::move(cInitEEval), std::move(sR));

  std::vector<std::unique_ptr<ff::ff_node>> crossoverMutationWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    crossoverMutationWorkers.push_back(std::unique_ptr<CrossoverMutation>(new CrossoverMutation(gen)));
  auto ECrossoverMutation = std::unique_ptr<EmitterCrossoverMutation>(new EmitterCrossoverMutation(crossoverProbability, mutationProbability, gen));
  ff::ff_Farm crossoverMutationFarm(std::move(crossoverMutationWorkers), std::move(ECrossoverMutation));


  ff::ff_Pipe<std::pair<std::vector<TId> *, double>> pipe_with_feedback(evaluateFarm, crossoverMutationFarm);
  pipe_with_feedback.wrap_around();
  pipe_with_feedback.isset_cleanup_nodes();
  ff::ff_Pipe<size_t> pipe(initializationFarm, pipe_with_feedback);
  pipe.isset_cleanup_nodes();
  if (pipe.run_and_wait_end() < 0)
    std::cout << " running myFarm" << std::endl;


}

template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::setRandomGraph(size_t nNodes) {
  graph_ = new undirectedGraph<TId, TValue>(nNodes);
  graph_->randomInit(seed);
}

template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::setGraph(graph<TId, TValue> *graph) {
  graph_ = graph;
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmFF<TId, TValue>::TSPGeneticAlgorithmFF():seed(0), multiplier(1), totalPopulation(500), nWorker(1) {
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmFF<TId, TValue>::TSPGeneticAlgorithmFF(int seed_):seed(seed_), multiplier(1), totalPopulation(500),
                                                                     nWorker(1) {
  gen.seed(seed);
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmFF<TId, TValue>::TSPGeneticAlgorithmFF(int seed_,
                                                          double crossoverProbability_,
                                                          double mutationProbability_)
    :seed(seed_),
     crossoverProbability(crossoverProbability_),
     mutationProbability(mutationProbability_),
     multiplier(1),
     totalPopulation(500),
     nWorker(1) {
  gen.seed(seed);
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::SetCrossoverProbability(double crossover_probability) {
  crossoverProbability = crossover_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::SetMutationProbability(double mutation_probability) {
  mutationProbability = mutation_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::SetMultiplier(int multiplier) {
  TSPGeneticAlgorithmFF::multiplier = multiplier;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::SetTotalPopulation(size_t total_population) {
  totalPopulation = total_population;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmFF<TId, TValue>::SetNWorker(int n_worker) {
  nWorker = n_worker;
}