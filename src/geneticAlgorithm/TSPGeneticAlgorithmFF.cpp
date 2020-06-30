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
  ff::ff_farm initializationFarm;
  std::vector<ff::ff_node *> initializerWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    initializerWorkers.push_back(new Initializer(graph_->getNodes(), gen));
  initializationFarm.add_workers(initializerWorkers);
  auto *emitterInitializer = new EmitterInitializer(totalPopulation, nWorker);
  initializationFarm.add_emitter(emitterInitializer);

  ff::ff_farm evaluateFarm;
  std::vector<ff::ff_node *> evaluateWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    evaluateWorkers.push_back(new Evaluate(graph_, gen));
  evaluateFarm.add_workers(evaluateWorkers);
  auto *cInitEEval = new collectorInitializerEmitterEvaluate(totalPopulation,
                                                             iteration,
                                                             totalPopulation * crossoverProbability
                                                                 + totalPopulation * mutationProbability);
  evaluateFarm.add_emitter(cInitEEval);

  auto *sR = new selectionReproductionNode(totalPopulation, multiplier, gen);
  evaluateFarm.add_collector(sR);
  //ff::ff_Pipe<size_t> pipe(initializationFarm, evaluateFarm);

  ff::ff_farm crossoverMutationFarm;
  std::vector<ff::ff_node *> crossoverMutationWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    crossoverMutationWorkers.push_back(new CrossoverMutation(gen));
  crossoverMutationFarm.add_workers(crossoverMutationWorkers);
  auto *ECrossoverMutation = new EmitterCrossoverMutation(crossoverProbability, mutationProbability, gen);
  crossoverMutationFarm.add_emitter(ECrossoverMutation);

  /*
  ff::ff_farm mutationFarm;
  std::vector<ff::ff_node *> mutationWorkers;
  for (size_t i = 0; i < nWorker; ++i)
    mutationWorkers.push_back(new Mutation(gen));
  mutationFarm.add_workers(mutationWorkers);
  auto *EMutation = new EmitterMutation(mutationProbability, gen);
  mutationFarm.add_emitter(EMutation);
*/
  ff::ff_Pipe<std::pair<std::vector<TId> *, double>> pipe_with_feedback(evaluateFarm, crossoverMutationFarm);
  pipe_with_feedback.wrap_around();
  ff::ff_Pipe<size_t> pipe(initializationFarm, pipe_with_feedback);

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