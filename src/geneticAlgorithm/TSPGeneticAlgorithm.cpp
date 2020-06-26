//
// Created by gs1010 on 31/05/20.
//
#define TIME

#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "TSPGeneticAlgorithm.h"
#include "../graph/undirectedGraph.h"
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::initializer() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  population.reserve(totalPopulation);
  //! Fill Population
  for (size_t i = 0; i < totalPopulation; i++) {
    std::vector<std::pair<TId, double>> chromosome_prob;
    for (auto &node : graph_->getNodes()) {
      chromosome_prob.push_back(std::make_pair(node, randomProbabilityGenerator()));
    }

    //! sort genes in each chromosome
    //TODO: aggiungere specifica se par o seq per funzioni libreria std
    std::sort(chromosome_prob.begin(),
              chromosome_prob.end(),
              [&](const std::pair<TId, double> &geneA, const std::pair<TId, double> &geneB) {
                return geneA.second < geneB.second;
              });

    std::vector<TId> chromosome;
    chromosome.reserve(chromosome_prob.size());
    for (auto &gene : chromosome_prob) {
      chromosome.push_back(gene.first);
    }
    chromosome.shrink_to_fit();
    population.push_back(chromosome);
  }
  population.shrink_to_fit();

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Initializer time (msecs): %ld\n", msec);
#endif

/*
for (std::vector<TId>& chromosome : population) {
  std::cout << std::endl<< "gene" << std::endl;

  for (TId& gene : chromosome) {
  std::cout<<gene<<std::endl;
  }
}
*/
}
template<typename TId, typename TValue>
double TSPGeneticAlgorithm<TId, TValue>::randomProbabilityGenerator() {
  return unif(gen);
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::selectionReproduction() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif


  //std::cout << std::endl << std::endl << std::endl;
  std::vector<double> randomNumber(multiplier * chromosomeEvals.size());
  //! Compute the avg value
  double avg = 0;
  for (auto chromosomeValue: chromosomeEvals) {
    avg += chromosomeValue.second;
  }
  avg /= chromosomeEvals.size();

  //! Compute the roulette wheel value
  //! Idea of the algorithm taken from "A genetic algorithm tutorial" of Darrell Whitley
  chromosomeEvals[0].second = (1 / chromosomeEvals[0].second) / (1 / avg);

  for (size_t i = 1; i < chromosomeEvals.size(); i++) {
    chromosomeEvals[i].second = chromosomeEvals[i - 1].second + (1 / chromosomeEvals[i].second) / (1 / avg);
  }
  //! Generation of the intermediatePopulation
  std::uniform_real_distribution<double> selectionDistribution{0, chromosomeEvals[chromosomeEvals.size() - 1].second};
  randomNumber[0] = selectionDistribution(gen);
  for (size_t i = 1; i < multiplier * chromosomeEvals.size(); i++) {
    randomNumber[i] = fmod(selectionDistribution(gen) + randomNumber[i - 1],
                           chromosomeEvals[chromosomeEvals.size() - 1].second);
  }

  std::sort(&randomNumber[0], &randomNumber[multiplier * chromosomeEvals.size()]);

  size_t chromosomeEvalsIndex = 0;
  //! Using Stochastic universal sampling (also known as "roulette wheel selection")
  for (size_t i = 0; i < multiplier * chromosomeEvals.size();) {
    if (chromosomeEvals[chromosomeEvalsIndex].second >= randomNumber[i]) {
      intermediatePopulation.push_back(population[chromosomeEvals[chromosomeEvalsIndex].first]);
      //std::cout << chromosomeEvals[chromosomeEvalsIndex].first << std::endl;
      i++;
    } else {
      chromosomeEvalsIndex++;
    }
  }

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("SelectionReproduction time (msecs): %ld\n", msec);
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::crossover() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  std::random_shuffle(intermediatePopulation.begin(), intermediatePopulation.end());
  // Possibile scegliere prima le coppie che dovranno essere usate nella fase di crossover, inserirle in un array secondario e modificarle in modo da
  // usare la vectorization
  std::uniform_int_distribution<size_t> crossoverDistribution{0, graph_->getNodesSize()};
  population.clear();
  size_t selected = 0;
  for (size_t i = 0; (i < intermediatePopulation.size() - 1) && (
      (intermediatePopulation.size() - 1) * crossoverProbability
          - selected) > 0; i++) {
    //! Scanning and selection algorithm 3.4 Prof. Ferragina notes of Algorithm Engineering course
    if (unif(gen)
        <= ((intermediatePopulation.size() - 1) * crossoverProbability - selected)
            / (intermediatePopulation.size() - i)) {
      selected++;
      size_t elemFirstChromosomeIndexA = crossoverDistribution(gen);
      size_t elemFirstChromosomeIndexB = crossoverDistribution(gen);

      if (elemFirstChromosomeIndexA > elemFirstChromosomeIndexB) {
        std::swap(elemFirstChromosomeIndexA, elemFirstChromosomeIndexB);
      }

      //! Adding to the crossover chromosome the first elemFirstChromosome genes from the chromosome A
      std::vector<TId> chromosomePop(&intermediatePopulation[i][elemFirstChromosomeIndexA],
                                     &intermediatePopulation[i][elemFirstChromosomeIndexB]);

/*
      for (size_t j = 0; j < elemFirstChromosome; j++) {
        chromosomePop.push_back(intermediatePopulation[i][j]);
      }
*/
      std::vector<TId> currentMissingIdA(&intermediatePopulation[i][0],
                                         &intermediatePopulation[i][elemFirstChromosomeIndexA]);

      std::vector<TId> currentMissingIdB(&intermediatePopulation[i][elemFirstChromosomeIndexB],
                                         &intermediatePopulation[i][intermediatePopulation[i].size()]);
      //! Adding to the crossover chromosome the remaining (chromosome.size() - elemFirstChromosome + 1) genes with
      //! the order that they appear in chromosome B
      for (size_t j = 0; j < intermediatePopulation[i + 1].size() /*&& currentMissingId.size()*/; j++) {
        auto itA = std::find(currentMissingIdA.begin(),
                            currentMissingIdA.end(),
                            intermediatePopulation[i + 1][j]);
        auto itB = std::find(currentMissingIdB.begin(),
                             currentMissingIdB.end(),
                             intermediatePopulation[i + 1][j]);

        if (itA != currentMissingIdA.end() || itB != currentMissingIdB.end()) {
          if(itA != currentMissingIdA.end()){
          chromosomePop.push_back(*itA);
          //! Possibile miglioria quando si usano molti nodi
          currentMissingIdA.erase(itA);
          }else{
            chromosomePop.push_back(*itB);
            //! Possibile miglioria quando si usano molti nodi
            currentMissingIdB.erase(itB);
          }
        }
      }

      chromosomePop.shrink_to_fit();
      population.push_back(chromosomePop);
    }
  }
  population.shrink_to_fit();

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Crossover time (msecs): %ld\n", msec);
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::mutation() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif
  std::uniform_int_distribution<size_t> mutationDistribution{0, graph_->getNodesSize() - 1};

  for (size_t i = 0; i < intermediatePopulation.size() - 1; i++) {

    double tmp = unif(gen);
    if (tmp <= mutationProbability) {

      size_t positionGenesA = mutationDistribution(gen);
      size_t positionGenesB = mutationDistribution(gen);

      std::vector<TId> chromosomePop(intermediatePopulation[i]);
      std::swap(chromosomePop[positionGenesA], chromosomePop[positionGenesB]);

      chromosomePop.shrink_to_fit();
      population.push_back(chromosomePop);
    }
  }
  population.shrink_to_fit();

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Mutation time (msecs): %ld\n", msec);
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::run(int iteration) {

  if (graph_ == nullptr) {
    setRandomGraph(500);
  }

  initializer();

  evaluate();
  //std::cout << chromosomeEvals[0].second << std::endl;

#ifdef VALUES
  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }
#endif


  for (int i = 0; i < iteration; i++) {
#ifdef TIME
    auto start = std::chrono::high_resolution_clock::now();
#endif
    selectionReproduction();
    crossover();
    mutation();
    adjustPopulation();
    intermediatePopulation.clear();
    chromosomeEvals.clear();
    rankedPopulation.clear();
    evaluate();
    //std::cout << chromosomeEvals[0].second << std::endl;
#ifdef TIME
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto msecs = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    printf("Intermediate time (msecs): %ld %d\n", msecs, i);
#endif

  }
#ifdef VALUES
  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::evaluate() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  size_t currentChromosomeIndex = 0;
  for (std::vector<TId> &chromosome : population) {
    double eval = 0;
    for (size_t i = 0; i < chromosome.size() - 1; i++) {
      eval += graph_->getValueEdge(chromosome[i], chromosome[i + 1]);
    }
    eval += graph_->getValueEdge(chromosome[chromosome.size() - 1], chromosome[0]);
    chromosomeEvals.push_back(std::make_pair(currentChromosomeIndex, eval / chromosome.size()));
    currentChromosomeIndex++;
  }

  //TODO: aggiungere specifica se par o seq per funzioni libreria std
  std::sort(chromosomeEvals.begin(),
            chromosomeEvals.end(),
            [&](const std::pair<TId, double> &chromosomeA, const std::pair<size_t, double> &chromosomeB) {
              return chromosomeA.second < chromosomeB.second;
            });

/*
  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }
*/
  for (auto &chromosome : chromosomeEvals) {
    rankedPopulation.push_back(population[chromosome.first]);
  }


#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Evaluate time (msecs): %ld\n", msec);
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::adjustPopulation() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  if (population.size() > rankedPopulation.size()) {
    std::random_shuffle(population.begin(), population.end());
    while (population.size() - rankedPopulation.size()) {
      population.pop_back();
    }
  } else {
    for (size_t i = 0; rankedPopulation.size() - population.size(); i++) {
      population.push_back(rankedPopulation[i]);
    }
  }

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("AdjustPopulation time (msecs): %ld\n", msec);
#endif

}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::setRandomGraph(size_t nNodes) {
  graph_ = new undirectedGraph<TId, TValue>(nNodes);
  graph_->randomInit(seed);
}

template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::setGraph(graph<TId, TValue> *graph) {
  graph_ = graph;
}
template<typename TId, typename TValue>
TSPGeneticAlgorithm<TId, TValue>::TSPGeneticAlgorithm():seed(0), multiplier(1), totalPopulation(500) {
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithm<TId, TValue>::TSPGeneticAlgorithm(int seed_):seed(seed_), multiplier(1), totalPopulation(500) {
  gen.seed(seed);
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithm<TId, TValue>::TSPGeneticAlgorithm(int seed_,
                                                      double crossoverProbability_,
                                                      double mutationProbability_)
    :seed(seed_),
     crossoverProbability(crossoverProbability_),
     mutationProbability(mutationProbability_),
     multiplier(1),
     totalPopulation(500) {
  gen.seed(seed);
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::SetCrossoverProbability(double crossover_probability) {
  crossoverProbability = crossover_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::SetMutationProbability(double mutation_probability) {
  mutationProbability = mutation_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::SetMultiplier(int multiplier) {
  TSPGeneticAlgorithm::multiplier = multiplier;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::SetTotalPopulation(size_t total_population) {
  totalPopulation = total_population;
}
