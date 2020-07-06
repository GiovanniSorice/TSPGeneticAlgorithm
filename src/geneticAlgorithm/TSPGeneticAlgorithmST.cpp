//
// Created by gs1010 on 26/06/20.
//

#include "TSPGeneticAlgorithmST.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <deque>
#include <thread>
#include "../graph/undirectedGraph.h"

//#define VALUES
//#define TIME

template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::initializer() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  std::deque<std::thread> coda;

  population.resize(totalPopulation);

  //! Fill Population
  for (auto range : populationRanges) {
    coda.emplace_back(std::thread([&](size_t start, size_t end) {
      for (size_t i = start; i < end; i++) {
        std::vector<std::pair<TId, double>> chromosome_prob;
        const std::vector<TId> &nodes = graph_->getNodes();

        chromosome_prob.reserve(nodes.size());

        for (auto &node : nodes) {
          chromosome_prob.emplace_back(node, randomProbabilityGenerator());
        }

        //! sort genes in each chromosome
        std::sort(chromosome_prob.begin(),
                  chromosome_prob.end(),
                  [&](const std::pair<TId, double> &geneA, const std::pair<TId, double> &geneB) {
                    return geneA.second < geneB.second;
                  });

        std::vector<TId> chromosome;
        chromosome.resize(chromosome_prob.size());
        for (size_t j = 0; j < chromosome_prob.size(); j++) {
          chromosome[j] = chromosome_prob[j].first;
        }
        chromosome.shrink_to_fit();
        population[i] = chromosome;
      }

    }, range.first, range.second));
  }

  for (auto &it : coda) {
    if (it.joinable())
      it.join();
  }

  coda.clear();

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
double TSPGeneticAlgorithmST<TId, TValue>::randomProbabilityGenerator() {
  return unif(gen);
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::selectionReproduction() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif


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
      intermediatePopulation.emplace_back(population[chromosomeEvals[chromosomeEvalsIndex].first]);
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
void TSPGeneticAlgorithmST<TId, TValue>::crossover() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  size_t totalCrossover = intermediatePopulation.size() * crossoverProbability;
  std::vector<std::pair<size_t, size_t>> ranges;
  size_t startRange = 0;
  size_t step = totalCrossover / nWorker;
  size_t remaining = totalCrossover - step * nWorker;
  size_t endRange = 0;

  while (remaining != 0) {
    startRange = endRange;
    endRange += step + 1;
    ranges.emplace_back(startRange, endRange);
    remaining--;
  }

  while (endRange < totalCrossover) {
    startRange = endRange;
    endRange += step;
    ranges.emplace_back(startRange, endRange);
  }

  std::random_shuffle(intermediatePopulation.begin(), intermediatePopulation.end());

  //! Selection of the chromosome used for the crossover (intermediatePopulation.size() - 2 because i exchange the value between i and i+1,
  //! so if i use only intermediatePopulation.size() - 1 the index will go out of bound and create some error)
  std::uniform_int_distribution<size_t> SelectionCrossoverDistribution{0, intermediatePopulation.size() - 2};


  std::vector<size_t> crossoverPosition(totalCrossover);

  for (size_t i = 0; i < totalCrossover; i++) {
    crossoverPosition[i] = SelectionCrossoverDistribution(gen);
  }


  std::uniform_int_distribution<size_t> crossoverDistribution{0, graph_->getNodesSize()};
  population.clear();
  population.resize(totalCrossover);

  std::deque<std::thread> coda;
  for (auto & range : ranges) {

    coda.emplace_back(std::thread([&](size_t start, size_t end) {

      for (size_t i = start; i < end; i++) {

        size_t elemFirstChromosomeIndexA = crossoverDistribution(gen);
        size_t elemFirstChromosomeIndexB = crossoverDistribution(gen);

        if (elemFirstChromosomeIndexA > elemFirstChromosomeIndexB) {
          std::swap(elemFirstChromosomeIndexA, elemFirstChromosomeIndexB);
        }

        //! Adding to the crossover chromosome the first elemFirstChromosome genes from the chromosome A
        std::vector<TId> chromosomePop(&intermediatePopulation[crossoverPosition[i]][elemFirstChromosomeIndexA],
                                       &intermediatePopulation[crossoverPosition[i]][elemFirstChromosomeIndexB]);

        auto interPopPtrBeginA= &intermediatePopulation[crossoverPosition[i]][0];
        auto interPopPtrEndA = &intermediatePopulation[crossoverPosition[i]][elemFirstChromosomeIndexA];

        auto interPopPtrBeginB = &intermediatePopulation[crossoverPosition[i]][elemFirstChromosomeIndexB];
        auto interPopPtrEndB = &intermediatePopulation[crossoverPosition[i]][intermediatePopulation[crossoverPosition[i]].size()];

        //! Adding to the crossover chromosome the remaining (chromosome.size() - elemFirstChromosome + 1) genes with
        //! the order that they appear in chromosome B
        for (size_t j = 0; j < intermediatePopulation[crossoverPosition[i] + 1].size() /*&& currentMissingId.size()*/;
             j++) {
          auto itA = std::find(interPopPtrBeginA,
                               interPopPtrEndA,
                               intermediatePopulation[crossoverPosition[i] + 1][j]);
          auto itB = std::find(interPopPtrBeginB,
                               interPopPtrEndB,
                               intermediatePopulation[crossoverPosition[i] + 1][j]);

          if (itA != interPopPtrEndA || itB != interPopPtrEndB) {
            if (itA != interPopPtrEndA) {
              chromosomePop.push_back(*itA);
            } else {
              chromosomePop.push_back(*itB);
            }
          }
        }

        chromosomePop.shrink_to_fit();
        population[i] = chromosomePop;
      }
    }, range.first, range.second));
  }


  for (auto &it : coda) {
    if (it.joinable())
      it.join();
  }

  coda.clear();

  population.shrink_to_fit();

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Crossover time (msecs): %ld\n", msec);
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::mutation() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif
  std::uniform_int_distribution<size_t> mutationDistribution{0, graph_->getNodesSize() - 1};

  for (size_t i = 0; i < intermediatePopulation.size() - 1; i++) {

    if (unif(gen) <= mutationProbability) {

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
void TSPGeneticAlgorithmST<TId, TValue>::run(int iteration) {

  if (graph_ == nullptr) {
    setRandomGraph(500);
  }
  setUpRanges();
  initializer();

#ifdef VALUES
  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }
#endif

  for (int i = 0; i < iteration; i++) {
#ifdef TIME
    auto start = std::chrono::high_resolution_clock::now();
#endif
    evaluate();
#ifdef VALUES
    std::cout <<"Best value: "<<chromosomeEvals[0].second << std::endl;
#endif
    selectionReproduction();
    crossover();
    mutation();
    adjustPopulation();
    intermediatePopulation.clear();
    chromosomeEvals.clear();
    rankedPopulation.clear();
#ifdef TIME
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto msecs = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    printf("Intermediate time (msecs): %ld %d\n", msecs, i);
#endif
  }
  evaluate();
  //std::cout << "Best value ST: " << chromosomeEvals[0].second << std::endl;

#ifdef VALUES
  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }
#endif
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::evaluate() {
#ifdef TIME
  auto start = std::chrono::high_resolution_clock::now();
#endif

  chromosomeEvals.resize(population.size());

  std::deque<std::thread> coda;

  for (auto range: populationRanges) {
    coda.emplace_back(std::thread([&](size_t start, size_t end) {
      for (size_t j = start; j < end; j++) {
        std::vector<TId> &chromosome = population[j];
        double eval = 0;
        for (size_t k = 0; k < chromosome.size() - 1; k++) {
          eval += graph_->getValueEdge(chromosome[k], chromosome[k + 1]);
        }
        eval += graph_->getValueEdge(chromosome[chromosome.size() - 1], chromosome[0]);
        chromosomeEvals[j] = std::make_pair(j, eval / chromosome.size());
      }
    }, range.first, range.second));
  }

  for (auto &it : coda) {
    if (it.joinable())
      it.join();
  }

  coda.clear();
  chromosomeEvals.shrink_to_fit();
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
  rankedPopulation.reserve(chromosomeEvals.size());
  for (auto &chromosome : chromosomeEvals) {
    rankedPopulation.push_back(population[chromosome.first]);
  }
  rankedPopulation.shrink_to_fit();

#ifdef TIME
  auto elapsed = std::chrono::high_resolution_clock::now() - start;
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  printf("Evaluate time (msecs): %ld\n", msec);
#endif
}

template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::adjustPopulation() {
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
void TSPGeneticAlgorithmST<TId, TValue>::setUpRanges() {

  size_t startRange = 0;
  size_t step = totalPopulation / nWorker;
  size_t remaining = totalPopulation - step * nWorker;
  size_t endRange = 0;

  while (remaining != 0) {
    startRange = endRange;
    endRange += step + 1;
    populationRanges.emplace_back(startRange, endRange);
    remaining--;
  }

  while (endRange < totalPopulation) {
    startRange = endRange;
    endRange += step;
    populationRanges.emplace_back(startRange, endRange);
  }

}

template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::setRandomGraph(size_t nNodes) {
  graph_ = new undirectedGraph<TId, TValue>(nNodes);
  graph_->randomInit(seed);
}

template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::setGraph(graph<TId, TValue> *graph) {
  graph_ = graph;
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmST<TId, TValue>::TSPGeneticAlgorithmST():seed(0), multiplier(1), totalPopulation(500), nWorker(1) {
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmST<TId, TValue>::TSPGeneticAlgorithmST(int seed_):seed(seed_), multiplier(1), totalPopulation(500),
                                                                     nWorker(1) {
  gen.seed(seed);
  crossoverProbability = unif(gen);
  mutationProbability = unif(gen);
}
template<typename TId, typename TValue>
TSPGeneticAlgorithmST<TId, TValue>::TSPGeneticAlgorithmST(int seed_,
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
void TSPGeneticAlgorithmST<TId, TValue>::SetCrossoverProbability(double crossover_probability) {
  crossoverProbability = crossover_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::SetMutationProbability(double mutation_probability) {
  mutationProbability = mutation_probability;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::SetMultiplier(int multiplier) {
  TSPGeneticAlgorithmST::multiplier = multiplier;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::SetTotalPopulation(size_t total_population) {
  totalPopulation = total_population;
}
template<typename TId, typename TValue>
void TSPGeneticAlgorithmST<TId, TValue>::SetNWorker(int n_worker) {
  nWorker = n_worker;
}
