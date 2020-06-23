//
// Created by gs1010 on 31/05/20.
//

#include <random>
#include <iostream>
#include <algorithm>
#include "TSPGeneticAlgorithm.h"
#include "undirectedGraph.h"
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::initializer(int totalPopulation) {
  graph_ = new undirectedGraph<TId, TValue>(20);
  graph_->randomInit();
  crossoverProbability = 0.2;
  elitePercentage = 0.1;
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
  std::cout << std::endl << std::endl << std::endl;

  int multiplier = 2;
  double randomNumber[multiplier * chromosomeEvals.size()];
  //! Compute the avg value
  double avg = 0;
  for (auto chromosomeValue: chromosomeEvals) {
    avg+= chromosomeValue.second;
  }
  avg /= chromosomeEvals.size();

  //! Compute the roulette wheel value
  //! Idea of the algorithm taken from "A genetic algorithm tutorial" of Darrell Whitley
  double tmp= (1 / chromosomeEvals[0].second) / (1 / avg);
  chromosomeEvals[0].second = (1 / chromosomeEvals[0].second) / (1 / avg);

  for (size_t i = 1; i < chromosomeEvals.size(); i++) {
    tmp= (1 / chromosomeEvals[i].second) / (1 / avg);
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
  for (size_t i = 0; i < multiplier * chromosomeEvals.size();) {
    if (chromosomeEvals[chromosomeEvalsIndex].second >= randomNumber[i]) {
      intermediatePopulation.push_back(population[chromosomeEvals[chromosomeEvalsIndex].first]);
      //std::cout << chromosomeEvals[chromosomeEvalsIndex].first << std::endl;
      i++;
    } else {
      chromosomeEvalsIndex++;
    }
  }



}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::crossover() {

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
      size_t elemFirstChromosome = crossoverDistribution(rd);
      std::vector<TId> chromosomePop;
      chromosomePop.reserve(intermediatePopulation[i].size());
      //! Adding to the crossover chromosome the first elemFirstChromosome genes from the chromosome A
      for (size_t j = 0; j < elemFirstChromosome; j++) {
        chromosomePop.push_back(intermediatePopulation[i][j]);
      }

      std::vector<TId> currentMissingId(&intermediatePopulation[i][elemFirstChromosome],
                                        &intermediatePopulation[i][intermediatePopulation[i].size()]);
      //! Adding to the crossover chromosome the remaining (chromosome.size() - elemFirstChromosome + 1) genes with
      //! the order that they appear in chromosome B
      for (size_t j = elemFirstChromosome; j < intermediatePopulation[i].size(); j++) {
        auto it = std::find_first_of(intermediatePopulation[i + 1].begin(),
                                     intermediatePopulation[i + 1].end(),
                                     currentMissingId.begin(),
                                     currentMissingId.end());
        chromosomePop.push_back(*it);

        currentMissingId.erase(std::find(currentMissingId.begin(), currentMissingId.end(), *it));
      }

      chromosomePop.shrink_to_fit();
      population.push_back(chromosomePop);
    }
  }
  population.shrink_to_fit();

}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::mutation() {

}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::run() {

}
template<typename TId, typename TValue>
void TSPGeneticAlgorithm<TId, TValue>::evaluate() {

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

  for (auto &chromosome : chromosomeEvals) {
    std::cout << chromosome.first << " " << chromosome.second << std::endl;
  }

}
