//
// Created by gs1010 on 27/06/20.
//

#ifndef TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMFF_H_
#define TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMFF_H_
#include <ff/ff.hpp>
#include "../graph/graph.h"
#include <ff/pipeline.hpp>
#include <ff/farm.hpp>
#include <random>
#include "geneticAlgorithm.h"
#include <iostream>
template<typename TId, typename TValue>
class TSPGeneticAlgorithmFF : public geneticAlgorithm {
 private:

  struct EmitterInitializer : ff::ff_monode_t<size_t> {

    explicit EmitterInitializer(const size_t totalPopulation_, const size_t nWorker_) :
        totalPopulation(totalPopulation_), nWorker(nWorker_) {}
    size_t *svc(size_t *) override {
      size_t totalSections = totalPopulation / nWorker;
      size_t remaining = totalPopulation - totalSections * nWorker;
      for (size_t j = 0; j < nWorker; ++j) {
        ff_send_out(new size_t(totalSections + (remaining ? 1 : 0)));
        if (remaining) {
          remaining--;
        }
      }
      return EOS;
    }

    size_t totalPopulation;
    size_t nWorker;
  };

  struct Initializer : ff::ff_node_t<size_t> {
    explicit Initializer(const std::vector<TId> &nodes_, std::mt19937 &gen_) : nodes(nodes_), gen(gen_) {}
    size_t *svc(size_t *generate) override {
      std::uniform_real_distribution<double> uniformDistribution{0, 1};
      for (size_t j = 0; j < *generate; ++j) {
        std::vector<std::pair<TId, double>> chromosome_prob;
        for (auto &node : nodes) {
          chromosome_prob.push_back(std::make_pair(node, uniformDistribution(gen)));
        }

        //! sort genes in each chromosome
        //TODO: aggiungere specifica se par o seq per funzioni libreria std
        std::sort(chromosome_prob.begin(),
                  chromosome_prob.end(),
                  [&](const std::pair<TId, double> &geneA, const std::pair<TId, double> &geneB) {
                    return geneA.second < geneB.second;
                  });

        auto *chromosome = new std::vector<TId>;
        chromosome->reserve(chromosome_prob.size());
        for (auto &gene : chromosome_prob) {
          chromosome->push_back(gene.first);
        }
        chromosome->shrink_to_fit();

        this->ff_send_out(chromosome);
      }
      delete generate;
      return EOS;
    }

    const std::vector<TId> &nodes;
    std::mt19937 &gen;
  };
  struct collectorInitializerEmitterEvaluate : ff::ff_minode_t<std::vector<TId>> {
    explicit collectorInitializerEmitterEvaluate(const size_t totalPopulation_, const size_t totalIterations_,
                                                 const size_t totalNewPopulation_) :
        totalPopulation(totalPopulation_), totalIterations(totalIterations_), currentIteration(0), totalNewPopulation(
        totalNewPopulation_) {}
    std::vector<TId> *svc(std::vector<TId> *in) override {
      population.push_back(in);

      if (currentIteration == totalIterations && population.size() == totalNewPopulation) {
        std::sort(rankedPopulation.begin(), rankedPopulation.end(),
                  [&](const std::pair<std::vector<TId> *, double> *chromosomeA,
                      const std::pair<std::vector<TId> *, double> *chromosomeB) {
                    //! I order backwards so I can take advantage of the pop_back later
                    return chromosomeA->second > chromosomeB->second;
                  });

        std::cout << "Best value: " << rankedPopulation.at(totalPopulation-1)->second << std::endl;

        if (population.size() > rankedPopulation.size()) {
          std::random_shuffle(population.begin(), population.end());
          while (population.size() - rankedPopulation.size()) {
            delete population.back();
            population.pop_back();
          }
        } else {
          while (totalPopulation - population.size()) {
            population.push_back(rankedPopulation.back()->first);
            delete rankedPopulation.back();
            rankedPopulation.pop_back();
          }
        }

        while (!rankedPopulation.empty()) {
          delete rankedPopulation.back()->first;
          delete rankedPopulation.back();
          rankedPopulation.pop_back();
        }
        return this->EOS;
      }


      //! Only for first iteration
      if (!currentIteration && population.size() == totalPopulation) {
        rankedPopulation.resize(totalPopulation);
        for (size_t i = 0; i < totalPopulation; i++) {
          rankedPopulation[i] = new std::pair<std::vector<TId> *, double>(population[i], 0);
          this->ff_send_out(rankedPopulation[i]);
        }
        population.clear();
        currentIteration++;
      }

      //! All other iterations
      if (currentIteration && population.size() == totalNewPopulation) {
        std::sort(rankedPopulation.begin(), rankedPopulation.end(),
                  [&](const std::pair<std::vector<TId> *, double> *chromosomeA,
                      const std::pair<std::vector<TId> *, double> *chromosomeB) {
                    //! I order backwards so I can take advantage of the pop_back later
                    return chromosomeA->second > chromosomeB->second;
                  });

        //std::cout<< "Best value: "<<rankedPopulation.at(totalPopulation - 1)->second<<std::endl;

        if (population.size() > rankedPopulation.size()) {
          std::random_shuffle(population.begin(), population.end());
          while (population.size() - rankedPopulation.size()) {
            delete population.back();
            population.pop_back();
          }
        } else {
          while (totalPopulation - population.size()) {
            population.push_back(rankedPopulation.back()->first);
            delete rankedPopulation.back();
            rankedPopulation.pop_back();
          }
        }

        while (!rankedPopulation.empty()) {
          delete rankedPopulation.back()->first;
          delete rankedPopulation.back();
          rankedPopulation.pop_back();
        }
        rankedPopulation.resize(totalPopulation);
        for (size_t i = 0; i < totalPopulation; i++) {
          rankedPopulation[i] = new std::pair<std::vector<TId> *, double>(population[i], 0);
          this->ff_send_out(rankedPopulation[i]);
        }
        population.clear();
        currentIteration++;
      }

      return this->GO_ON;
    }
    std::vector<std::vector<TId> *> population;
    std::vector<std::pair<std::vector<TId> *, double> *> rankedPopulation;
    size_t totalPopulation;
    size_t totalIterations;
    size_t currentIteration;
    size_t totalNewPopulation;
  };


  struct Evaluate : ff::ff_node_t<std::pair<std::vector<TId> *, double>> {
    explicit Evaluate(const graph<TId, TValue> *graphPtr, std::mt19937 &generator): graph_(graphPtr), gen(generator) {}
    std::pair<std::vector<TId> *, double> *svc(std::pair<std::vector<TId> *, double> *chromosome) override {

      double eval = 0;
      for (size_t i = 0; i < chromosome->first->size() - 1; i++) {
        eval += graph_->getValueEdge(chromosome->first->at(i), chromosome->first->at(i + 1));
      }
      eval += graph_->getValueEdge(chromosome->first->at(chromosome->first->size() - 1), chromosome->first->at(0));
      chromosome->second = eval/ chromosome->first->size();
      this->ff_send_out(chromosome);

      return this->GO_ON;
    }
    const graph<TId, TValue> *graph_;
    std::mt19937 &gen;
  };

  struct selectionReproductionNode : ff::ff_minode_t<std::pair<std::vector<TId> *, double>> {
    explicit selectionReproductionNode(const size_t totalPopulation_, const int multiplier_, std::mt19937 &generator) :
        totalPopulation(totalPopulation_), multiplier(multiplier_), gen(generator) {}
    std::pair<std::vector<TId> *, double> *svc(std::pair<std::vector<TId> *, double> *in) override {
      evaluatedPopulation.push_back(in);
      if (evaluatedPopulation.size() == totalPopulation) {
        intermediatePopulation = new std::vector<std::vector<TId> *>();
        intermediatePopulation->reserve(multiplier * totalPopulation);
        std::vector<double> randomNumber(multiplier *evaluatedPopulation.size());
        //! Compute the avg value
        double avg = 0;
        for (auto chromosomeValue: evaluatedPopulation) {
          avg += chromosomeValue->second;
        }
        avg /= evaluatedPopulation.size();

        //! Compute the roulette wheel value
        //! Idea of the algorithm taken from "A genetic algorithm tutorial" of Darrell Whitley
        std::vector<double> rouletteWheelValue(evaluatedPopulation.size());
        rouletteWheelValue[0] = (1 / evaluatedPopulation[0]->second) / (1 / avg);

        for (size_t i = 1; i < evaluatedPopulation.size(); i++) {
          rouletteWheelValue[i] =
              rouletteWheelValue[i - 1] + (1 / evaluatedPopulation[i]->second) / (1 / avg);
        }
        //! Generation of the intermediatePopulation
        std::uniform_real_distribution<double>
            selectionDistribution{0, rouletteWheelValue[rouletteWheelValue.size() - 1]};
        randomNumber[0] = selectionDistribution(gen);
        for (size_t i = 1; i < multiplier * rouletteWheelValue.size(); i++) {
          randomNumber[i] = fmod(selectionDistribution(gen) + randomNumber[i - 1],
                                 rouletteWheelValue[evaluatedPopulation.size() - 1]);
        }

        std::sort(&randomNumber[0], &randomNumber[multiplier * evaluatedPopulation.size()]);

        size_t chromosomeEvalsIndex = 0;
        //! Using Stochastic universal sampling (also known as "roulette wheel selection")
        for (size_t i = 0; i < multiplier * evaluatedPopulation.size();) {
          if (rouletteWheelValue[chromosomeEvalsIndex] >= randomNumber[i]) {
            intermediatePopulation->push_back(evaluatedPopulation[chromosomeEvalsIndex]->first);
            //std::cout << chromosomeEvals[chromosomeEvalsIndex].first << std::endl;
            i++;
          } else {
            chromosomeEvalsIndex++;
          }
        }
        this->ff_send_out(intermediatePopulation);
        evaluatedPopulation.clear();
      }

      return this->GO_ON;
    }
    std::vector<std::pair<std::vector<TId> *, double> *> evaluatedPopulation;
    std::vector<std::vector<TId> *> *intermediatePopulation;

    size_t totalPopulation;
    int multiplier;
    std::mt19937 &gen;

  };

  struct EmitterCrossoverMutation : ff::ff_monode_t<std::vector<std::vector<TId> *>> {

    explicit EmitterCrossoverMutation(const double crossoverProbability_,
                                      const double mutationProbability_, std::mt19937 &generator) :
        crossoverProbability(crossoverProbability_), mutationProbability(mutationProbability_), gen(generator) {}
    std::vector<std::vector<TId> *> *svc(std::vector<std::vector<TId> *> *intermediatePopulation) override {
      //!Crossover emitter
      size_t totalCrossover = intermediatePopulation->size() * crossoverProbability;
      std::random_shuffle(intermediatePopulation->begin(), intermediatePopulation->end());

      //! Selection of the chromosome used for the crossover (intermediatePopulation.size() - 2 because i exchange the value between i and i+1,
      //! so if i use only intermediatePopulation.size() - 1 the index will go out of bound and create some error)
      std::uniform_int_distribution<size_t> SelectionCrossoverDistribution{0, intermediatePopulation->size() - 2};

      std::vector<size_t> crossoverPosition(totalCrossover);

      for (size_t i = 0; i < totalCrossover; i++) {
        crossoverPosition[i] = SelectionCrossoverDistribution(gen);
      }

      for (auto index: crossoverPosition) {
        this->ff_send_out(new std::pair<std::vector<TId> *, std::vector<TId> *>(intermediatePopulation->at(index),
                                                                                intermediatePopulation->at(index + 1)));
      }

      //!Mutation emitter
      size_t totalMutation = intermediatePopulation->size() * mutationProbability;
      std::random_shuffle(intermediatePopulation->begin(), intermediatePopulation->end());

      std::uniform_int_distribution<size_t> SelectionMutationDistribution{0, intermediatePopulation->size() - 1};

      std::vector<size_t> mutationPosition(totalMutation);

      for (size_t i = 0; i < totalMutation; i++) {
        mutationPosition[i] = SelectionMutationDistribution(gen);
      }

      for (auto index: mutationPosition) {
        this->ff_send_out(new std::pair<std::vector<TId> *, std::vector<TId> *>(intermediatePopulation->at(index),
                                                                                nullptr));
      }
      delete intermediatePopulation;
      return this->GO_ON;
    }

    double crossoverProbability;
    double mutationProbability;
    std::mt19937 &gen;
  };

  struct CrossoverMutation : ff::ff_node_t<std::pair<std::vector<TId> *, std::vector<TId> *>, std::vector<TId>> {
    explicit CrossoverMutation(std::mt19937 &generator) : gen(generator) {}
    std::vector<TId> *svc(std::pair<std::vector<TId> *,
                                    std::vector<TId> *> *crossoverPair) override {
      if (crossoverPair->second != nullptr) {
        //!Crossover

        std::uniform_int_distribution<size_t> crossoverDistribution{0, crossoverPair->first->size()};

        size_t elemFirstChromosomeIndexA = crossoverDistribution(gen);
        size_t elemFirstChromosomeIndexB = crossoverDistribution(gen);

        if (elemFirstChromosomeIndexA > elemFirstChromosomeIndexB) {
          std::swap(elemFirstChromosomeIndexA, elemFirstChromosomeIndexB);
        }

        //! Adding to the crossover chromosome the first elemFirstChromosome genes from the chromosome A
        auto *chromosomePop = new std::vector<TId>
            (crossoverPair->first->begin() + elemFirstChromosomeIndexA,
             crossoverPair->first->begin() + elemFirstChromosomeIndexB);

        auto interPopPtrBeginA = crossoverPair->first->begin();
        auto interPopPtrEndA = crossoverPair->first->begin() + elemFirstChromosomeIndexA;

        auto interPopPtrBeginB = crossoverPair->first->begin() + elemFirstChromosomeIndexB;
        auto interPopPtrEndB =
            crossoverPair->first->end();

        //! Adding to the crossover chromosome the remaining (chromosome.size() - elemFirstChromosome + 1) genes with
        //! the order that they appear in chromosome B
        for (size_t j = 0; j < crossoverPair->first->size() /*&& currentMissingId.size()*/;
             j++) {
          auto itA = std::find(interPopPtrBeginA,
                               interPopPtrEndA,
                               crossoverPair->second->at(j));
          auto itB = std::find(interPopPtrBeginB,
                               interPopPtrEndB,
                               crossoverPair->second->at(j));

          if (itA != interPopPtrEndA || itB != interPopPtrEndB) {
            if (itA != interPopPtrEndA) {
              chromosomePop->push_back(*itA);
            } else {
              chromosomePop->push_back(*itB);
            }
          }
        }

        //TODO: Verificare correttezza delete
        delete crossoverPair;
        this->ff_send_out(chromosomePop);
      } else {
        //!Mutation
        std::vector<TId> *chromosome = crossoverPair->first;
        std::uniform_int_distribution<size_t> mutationDistribution{0, chromosome->size() - 1};

        size_t positionGenesA = mutationDistribution(gen);
        size_t positionGenesB = mutationDistribution(gen);

        std::vector<TId> *chromosomePop = new std::vector<TId>(*chromosome);
        std::swap(chromosomePop->at(positionGenesA), chromosomePop->at(positionGenesB));
        chromosomePop->shrink_to_fit();
        delete crossoverPair;
        this->ff_send_out(chromosomePop);
      }

      return this->GO_ON;
    }

    std::mt19937 &gen;
  };

  struct EmitterMutation : ff::ff_monode_t<std::vector<std::vector<TId> *>> {

    explicit EmitterMutation(const double mutationProbability_, std::mt19937 &generator) :
        mutationProbability(mutationProbability_), gen(generator) {}
    std::vector<std::vector<TId> *> *svc(std::vector<std::vector<TId> *> *intermediatePopulation) override {
      size_t totalMutation = intermediatePopulation->size() * mutationProbability;
      std::random_shuffle(intermediatePopulation->begin(), intermediatePopulation->end());

      std::uniform_int_distribution<size_t> SelectionMutationDistribution{0, intermediatePopulation->size() - 1};

      std::vector<size_t> mutationPosition(totalMutation);

      for (size_t i = 0; i < totalMutation; i++) {
        mutationPosition[i] = SelectionMutationDistribution(gen);
      }

      for (auto index: mutationPosition) {
        this->ff_send_out(intermediatePopulation->at(index));
      }
      return this->GO_ON;
    }

    double mutationProbability;
    std::mt19937 &gen;
  };

  struct Mutation : ff::ff_node_t<std::vector<TId>> {
    explicit Mutation(std::mt19937 &generator) : gen(generator) {}
    std::vector<TId> *svc(std::vector<TId> *chromosome) override {
      std::uniform_int_distribution<size_t> mutationDistribution{0, chromosome->size() - 1};

      size_t positionGenesA = mutationDistribution(gen);
      size_t positionGenesB = mutationDistribution(gen);

      std::vector<TId> *chromosomePop = new std::vector<TId>(*chromosome);
      std::swap(chromosomePop[positionGenesA], chromosomePop[positionGenesB]);
      chromosomePop->shrink_to_fit();
      this->ff_send_out(chromosomePop);

      return this->GO_ON;
    }

    std::mt19937 &gen;
  };

  graph<TId, TValue> *graph_;
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
 public:
  void initializer() override {};
  void evaluate() override {};
  void selectionReproduction() override {};
  void crossover() override {};
  void mutation() override {};
  void run(int iteration) override;
  explicit TSPGeneticAlgorithmFF();
  explicit TSPGeneticAlgorithmFF(int seed_);
  explicit TSPGeneticAlgorithmFF(int seed_, double crossoverProbability_, double mutationProbability_);
  void setGraph(graph<TId, TValue> *graph);
  void setRandomGraph(size_t nNodes);
  void SetCrossoverProbability(double crossover_probability);
  void SetMutationProbability(double mutation_probability);
  void SetMultiplier(int multiplier);
  void SetTotalPopulation(size_t total_population);
  void SetNWorker(int n_worker);

};

template
class TSPGeneticAlgorithmFF<int, int>;
template
class TSPGeneticAlgorithmFF<int, double>;
template
class TSPGeneticAlgorithmFF<unsigned long long, double>;

#endif //TSPGENETICALGORITHM_SRC_GENETICALGORITHM_TSPGENETICALGORITHMFF_H_
