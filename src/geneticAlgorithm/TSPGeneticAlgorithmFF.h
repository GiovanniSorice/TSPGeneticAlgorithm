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
    void svc_end() { std::cout << "fine EmitterInitializer\n"; }

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
      return this->EOS;
    }

    const std::vector<TId> &nodes;
    std::mt19937 &gen;
  };
  struct collectorInitializerEmitterEvaluate : ff::ff_minode_t<std::vector<TId>> {
    explicit collectorInitializerEmitterEvaluate(const size_t totalPopulation_) :
        totalPopulation(totalPopulation_) {}
    std::vector<TId> *svc(std::vector<TId> *in) override {
      population.push_back(in);
      if (population.size() == totalPopulation) {
        rankedPopulation.resize(totalPopulation);
        for (size_t i = 0; i < totalPopulation; i++) {
          rankedPopulation[i] = new std::pair<std::vector<TId> *, double>(population[i], 0);
          this->ff_send_out(rankedPopulation[i]);
        }
      }
      return this->GO_ON;
    }
    void svc_end() { std::cout << "fine receiver\n"; }
    std::vector<std::vector<TId> *> population;
    std::vector<std::pair<std::vector<TId> *, double> *> rankedPopulation;
    size_t totalPopulation;
  };


  struct Evaluate : ff::ff_node_t<std::pair<std::vector<TId> *, double>> {
    explicit Evaluate(const graph<TId, TValue> *graphPtr, std::mt19937 &generator): graph_(graphPtr), gen(generator) {}
    std::pair<std::vector<TId> *, double> *svc(std::pair<std::vector<TId> *, double> *chromosome) override {

      double eval = 0;
        for (size_t i = 0; i < chromosome->first->size() - 1; i++) {
          eval += graph_->getValueEdge(chromosome->first->at(i), chromosome->first->at(i+1));
        }
        eval += graph_->getValueEdge(chromosome->first->at(chromosome->first->size() - 1), chromosome->first->at(0));
      chromosome->second=eval;
      this->ff_send_out(chromosome);

      return this->EOS;
    }

    const graph<TId,TValue>* graph_;
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
