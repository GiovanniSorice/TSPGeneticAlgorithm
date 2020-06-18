//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#define TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
#include "graph.h"
#include <unordered_map>
#include <vector>
#include "geneticAlgorithm.h"
template<class T, class TId, class TValue>
class TSPGeneticAlgorithm : public geneticAlgorithm<T> {
 private:
  graph<T, TId, TValue>* graph_;
 public:
  virtual void initializer();
  virtual void evaluate();
  virtual void selectionReproduction();
  virtual void crossover();
  virtual void mutation();
  virtual void run();

};

#endif //TSPGENETICALGORITHM__TSPGENETICALGORITHM_H_
