//
// Created by gs1010 on 31/05/20.
//

#include "TSPGeneticAlgorithm.h"
template<class T, class TId, class TValue>
void TSPGeneticAlgorithm<T, TId, TValue>::initializer() {
  graph_ = new graph<T, TId, TValue>(20);
  graph_->randomInit();

}
