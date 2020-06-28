//
// Created by gs1010 on 31/05/20.
//

#include <algorithm>
#include <random>
#include <iostream>
#include "undirectedGraph.h"

template<typename TId, typename TValue>
undirectedGraph<TId, TValue>::undirectedGraph(unsigned long long nNode):
    graph< TId, TValue>(nNode) {
  for (size_t i = 0; i < nNode; i++) {
    this->nodes[i] = i;
  }
  randomInit();
}
template<typename TId, typename TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_):
    graph<TId, TValue>(nodes_) {}

template<typename TId, typename TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_,
                                              edgeType<TId, TValue> edges_):
    graph<TId, TValue>(nodes_) {
  edges = edges_;
}
/** Random initialization of the edges
 */
template<typename TId, typename TValue>
void undirectedGraph<TId, TValue>::randomInit(int seed) {
  //edgeType<TId, TValue> tmpEdges(this->nodes.size());
    std::random_device rd;
    std::mt19937 gen(rd());
  if (seed) {
    gen.seed(seed);
  }
  std::uniform_real_distribution<double> unif(0, 100);
  this->edges.clear();
  for (size_t i = 0; i < this->nodes.size(); i++) {
    //! Initialization of the nodes[i] edges unordered_map
    this->edges.insert(std::make_pair(this->nodes[i], std::unordered_map<TId, TValue>()));
    for (size_t j = i + 1; j < this->nodes.size(); j++) {
      this->edges[this->nodes[i]].emplace(this->nodes[j], unif(gen));
    }
  }
/*
  for (auto &x: this->edges) {
    std::cout << x.first << ": " << std::endl;
    for (auto &y: x.second) {
      std::cout << "    " << y.first << ": " << y.second << std::endl;
    }
  }
*/
}
/** Get the value of the edge from node keyFrom to node keyTo
 * @param keyFrom: key of the node where the edge starts
 * @param keyTo: key of the node where the edge arrives
 */
template<typename TId, typename TValue>
TValue undirectedGraph<TId, TValue>::getValueEdge(TId keyFrom, TId keyTo) const {

  auto nodeEdges = this->edges.find(keyFrom);
  if (nodeEdges == this->edges.end()) {
    throw "node not in the list";
  }

  auto edge = nodeEdges->second.find(keyTo);
  if (edge == nodeEdges->second.end()) {
    return this->edges.find(keyTo)->second.find(keyFrom)->second;
  } else {
    return edge->second;
  }
}
