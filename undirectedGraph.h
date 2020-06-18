//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#define TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#include <unordered_map>
#include <vector>
#include "graph.h"
template<class TId = unsigned long long, class TValue = unsigned long long>
class undirectedGraph : public graph<std::vector<std::unordered_map<TId, TValue>>, TId, TValue> {
 public:
  ~undirectedGraph() override = default;
  explicit undirectedGraph() = default;
  explicit undirectedGraph(unsigned long long nNode);
  explicit undirectedGraph(std::vector<TId> nodes_);
  explicit undirectedGraph(std::vector<TId> nodes_, std::vector<std::unordered_map<TId, TValue>> edge_);
  void randomInit() override;
  TValue getValueEdge() override;

};

#endif //TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
