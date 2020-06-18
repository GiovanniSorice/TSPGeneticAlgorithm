//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#define TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#include <unordered_map>
#include <vector>
#include "graph.h"
template<class TId = unsigned long long, class TValue = unsigned long long>
using edgeType = std::unordered_map<TId,std::unordered_map<TId, TValue>>;

template<class TId = unsigned long long, class TValue = unsigned long long>
class undirectedGraph : public graph<edgeType<TId,TValue>, TId, TValue> {

 public:
  ~undirectedGraph() override = default;
  explicit undirectedGraph() = default;
  explicit undirectedGraph(unsigned long long nNode);
  explicit undirectedGraph(std::vector<TId> nodes_);
  explicit undirectedGraph(std::vector<TId> nodes_, edgeType<TId, TValue> edge_);
  void randomInit() override;
  TValue getValueEdge(TId keyFrom, TId keyTo) override;

};

template
class undirectedGraph<int, int>;
template
class undirectedGraph<unsigned long long, unsigned long long>;

#endif //TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
