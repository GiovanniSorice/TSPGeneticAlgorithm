//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#define TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
#include <unordered_map>
#include <vector>
#include "graph.h"
template<typename TId = unsigned long long, typename TValue = double>
using edgeType = std::unordered_map<TId,std::unordered_map<TId, TValue>>;

template<typename TId = unsigned long long, typename TValue = unsigned long long>
class undirectedGraph : public graph<TId, TValue> {
 protected:
  edgeType<TId, TValue> edges;

 public:
  ~undirectedGraph() override = default;
  explicit undirectedGraph() = default;
  explicit undirectedGraph(unsigned long long nNode);
  explicit undirectedGraph(std::vector<TId> nodes_);
  explicit undirectedGraph(std::vector<TId> nodes_, edgeType<TId, TValue> edges_);
  void randomInit(int seed = 0) override;
  TValue getValueEdge(TId keyFrom, TId keyTo)const override;
};

template
class undirectedGraph<int, int>;
template
class undirectedGraph<int, double>;
template
class undirectedGraph<unsigned long long, double>;

#endif //TSPGENETICALGORITHM__UNDIRECTEDGRAPH_H_
