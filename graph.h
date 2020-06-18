//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__GRAPH_H_
#define TSPGENETICALGORITHM__GRAPH_H_
#include <unordered_map>
#include <vector>
template<class T, class TId = unsigned long long, class TValue = unsigned long long>
class graph {
 protected:
  T edge;
  std::vector<TId> nodes;
 public:
  virtual ~graph() = default;
  explicit graph() = default;
  virtual void randomInit() = 0;
  virtual TValue getValueEdge(TId keyFrom, TId keyTo) = 0;
  explicit graph(unsigned long long nNode) : nodes(nNode) {};
  inline explicit graph(std::vector<TId> nodes_){
    nodes=nodes_;
  };
  inline explicit graph(std::vector<TId> nodes_, T edge_) {
    nodes = nodes_;
    edge=edge_;
  };
};
#endif //TSPGENETICALGORITHM__GRAPH_H_
