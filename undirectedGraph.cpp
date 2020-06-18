//
// Created by gs1010 on 31/05/20.
//

#include "undirectedGraph.h"
template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(unsigned long long nNode):
    graph<std::vector<std::unordered_map<TId, TValue>>, TId, TValue>(nNode) {
  randomInit();
}
template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_):
    graph<std::vector<std::unordered_map<TId, TValue>>, TId, TValue>(nodes_) {}

template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_,
                                              std::vector<std::unordered_map<TId, TValue>> edge_):
    graph<std::vector<std::unordered_map<TId, TValue>>, TId, TValue>(nodes_, edge_) {}

template<class TId, class TValue>
void undirectedGraph<TId, TValue>::randomInit() {
  std::vector<std::unordered_map<TId, TValue>> tmp(this->nodes.size());
  for(int i=0; i< this->nodes.size(); i++){
    for(int j=i+1;j< this->nodes.size(); j++){
      //tmp[i];
    }
  }
}
template<class TId, class TValue>
TValue undirectedGraph<TId, TValue>::getValueEdge() {
  return nullptr;
}
