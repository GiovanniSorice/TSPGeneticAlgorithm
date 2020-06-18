//
// Created by gs1010 on 31/05/20.
//

#include <algorithm>
#include "undirectedGraph.h"
template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(unsigned long long nNode):
    graph<edgeType<TId, TValue>, TId, TValue>(nNode) {
  randomInit();
}
template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_):
    graph<edgeType<TId, TValue>, TId, TValue>(nodes_) {}

template<class TId, class TValue>
undirectedGraph<TId, TValue>::undirectedGraph(std::vector<TId> nodes_,
                                              edgeType<TId, TValue> edge_):
    graph<edgeType<TId, TValue>, TId, TValue>(nodes_, edge_) {}

template<class TId, class TValue>
void undirectedGraph<TId, TValue>::randomInit() {
  edgeType<TId, TValue> tmp(this->nodes.size());
  for (int i = 0; i < this->nodes.size(); i++) {
    for (int j = i + 1; j < this->nodes.size(); j++) {
      //tmp[i];
    }
  }
}
template<class TId, class TValue>
TValue undirectedGraph<TId, TValue>::getValueEdge(TId keyFrom, TId keyTo) {
  //TODO: Aggiungere execution policy

  auto nodeEdges = this->edge.find(keyFrom);
  if (nodeEdges == this->edge.end()) {
    nodeEdges = this->edge.find(keyTo);
    return nodeEdges->second.find(keyFrom)->second;
  }else{
    return nodeEdges->second.find(keyTo)->second;
  }
}
