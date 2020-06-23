#include <iostream>
#include "undirectedGraph.h"
#include "TSPGeneticAlgorithm.h"

int main() {
  std::cout << "Hello, World!" << std::endl;
  undirectedGraph<int, double> g(5);
  std::cout <<g.getValueEdge(4,2)<<std::endl;
  std::cout << g.getValueEdge(2, 4) << std::endl;

  std::cout << "Hello, World!" << std::endl;

  std::unordered_map<std::string, std::unordered_map<std::string, int>> x(5);
  x["s"];

  TSPGeneticAlgorithm<int, double> tsp;
  tsp.initializer(50);
  tsp.evaluate();
  tsp.selectionReproduction();
}