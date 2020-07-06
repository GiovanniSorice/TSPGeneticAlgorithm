#include <iostream>
#include <chrono>
#include "src/graph/undirectedGraph.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithmST.h"

int main(int argc, char *argv[]) {

  if (argc < 4) {
    std::cout << "Usage is " << argv[0]
              << " numberOfNode population nIteration [nWorker] [seed] [crossoverProbability] [mutationProbability]" << std::endl;
    return (-1);
  }

  const size_t numberOfNode = atoi(argv[1]);
  const size_t nPopulation = atoi(argv[2]);

  if (nPopulation < 2) {
    std::cout << "Population must be greater than or equal to 2"
              << std::endl;
    return (-1);
  }

  const int nIteration = atoi(argv[3]);

  int nWorker = 1;
  int seed = 0;
  double crossoverProbability = 0.2;
  double mutationProbability = 0.1;

  if (argc > 4) {
    nWorker = atoi(argv[4]);
  }

  if (argc > 5) {
    seed = atoi(argv[5]);
  }

  if (argc > 6) {
    crossoverProbability = std::stod(argv[6]);
  }
  if (argc > 7) {
    mutationProbability = std::stod(argv[7]);
  }

  /*
  undirectedGraph<int, double> g(200);
  g.randomInit(5);
  */

  TSPGeneticAlgorithmST<int, double> tspST(seed, crossoverProbability, mutationProbability);
  tspST.SetMultiplier(1);
  tspST.SetTotalPopulation(nPopulation);
  tspST.setRandomGraph(numberOfNode);
  tspST.SetNWorker(nWorker);
  //tsp.setGraph(&g);
  auto startST = std::chrono::high_resolution_clock::now();
  tspST.run(nIteration);
  auto elapsedST = std::chrono::high_resolution_clock::now() - startST;
  auto msecST = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedST).count();
  printf("%d %ld %ld %ld %d\n", nWorker, msecST, numberOfNode, nPopulation, nIteration);

}