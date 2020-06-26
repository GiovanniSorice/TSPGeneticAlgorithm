#include <iostream>
#include <chrono>
#include "src/graph/undirectedGraph.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithm.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithmOMP.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithmST.h"

int main(int argc, char *argv[]) {

  if (argc < 4) {
    std::cout << "Usage is " << argv[0]
              << " numberOfNode population nIteration [nWorker] [seed] [crossoverProbability] [mutationProbability]" << std::endl;
    return (-1);
  }

  const size_t numberOfNode = atoi(argv[1]);
  const size_t population = atoi(argv[2]);
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

/*
  TSPGeneticAlgorithm<int, double> tspSeq(seed, crossoverProbability, mutationProbability);
  tspSeq.SetMultiplier(1);
  tspSeq.SetTotalPopulation(population);
  //auto startGraph = std::chrono::high_resolution_clock::now();
  tspSeq.setRandomGraph(numberOfNode);
  //auto elapsedGraph = std::chrono::high_resolution_clock::now() - startGraph;
  //auto msecGraph = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedGraph).count();
  //printf("setRandomGraph time (msecs): %ld\n", msecGraph);

  auto startSeq = std::chrono::high_resolution_clock::now();
  tspSeq.run(nIteration);
  auto elapsedSeq = std::chrono::high_resolution_clock::now() - startSeq;
  auto msecSeq = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedSeq).count();
  printf("Final Seq time (msecs): %ld\n", msecSeq);
*/


  std::cout<<std::endl;

  TSPGeneticAlgorithmST<int, double> tspST(seed, crossoverProbability, mutationProbability);
  tspST.SetMultiplier(1);
  tspST.SetTotalPopulation(population);
  tspST.setRandomGraph(numberOfNode);
  tspST.SetNWorker(nWorker);
  //tsp.setGraph(&g);
  auto startST = std::chrono::high_resolution_clock::now();
  tspST.run(nIteration);
  auto elapsedST = std::chrono::high_resolution_clock::now() - startST;
  auto msecST = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedST).count();
  printf("Final ST time (msecs): %ld\n", msecST);


  /*
  TSPGeneticAlgorithmOMP<int, double> tspOMP(seed, crossoverProbability, mutationProbability);
  tspOMP.SetMultiplier(1);
  tspOMP.SetTotalPopulation(population);
  tspOMP.setRandomGraph(numberOfNode);
  tspOMP.SetNWorker(nWorker);
  //tsp.setGraph(&g);
  auto startOMP = std::chrono::high_resolution_clock::now();
  tspOMP.run(nIteration);
  auto elapsedOMP = std::chrono::high_resolution_clock::now() - startOMP;
  auto msecOMP = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedOMP).count();
  printf("Final OMP time (msecs): %ld\n", msecOMP);
  */
}