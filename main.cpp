#include <iostream>
#include <chrono>
#include "src/graph/undirectedGraph.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithm.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithmST.h"
#include "src/geneticAlgorithm/TSPGeneticAlgorithmFF.h"

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


  TSPGeneticAlgorithm<int, double> tspSeq(seed, crossoverProbability, mutationProbability);
  tspSeq.SetMultiplier(1);
  tspSeq.SetTotalPopulation(nPopulation);
  tspSeq.setRandomGraph(numberOfNode);
  auto startSeq = std::chrono::high_resolution_clock::now();
  tspSeq.run(nIteration);
  auto elapsedSeq = std::chrono::high_resolution_clock::now() - startSeq;
  auto msecSeq = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedSeq).count();
  printf("Final Seq time (msecs): %ld\n", msecSeq);



  std::cout<<std::endl;

  TSPGeneticAlgorithmST<int, double> tspST(seed, crossoverProbability, mutationProbability);
  tspST.SetMultiplier(1);
  tspST.SetTotalPopulation(nPopulation);
  tspST.setRandomGraph(numberOfNode);
  tspST.SetNWorker(nWorker);
  auto startST = std::chrono::high_resolution_clock::now();
  tspST.run(nIteration);
  auto elapsedST = std::chrono::high_resolution_clock::now() - startST;
  auto msecST = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedST).count();
  printf("Final ST time (msecs): %ld\n", msecST);

  std::cout << std::endl;

  TSPGeneticAlgorithmFF<int, double> tspFF(seed, crossoverProbability, mutationProbability);
  tspFF.SetMultiplier(1);
  tspFF.SetTotalPopulation(nPopulation);
  tspFF.setRandomGraph(numberOfNode);
  tspFF.SetNWorker(nWorker);
  auto startFF = std::chrono::high_resolution_clock::now();
  tspFF.run(nIteration);
  auto elapsedFF = std::chrono::high_resolution_clock::now() - startFF;
  auto msecFF = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedFF).count();
  printf("Final FF time (msecs): %ld\n", msecFF);

}