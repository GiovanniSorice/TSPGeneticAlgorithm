//
// Created by gs1010 on 31/05/20.
//

#ifndef TSPGENETICALGORITHM__GENETICALGORITHM_H_
#define TSPGENETICALGORITHM__GENETICALGORITHM_H_

class geneticAlgorithm {
 public:
  //TODO: Lasciare solo run come public, quindi togliere il resto delle funzioni
  virtual void initializer() = 0;
  virtual void evaluate() = 0;
  virtual void selectionReproduction() = 0;
  virtual void crossover() = 0;
  virtual void mutation() = 0;
  virtual void run(int iteration) = 0;

};

#endif //TSPGENETICALGORITHM__GENETICALGORITHM_H_
