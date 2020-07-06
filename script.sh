#!/bin/sh
clear
echo "Running first test."

for nodes in 500 1000 2000; do
  echo "nw time nodes chromosome iteration" >>resultsSeq${nodes}.txt
  echo nodes $nodes
  for chromosome in 500 5000 20000; do
    build/mainSeq ${nodes} ${chromosome} 10 1 123>>resultsSeq${nodes}.txt
  done
done
echo "Fine Seq."

for nodes in 2000; do
  echo "nw time nodes chromosome iteration" >>resultsST${nodes}${chromosome}.txt
  echo nodes $nodes
  for chromosome in 500 5000 20000; do
    for w in {1..250}; do
      echo $w
      build/mainST ${nodes} ${chromosome} 10 $w 123 >>resultsST${nodes}${chromosome}.txt
    done
  done
done
echo "Fine ST."

for nodes in 2000; do
  echo "nw time nodes chromosome iteration" >>resultsFF${nodes}${chromosome}.txt
  echo nodes $nodes
  for chromosome in 500 5000 20000; do
    for w in {1..250}; do
      echo $w
      build/mainFF ${nodes} ${chromosome} 10 $w 123 >>resultsFF${nodes}${chromosome}.txt
    done
  done
done
echo "Fine FF."
