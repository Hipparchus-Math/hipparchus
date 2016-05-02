# Genetic Algorithms

## Overview

The genetics package provides a framework and implementations for
genetic algorithms.


## GA Framework

[GeneticAlgorithm](../apidocs/org/hipparchus/genetics/GeneticAlgorithm.html)
provides an execution framework for Genetic Algorithms (GA).
[Populations](../apidocs/org/hipparchus/genetics/Population.html),
consisting of [Chromosomes](../apidocs/org/hipparchus/genetics/Chromosome.html)
are evolved by the `GeneticAlgorithm` until a
[StoppingCondition](../apidocs/org/hipparchus/genetics/StoppingCondition.html)
is reached. Evolution is determined by
[SelectionPolicy](../apidocs/org/hipparchus/genetics/SelectionPolicy.html),
[MutationPolicy](../apidocs/org/hipparchus/genetics/MutationPolicy.html)
and [Fitness](../apidocs/org/hipparchus/genetics/Fitness.html).

The GA itself is implemented by the `evolve` method of the
`GeneticAlgorithm` class, which looks like this:

    public Population evolve(Population initial, StoppingCondition condition) {
        Population current = initial;
        while (!condition.isSatisfied(current)) {
            current = nextGeneration(current);
        }
        return current;
    }

The `nextGeneration` method implements the following algorithm:

1. Get nextGeneration population to fill from `current` generation, using its nextGeneration method
  * Apply configured `SelectionPolicy` to select a pair of parents from `current`
  * With probability = [getCrossoverRate\(\)](../apidocs/org/hipparchus/genetics/GeneticAlgorithm.html#getCrossoverRate--), apply configured `CrossoverPolicy` to parents
  * With probability = [getMutationRate\(\)](../apidocs/org/hipparchus/genetics/GeneticAlgorithm.html#getMutationRate--), apply configured `MutationPolicy` to each of the offspring
  * Add offspring individually to nextGeneration, space permitting
1. Return nextGeneration


## Implementation

Here is an example GA execution:

    // initialize a new genetic algorithm
    GeneticAlgorithm ga = new GeneticAlgorithm(
        new OnePointCrossover&lt;Integer&gt;(),
        1,
        new RandomKeyMutation(),
        0.10,
        new TournamentSelection(TOURNAMENT_ARITY)
    );
            
    // initial population
    Population initial = getInitialPopulation();
            
    // stopping condition
    StoppingCondition stopCond = new FixedGenerationCount(NUM_GENERATIONS);
            
    // run the algorithm
    Population finalPopulation = ga.evolve(initial, stopCond);
            
    // best chromosome from the final population
    Chromosome bestFinal = finalPopulation.getFittestChromosome();
    
The arguments to the `GeneticAlgorithm` constructor above are:

| Parameter | value in example | meaning |
| --- | --- | --- |
| crossoverPolicy | [OnePointCrossover](../apidocs/org/hipparchus/genetics/OnePointCrossover.html) | A random crossover point is selected and the first part from each parent is copied to the corresponding child, and the second parts are copied crosswise. |
| crossoverRate | 1 | Always apply crossover |
| mutationPolicy | [RandomKeyMutation](../apidocs/org/hipparchus/genetics/RandomKeyMutation.html) | Changes a randomly chosen element of the array representation to a random value uniformly distributed in [0,1]. |
| mutationRate | .1 | Apply mutation with probability 0.1 - that is, 10% of the time. |
| selectionPolicy | [TournamentSelection](../apidocs/org/hipparchus/genetics/TournamentSelection.html) | Each of the two selected chromosomes is selected based on an n-ary tournament -- this is done by drawing n random chromosomes without replacement from the population, and then selecting the fittest chromosome among them. |

The algorithm starts with an `initial` population of `Chromosomes.` and executes until
the specified [StoppingCondition](../apidocs/org/hipparchus/genetics/StoppingCondition.html)
is reached.  In the example above, a
[FixedGenerationCount](../apidocs/org/hipparchus/genetics/FixedGenerationCount.html)
stopping condition is used, which means the algorithm proceeds through a fixed number of generations.
