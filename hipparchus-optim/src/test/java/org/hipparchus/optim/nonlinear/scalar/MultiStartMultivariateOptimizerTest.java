/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.optim.nonlinear.scalar;

import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.PointValuePair;
import org.hipparchus.optim.SimpleValueChecker;
import org.hipparchus.optim.nonlinear.scalar.gradient.CircleScalar;
import org.hipparchus.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer;
import org.hipparchus.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.hipparchus.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.hipparchus.random.GaussianRandomGenerator;
import org.hipparchus.random.JDKRandomGenerator;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.UncorrelatedRandomVectorGenerator;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.random.HaltonSequenceGenerator;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class MultiStartMultivariateOptimizerTest {
    @Test
    void testCircleFitting() {
        CircleScalar circle = new CircleScalar();
        circle.addPoint( 30.0,  68.0);
        circle.addPoint( 50.0,  -6.0);
        circle.addPoint(110.0, -20.0);
        circle.addPoint( 35.0,  15.0);
        circle.addPoint( 45.0,  97.0);
        // TODO: the wrapper around NonLinearConjugateGradientOptimizer is a temporary hack for
        // version 3.1 of the library. It should be removed when NonLinearConjugateGradientOptimizer
        // will officially be declared as implementing MultivariateDifferentiableOptimizer
        GradientMultivariateOptimizer underlying
            = new NonLinearConjugateGradientOptimizer(NonLinearConjugateGradientOptimizer.Formula.POLAK_RIBIERE,
                                                      new SimpleValueChecker(1e-10, 1e-10));
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(753289573253L);
        RandomVectorGenerator generator
            = new UncorrelatedRandomVectorGenerator(new double[] { 50, 50 },
                                                    new double[] { 10, 10 },
                                                    new GaussianRandomGenerator(g));
        int nbStarts = 10;
        MultiStartMultivariateOptimizer optimizer
            = new MultiStartMultivariateOptimizer(underlying, nbStarts, generator);
        PointValuePair optimum
            = optimizer.optimize(new MaxEval(1000),
                                 circle.getObjectiveFunction(),
                                 circle.getObjectiveFunctionGradient(),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 98.680, 47.345 }));
        assertEquals(1000, optimizer.getMaxEvaluations());
        PointValuePair[] optima = optimizer.getOptima();
        assertEquals(nbStarts, optima.length);
        for (PointValuePair o : optima) {
            // we check the results of all intermediate restarts here (there are 10 such results)
            Vector2D center = new Vector2D(o.getPointRef()[0], o.getPointRef()[1]);
            assertTrue(69.9592 < circle.getRadius(center));
            assertTrue(69.9602 > circle.getRadius(center));
            assertTrue(96.0745 < center.getX());
            assertTrue(96.0762 > center.getX());
            assertTrue(48.1344 < center.getY());
            assertTrue(48.1354 > center.getY());
        }

        assertTrue(optimizer.getEvaluations() > 850);
        assertTrue(optimizer.getEvaluations() < 900);

        assertEquals(3.1267527, optimum.getValue(), 1e-8);
    }

    @Test
    void testRosenbrock() {
        Rosenbrock rosenbrock = new Rosenbrock();
        SimplexOptimizer underlying
            = new SimplexOptimizer(new SimpleValueChecker(-1, 1e-3));
        NelderMeadSimplex simplex = new NelderMeadSimplex(new double[][] {
                { -1.2,  1.0 },
                { 0.9, 1.2 } ,
                {  3.5, -2.3 }
            });
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(16069223052L);
        RandomVectorGenerator generator
            = new UncorrelatedRandomVectorGenerator(2, new GaussianRandomGenerator(g));
        int nbStarts = 10;
        MultiStartMultivariateOptimizer optimizer
            = new MultiStartMultivariateOptimizer(underlying, nbStarts, generator);
        PointValuePair optimum
            = optimizer.optimize(new MaxEval(1100),
                                 new ObjectiveFunction(rosenbrock),
                                 GoalType.MINIMIZE,
                                 simplex,
                                 new InitialGuess(new double[] { -1.2, 1.0 }));
        assertEquals(nbStarts, optimizer.getOptima().length);

        assertEquals(rosenbrock.getCount(), optimizer.getEvaluations());
        assertTrue(optimizer.getEvaluations() > 900);
        assertTrue(optimizer.getEvaluations() < 1200);
        assertTrue(optimum.getValue() < 5e-5);
    }

    @Test
    void testMissingMaxEval() {
        SimplexOptimizer underlying = new SimplexOptimizer(new SimpleValueChecker(1e-10, 1e-10));
        RandomVectorGenerator generator = new UncorrelatedRandomVectorGenerator(2, new GaussianRandomGenerator(new JDKRandomGenerator()));
        MultiStartMultivariateOptimizer optimizer = new MultiStartMultivariateOptimizer(underlying, 10, generator);
        assertThrows(MathIllegalStateException.class, () ->
            optimizer.optimize(new ObjectiveFunction(x -> 0), GoalType.MINIMIZE, new InitialGuess(new double[] { 0, 0 }))
        );
    }

    @Test
    void testMissingInitialGuess() {
        SimplexOptimizer underlying = new SimplexOptimizer(new SimpleValueChecker(1e-10, 1e-10));
        RandomVectorGenerator generator = new UncorrelatedRandomVectorGenerator(2, new GaussianRandomGenerator(new JDKRandomGenerator()));
        MultiStartMultivariateOptimizer optimizer = new MultiStartMultivariateOptimizer(underlying, 10, generator);
        assertThrows(MathIllegalStateException.class, () ->
            optimizer.optimize(new MaxEval(100), new ObjectiveFunction(x -> 0), GoalType.MINIMIZE)
        );
    }

    @Test
    void testAllRunsFail() {
        JDKRandomGenerator g = new JDKRandomGenerator();
        g.setSeed(16069223052L);
        RandomVectorGenerator generator = new UncorrelatedRandomVectorGenerator(2, new GaussianRandomGenerator(g));
        
        // Optimizer that always fails
        MultivariateOptimizer underlying = new MultivariateOptimizer(new SimpleValueChecker(1e-10, 1e-10)) {
            @Override
            protected PointValuePair doOptimize() {
                throw new RuntimeException("failed");
            }
        };
        
        MultiStartMultivariateOptimizer optimizer = new MultiStartMultivariateOptimizer(underlying, 10, generator);
        RuntimeException ex = assertThrows(RuntimeException.class, () ->
            optimizer.optimize(new MaxEval(100), new ObjectiveFunction(x -> 0), GoalType.MINIMIZE, new InitialGuess(new double[] { 0, 0 }))
        );
        assertEquals("failed", ex.getMessage());
    }

    @Test
    void testVectorGeneratorRetry() {
        // Use an optimizer that supports bounds, like BOBYQA
        MultivariateOptimizer underlying = new org.hipparchus.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer(5);
        // Generator that returns values outside bounds
        RandomVectorGenerator generator = new RandomVectorGenerator() {
            private int count = 0;
            @Override
            public double[] nextVector() {
                return (count++ == 0) ? new double[] { 10, 10 } : new double[] { 0.5, 0.5 };
            }
        };
        MultiStartMultivariateOptimizer optimizer = new MultiStartMultivariateOptimizer(underlying, 2, generator);
        optimizer.optimize(new MaxEval(100), new ObjectiveFunction(x -> (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5)), GoalType.MINIMIZE,
                           new InitialGuess(new double[] { 0.5, 0.5 }),
                           new SimpleBounds(new double[] { 0, 0 }, new double[] { 1, 1 }));
        assertEquals(2, optimizer.getOptima().length);
    }

    @Test
    void testVectorGeneratorExhausted() {
        MultivariateOptimizer underlying = new org.hipparchus.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer(5);
        // Generator that always returns null
        RandomVectorGenerator generator = () -> null;
        MultiStartMultivariateOptimizer optimizer = new MultiStartMultivariateOptimizer(underlying, 2, generator);
        assertThrows(MathIllegalStateException.class, () ->
            optimizer.optimize(new MaxEval(100), new ObjectiveFunction(x -> 0), GoalType.MINIMIZE,
                               new InitialGuess(new double[] { 0.5, 0.5 }))
        );
    }

    private static class Rosenbrock implements MultivariateFunction {
        private int count;

        public Rosenbrock() {
            count = 0;
        }

        public double value(double[] x) {
            ++count;
            double a = x[1] - x[0] * x[0];
            double b = 1 - x[0];
            return 100 * a * a + b * b;
        }

        public int getCount() {
            return count;
        }
    }
}
