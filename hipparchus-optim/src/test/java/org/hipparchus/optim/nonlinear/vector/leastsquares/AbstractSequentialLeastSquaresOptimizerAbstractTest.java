/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.not;
import static org.hamcrest.CoreMatchers.sameInstance;

import java.io.IOException;
import java.util.Arrays;

import org.hamcrest.MatcherAssert;
import org.hipparchus.analysis.MultivariateMatrixFunction;
import org.hipparchus.analysis.MultivariateVectorFunction;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.BlockRealMatrix;
import org.hipparchus.linear.DiagonalMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.SimpleVectorValueChecker;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresOptimizer.Optimum;
import org.hipparchus.optim.nonlinear.vector.leastsquares.LeastSquaresProblem.Evaluation;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
import org.junit.Assert;
import org.junit.Test;

/**
 * Some of the unit tests are re-implementations of the MINPACK <a
 * href="http://www.netlib.org/minpack/ex/file17">file17</a> and <a
 * href="http://www.netlib.org/minpack/ex/file22">file22</a> test files. The
 * redistribution policy for MINPACK is available <a href="http://www.netlib.org/minpack/disclaimer">here</a>.
 * <p/>
 * <T> Concrete implementation of an optimizer.
 *
 */
public abstract class AbstractSequentialLeastSquaresOptimizerAbstractTest {

    /** default absolute tolerance of comparisons */
    public static final double TOl = 1e-10;
    
    /**
     * The subject under test.
     */
    protected SequentialGaussNewtonOptimizer optimizer;

    public LeastSquaresBuilder base() {
        return new LeastSquaresBuilder()
                .checkerPair(new SimpleVectorValueChecker(1e-6, 1e-6))
                .maxEvaluations(100)
                .maxIterations(getMaxIterations());
    }

    public LeastSquaresBuilder builder(CircleVectorial c) {
        final double[] weights = new double[c.getN()];
        Arrays.fill(weights, 1.0);
        return base()
                .model(c.getModelFunction(), c.getModelFunctionJacobian())
                .target(new double[c.getN()])
                .weight(new DiagonalMatrix(weights));
    }

    public LeastSquaresBuilder builder(StatisticalReferenceDataset dataset) {
        StatisticalReferenceDataset.LeastSquaresProblem problem
                = dataset.getLeastSquaresProblem();
        final double[] weights = new double[dataset.getNumObservations()];
        Arrays.fill(weights, 1.0);
        return base()
                .model(problem.getModelFunction(), problem.getModelFunctionJacobian())
                .target(dataset.getData()[1])
                .weight(new DiagonalMatrix(weights))
                .start(dataset.getStartingPoint(0));
    }

    public void fail(LeastSquaresOptimizer optimizer) {
        Assert.fail("Expected Exception from: " + optimizer.toString());
    }

    /**
     * Check the value of a vector.
     * @param tolerance the absolute tolerance of comparisons
     * @param actual the vector to test
     * @param expected the expected values
     */
    public void assertEquals(double tolerance, RealVector actual, double... expected){
        for (int i = 0; i < expected.length; i++) {
            Assert.assertEquals(expected[i], actual.getEntry(i), tolerance);
        }
        Assert.assertEquals(expected.length, actual.getDimension());
    }

    /**
     * @return the default number of allowed iterations (which will be used when not
     *         specified otherwise).
     */
    public abstract int getMaxIterations();

    /**
     * Set an instance of the optimizer under test.
     * @param evaluation previous evaluation
     */
    public abstract void defineOptimizer(Evaluation evaluation);

    @Test
    public void testGetIterations() {
        LeastSquaresProblem lsp = base()
                .target(new double[]{1})
                .weight(new DiagonalMatrix(new double[]{1}))
                .start(new double[]{3})
                .model(new MultivariateJacobianFunction() {
                    public Pair<RealVector, RealMatrix> value(final RealVector point) {
                        return new Pair<RealVector, RealMatrix>(
                                new ArrayRealVector(
                                        new double[]{
                                                FastMath.pow(point.getEntry(0), 4)
                                        },
                                        false),
                                new Array2DRowRealMatrix(
                                        new double[][]{
                                                {0.25 * FastMath.pow(point.getEntry(0), 3)}
                                        },
                                        false)
                        );
                    }
                })
                .build();

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(lsp);

        //TODO more specific test? could pass with 'return 1;'
        Assert.assertTrue(optimum.getIterations() > 0);
    }

    @Test
    public void testTrivial() {
        LinearProblem problem
                = new LinearProblem(new double[][]{{2}},
                new double[]{3});
        LeastSquaresProblem ls = problem.getBuilder().build();

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(ls);

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 1.5);
        Assert.assertEquals(0.0, optimum.getResiduals().getEntry(0), TOl);
    }

    @Test
    public void testQRColumnsPermutation() {
        
        LinearProblem problem
                = new LinearProblem(new double[][]{{1, -1}, {0, 2}, {1, -2}},
                new double[]{4, 6, 1});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(problem.getBuilder().build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 7, 3);
        assertEquals(TOl, optimum.getResiduals(), 0, 0, 0);
    }

    @Test
    public void testNoDependency() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {2, 0, 0, 0, 0, 0},
                {0, 2, 0, 0, 0, 0},
                {0, 0, 2, 0, 0, 0},
                {0, 0, 0, 2, 0, 0},
                {0, 0, 0, 0, 2, 0},
                {0, 0, 0, 0, 0, 2}
        }, new double[]{0, 1.1, 2.2, 3.3, 4.4, 5.5});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(problem.getBuilder().build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        for (int i = 0; i < problem.target.length; ++i) {
            Assert.assertEquals(0.55 * i, optimum.getPoint().getEntry(i), TOl);
        }
    }

    @Test
    public void testOneSet() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 0, 0},
                {-1, 1, 0},
                {0, -1, 1}
        }, new double[]{1, 1, 1});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(problem.getBuilder().build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 1, 2, 3);
    }

    @Test
    public void testTwoSets() {
        double epsilon = 1e-7;
        LinearProblem problem = new LinearProblem(new double[][]{
                {2, 1, 0, 4, 0, 0},
                {-4, -2, 3, -7, 0, 0},
                {4, 1, -2, 8, 0, 0},
                {0, -3, -12, -1, 0, 0},
                {0, 0, 0, 0, epsilon, 1},
                {0, 0, 0, 0, 1, 1}
        }, new double[]{2, -9, 2, 2, 1 + epsilon * epsilon, 2});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(problem.getBuilder().build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 3, 4, -1, -2, 1 + epsilon, 1 - epsilon);
    }

    @Test
    public void testNonInvertible() throws Exception {
        try {
            LinearProblem problem = new LinearProblem(new double[][]{
                {1, 2, -3},
                {2, 1, 3},
                {-3, 0, -9}
            }, new double[]{1, 1, 1});
            
            defineOptimizer(null);

            optimizer.optimize(problem.getBuilder().build());

            fail(optimizer);
        } catch (MathIllegalArgumentException miae) {
            // expected
        } catch (MathIllegalStateException mise) {
            // expected
        }
    }

    @Test
    public void testIllConditioned() {
        LinearProblem problem1 = new LinearProblem(new double[][]{
                {10, 7, 8, 7},
                {7, 5, 6, 5},
                {8, 6, 10, 9},
                {7, 5, 9, 10}
        }, new double[]{32, 23, 33, 31});
        final double[] start = {0, 1, 2, 3};

        defineOptimizer(null);
        Optimum optimum = optimizer
                .optimize(problem1.getBuilder().start(start).build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 1, 1, 1, 1);

        LinearProblem problem2 = new LinearProblem(new double[][]{
                {10.00, 7.00, 8.10, 7.20},
                {7.08, 5.04, 6.00, 5.00},
                {8.00, 5.98, 9.89, 9.00},
                {6.99, 4.99, 9.00, 9.98}
        }, new double[]{32, 23, 33, 31});

        optimum = optimizer.optimize(problem2.getBuilder().start(start).build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(1e-8, optimum.getPoint(), -81, 137, -34, 22);
    }

    @Test
    public void testMoreEstimatedParametersSimple() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {3, 2, 0, 0},
                {0, 1, -1, 1},
                {2, 0, 1, 0}
        }, new double[]{7, 3, 5});

        defineOptimizer(null);
        Optimum optimum = optimizer
                .optimize(problem.getBuilder().start(new double[]{7, 6, 5, 4}).build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
    }

    @Test
    public void testMoreEstimatedParametersUnsorted() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 1, 0, 0, 0, 0},
                {0, 0, 1, 1, 1, 0},
                {0, 0, 0, 0, 1, -1},
                {0, 0, -1, 1, 0, 1},
                {0, 0, 0, -1, 1, 0}
        }, new double[]{3, 12, -1, 7, 1});

        defineOptimizer(null);
        Optimum optimum = optimizer.optimize(
                problem.getBuilder().start(new double[]{2, 2, 2, 2, 2, 2}).build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        RealVector point = optimum.getPoint();
        //the first two elements are under constrained
        //check first two elements obey the constraint: sum to 3
        Assert.assertEquals(3, point.getEntry(0) + point.getEntry(1), TOl);
        //#constrains = #states fro the last 4 elements
        assertEquals(TOl, point.getSubVector(2, 4), 3, 4, 5, 6);
    }

    @Test
    public void testRedundantEquations() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 1},
                {1, -1},
                {1, 3}
        }, new double[]{3, 1, 5});

        defineOptimizer(null);
        Optimum optimum = optimizer
                .optimize(problem.getBuilder().start(new double[]{1, 1}).build());

        Assert.assertEquals(0, optimum.getRMS(), TOl);
        assertEquals(TOl, optimum.getPoint(), 2, 1);
    }

    @Test
    public void testInconsistentEquations() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 1},
                {1, -1},
                {1, 3}
        }, new double[]{3, 1, 4});

        defineOptimizer(null);
        Optimum optimum = optimizer
                .optimize(problem.getBuilder().start(new double[]{1, 1}).build());

        //TODO what is this actually testing?
        Assert.assertTrue(optimum.getRMS() > 0.1);
    }

    @Test
    public void testInconsistentSizes1() {
        try {
            LinearProblem problem
                    = new LinearProblem(new double[][]{{1, 0},
                    {0, 1}},
                    new double[]{-1, 1});

            defineOptimizer(null);
            //TODO why is this part here? hasn't it been tested already?
            Optimum optimum = optimizer.optimize(problem.getBuilder().build());

            Assert.assertEquals(0, optimum.getRMS(), TOl);
            assertEquals(TOl, optimum.getPoint(), -1, 1);

            //TODO move to builder test
            optimizer.optimize(
                    problem.getBuilder().weight(new DiagonalMatrix(new double[]{1})).build());

            fail(optimizer);
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testInconsistentSizes2() {
        try {
            LinearProblem problem
                    = new LinearProblem(new double[][]{{1, 0}, {0, 1}},
                    new double[]{-1, 1});

            defineOptimizer(null);
            Optimum optimum = optimizer.optimize(problem.getBuilder().build());

            Assert.assertEquals(0, optimum.getRMS(), TOl);
            assertEquals(TOl, optimum.getPoint(), -1, 1);

            //TODO move to builder test
            optimizer.optimize(
                    problem.getBuilder()
                            .target(new double[]{1})
                            .weight(new DiagonalMatrix(new double[]{1}))
                            .build()
            );

            fail(optimizer);
        } catch (MathIllegalArgumentException e) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
        }
    }

    @Test
    public void testSequential() throws Exception {

        CircleVectorial circleAll    = new CircleVectorial();
        CircleVectorial circleFirst  = new CircleVectorial();
        CircleVectorial circleSecond = new CircleVectorial();
        circleAll.addPoint( 30.0,  68.0);
        circleFirst.addPoint( 30.0,  68.0);
        circleAll.addPoint( 50.0,  -6.0);
        circleFirst.addPoint( 50.0,  -6.0);
        circleAll.addPoint(110.0, -20.0);
        circleFirst.addPoint(110.0, -20.0);
        circleAll.addPoint( 35.0,  15.0);
        circleSecond.addPoint( 35.0,  15.0);
        circleAll.addPoint( 45.0,  97.0);
        circleSecond.addPoint( 45.0,  97.0);

        // first use the 5 observations in one run
        defineOptimizer(null);
        Optimum oneRun = optimizer.optimize(builder(circleAll).
                                            checkerPair(new SimpleVectorValueChecker(1e-3, 1e-3)).
                                            start(new double[] { 98.680, 47.345 }).
                                            build());
        Assert.assertEquals(2, oneRun.getJacobian().getColumnDimension());
        Vector2D oneShotCenter = new Vector2D(oneRun.getPoint().getEntry(0), oneRun.getPoint().getEntry(1));
        Assert.assertEquals(96.075901, oneShotCenter.getX(),               1.0e-6);
        Assert.assertEquals(48.135169, oneShotCenter.getY(),               1.0e-6);
        Assert.assertEquals(69.960161, circleAll.getRadius(oneShotCenter), 1.0e-6);

        // then split the observations in two sets,
        // the first one using 3 observations and the second one using 2 observations
        defineOptimizer(null);
        Optimum firstRun = optimizer.optimize(builder(circleFirst).
                                              checkerPair(new SimpleVectorValueChecker(1e-3, 1e-3)).
                                              start(new double[] { 98.680, 47.345 }).
                                              build());
        Assert.assertEquals(2, firstRun.getJacobian().getColumnDimension());
        Vector2D firstRunCenter = new Vector2D(firstRun.getPoint().getEntry(0), firstRun.getPoint().getEntry(1));
        Assert.assertEquals(93.650000, firstRunCenter.getX(),               1.0e-6);
        Assert.assertEquals(45.500000, firstRunCenter.getY(),               1.0e-6);
        Assert.assertEquals(67.896265, circleAll.getRadius(firstRunCenter), 1.0e-6);

        // for the second run, we start from the state and covariance only,
        // instead of using the evaluation "firstRun" that we have
        // (we could have used firstRun directly, but this is for testing this feature)
        optimizer = optimizer.withAPrioriData(firstRun.getPoint(), firstRun.getCovariances(1.0e-8));
        Optimum  secondRun       = optimizer.optimize(builder(circleSecond).
                                                      checkerPair(new SimpleVectorValueChecker(1e-3, 1e-3)).
                                                      start(new double[] { firstRunCenter.getX(), firstRunCenter.getY() }).
                                                      build());
        Assert.assertEquals(2, secondRun.getJacobian().getColumnDimension());
        Vector2D secondRunCenter = new Vector2D(secondRun.getPoint().getEntry(0), secondRun.getPoint().getEntry(1));
        Assert.assertEquals(97.070437, secondRunCenter.getX(),               1.0e-6);
        Assert.assertEquals(49.039898, secondRunCenter.getY(),               1.0e-6);
        Assert.assertEquals(70.789016, circleAll.getRadius(secondRunCenter), 1.0e-6);
        
    }

    protected void doTestStRD(final StatisticalReferenceDataset dataset,
                              final double errParams, final double errParamsSd) {

        defineOptimizer(null);
        final Optimum optimum = optimizer.optimize(builder(dataset).build());

        final RealVector actual = optimum.getPoint();
        for (int i = 0; i < actual.getDimension(); i++) {
            double expected = dataset.getParameter(i);
            double delta = FastMath.abs(errParams * expected);
            Assert.assertEquals(dataset.getName() + ", param #" + i,
                    expected, actual.getEntry(i), delta);
        }
    }

    @Test
    public void testKirby2() throws IOException {
        doTestStRD(StatisticalReferenceDatasetFactory.createKirby2(), 1E-7, 1E-7);
    }

    @Test
    public void testHahn1() throws IOException {
        doTestStRD(StatisticalReferenceDatasetFactory.createHahn1(), 1E-7, 1E-4);
    }

    @Test
    public void testPointCopy() {
        LinearProblem problem = new LinearProblem(new double[][]{
                {1, 0, 0},
                {-1, 1, 0},
                {0, -1, 1}
        }, new double[]{1, 1, 1});
        //mutable boolean
        final boolean[] checked = {false};

        final LeastSquaresBuilder builder = problem.getBuilder()
                .checker(new ConvergenceChecker<Evaluation>() {
                    public boolean converged(int iteration, Evaluation previous, Evaluation current) {
                        MatcherAssert.assertThat(
                                previous.getPoint(),
                                not(sameInstance(current.getPoint())));
                        Assert.assertArrayEquals(new double[3], previous.getPoint().toArray(), 0);
                        Assert.assertArrayEquals(new double[] {1, 2, 3}, current.getPoint().toArray(), TOl);
                        checked[0] = true;
                        return true;
                    }
                });
        defineOptimizer(null);
        optimizer.optimize(builder.build());

        MatcherAssert.assertThat(checked[0], is(true));
    }
    
    @Test
    public void testPointDifferentDim() {
        LinearProblem problem
        = new LinearProblem(new double[][]{{2}},
        new double[]{3});
        
        LeastSquaresProblem lsp = problem.getBuilder().build();
        defineOptimizer(new AbstractEvaluation(2){
            public RealMatrix getJacobian() {
                return MatrixUtils.createRealMatrix(2, 2);
            }
            public RealVector getPoint() {
                return MatrixUtils.createRealVector(2);
            }
            public RealVector getResiduals() {
                return MatrixUtils.createRealVector(2);
            }
        });
        try {
            optimizer.optimize(lsp);
            fail(optimizer);
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, mise.getSpecifier());
        }
    }

    class LinearProblem {
        private final RealMatrix factors;
        private final double[] target;

        public LinearProblem(double[][] factors, double[] target) {
            this.factors = new BlockRealMatrix(factors);
            this.target = target;
        }

        public double[] getTarget() {
            return target;
        }

        public MultivariateVectorFunction getModelFunction() {
            return new MultivariateVectorFunction() {
                public double[] value(double[] params) {
                    return factors.operate(params);
                }
            };
        }

        public MultivariateMatrixFunction getModelFunctionJacobian() {
            return new MultivariateMatrixFunction() {
                public double[][] value(double[] params) {
                    return factors.getData();
                }
            };
        }

        public LeastSquaresBuilder getBuilder() {
            final double[] weights = new double[target.length];
            Arrays.fill(weights, 1.0);
            return base()
                    .model(getModelFunction(), getModelFunctionJacobian())
                    .target(target)
                    .weight(new DiagonalMatrix(weights))
                    .start(new double[factors.getColumnDimension()]);
        }
    }

}
