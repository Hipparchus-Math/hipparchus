/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import java.util.Arrays;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code AIRPORT}.
 *
 * <p>The problem selects one point for each of 42 Brazilian cities and
 * minimizes the sum of the squared pairwise distances. Each selected point
 * must remain inside the disk associated with its city. The original CUTEst dimensions are retained: 84 variables and
 * 42 nonlinear inequality constraints.</p>
 */
public class AIRPORTTest {

    /** Number of cities. */
    private static final int CITIES = 42;

    /** Two coordinates for every city. */
    private static final int N = 2 * CITIES;

    /** Number of nonlinear inequalities. */
    private static final int M = CITIES;

    /** Standard CUTEst reference objective. */
    private static final double EXPECTED_OBJECTIVE =
            4.7952695811e+04;

    /** Squared-distance upper bounds from AIRPORT.SIF. */
    private static final double[] RADII = {
            0.09, 0.3, 0.09, 0.45, 0.5, 0.04,
            0.1, 0.02, 0.02, 0.07, 0.4, 0.045,
            0.05, 0.056, 0.36, 0.08, 0.07, 0.36,
            0.67, 0.38, 0.37, 0.05, 0.4, 0.66,
            0.05, 0.07, 0.08, 0.3, 0.31, 0.49,
            0.09, 0.46, 0.12, 0.07, 0.07, 0.09,
            0.05, 0.13, 0.16, 0.46, 0.25, 0.1
    };

    /** X coordinates of the city centers. */
    private static final double[] CENTER_X = {
            -6.3, -7.8, -9, -7.2, -5.7, -1.9,
            -3.5, -0.5, 1.4, 4, 2.1, 5.5,
            5.7, 5.7, 3.8, 5.3, 4.7, 3.3,
            0, -1, -0.4, 4.2, 3.2, 1.7,
            3.3, 2, 0.7, 0.1, -0.1, -3.5,
            -4, -2.7, -0.5, -2.9, -1.2, -0.4,
            -0.1, -1, -1.7, -2.1, -1.8, 0
    };

    /** Y coordinates of the city centers. */
    private static final double[] CENTER_Y = {
            8, 5.1, 2, 2.6, 5.5, 7.1,
            5.9, 6.6, 6.1, 5.6, 4.9, 4.7,
            4.3, 3.6, 4.1, 3, 2.4, 3,
            4.7, 3.4, 2.3, 1.5, 0.5, -1.7,
            -2, -3.1, -3.5, -2.4, -1.3, 0,
            -1.7, -2.1, -0.4, -2.9, -3.4, -4.3,
            -5.2, -6.5, -7.5, -6.4, -5.1, 0
    };

    /**
     * Sum of squared pairwise distances:
     *
     * <pre>
     * f(x, y) = sum(i &lt; j)
     *           [ (x_i - x_j)^2 + (y_i - y_j)^2 ].
     * </pre>
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {
            double value = 0.0;

            for (int i = 0; i < CITIES - 1; ++i) {
                final double xi =
                        point.getEntry(xIndex(i));

                final double yi =
                        point.getEntry(yIndex(i));

                for (int j = i + 1; j < CITIES; ++j) {
                    final double dx =
                            xi - point.getEntry(xIndex(j));

                    final double dy =
                            yi - point.getEntry(yIndex(j));

                    value += dx * dx + dy * dy;
                }
            }

            return value;
        }

        @Override
        public RealVector gradient(final RealVector point) {
            final RealVector gradient =
                    new ArrayRealVector(N);

            for (int i = 0; i < CITIES - 1; ++i) {
                final double xi =
                        point.getEntry(xIndex(i));

                final double yi =
                        point.getEntry(yIndex(i));

                for (int j = i + 1; j < CITIES; ++j) {
                    final double dx =
                            xi - point.getEntry(xIndex(j));

                    final double dy =
                            yi - point.getEntry(yIndex(j));

                    final double gx = 2.0 * dx;
                    final double gy = 2.0 * dy;

                    gradient.addToEntry(xIndex(i), gx);
                    gradient.addToEntry(yIndex(i), gy);

                    gradient.addToEntry(xIndex(j), -gx);
                    gradient.addToEntry(yIndex(j), -gy);
                }
            }

            return gradient;
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            final RealMatrix hessian =
                    new OpenMapRealMatrix(N, N);

            for (int i = 0; i < CITIES - 1; ++i) {
                for (int j = i + 1; j < CITIES; ++j) {
                    addPairBlock(
                            hessian,
                            xIndex(i),
                            xIndex(j));

                    addPairBlock(
                            hessian,
                            yIndex(i),
                            yIndex(j));
                }
            }

            return hessian;
        }

        private static void addPairBlock(final RealMatrix hessian,
                                         final int first,
                                         final int second) {
            hessian.addToEntry(first, first, 2.0);
            hessian.addToEntry(second, second, 2.0);
            hessian.addToEntry(first, second, -2.0);
            hessian.addToEntry(second, first, -2.0);
        }
    }

    /**
     * Quadratic disk constraints:
     *
     * <pre>
     * (x_i - cx_i)^2 + (y_i - cy_i)^2 <= r_i,
     * i = 0,...,41.
     * </pre>
     *
     * <p>Hipparchus represents inequalities as c(x) >= 0, hence
     * c_i(x) = r_i - (x_i - cx_i)^2 - (y_i - cy_i)^2.</p>
     */
    private static final class AirportInequalities
            extends InequalityConstraint {

        AirportInequalities() {
            super(new ArrayRealVector(M));
        }

        @Override
        public int dim() {
            return N;
        }

        @Override
        public RealVector value(final RealVector point) {
            final RealVector constraints =
                    new ArrayRealVector(M);

            for (int i = 0; i < CITIES; ++i) {
                final double dx =
                        point.getEntry(xIndex(i)) -
                        CENTER_X[i];

                final double dy =
                        point.getEntry(yIndex(i)) -
                        CENTER_Y[i];

                constraints.setEntry(
                        i,
                        RADII[i] - dx * dx - dy * dy);
            }

            return constraints;
        }

        @Override
        public RealMatrix jacobian(final RealVector point) {
            final RealMatrix jacobian =
                    new OpenMapRealMatrix(M, N);

            for (int i = 0; i < CITIES; ++i) {
                jacobian.setEntry(
                        i,
                        xIndex(i),
                        -2.0 *
                        (point.getEntry(xIndex(i)) -
                         CENTER_X[i]));

                jacobian.setEntry(
                        i,
                        yIndex(i),
                        -2.0 *
                        (point.getEntry(yIndex(i)) -
                         CENTER_Y[i]));
            }

            return jacobian;
        }
    }

    /**
     * AIRPORT.SIF does not define an explicit start point, so CUTEst uses
     * the default zero vector.
     */
    private static double[] initialPoint() {
        return new double[N];
    }

    /** Original variable order: X(1), Y(1), X(2), Y(2), ... */
    private static int xIndex(final int city) {
        return 2 * city;
    }

    private static int yIndex(final int city) {
        return 2 * city + 1;
    }

    private static SimpleBounds bounds() {
        final double[] lower = new double[N];
        final double[] upper = new double[N];

        Arrays.fill(lower, -10.0);
        Arrays.fill(upper, 10.0);

        return new SimpleBounds(lower, upper);
    }

    @Test
    public void testAIRPORT() {
        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(initialPoint()),
                        new ObjectiveFunction(new Objective()),
                        new AirportInequalities(),
                        bounds(),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}