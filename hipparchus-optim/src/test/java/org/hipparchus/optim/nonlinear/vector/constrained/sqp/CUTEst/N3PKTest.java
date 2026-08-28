/*
 * Licensed to the Hipparchus project under one or more contributor
 * license agreements. See the NOTICE file distributed with this work
 * for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.CUTEst;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxIter;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.junit.jupiter.api.Test;

/** CUTEst problem {@code 3PK}. */
public class N3PKTest {

    private static final int N = 30;
    private static final int M = 42;
    private static final double EXPECTED_OBJECTIVE = 1.720119;

    private static final RealMatrix DESIGN = new OpenMapRealMatrix(M, N);
    private static final double[] TARGET = new double[M];
    private static final double[] SCALE = new double[M];

    private static final double[] START = {
        0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 0.5, 0.5, 0.5, 0.5,
        0.5, 0.5, 0.5, 0.5, 0.5,
        100.0, 140.0, 120.0, 20.0, 20.0,
        200.0, 180.0, 20.0, 600.0, 40.0,
        50.0, 30.0, 70.0, 150.0, 20.0
    };

    private static final double[] LOWER = constantArray(0.0);
    private static final double[] UPPER =
            constantArray(Double.POSITIVE_INFINITY);

    static {
        initializeModel();
    }

    private static final class Objective extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            double objective = 0.0;

            for (int row = 0; row < M; ++row) {
                final double residual =
                        rowProduct(row, point) - TARGET[row];
                objective += residual * residual / SCALE[row];
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            for (int row = 0; row < M; ++row) {
                final double residual =
                        rowProduct(row, point) - TARGET[row];
                final double factor =
                        2.0 * residual / SCALE[row];

                for (int column = 0; column < N; ++column) {
                    final double coefficient =
                            DESIGN.getEntry(row, column);

                    if (coefficient != 0.0) {
                        gradient[column] +=
                                factor * coefficient;
                    }
                }
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            return new OpenMapRealMatrix(N, N);
        }
    }

    private static void initializeModel() {

        int row = 0;

        final double[] prior = {
            0.010000, 0.007143, 0.008333, 0.050000, 0.050000,
            0.005000, 0.005556, 0.050000, 0.001667, 0.025000,
            0.020000, 0.033333, 0.014286, 0.006667, 0.050000
        };

        for (int i = 0; i < prior.length; ++i) {
            setRow(row++, 1.0, 1.0,
                   new int[] {15 + i},
                   new double[] {prior[i]});
        }

        for (int block = 0; block < 3; ++block) {
            final int first = 5 * block;
            setRow(row++, 1.0, 1.0,
                   new int[] {
                       first, first + 1, first + 2,
                       first + 3, first + 4
                   },
                   new double[] {
                       0.4, 0.4, 0.4, 0.4, 0.4
                   });
        }

        final int countingBase = row;

        for (int k = 0; k < 9; ++k) {
            TARGET[countingBase + k] = 1.0;
            SCALE[countingBase + k] = 1.0e-4;
        }

        add(countingBase + 7, 0, 200.0 / 1450.0);
        add(countingBase + 4, 0, 200.0 / 260.0);
        add(countingBase + 2, 0, 200.0 / 1915.0);

        add(countingBase + 8, 1, 480.0 / 990.0);
        add(countingBase + 7, 1, 480.0 / 1450.0);
        add(countingBase + 6, 1, 480.0 / 670.0);
        add(countingBase + 2, 1, 480.0 / 1915.0);

        add(countingBase + 2, 2, 120.0 / 1915.0);

        add(countingBase + 7, 3, 360.0 / 1450.0);
        add(countingBase + 2, 3, 360.0 / 1915.0);

        add(countingBase + 8, 4, 560.0 / 990.0);
        add(countingBase + 7, 4, 560.0 / 1450.0);
        add(countingBase + 2, 4, 560.0 / 1915.0);

        add(countingBase + 0, 5, 240.0 / 910.0);

        add(countingBase + 8, 6, 400.0 / 990.0);
        add(countingBase + 7, 6, 400.0 / 1450.0);
        add(countingBase + 6, 6, 400.0 / 670.0);
        add(countingBase + 2, 6, 400.0 / 1915.0);
        add(countingBase + 0, 6, 400.0 / 910.0);

        add(countingBase + 2, 7, 420.0 / 1915.0);
        add(countingBase + 0, 7, 420.0 / 910.0);

        add(countingBase + 7, 8, 180.0 / 1450.0);
        add(countingBase + 2, 8, 180.0 / 1915.0);
        add(countingBase + 0, 8, 180.0 / 910.0);

        add(countingBase + 8, 9, 320.0 / 990.0);
        add(countingBase + 7, 9, 320.0 / 1450.0);
        add(countingBase + 2, 9, 320.0 / 1915.0);
        add(countingBase + 0, 9, 320.0 / 910.0);

        add(countingBase + 1, 10, 20.0 / 175.0);
        add(countingBase + 0, 10, 20.0 / 910.0);

        add(countingBase + 1, 11, 60.0 / 175.0);

        add(countingBase + 2, 12, 40.0 / 1915.0);
        add(countingBase + 1, 12, 40.0 / 175.0);
        add(countingBase + 0, 12, 40.0 / 910.0);

        add(countingBase + 5, 13, 120.0 / 80.0);

        add(countingBase + 8, 14, 20.0 / 990.0);
        add(countingBase + 5, 14, 20.0 / 80.0);

        add(countingBase + 7, 15, 1.0 / 1450.0);
        add(countingBase + 3, 15, 1.0 / 450.0);

        add(countingBase + 7, 16, 1.0 / 1450.0);
        add(countingBase + 4, 16, 1.0 / 260.0);

        add(countingBase + 8, 17, 1.0 / 990.0);
        add(countingBase + 7, 17, 1.0 / 1450.0);
        add(countingBase + 6, 17, 1.0 / 670.0);

        add(countingBase + 7, 18, 1.0 / 1450.0);

        add(countingBase + 8, 19, 1.0 / 990.0);
        add(countingBase + 7, 19, 1.0 / 1450.0);

        add(countingBase + 3, 20, 1.0 / 450.0);
        add(countingBase + 4, 21, 1.0 / 260.0);

        add(countingBase + 8, 22, 1.0 / 990.0);
        add(countingBase + 6, 22, 1.0 / 670.0);

        add(countingBase + 3, 23, 1.0 / 450.0);
        add(countingBase + 2, 23, 1.0 / 1915.0);

        add(countingBase + 8, 24, 1.0 / 990.0);

        add(countingBase + 6, 25, 1.0 / 670.0);
        add(countingBase + 1, 25, 1.0 / 175.0);
        add(countingBase + 0, 25, 1.0 / 910.0);

        add(countingBase + 6, 26, 1.0 / 670.0);
        add(countingBase + 1, 26, 1.0 / 175.0);

        add(countingBase + 6, 27, 1.0 / 670.0);

        add(countingBase + 6, 28, 1.0 / 670.0);
        add(countingBase + 2, 28, 1.0 / 1915.0);
        add(countingBase + 1, 28, 1.0 / 175.0);
        add(countingBase + 0, 28, 1.0 / 910.0);

        add(countingBase + 6, 29, 1.0 / 670.0);
        add(countingBase + 5, 29, 1.0 / 80.0);

        row += 9;

        for (int block = 0; block < 3; ++block) {
            final int first = 5 * block;

            for (int selected = 0; selected < 5; ++selected) {
                final int[] indices = new int[5];
                final double[] coefficients = new double[5];

                for (int j = 0; j < 5; ++j) {
                    indices[j] = first + j;
                    coefficients[j] =
                            j == selected ? -0.8 : 0.2;
                }

                setRow(row++, 0.0, 0.5,
                       indices, coefficients);
            }
        }
    }

    private static void setRow(final int row,
                               final double target,
                               final double scale,
                               final int[] indices,
                               final double[] coefficients) {

        TARGET[row] = target;
        SCALE[row] = scale;

        for (int i = 0; i < indices.length; ++i) {
            DESIGN.setEntry(row, indices[i], coefficients[i]);
        }
    }

    private static void add(final int row,
                            final int column,
                            final double value) {
        DESIGN.addToEntry(row, column, value);
    }

    private static double rowProduct(final int row,
                                     final RealVector point) {

        double value = 0.0;

        for (int column = 0; column < N; ++column) {
            final double coefficient =
                    DESIGN.getEntry(row, column);

            if (coefficient != 0.0) {
                value += coefficient *
                         point.getEntry(column);
            }
        }

        return value;
    }

    private static double[] constantArray(final double value) {

        final double[] array = new double[N];

        for (int i = 0; i < N; ++i) {
            array[i] = value;
        }

        return array;
    }

    @Test
    public void test3PK() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newExternalOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(5000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}