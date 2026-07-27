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
import org.hipparchus.linear.MatrixUtils;
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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/**
 * CUTEst problem {@code NONMSQRT}.
 *
 * <p>This is the "non-matrix square-root problem", obtained from an error
 * in the formulation of a matrix square-root problem by Nocedal and Liu.
 * The original CUTEst instance uses a 3 by 3 matrix, hence 9 variables.</p>
 *
 * <p>The objective is</p>
 *
 * <pre>
 * f(X) = sum(i,j) [X(i,j) sum(k) X(i,k) - A(i,j)]^2,
 * </pre>
 *
 * <p>where {@code A = B B} and {@code B(i,j) = sin(k^2)} in row-major
 * order, except for {@code B(3,1) = 0}.</p>
 */
public class NONMSQRTTest {

    /** Matrix dimension used by the original CUTEst instance. */
    private static final int P = 3;

    /** Number of optimization variables. */
    private static final int N = P * P;

    /** Expected objective reported by the SIF for P = 3. */
    private static final double EXPECTED_OBJECTIVE = 0.7518002333474424;

    /** Matrix B defined by the SIF. */
    private static final RealMatrix B = buildB();

    /** Constant matrix A = B B. */
    private static final RealMatrix A = B.multiply(B);

    /** Official starting point. */
    private static final double[] START = buildStartPoint();

    /** All variables are free. */
    private static final double[] LOWER = buildBounds(
            Double.NEGATIVE_INFINITY);

    /** All variables are free. */
    private static final double[] UPPER = buildBounds(
            Double.POSITIVE_INFINITY);

    /**
     * NONMSQRT objective.
     */
    private static final class Objective
            extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return N;
        }

        @Override
        public double value(final RealVector point) {

            double objective = 0.0;

            for (int i = 0; i < P; ++i) {

                double rowSum = 0.0;
                for (int k = 0; k < P; ++k) {
                    rowSum += point.getEntry(index(i, k));
                }

                for (int j = 0; j < P; ++j) {
                    final double residual =
                            point.getEntry(index(i, j)) *
                            rowSum -
                            A.getEntry(i, j);

                    objective += residual * residual;
                }
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            for (int i = 0; i < P; ++i) {

                double rowSum = 0.0;
                for (int k = 0; k < P; ++k) {
                    rowSum += point.getEntry(index(i, k));
                }

                final double[] residuals =
                        new double[P];

                double residualDotRow = 0.0;

                for (int j = 0; j < P; ++j) {
                    final double xij =
                            point.getEntry(index(i, j));

                    final double residual =
                            xij * rowSum -
                            A.getEntry(i, j);

                    residuals[j] = residual;
                    residualDotRow += residual * xij;
                }

                /*
                 * For
                 *
                 * r(i,j) = X(i,j) s(i) - A(i,j),
                 * s(i)   = sum(k) X(i,k),
                 *
                 * d r(i,j) / d X(i,l)
                 *     = delta(j,l) s(i) + X(i,j).
                 */
                for (int j = 0; j < P; ++j) {
                    gradient[index(i, j)] =
                            2.0 *
                            (rowSum * residuals[j] +
                             residualDotRow);
                }
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            /*
             * SQPOptimizerS2 uses its BFGS approximation in this test.
             * The exact Hessian is therefore intentionally not supplied.
             */
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Build matrix B exactly as specified by the SIF.
     *
     * @return matrix B
     */
    private static RealMatrix buildB() {

        final RealMatrix matrix =
                MatrixUtils.createRealMatrix(P, P);

        int k = 0;

        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < P; ++j) {
                ++k;
                matrix.setEntry(
                        i,
                        j,
                        FastMath.sin((double) k * k));
            }
        }

        /*
         * The SIF explicitly overwrites B(3,1).
         * Indices are zero-based in Java.
         */
        matrix.setEntry(2, 0, 0.0);

        return matrix;
    }

    /**
     * Build the official starting point.
     *
     * <p>For ordinary entries, {@code X(i,j) = 0.2 B(i,j)}. The SIF
     * explicitly overrides {@code X(3,1)} with
     * {@code -0.8 sin(7^2)}.</p>
     *
     * @return starting point in row-major order
     */
    private static double[] buildStartPoint() {

        final double[] start = new double[N];

        for (int i = 0; i < P; ++i) {
            for (int j = 0; j < P; ++j) {
                start[index(i, j)] =
                        0.2 * B.getEntry(i, j);
            }
        }

        start[index(2, 0)] =
                -0.8 * FastMath.sin(49.0);

        return start;
    }

    /**
     * Build a constant bound vector.
     *
     * @param value bound value
     * @return bound vector
     */
    private static double[] buildBounds(final double value) {

        final double[] bounds = new double[N];

        for (int i = 0; i < N; ++i) {
            bounds[i] = value;
        }

        return bounds;
    }

    /**
     * Convert a matrix position to the row-major optimization-vector index.
     *
     * @param row matrix row
     * @param column matrix column
     * @return vector index
     */
    private static int index(final int row,
                             final int column) {
        return row * P + column;
    }

    @Test
    public void testNONMSQRT() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newExternalOption();
        option.setMaxteration(1000);

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