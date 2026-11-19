/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
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

import java.util.Arrays;

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

/**
 * CUTEst problem {@code QR3DLS}.
 *
 * <p>The problem computes a QR factorization of a tridiagonal matrix
 * {@code A}. The unknowns are the entries of a dense orthogonal matrix
 * {@code Q} and the entries of an upper-triangular matrix {@code R}.</p>
 *
 * <p>The default SIF instance uses matrix order {@code 5}, hence it contains
 * 25 entries of {@code Q} and 15 stored entries of {@code R}, for a total of
 * 40 variables. The objective is the sum of squares of 15 orthogonality
 * residuals and 25 factorization residuals.</p>
 */
public class QR3DLSTest {

    /** Matrix order selected by the SIF test instance. */
    private static final int MATRIX_ORDER = 5;

    /** Number of dense Q variables. */
    private static final int Q_VARIABLES =
            MATRIX_ORDER * MATRIX_ORDER;

    /** Number of upper-triangular R variables. */
    private static final int R_VARIABLES =
            MATRIX_ORDER * (MATRIX_ORDER + 1) / 2;

    /** Total number of optimization variables. */
    private static final int N =
            Q_VARIABLES + R_VARIABLES;

    /** Value reported by the SIF LO SOLTN field. */
    private static final double EXPECTED_OBJECTIVE = 0.0;

    /** Matrix appearing in A = QR. */
    private static final RealMatrix A = buildMatrixA();

    /** Official CUTEst starting point. */
    private static final double[] START = buildStartPoint();

    /** Variable lower bounds. */
    private static final double[] LOWER = buildLowerBounds();

    /** Variable upper bounds. */
    private static final double[] UPPER = buildUpperBounds();

    /**
     * QR least-squares objective.
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

            /*
             * Orthogonality residuals:
             *
             *     sum_k Q(i,k) Q(j,k) - delta(i,j)
             *
             * Only the upper triangle is included, exactly as in the SIF.
             */
            for (int i = 0; i < MATRIX_ORDER; ++i) {
                for (int j = i; j < MATRIX_ORDER; ++j) {

                    double scalarProduct = 0.0;

                    for (int k = 0; k < MATRIX_ORDER; ++k) {
                        scalarProduct +=
                                point.getEntry(qIndex(i, k)) *
                                point.getEntry(qIndex(j, k));
                    }

                    final double residual =
                            scalarProduct -
                            (i == j ? 1.0 : 0.0);

                    objective += residual * residual;
                }
            }

            /*
             * Factorization residuals:
             *
             *     sum_{k=0}^{j} Q(i,k) R(k,j) - A(i,j)
             *
             * R(k,j) exists only for k <= j.
             */
            for (int i = 0; i < MATRIX_ORDER; ++i) {
                for (int j = 0; j < MATRIX_ORDER; ++j) {

                    double product = 0.0;

                    for (int k = 0; k <= j; ++k) {
                        product +=
                                point.getEntry(qIndex(i, k)) *
                                point.getEntry(rIndex(k, j));
                    }

                    final double residual =
                            product - A.getEntry(i, j);

                    objective += residual * residual;
                }
            }

            return objective;
        }

        @Override
        public RealVector gradient(final RealVector point) {

            final double[] gradient = new double[N];

            /*
             * Gradient of the orthogonality least-squares terms.
             */
            for (int i = 0; i < MATRIX_ORDER; ++i) {
                for (int j = i; j < MATRIX_ORDER; ++j) {

                    double scalarProduct = 0.0;

                    for (int k = 0; k < MATRIX_ORDER; ++k) {
                        scalarProduct +=
                                point.getEntry(qIndex(i, k)) *
                                point.getEntry(qIndex(j, k));
                    }

                    final double residual =
                            scalarProduct -
                            (i == j ? 1.0 : 0.0);

                    if (i == j) {
                        for (int k = 0; k < MATRIX_ORDER; ++k) {
                            gradient[qIndex(i, k)] +=
                                    4.0 * residual *
                                    point.getEntry(qIndex(i, k));
                        }
                    } else {
                        for (int k = 0; k < MATRIX_ORDER; ++k) {
                            gradient[qIndex(i, k)] +=
                                    2.0 * residual *
                                    point.getEntry(qIndex(j, k));

                            gradient[qIndex(j, k)] +=
                                    2.0 * residual *
                                    point.getEntry(qIndex(i, k));
                        }
                    }
                }
            }

            /*
             * Gradient of the factorization least-squares terms.
             */
            for (int i = 0; i < MATRIX_ORDER; ++i) {
                for (int j = 0; j < MATRIX_ORDER; ++j) {

                    double product = 0.0;

                    for (int k = 0; k <= j; ++k) {
                        product +=
                                point.getEntry(qIndex(i, k)) *
                                point.getEntry(rIndex(k, j));
                    }

                    final double residual =
                            product - A.getEntry(i, j);

                    for (int k = 0; k <= j; ++k) {

                        final int q = qIndex(i, k);
                        final int r = rIndex(k, j);

                        gradient[q] +=
                                2.0 * residual *
                                point.getEntry(r);

                        gradient[r] +=
                                2.0 * residual *
                                point.getEntry(q);
                    }
                }
            }

            return new ArrayRealVector(gradient, false);
        }

        @Override
        public RealMatrix hessian(final RealVector point) {
            /*
             * SQPOptimizerS2 uses its BFGS approximation in this test.
             */
            return new OpenMapRealMatrix(N, N);
        }
    }

    /**
     * Build the matrix A exactly as specified by the SIF instance.
     *
     * <p>For the last diagonal entry the SIF explicitly assigns
     * {@code 2 * M}, not {@code 2 * M / M}. This implementation deliberately
     * preserves that literal definition.</p>
     *
     * @return matrix A
     */
    private static RealMatrix buildMatrixA() {

        final RealMatrix matrix =
                new OpenMapRealMatrix(
                        MATRIX_ORDER, MATRIX_ORDER);

        final double order = MATRIX_ORDER;

        matrix.setEntry(0, 0, 2.0 / order);
        matrix.setEntry(0, 1, 0.0);

        for (int i = 1; i < MATRIX_ORDER - 1; ++i) {

            final double offDiagonal =
                    (1.0 - (i + 1.0)) / order;

            matrix.setEntry(i, i - 1, offDiagonal);
            matrix.setEntry(
                    i, i,
                    2.0 * (i + 1.0) / order);
            matrix.setEntry(i, i + 1, offDiagonal);
        }

        matrix.setEntry(
                MATRIX_ORDER - 1,
                MATRIX_ORDER - 2,
                (1.0 - MATRIX_ORDER) / order);

        matrix.setEntry(
                MATRIX_ORDER - 1,
                MATRIX_ORDER - 1,
                2.0 * MATRIX_ORDER);

        return matrix;
    }

    /**
     * Build the CUTEst starting point.
     *
     * <p>{@code Q} starts at the identity. {@code R} starts from the upper
     * triangular part of {@code A}.</p>
     *
     * @return starting point
     */
    private static double[] buildStartPoint() {

        final double[] start = new double[N];

        for (int i = 0; i < MATRIX_ORDER; ++i) {
            start[qIndex(i, i)] = 1.0;
        }

        for (int i = 0; i < MATRIX_ORDER - 1; ++i) {
            start[rIndex(i, i)] = A.getEntry(i, i);
            start[rIndex(i, i + 1)] = A.getEntry(i, i + 1);
        }

        start[rIndex(MATRIX_ORDER - 1,
                     MATRIX_ORDER - 1)] =
                A.getEntry(MATRIX_ORDER - 1, MATRIX_ORDER - 1);

        return start;
    }

    /**
     * Build lower variable bounds.
     *
     * <p>All variables are free except the diagonal entries of {@code R},
     * which are nonnegative.</p>
     *
     * @return lower bounds
     */
    private static double[] buildLowerBounds() {

        final double[] lower = new double[N];
        Arrays.fill(lower, Double.NEGATIVE_INFINITY);

        for (int i = 0; i < MATRIX_ORDER; ++i) {
            lower[rIndex(i, i)] = 0.0;
        }

        return lower;
    }

    /**
     * Build upper variable bounds.
     *
     * @return upper bounds
     */
    private static double[] buildUpperBounds() {

        final double[] upper = new double[N];
        Arrays.fill(upper, Double.POSITIVE_INFINITY);
        return upper;
    }

    /**
     * Get the row-major index of Q(i,j).
     *
     * @param row row index
     * @param column column index
     * @return vector index
     */
    private static int qIndex(final int row,
                              final int column) {
        return row * MATRIX_ORDER + column;
    }

    /**
     * Get the packed upper-triangular index of R(i,j).
     *
     * @param row row index
     * @param column column index
     * @return vector index
     */
    private static int rIndex(final int row,
                              final int column) {

        if (row < 0 ||
            row >= MATRIX_ORDER ||
            column < row ||
            column >= MATRIX_ORDER) {
            throw new IllegalArgumentException(
                    "Invalid upper-triangular R index (" +
                    row + ", " + column + ").");
        }

        int index = Q_VARIABLES;

        for (int i = 0; i < row; ++i) {
            index += MATRIX_ORDER - i;
        }

        return index + column - row;
    }

    @Test
    public void testQR3DLS() {

        final SQPOptimizerS2 optimizer =
                CUTEstProblemUtils.newOptimizer();

        final SQPOption option =
                CUTEstProblemUtils.newForwardDifferenceOption();

        final LagrangeSolution solution =
                optimizer.optimize(
                        new MaxIter(10000),
                        new InitialGuess(START),
                        new ObjectiveFunction(new Objective()),
                        new SimpleBounds(LOWER, UPPER),
                        option);

        CUTEstProblemUtils.assertExpectedObjective(
                EXPECTED_OBJECTIVE,
                solution);
    }
}