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
 * This is not the original file distributed by the Apache Software Foundation.
 * It has been extended by the Hipparchus project.
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.DecompositionSolver;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;

/**
 * Calculates a Cholesky decomposition with adaptive diagonal regularization.
 * <p>
 * Given a symmetric matrix {@code A}, this class attempts to compute a lower
 * triangular matrix {@code L} such that:
 * </p>
 * <pre>
 * A + shift * I = L L^T
 * </pre>
 * <p>
 * where {@code shift >= 0} is chosen automatically when needed to enforce
 * strict positive definiteness.
 * </p>
 * <p>
 * The decomposition is therefore performed on the regularized matrix, which can
 * be retrieved using {@link #getRegularizedMatrix()}.
 * </p>
 *
 * @see org.hipparchus.linear.CholeskyDecomposition
 * @see org.hipparchus.linear.SemiDefinitePositiveCholeskyDecomposition
 */
public class RegularizedCholeskyDecomposition {

    /**
     * Default threshold above which off-diagonal elements are considered too
     * different and the matrix is not symmetric.
     */
    public static final double DEFAULT_RELATIVE_SYMMETRY_THRESHOLD = 1.0e-15;

    /**
     * Default threshold below which diagonal pivots are considered too small.
     */
    public static final double DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD = 1.0e-14;

    /**
     * Default growth factor for the diagonal shift.
     */
    public static final double DEFAULT_DIAGONAL_GROWTH_FACTOR = 2.0;

    /**
     * Default maximum number of regularization attempts.
     */
    public static final int DEFAULT_MAX_REGULARIZATION_STEPS = 50;

    /**
     * Upper-triangular row-oriented storage for L^T.
     */
    private final double[][] lTData;

    /**
     * Cached value of L.
     */
    private RealMatrix cachedL;

    /**
     * Cached value of LT.
     */
    private RealMatrix cachedLT;

    /**
     * Regularized matrix actually factorized.
     */
    private final RealMatrix regularizedMatrix;

    /**
     * Final diagonal shift added to the input matrix.
     */
    private final double diagonalShift;

    /**
     * Build the decomposition using default thresholds.
     *
     * @param matrix matrix to decompose
     */
    public RegularizedCholeskyDecomposition(final RealMatrix matrix) {
        this(matrix,
             DEFAULT_RELATIVE_SYMMETRY_THRESHOLD,
             DEFAULT_ABSOLUTE_POSITIVITY_THRESHOLD,
             DEFAULT_DIAGONAL_GROWTH_FACTOR,
             DEFAULT_MAX_REGULARIZATION_STEPS);
    }

    /**
     * Build the decomposition using custom symmetry and positivity thresholds.
     *
     * @param matrix matrix to decompose
     * @param relativeSymmetryThreshold symmetry threshold
     * @param absolutePositivityThreshold positivity threshold
     */
    public RegularizedCholeskyDecomposition(final RealMatrix matrix,
                                            final double relativeSymmetryThreshold,
                                            final double absolutePositivityThreshold) {
        this(matrix,
             relativeSymmetryThreshold,
             absolutePositivityThreshold,
             DEFAULT_DIAGONAL_GROWTH_FACTOR,
             DEFAULT_MAX_REGULARIZATION_STEPS);
    }

    /**
     * Build the decomposition using fully custom settings.
     *
     * @param matrix matrix to decompose
     * @param relativeSymmetryThreshold symmetry threshold
     * @param absolutePositivityThreshold positivity threshold
     * @param diagonalGrowthFactor growth factor applied to the diagonal shift
     * @param maxRegularizationSteps maximum number of regularization attempts
     */
    public RegularizedCholeskyDecomposition(final RealMatrix matrix,
                                            final double relativeSymmetryThreshold,
                                            final double absolutePositivityThreshold,
                                            final double diagonalGrowthFactor,
                                            final int maxRegularizationSteps) {
        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(),
                                                   matrix.getColumnDimension());
        }

        final double[][] original = copyAndCheckSymmetry(matrix, relativeSymmetryThreshold);
        final FactorizationResult result = factorizeWithRegularization(original,
                                                                       absolutePositivityThreshold,
                                                                       diagonalGrowthFactor,
                                                                       maxRegularizationSteps);

        this.lTData = result.lTData;
        this.regularizedMatrix = new Array2DRowRealMatrix(result.regularizedMatrix, false);
        this.diagonalShift = result.diagonalShift;
        this.cachedL = null;
        this.cachedLT = null;
    }

    /**
     * Get the final diagonal shift added to the original matrix.
     *
     * @return diagonal shift
     */
    public double getDiagonalShift() {
        return diagonalShift;
    }

    /**
     * Get the regularized matrix actually factorized.
     *
     * @return regularized matrix
     */
    public RealMatrix getRegularizedMatrix() {
        return regularizedMatrix;
    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>
     * L is a lower-triangular matrix.
     * </p>
     *
     * @return the L matrix
     */
    public RealMatrix getL() {
        if (cachedL == null) {
            cachedL = getLT().transpose();
        }
        return cachedL;
    }

    /**
     * Returns the transpose of the matrix L of the decomposition.
     * <p>
     * L<sup>T</sup> is an upper-triangular matrix.
     * </p>
     *
     * @return the transpose of the matrix L
     */
    public RealMatrix getLT() {
        if (cachedLT == null) {
            cachedLT = MatrixUtils.createRealMatrix(lTData);
        }
        return cachedLT;
    }

    /**
     * Return the determinant of the regularized matrix.
     *
     * @return determinant of the regularized matrix
     */
    public double getDeterminant() {
        double determinant = 1.0;
        for (int i = 0; i < lTData.length; ++i) {
            final double d = lTData[i][i];
            determinant *= d * d;
        }
        return determinant;
    }

    /**
     * Get a solver for finding the solution of the regularized system.
     *
     * @return a solver
     */
    public DecompositionSolver getSolver() {
        return new Solver();
    }

    /**
     * Copy matrix data and check symmetry.
     *
     * @param matrix input matrix
     * @param relativeSymmetryThreshold symmetry threshold
     * @return symmetric copy of the matrix data
     */
    private static double[][] copyAndCheckSymmetry(final RealMatrix matrix,
                                                   final double relativeSymmetryThreshold) {
        final int order = matrix.getRowDimension();
        final double[][] data = matrix.getData();

        for (int i = 0; i < order; ++i) {
            final double[] rowI = data[i];
            for (int j = i + 1; j < order; ++j) {
                final double[] rowJ = data[j];
                final double aij = rowI[j];
                final double aji = rowJ[i];
                final double maxDelta =
                        relativeSymmetryThreshold * FastMath.max(FastMath.abs(aij), FastMath.abs(aji));

                if (FastMath.abs(aij - aji) > maxDelta) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SYMMETRIC_MATRIX,
                                                           i, j, relativeSymmetryThreshold);
                }

                rowJ[i] = aij;
            }
        }

        return data;
    }

    /**
     * Perform factorization with adaptive diagonal regularization.
     *
     * @param original original symmetric matrix
     * @param positivityThreshold positivity threshold
     * @param diagonalGrowthFactor growth factor for the shift
     * @param maxSteps maximum number of regularization attempts
     * @return factorization result
     */
    private static FactorizationResult factorizeWithRegularization(final double[][] original,
                                                                   final double positivityThreshold,
                                                                   final double diagonalGrowthFactor,
                                                                   final int maxSteps) {
        final int n = original.length;
        final double[] originalDiag = new double[n];

        for (int i = 0; i < n; ++i) {
            originalDiag[i] = original[i][i];
        }

        double shift = 0.0;

        // Initial estimate of the shift.
        for (int i = 0; i < n; ++i) {
            shift = FastMath.max(shift, positivityThreshold - originalDiag[i]);
            for (int j = i + 1; j < n; ++j) {
                final double aij = original[i][j];
                final double ajj = originalDiag[j];
                double candidate = -FastMath.min(originalDiag[i], ajj);
                final double scale = FastMath.abs(originalDiag[i] - ajj) + FastMath.abs(aij);
                if (scale > 0.0) {
                    candidate += (aij * aij) / scale;
                }
                shift = FastMath.max(shift, candidate);
            }
        }

        // If no shift appears necessary, try the unmodified matrix first.
        if (shift <= 0.0) {
            final AttemptResult plain = attemptFactorization(original, 0.0, positivityThreshold);
            if (plain.success) {
                return new FactorizationResult(plain.lTData, plain.regularizedMatrix, 0.0);
            }
            shift = plain.additionalShift;
        }

        double currentShift = shift;
        for (int step = 0; step < maxSteps; ++step) {
            currentShift *= diagonalGrowthFactor;

            final AttemptResult attempt = attemptFactorization(original, currentShift, positivityThreshold);
            if (attempt.success) {
                return new FactorizationResult(attempt.lTData,
                                               attempt.regularizedMatrix,
                                               currentShift);
            }

            currentShift += attempt.additionalShift;
        }

        throw new MathIllegalArgumentException(LocalizedCoreFormats.NOT_POSITIVE_DEFINITE_MATRIX);
    }

    /**
     * Attempt factorization of {@code original + shift * I}.
     *
     * @param original original matrix
     * @param shift diagonal shift
     * @param positivityThreshold positivity threshold
     * @return attempt result
     */
    private static AttemptResult attemptFactorization(final double[][] original,
                                                      final double shift,
                                                      final double positivityThreshold) {
        final int n = original.length;
        final double[][] r = new double[n][n];
        final double[][] regularized = new double[n][n];

        for (int i = 0; i < n; ++i) {
            System.arraycopy(original[i], 0, regularized[i], 0, n);
            regularized[i][i] += shift;
        }

        for (int j = 0; j < n; ++j) {
            double pivotCandidate = 0.0;

            for (int i = 0; i <= j; ++i) {
                pivotCandidate = regularized[i][j];

                for (int k = 0; k < i; ++k) {
                    pivotCandidate -= r[k][i] * r[k][j];
                }

                if (i < j) {
                    if (FastMath.abs(r[i][i]) <= positivityThreshold) {
                        final double additional = computeAdditionalShift(r, i, pivotCandidate, positivityThreshold);
                        return AttemptResult.failure(additional);
                    }
                    r[i][j] = pivotCandidate / r[i][i];
                } else {
                    if (pivotCandidate < positivityThreshold) {
                        final double additional = computeAdditionalShift(r, j, pivotCandidate, positivityThreshold);
                        return AttemptResult.failure(additional);
                    }
                    r[j][j] = FastMath.sqrt(pivotCandidate);
                }
            }
        }

        return AttemptResult.success(r, regularized);
    }

    /**
     * Compute the extra diagonal correction when a pivot fails.
     *
     * @param r current upper-triangular factor
     * @param pivotIndex failed pivot index
     * @param pivotCandidate failed pivot value
     * @param positivityThreshold positivity threshold
     * @return additional shift
     */
    private static double computeAdditionalShift(final double[][] r,
                                                 final int pivotIndex,
                                                 final double pivotCandidate,
                                                 final double positivityThreshold) {
        final double sumx = computeRecurrenceNorm(r, pivotIndex);
        if (sumx <= 0.0 || Double.isNaN(sumx) || Double.isInfinite(sumx)) {
            return positivityThreshold - pivotCandidate;
        }
        return positivityThreshold - pivotCandidate / sumx;
    }

    /**
     * Compute the norm-like recurrence used during pivot correction.
     *
     * @param r current upper-triangular factor
     * @param pivotIndex pivot index
     * @return recurrence norm
     */
    private static double computeRecurrenceNorm(final double[][] r, final int pivotIndex) {
        if (pivotIndex == 0) {
            return 1.0;
        }

        final double[] work = new double[pivotIndex + 1];
        work[pivotIndex] = 1.0;
        double sumx = 1.0;

        int k = pivotIndex;
        while (k > 0) {
            double sum = 0.0;
            for (int i = k; i <= pivotIndex; ++i) {
                sum -= r[k - 1][i] * work[i];
            }

            final double diag = r[k - 1][k - 1];
            if (FastMath.abs(diag) < FastMath.ulp(1.0)) {
                return Double.POSITIVE_INFINITY;
            }

            work[k - 1] = sum / diag;
            sumx += work[k - 1] * work[k - 1];
            --k;
        }

        return sumx;
    }

    /**
     * Internal factorization result.
     */
    private static final class FactorizationResult {

        /** LT matrix data. */
        private final double[][] lTData;

        /** Regularized matrix data. */
        private final double[][] regularizedMatrix;

        /** Final diagonal shift. */
        private final double diagonalShift;

        FactorizationResult(final double[][] lTData,
                            final double[][] regularizedMatrix,
                            final double diagonalShift) {
            this.lTData = lTData;
            this.regularizedMatrix = regularizedMatrix;
            this.diagonalShift = diagonalShift;
        }
    }

    /**
     * Result of one factorization attempt.
     */
    private static final class AttemptResult {

        /** Whether the attempt succeeded. */
        private final boolean success;

        /** LT matrix data if successful. */
        private final double[][] lTData;

        /** Regularized matrix data if successful. */
        private final double[][] regularizedMatrix;

        /** Additional shift suggested after failure. */
        private final double additionalShift;

        private AttemptResult(final boolean success,
                              final double[][] lTData,
                              final double[][] regularizedMatrix,
                              final double additionalShift) {
            this.success = success;
            this.lTData = lTData;
            this.regularizedMatrix = regularizedMatrix;
            this.additionalShift = additionalShift;
        }

        static AttemptResult success(final double[][] lTData,
                                     final double[][] regularizedMatrix) {
            return new AttemptResult(true, lTData, regularizedMatrix, 0.0);
        }

        static AttemptResult failure(final double additionalShift) {
            return new AttemptResult(false, null, null, additionalShift);
        }
    }

    /**
     * Specialized solver.
     */
    private class Solver implements DecompositionSolver {

        /** {@inheritDoc} */
        @Override
        public boolean isNonSingular() {
            return true;
        }

        /** {@inheritDoc} */
        @Override
        public RealVector solve(final RealVector b) {
            final int m = lTData.length;
            if (b.getDimension() != m) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                       b.getDimension(), m);
            }

            final double[] x = b.toArray();

            // Solve L Y = b
            for (int j = 0; j < m; j++) {
                final double[] lJ = lTData[j];
                x[j] /= lJ[j];
                final double xJ = x[j];
                for (int i = j + 1; i < m; i++) {
                    x[i] -= xJ * lJ[i];
                }
            }

            // Solve L^T X = Y
            for (int j = m - 1; j >= 0; j--) {
                x[j] /= lTData[j][j];
                final double xJ = x[j];
                for (int i = 0; i < j; i++) {
                    x[i] -= xJ * lTData[i][j];
                }
            }

            return new ArrayRealVector(x, false);
        }

        /** {@inheritDoc} */
        @Override
        public RealMatrix solve(final RealMatrix b) {
            final int m = lTData.length;
            if (b.getRowDimension() != m) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                       b.getRowDimension(), m);
            }

            final int nColB = b.getColumnDimension();
            final double[][] x = b.getData();

            // Solve L Y = b
            for (int j = 0; j < m; j++) {
                final double[] lJ = lTData[j];
                final double lJJ = lJ[j];
                final double[] xJ = x[j];
                for (int k = 0; k < nColB; ++k) {
                    xJ[k] /= lJJ;
                }
                for (int i = j + 1; i < m; i++) {
                    final double[] xI = x[i];
                    final double lJI = lJ[i];
                    for (int k = 0; k < nColB; ++k) {
                        xI[k] -= xJ[k] * lJI;
                    }
                }
            }

            // Solve L^T X = Y
            for (int j = m - 1; j >= 0; j--) {
                final double lJJ = lTData[j][j];
                final double[] xJ = x[j];
                for (int k = 0; k < nColB; ++k) {
                    xJ[k] /= lJJ;
                }
                for (int i = 0; i < j; i++) {
                    final double[] xI = x[i];
                    final double lIJ = lTData[i][j];
                    for (int k = 0; k < nColB; ++k) {
                        xI[k] -= xJ[k] * lIJ;
                    }
                }
            }

            return new Array2DRowRealMatrix(x, false);
        }

        /**
         * Get the inverse of the regularized matrix.
         *
         * @return inverse matrix
         */
        @Override
        public RealMatrix getInverse() {
            return solve(MatrixUtils.createRealIdentityMatrix(lTData.length));
        }

        /** {@inheritDoc} */
        @Override
        public int getRowDimension() {
            return lTData.length;
        }

        /** {@inheritDoc} */
        @Override
        public int getColumnDimension() {
            return lTData[0].length;
        }

        /**
     * Get the final diagonal shift added to the original matrix.
     *
     * @return diagonal shift
     */
    public double getDiagonalShift() {
        return diagonalShift;
    }
    }
}
